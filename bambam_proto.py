import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass
import pysam
import time


# DNA base encoding
A = 0b00  # 0
C = 0b01  # 1
G = 0b10  # 2
T = 0b11  # 3


@dataclass(frozen=True, order=True)
class Location:
    chrom: str
    start: int
    end: int


def encode_kmer(seq: str) -> int | None:
    """Encodes a DNA sequence into a 2-bit integer representation."""
    kmer = 0
    for char in seq.upper():
        kmer <<= 2
        if char == 'A':
            kmer |= A
        elif char == 'C':
            kmer |= C
        elif char == 'G':
            kmer |= G
        elif char == 'T':
            kmer |= T
        else:
            return None  # Equivalent to returning Option::None for 'N'
    return kmer


def build_rare_kmers(ref_file: str, kmer_len: int, min_count: int, max_count: int) -> dict[
    str, list[tuple[int, int, int]]]:
    """
    Builds the rare k-mer library directly from the reference FASTA using a rolling bitmask.
    Pass 1: Count frequencies.
    Pass 2: Record coordinates of k-mers within the [min_count, max_count] range.
    """

    # Quick lookup for bases
    base_to_int = {'A': A, 'C': C, 'G': G, 'T': T}

    # this mask keeps our integer exactly at the size of our k-mer
    # e.g., for k=3, 2*k = 6 bits. 1 << 6 is 1000000 (binary). minus 1 is 111111.
    mask = (1 << (2 * kmer_len)) - 1

    kmer_counts = defaultdict(int)

    print("Pass 1: Counting k-mers in the reference genome (this MAY take a moment)...")
    with pysam.FastaFile(ref_file) as fasta:
        for chrom in fasta.references:
            seq = fasta.fetch(chrom).upper()
            kmer_val = 0
            valid_bases = 0

            for char in seq:
                if char in base_to_int:
                    # Shift left by 2, apply mask to trim overflow, add new base
                    kmer_val = ((kmer_val << 2) & mask) | base_to_int[char]
                    valid_bases += 1

                    # Once we've seen enough valid consecutive bases, start counting
                    if valid_bases >= kmer_len:
                        kmer_counts[kmer_val] += 1
                else:
                    # We hit an 'N' or unknown base. Reset the counter.
                    kmer_val = 0
                    valid_bases = 0

    print(f"Total unique k-mers found: {len(kmer_counts)}")

    # Filter to keep only the "rare" ones based on min and max
    valid_kmers = {k for k, v in kmer_counts.items() if min_count <= v <= max_count}
    print(f"K-mers within frequency range ({min_count}-{max_count}): {len(valid_kmers)}")

    # Free up memory before Pass 2
    del kmer_counts

    expected_kmers = defaultdict(list)

    print("Pass 2: Locating valid rare k-mers in the reference...")
    with pysam.FastaFile(ref_file) as fasta:
        for chrom in fasta.references:
            seq = fasta.fetch(chrom).upper()
            kmer_val = 0
            valid_bases = 0

            for i, char in enumerate(seq):
                if char in base_to_int:
                    kmer_val = ((kmer_val << 2) & mask) | base_to_int[char]
                    valid_bases += 1

                    if valid_bases >= kmer_len:
                        if kmer_val in valid_kmers:
                            # Calculate exactly where this k-mer started
                            start = i - kmer_len + 1
                            end = i + 1  # Exclusive end for standard BED/0-based indexing
                            expected_kmers[chrom].append((start, end, kmer_val))
                else:
                    kmer_val = 0
                    valid_bases = 0

    return dict(expected_kmers)


def decode_kmer(val: int, length: int) -> str:
    """Decodes a 2-bit integer back into a DNA string."""
    chars = []
    for _ in range(length):
        base_val = val & 0b11
        if base_val == A:
            chars.append('A')
        elif base_val == C:
            chars.append('C')
        elif base_val == G:
            chars.append('G')
        elif base_val == T:
            chars.append('T')
        val >>= 2
    return "".join(reversed(chars))


def edit_distance(s1: str, s2: str) -> int:
    """Calculates the Levenshtein distance between two strings."""
    if len(s1) < len(s2):
        return edit_distance(s2, s1)
    if len(s2) == 0:
        return len(s1)
    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]


def apply_dynamic_threshold_tolerance(expected_kmers: dict[str, list[tuple[int, int, int]]],
                            base_tolerance: int,
                            distance_threshold: int = 5000) -> dict[str, list[tuple[int, int, int, int]]]:
    """
    Evaluates the density of k-mers and assigns a specific tolerance to each.
    Returns: dict of chrom -> list of (start, end, kmer_val, assigned_tolerance)
    """
    dynamic_kmers = defaultdict(list)

    for chrom, kmers in expected_kmers.items():
        n = len(kmers)
        for i in range(n):
            start, end, kmer_val = kmers[i]

            # Calculate distance to left and right neighbors (default to 0 if at the very edge)
            dist_left = (start - kmers[i - 1][1]) if i > 0 else 0
            dist_right = (kmers[i + 1][0] - end) if i < n - 1 else 0

            # Use the maximum distance to a neighbor to determine isolation
            isolation_distance = dist_left + dist_right

            if isolation_distance > distance_threshold:
                assigned_tolerance = 0
            else:
                assigned_tolerance = base_tolerance

            dynamic_kmers[chrom].append((start, end, kmer_val, assigned_tolerance))

    return dict(dynamic_kmers)


def filter_bam(bam_file: str, output_bam: str, expected_kmers: dict[str, list[tuple[int, int, int, int]]], kmer_len: int,
               min_pct: float, min_count: int) -> None:
    """
    Filters a BAM file based on absolute positional compliance of rare k-mers.
    Includes a dual-path check to tolerate indels/SNPs up to a specified edit distance.
    """
    import pysam
    import bisect

    zero_kmers = 0
    unmapped = 0

    with pysam.AlignmentFile(bam_file, "rb") as bam_in:
        with pysam.AlignmentFile(output_bam, "wb", header=bam_in.header) as bam_out:

            for read in bam_in.fetch(until_eof=True):
                if read.is_unmapped or read.query_sequence is None:
                    unmapped += 1
                    continue

                chrom = read.reference_name
                if chrom not in expected_kmers:
                    continue

                ref_start = read.reference_start
                ref_end = read.reference_end
                chrom_kmers = expected_kmers[chrom]

                # 1. BINARY SEARCH: Find expected k-mers
                idx = bisect.bisect_left(chrom_kmers, (ref_start, 0, 0, 0))

                kmers_in_range = []
                while idx < len(chrom_kmers):
                    k_start, k_end, k_val, local_tolerance = chrom_kmers[idx]

                    if k_start >= ref_end:
                        break

                    if k_end <= ref_end:
                        kmers_in_range.append((k_start, k_end, k_val, local_tolerance))

                    idx += 1

                if len(kmers_in_range) == 0:
                    zero_kmers += 1
                # Condition 1: Minimum count threshold
                if len(kmers_in_range) < min_count:
                    continue

                # 2. COMPLIANCE CHECK
                aligned_pairs = read.get_aligned_pairs(matches_only=True)
                ref_to_query = {ref_pos: query_pos for query_pos, ref_pos in aligned_pairs}

                valid_kmer_count = 0

                for k_start, k_end, expected_val, local_tolerance in kmers_in_range:
                    if k_start not in ref_to_query or (k_end - 1) not in ref_to_query:
                        continue  # Indel exactly at the boundary. Skip.

                    q_start = ref_to_query[k_start]
                    q_end_inclusive = ref_to_query[k_end - 1]

                    # Check B: Is the length difference within our allowed tolerance?
                    actual_len = (q_end_inclusive - q_start) + 1
                    length_diff = abs(actual_len - kmer_len)

                    if length_diff > local_tolerance:
                        continue  # Indel changed the length too much. Skip.

                    # Check C: Perfect or Tolerated
                    read_seq_slice = read.query_sequence[q_start: q_end_inclusive + 1]

                    # Path 1: The Fast Path (Exact Length & Perfect Match)
                    if length_diff == 0:
                        actual_val = encode_kmer(read_seq_slice)
                        if actual_val == expected_val:
                            valid_kmer_count += 1
                            continue  # Match found

                    # Path 2: The Fallback Path (Indels or SNPs)
                    if local_tolerance > 0:
                        expected_seq = decode_kmer(expected_val, kmer_len)
                        if edit_distance(expected_seq, read_seq_slice) <= local_tolerance:
                            valid_kmer_count += 1

                total_expected = len(kmers_in_range)
                survival_rate = (valid_kmer_count / total_expected) * 100 if total_expected > 0 else 0

                # Condition 2: The Percentage Threshold
                if total_expected >= min_count and (survival_rate >= min_pct or min_count == 0):
                    bam_out.write(read)

    print(f"Number of reads with zero rare kmers: {zero_kmers}")
    print(f"Number of unmapped reads: {unmapped}")

def main():
    start_time = time.time()
    parser = argparse.ArgumentParser(description="Filter BAM files based on rare k-mers.")

    # Positional arguments
    parser.add_argument("input_bam_file", help="The path to the input BAM file to be filtered.")
    parser.add_argument("ref_file", help="The file containing the rare k-mers as well as the rest of the genome")

    # Optional arguments
    parser.add_argument("-o", "--output", default="./out.bam",
                        help="The path where the filtered output BAM file will be written.")
    parser.add_argument("--min", type=int, default=5, help="Minimal k-mer count to be considered for filtering.")
    parser.add_argument("--max", type=int, default=10, help="Maximal k-mer count to be considered for filtering.")
    parser.add_argument("-i", "--inf", type=int, default=1,
                        help="Minimal count of k-mers in a read for the read to be written in the output file.")
    parser.add_argument("-p", "--pct", type=float, default=100.0, help="Minimum percentage of intact expected rare k-mers required to keep a read.")
    parser.add_argument("-l", "--len", type=int, default=31, help="Length of the unique k-mer (higher yields more values)")
    parser.add_argument("-t", "--tolerance", type=int, default=0,
                        help="Maximum allowed edit distance (or base tolerance if --dyn-tol is used) for a k-mer to still be considered a match.")
    parser.add_argument("--dyn-tol", action="store_true",
                        help="Enable dynamic, density-aware tolerance based on k-mer isolation.")
    parser.add_argument("--bed", type=str, help="Optional: Path to output a BED file of the rare k-mer coordinates.")

    args = parser.parse_args()

    print("Building the k-mer dictionary...")
    kmers = build_rare_kmers(args.ref_file, args.len, args.min, args.max)
    print("Loading coordinates and encoding kmers...")

    if args.dyn_tol:
        print("Applying dynamic threshold tolerance...")
        kmers = apply_dynamic_threshold_tolerance(kmers, args.tolerance, distance_threshold=5000)
    else:
        print(f"Applying static tolerance of {args.tolerance} to all k-mers...")
        # Convert 3-tuples to 4-tuples using the static tolerance
        for chrom in kmers:
            kmers[chrom] = [(start, end, val, args.tolerance) for start, end, val in kmers[chrom]]

    print(f"Successfully loaded {sum(len(v) for v in kmers.values())} unique rare k-mers.")

    # Output bed file in case one wants to know where the rare k-mers are located

    if args.bed:
        print(f"Exporting k-mer coordinates to {args.bed}...")
        with open(args.bed, 'w') as bed_out:
            for chrom, coords_list in kmers.items():
                for start, end, kmer_val, _ in coords_list:
                    # BED format is: chrom \t start \t end \t name
                    kmer_seq = decode_kmer(kmer_val, args.len)
                    # using the kmer itself as the name
                    bed_out.write(f"{chrom}\t{start}\t{end}\t{kmer_seq}\n")
        print("BED export complete.")

    print(f"Processing BAM file: {args.input_bam_file}...")
    try:
        # If the user types --min 50, we pass 50.0
        filter_bam(args.input_bam_file, args.output, kmers, args.len, args.pct, args.inf)
        print(f"BAM filtering complete. Saved to {args.output}")
    except Exception as e:
        print(e, file=sys.stderr)

    end_time = time.time()
    print(f"Took {end_time - start_time} seconds to run. ")


if __name__ == "__main__":
    main()