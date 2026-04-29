import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass
import pysam
import bisect
import functools
import operator


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


def load_bed(bed_file: str) -> list[Location]:
    """Loads rare k-mer coordinates from a BED file."""
    coords = []
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split('\t')
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    coords.append(Location(chrom, start, end))
    except Exception as e:
        print(f"Mistake reading BED file: {e}", file=sys.stderr)
        sys.exit(1)

    coords.sort()  # Sorts by chrom, then start, then end
    return coords


def load_expected_kmers(bed_file: str, ref_file: str) -> dict[str, list[tuple[int, int, int]]]:
    """Builds a dictionary of chrom -> sorted list of (start, end, kmer_val)"""
    coords = load_bed(bed_file)
    expected_kmers = defaultdict(list)

    with pysam.FastaFile(ref_file) as fasta:
        for loc in coords:
            if loc.chrom not in fasta.references:
                continue
            seq_slice = fasta.fetch(loc.chrom, loc.start, loc.end)
            kmer_val = encode_kmer(seq_slice)

            if kmer_val is not None:
                expected_kmers[loc.chrom].append((loc.start, loc.end, kmer_val))

    # The BED file was already sorted, so these lists are naturally sorted by start position!
    return dict(expected_kmers)


def filter_bam(bam_file: str, output_bam: str, expected_kmers: dict[str, list[tuple[int, int, int]]], kmer_len: int, min_pct: float):
    """
    Filters a BAM file based on absolute positional compliance of rare k-mers.
    Uses a relaxed threshold (>= 1 match) to account for Oxford Nanopore indel rates.
    """
    import pysam
    import bisect

    with pysam.AlignmentFile(bam_file, "rb") as bam_in:
        # Open the output BAM file, copying the header from the input
        with pysam.AlignmentFile(output_bam, "wb", header=bam_in.header) as bam_out:

            for read in bam_in.fetch(until_eof=True):
                # Throw away unmapped reads or reads with no sequence
                if read.is_unmapped or read.query_sequence is None:
                    continue

                chrom = read.reference_name
                if chrom not in expected_kmers:
                    continue  # No rare k-mers on this chromosome at all -> throw away

                ref_start = read.reference_start
                ref_end = read.reference_end
                chrom_kmers = expected_kmers[chrom]

                # 1. BINARY SEARCH: Find expected k-mers that fall in this read's mapping area
                idx = bisect.bisect_left(chrom_kmers, (ref_start, 0, 0))

                kmers_in_range = []
                while idx < len(chrom_kmers):
                    k_start, k_end, k_val = chrom_kmers[idx]

                    if k_start >= ref_end:
                        break  # We have passed the end of the read

                    if k_end <= ref_end:
                        kmers_in_range.append((k_start, k_end, k_val))

                    idx += 1

                # Condition 1: If there are NO rare k-mers supposed to be here -> throw away
                if not kmers_in_range:
                    continue

                # 2. COMPLIANCE CHECK: Map reference coordinates to read coordinates
                aligned_pairs = read.get_aligned_pairs(matches_only=True)
                ref_to_query = {ref_pos: query_pos for query_pos, ref_pos in aligned_pairs}

                valid_kmer_count = 0

                for k_start, k_end, expected_val in kmers_in_range:
                    # Check A: Does the alignment cleanly cover the start and end of the k-mer?
                    if k_start not in ref_to_query or (k_end - 1) not in ref_to_query:
                        continue  # Indel at the boundary. Skip this k-mer.

                    q_start = ref_to_query[k_start]
                    q_end_inclusive = ref_to_query[k_end - 1]

                    # Check B: Did an insertion/deletion happen *inside* the k-mer?
                    if (q_end_inclusive - q_start) != (kmer_len - 1):
                        continue  # Indel ruined the internal length. Skip this k-mer.

                    # Check C: Does the actual sequence match our expected bits perfectly?
                    read_seq_slice = read.query_sequence[q_start: q_end_inclusive + 1]
                    actual_val = encode_kmer(read_seq_slice)

                    if actual_val == expected_val:
                        valid_kmer_count += 1  # We found a perfect match!

                total_expected = len(kmers_in_range)

                # Let's print the stats for the first few reads to see what's happening
                survival_rate = (valid_kmer_count / total_expected) * 100 if total_expected > 0 else 0

                # Condition 2: The Percentage Threshold
                # Let's demand that at least 5% of the expected k-mers are perfectly intact
                if total_expected > 0 and survival_rate >= min_pct:
                    bam_out.write(read)

def main():
    parser = argparse.ArgumentParser(description="Filter BAM files based on rare k-mers.")

    # Positional arguments
    parser.add_argument("input_bam_file", help="The path to the input BAM file to be filtered.")
    parser.add_argument("ref_file", help="The file containing the rare k-mers as well as the rest of the genome")
    parser.add_argument("bed_file", help="The file containing the locations of the rare k-mers in the reference")

    # Optional arguments
    parser.add_argument("-o", "--output", default="./out.bam",
                        help="The path where the filtered output BAM file will be written.")
    parser.add_argument("--min", type=int, default=5, help="Minimal k-mer count to be considered for filtering.")
    parser.add_argument("--max", type=int, default=10, help="Maximal k-mer count to be considered for filtering.")
    parser.add_argument("-p", "--pct", type=float, default=5.0, help="Minimum percentage of intact expected rare k-mers required to keep a read.")
    parser.add_argument("--len", type=int, default=31, help="Length of the unique k-mer (lower yields more values)")
    parser.add_argument("-w", "--weight", type=int, default=10,
                        help="Significance weight of each k-mer for filtering criteria.")

    args = parser.parse_args()

    print("Loading coordinates and encoding kmers...")
    kmers = load_expected_kmers(args.bed_file, args.ref_file)

    print(f"Successfully loaded {functools.reduce(operator.add, map(len, kmers.values()), 0)} unique rare k-mers.")

    print(f"Processing BAM file: {args.input_bam_file}...")
    try:
        # If the user types --min 50, we pass 50.0
        filter_bam(args.input_bam_file, args.output, kmers, args.len, args.pct)
        print(f"BAM filtering complete. Saved to {args.output}")
    except Exception as e:
        print(e, file=sys.stderr)


if __name__ == "__main__":
    main()