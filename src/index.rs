use anyhow::{Context, Result};
use noodles::fasta;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::BufReader;

use crate::kmer::encode_base;
use crate::types::{KmerLibrary, KmerVal, RareKmer};

/// builds the rare k-mer library from the reference FASTA using a rolling bitmask
pub fn build_rare_kmers(
    ref_file: &str,
    kmer_len: usize,
    min_count: u32,
    max_count: u32,
    base_tolerance: usize,
) -> Result<KmerLibrary> {
    // prevent overflow if asked for a k-mer of exactly 32
    let mask: KmerVal = if kmer_len == 32 {
        KmerVal::MAX
    } else {
        ((1 as KmerVal) << (2 * kmer_len)) - 1
    };

    println!("Pass 1: Counting k-mers in the reference genome...");
    let mut kmer_counts: FxHashMap<KmerVal, u32> = FxHashMap::default();

    // pass 1: count all k-mers
    let mut reader = File::open(ref_file)
        .map(BufReader::new)
        .map(fasta::io::Reader::new)?;

    for result in reader.records() {
        let record = result.context("Failed to read FASTA record")?;
        let mut kmer_val: KmerVal = 0;
        let mut valid_bases: usize = 0;

        for &byte in record.sequence().as_ref() {
            if let Some(base_bits) = encode_base(byte) {
                kmer_val = ((kmer_val << 2) & mask) | base_bits;
                valid_bases += 1;

                if valid_bases >= kmer_len {
                    *kmer_counts.entry(kmer_val).or_insert(0) += 1;
                }
            } else {
                // hit an 'N' or invalid base, reset the rolling window.
                kmer_val = 0;
                valid_bases = 0;
            }
        }
    }

    println!("Total unique k-mers found: {}", kmer_counts.len());

    // filter into a much smaller HashSet containing only the valid rare k-mers
    let mut valid_kmers: FxHashSet<KmerVal> = FxHashSet::default();
    for (kmer, count) in kmer_counts.into_iter() {
        if count >= min_count && count <= max_count {
            valid_kmers.insert(kmer);
        }
    }

    println!("K-mers within frequency range ({}-{}): {}", min_count, max_count, valid_kmers.len());

    println!("Pass 2: Locating valid rare k-mers in the reference...");
    let mut expected_kmers: KmerLibrary = FxHashMap::default();

    // re-open the reader for pass 2 to avoid storing the reference in RAM
    let mut reader = File::open(ref_file)
        .map(BufReader::new)
        .map(fasta::io::Reader::new)?;

    for result in reader.records() {
        let record = result.context("Failed to read FASTA record")?;

        // noodles returns the sequence name as a byte slice, convert it to a string for our HashMap key
        let chrom = String::from_utf8_lossy(record.name()).into_owned();
        let mut chrom_kmers = Vec::new();

        let mut kmer_val: KmerVal = 0;
        let mut valid_bases: usize = 0;

        for (i, &byte) in record.sequence().as_ref().iter().enumerate() {
            if let Some(base_bits) = encode_base(byte) {
                kmer_val = ((kmer_val << 2) & mask) | base_bits;
                valid_bases += 1;

                if valid_bases >= kmer_len && valid_kmers.contains(&kmer_val) {
                    let start = i + 1 - kmer_len;
                    let end = i + 1; // Exclusive end

                    chrom_kmers.push(RareKmer {
                        start,
                        end,
                        val: kmer_val,
                        local_tolerance: base_tolerance, // Initialize with baseline tolerance
                    });
                }
            } else {
                kmer_val = 0;
                valid_bases = 0;
            }
        }

        if !chrom_kmers.is_empty() {
            expected_kmers.insert(chrom, chrom_kmers);
        }
    }

    Ok(expected_kmers)
}


/// evaluates the density of k-mers and assigns a specific tolerance to each.
/// instead of recreating the dictionary, we mutate the existing one in-place for performance.
pub fn apply_dynamic_threshold_tolerance(
    mut library: KmerLibrary,
    base_tolerance: usize,
    distance_threshold: usize,
) -> KmerLibrary {
    for kmers in library.values_mut() {
        let n = kmers.len();
        if n == 0 {
            continue;
        }

        for i in 0..n {
            // saturating_sub prevents integer underflow (crashing) if the distances are weird
            // TODO: check
            let dist_left = if i > 0 {
                kmers[i].start.saturating_sub(kmers[i - 1].end)
            } else {
                0
            };

            let dist_right = if i < n - 1 {
                kmers[i + 1].start.saturating_sub(kmers[i].end)
            } else {
                0
            };

            let isolation_distance = dist_left + dist_right;

            if isolation_distance > distance_threshold {
                kmers[i].local_tolerance = 0;
            } else {
                kmers[i].local_tolerance = base_tolerance;
            }
        }
    }

    library
}