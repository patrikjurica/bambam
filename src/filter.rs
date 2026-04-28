use anyhow::{Context, Result};
use noodles::bam;
use noodles::bgzf;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::Record;
use rustc_hash::FxHashMap;
use std::fs::File;

use crate::distance::edit_distance;
use crate::kmer::{decode_kmer, encode_kmer};
use crate::types::KmerLibrary;


// TODO: get back to this


/// filters a BAM file based on absolute positional compliance of rare k-mers.
pub fn filter_bam(
    bam_in_path: &str,
    bam_out_path: &str,
    expected_kmers: &KmerLibrary,
    kmer_len: usize,
    min_pct: f64,
    min_count: usize,
) -> Result<()> {
    let mut reader = bam::io::Reader::new(
        bgzf::Reader::new(File::open(bam_in_path).context("Failed to open input BAM")?)
    );
    let header = reader.read_header().context("Failed to read BAM header")?;

    let mut writer = bam::io::Writer::new(
        bgzf::Writer::new(File::create(bam_out_path).context("Failed to create output BAM")?)
    );
    writer.write_header(&header).context("Failed to write BAM header")?;

    let mut zero_kmers = 0;
    let mut unmapped = 0;

    for result in reader.records() {
        let record = result.context("Failed to read BAM record")?;

        if record.flags().is_unmapped() || record.sequence().is_empty() {
            unmapped += 1;
            continue;
        }

        // extract the chromosome name using the reference ID
        let ref_id = match record.reference_sequence_id() {
            Some(Ok(id)) => id,
            _ => continue,
        };

        let chrom_name = match header.reference_sequences().get_index(ref_id) {
            Some((name_key, _metadata_value)) => {
                // the name is the key (.0)
                String::from_utf8_lossy(name_key.as_ref()).into_owned()
            }
            None => continue, // skip the read if the header is missing this ID
        };

        let chrom_kmers = match expected_kmers.get(&chrom_name) {
            Some(kmers) => kmers,
            None => continue,
        };

        // convert 1-based noodles positions to 0-based usize
        let ref_start = match record.alignment_start() {
            Some(Ok(pos)) => usize::from(pos) - 1,
            _ => continue, // Skip if unmapped or corrupted
        };

        let ref_end = match record.alignment_end() {
            Some(Ok(pos)) => usize::from(pos) - 1,
            _ => continue,
        };

        // 1. BINARY SEARCH
        // find the index where our read starts in the k-mer library
        let mut idx = match chrom_kmers.binary_search_by_key(&ref_start, |k| k.start) {
            Ok(i) => i,
            Err(i) => i,
        };

        let mut kmers_in_range = Vec::new();
        while idx < chrom_kmers.len() {
            let kmer = &chrom_kmers[idx];
            if kmer.start >= ref_end {
                break;
            }
            if kmer.end <= ref_end {
                kmers_in_range.push(kmer);
            }
            idx += 1;
        }

        if kmers_in_range.is_empty() {
            zero_kmers += 1;
        }

        // condition 1: minimum count threshold
        if kmers_in_range.len() < min_count {
            continue;
        }

        // 2. BUILD CIGAR MAPPING
        // (replace pysam.get_aligned_pairs)
        // maps 0-based reference coordinates to 0-based query sequence coordinates
        let mut ref_to_query: FxHashMap<usize, usize> = FxHashMap::default();
        let mut curr_ref = ref_start;
        let mut curr_query = 0;

        for op_result in record.cigar().iter() {
            let op = op_result?;
            let len = op.len();

            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for _ in 0..len {
                        ref_to_query.insert(curr_ref, curr_query);
                        curr_ref += 1;
                        curr_query += 1;
                    }
                }
                Kind::Insertion | Kind::SoftClip => {
                    curr_query += len; // read moves forward, reference stays still
                }
                Kind::Deletion | Kind::Skip => {
                    curr_ref += len; // reference moves forward, read stays still
                }
                Kind::HardClip | Kind::Pad => {} // doesn't consume sequence or reference
            }
        }

        // 3. COMPLIANCE CHECK
        let mut valid_kmer_count = 0;
        let read_seq = record.sequence();

        for kmer in &kmers_in_range {
            // check boundaries
            let q_start = match ref_to_query.get(&kmer.start) {
                Some(&pos) => pos,
                None => continue,
            };
            let q_end_inclusive = match ref_to_query.get(&(kmer.end - 1)) {
                Some(&pos) => pos,
                None => continue,
            };

            let actual_len = (q_end_inclusive - q_start) + 1;
            let length_diff = actual_len.abs_diff(kmer_len);

            if length_diff > kmer.local_tolerance {
                continue;
            }

            // extract the byte slice of the sequence
            // noodles sequence extraction requires mapping via `.iter().map()`
            let slice: Vec<u8> = read_seq.iter()
                .skip(q_start)
                .take(actual_len)
                .map(u8::from) // Convert noodles Base to u8
                .collect();

            // path 1: direct
            if length_diff == 0 {
                if let Some(actual_val) = encode_kmer(&slice) {
                    if actual_val == kmer.val {
                        valid_kmer_count += 1;
                        continue;
                    }
                }
            }

            // path 2: distance
            if kmer.local_tolerance > 0 {
                let expected_seq = decode_kmer(kmer.val, kmer_len);
                if edit_distance(expected_seq.as_bytes(), &slice) <= kmer.local_tolerance {
                    valid_kmer_count += 1;
                }
            }
        }

        // 4. SURVIVAL RATE CALCULATION
        let total_expected = kmers_in_range.len();
        let survival_rate = if total_expected > 0 {
            (valid_kmer_count as f64 / total_expected as f64) * 100.0
        } else {
            0.0
        };

        if total_expected >= min_count && (survival_rate >= min_pct || min_count == 0) {
            writer.write_record(&header, &record).context("Failed to write record")?;
        }
    }

    println!("Number of reads with zero rare kmers: {}", zero_kmers);
    println!("Number of unmapped reads: {}", unmapped);

    Ok(())
}