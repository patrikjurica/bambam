use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::Record; // Required trait for .alignment_start(), etc.

use crate::distance::edit_distance;
use crate::kmer::{decode_kmer, encode_kmer};
use crate::types::RareKmer;

/// Evaluates a single sequence against the k-mer library.
/// Returns a tuple: (passes_threshold, has_zero_kmers)
#[allow(clippy::too_many_arguments)]
pub(crate) fn evaluate_alignment(
    record: &bam::Record,
    active_seq: &[u8],
    chrom_kmers: &[RareKmer],
    ref_to_query_buffer: &mut Vec<usize>,
    kmer_len: usize,
    min_pct: f64,
    min_count: usize,
    use_base_seq: bool,
) -> (bool, bool) {
    let ref_start = match record.alignment_start() {
        Some(Ok(pos)) => usize::from(pos) - 1,
        _ => return (false, false),
    };

    let ref_end = match record.alignment_end() {
        Some(Ok(pos)) => usize::from(pos) - 1,
        _ => return (false, false),
    };

    let mut idx = match chrom_kmers.binary_search_by_key(&ref_start, |k| k.start) {
        Ok(i) => i,
        Err(i) => i,
    };

    let mut kmers_in_range = Vec::new();
    while idx < chrom_kmers.len() {
        let kmer = &chrom_kmers[idx];
        if kmer.start >= ref_end { break; }
        if kmer.end <= ref_end { kmers_in_range.push(kmer); }
        idx += 1;
    }

    let has_zero = kmers_in_range.is_empty();

    if kmers_in_range.len() < min_count {
        return (false, has_zero);
    }

    let ref_len = ref_end.saturating_sub(ref_start);
    ref_to_query_buffer.clear();
    ref_to_query_buffer.resize(ref_len + 1, usize::MAX);

    let mut curr_ref = ref_start;
    let mut curr_query = 0;

    // Build CIGAR Map
    if let Ok(cigar) = record.cigar().iter().collect::<Result<Vec<_>, _>>() {
        for op in cigar {
            let len = op.len();
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for _ in 0..len {
                        if curr_ref >= ref_start && curr_ref <= ref_end {
                            ref_to_query_buffer[curr_ref - ref_start] = curr_query;
                        }
                        curr_ref += 1;
                        curr_query += 1;
                    }
                }
                Kind::Insertion | Kind::SoftClip => { curr_query += len; }
                Kind::Deletion | Kind::Skip => { curr_ref += len; }
                Kind::HardClip => { if use_base_seq { curr_query += len; } }
                Kind::Pad => {}
            }
        }
    }

    let mut valid_kmer_count = 0;

    // Compliance Check
    for kmer in &kmers_in_range {
        let mut start_idx = kmer.start.saturating_sub(ref_start);
        let mut end_idx = (kmer.end - 1).saturating_sub(ref_start);

        if start_idx >= ref_to_query_buffer.len() || end_idx >= ref_to_query_buffer.len() {
            continue;
        }

        let mut q_start = ref_to_query_buffer[start_idx];
        while q_start == usize::MAX && start_idx < end_idx {
            start_idx += 1;
            q_start = ref_to_query_buffer[start_idx];
        }

        let mut q_end_inclusive = ref_to_query_buffer[end_idx];
        while q_end_inclusive == usize::MAX && end_idx > start_idx {
            end_idx -= 1;
            q_end_inclusive = ref_to_query_buffer[end_idx];
        }

        if q_start == usize::MAX || q_end_inclusive == usize::MAX {
            continue;
        }

        if q_start >= active_seq.len() || q_end_inclusive >= active_seq.len() {
            continue; // strict bounds safety check
        }

        let actual_len = (q_end_inclusive - q_start) + 1;
        let length_diff = actual_len.abs_diff(kmer_len);

        if length_diff > kmer.local_tolerance {
            continue;
        }

        let slice = &active_seq[q_start..=q_end_inclusive];

        if length_diff == 0 {
            if let Some(actual_val) = encode_kmer(slice) {
                if actual_val == kmer.val {
                    valid_kmer_count += 1;
                    continue;
                }
            }
        }

        if kmer.local_tolerance > 0 {
            let expected_seq = decode_kmer(kmer.val, kmer_len);
            if edit_distance(expected_seq.as_bytes(), slice) <= kmer.local_tolerance {
                valid_kmer_count += 1;
            }
        }
    }

    let total_expected = kmers_in_range.len();
    let survival_rate = if total_expected > 0 {
        (valid_kmer_count as f64 / total_expected as f64) * 100.0
    } else {
        0.0
    };

    let passes = total_expected >= min_count && (survival_rate >= min_pct || min_count == 0);
    (passes, has_zero)
}
