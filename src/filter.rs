use anyhow::{Context, Result};
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::Record;
use std::fs::File;
use std::io::{BufReader, BufWriter};

use crate::distance::edit_distance;
use crate::kmer::{decode_kmer, encode_kmer};
use crate::types::{KmerLibrary, RareKmer};

/// Filters a BAM file based on absolute positional compliance of rare k-mers.
/// Acts as a router: if `primary_only` is true, it uses a high-speed sequential stream.
/// If false, it groups alignments by QNAME (Requires `samtools sort -n`).
pub fn filter_bam(
    bam_in_path: &str,
    bam_out_path: &str,
    expected_kmers: &KmerLibrary,
    kmer_len: usize,
    min_pct: f64,
    min_count: usize,
    primary_only: bool,
) -> Result<()> {
    let mut reader = bam::io::Reader::new(BufReader::with_capacity(
        128 * 1024,
        File::open(bam_in_path).context("Failed to open input BAM")?,
    ));
    let header = reader.read_header().context("Failed to read BAM header")?;

    let mut writer = bam::io::Writer::new(BufWriter::with_capacity(
        128 * 1024,
        File::create(bam_out_path).context("Failed to create output BAM")?,
    ));
    writer.write_header(&header).context("Failed to write BAM header")?;

    if primary_only {
        println!("Running in PRIMARY ONLY mode. Processing sequentially...");
        process_primary_stream(
            &mut reader, &mut writer, &header, expected_kmers, kmer_len, min_pct, min_count
        )?;
    } else {
        println!("Running in GROUPED mode. Expecting name-sorted BAM...");
        process_grouped_stream(
            &mut reader, &mut writer, &header, expected_kmers, kmer_len, min_pct, min_count
        )?;
    }

    Ok(())
}

// ==============================================================================
// STREAM 1: HIGH-SPEED PRIMARY ONLY
// ==============================================================================

fn process_primary_stream<R: std::io::Read, W: std::io::Write>(
    reader: &mut bam::io::Reader<R>,
    writer: &mut bam::io::Writer<W>,
    header: &noodles::sam::Header,
    expected_kmers: &KmerLibrary,
    kmer_len: usize,
    min_pct: f64,
    min_count: usize,
) -> Result<()> {
    let mut zero_kmers = 0;
    let mut unmapped = 0;

    let mut ref_to_query_buffer: Vec<usize> = Vec::with_capacity(100_000);
    let mut decoded_seq_buffer: Vec<u8> = Vec::with_capacity(100_000);

    for result in reader.records() {
        let record = result.context("Failed to read BAM record")?;

        if record.flags().is_unmapped() || record.sequence().is_empty() {
            unmapped += 1;
            continue;
        }

        let chrom_name = match extract_chrom_name(&record, header) {
            Some(name) => name,
            None => continue,
        };

        let chrom_kmers = match expected_kmers.get(&chrom_name) {
            Some(kmers) => kmers,
            None => continue,
        };

        decoded_seq_buffer.clear();
        decoded_seq_buffer.extend(record.sequence().iter().map(u8::from));

        let (passes, has_zero) = evaluate_alignment(
            &record,
            &decoded_seq_buffer,
            chrom_kmers,
            &mut ref_to_query_buffer,
            kmer_len,
            min_pct,
            min_count,
            false, // never using borrowed base sequences here
        );

        if has_zero { zero_kmers += 1; }
        if passes {
            writer.write_record(header, &record).context("Failed to write record")?;
        }
    }

    println!("Number of alignments with zero rare kmers: {}", zero_kmers);
    println!("Number of unmapped/skipped alignments: {}", unmapped);
    Ok(())
}

// ==============================================================================
// STREAM 2: GROUPED SECONDARY ALIGNMENTS (NAME-SORTED)
// ==============================================================================

fn process_grouped_stream<R: std::io::Read, W: std::io::Write>(
    reader: &mut bam::io::Reader<R>,
    writer: &mut bam::io::Writer<W>,
    header: &noodles::sam::Header,
    expected_kmers: &KmerLibrary,
    kmer_len: usize,
    min_pct: f64,
    min_count: usize,
) -> Result<()> {
    let mut zero_kmers = 0;
    let mut unmapped = 0;

    let mut ref_to_query_buffer = Vec::with_capacity(100_000);
    let mut base_seq_buffer = Vec::with_capacity(100_000);
    let mut rc_seq_buffer = Vec::with_capacity(100_000);
    let mut local_seq_buffer = Vec::with_capacity(100_000);

    let mut read_group = Vec::new();
    let mut current_qname: Option<Vec<u8>> = None;

    for result in reader.records() {
        let record = result.context("Failed to read BAM record")?;

        let qname = match record.name() {
            Some(name) => {
                let bytes: &[u8] = name.as_ref();
                bytes.to_vec()
            }
            None => continue,
        };

        if Some(&qname) != current_qname.as_ref() {
            if !read_group.is_empty() {
                process_read_group(
                    &read_group, writer, header, expected_kmers, kmer_len, min_pct, min_count,
                    &mut zero_kmers, &mut unmapped, &mut ref_to_query_buffer,
                    &mut base_seq_buffer, &mut rc_seq_buffer, &mut local_seq_buffer,
                )?;
            }
            read_group.clear();
            current_qname = Some(qname);
        }
        read_group.push(record);
    }

    if !read_group.is_empty() {
        process_read_group(
            &read_group, writer, header, expected_kmers, kmer_len, min_pct, min_count,
            &mut zero_kmers, &mut unmapped, &mut ref_to_query_buffer,
            &mut base_seq_buffer, &mut rc_seq_buffer, &mut local_seq_buffer,
        )?;
    }

    println!("Number of alignments with zero rare kmers: {}", zero_kmers);
    println!("Number of unmapped/skipped alignments: {}", unmapped);
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn process_read_group<W: std::io::Write>(
    group: &[bam::Record],
    writer: &mut bam::io::Writer<W>,
    header: &noodles::sam::Header,
    expected_kmers: &KmerLibrary,
    kmer_len: usize,
    min_pct: f64,
    min_count: usize,
    zero_kmers: &mut usize,
    unmapped: &mut usize,
    ref_to_query_buffer: &mut Vec<usize>,
    base_seq_buffer: &mut Vec<u8>,
    rc_seq_buffer: &mut Vec<u8>,
    local_seq_buffer: &mut Vec<u8>,
) -> Result<()> {
    base_seq_buffer.clear();
    let mut base_is_rc = false;

    for rec in group {
        let seq = rec.sequence();
        if seq.len() > base_seq_buffer.len() {
            base_seq_buffer.clear();
            base_seq_buffer.extend(seq.iter().map(u8::from));
            base_is_rc = rec.flags().is_reverse_complemented();
        }
    }

    if base_seq_buffer.is_empty() {
        *unmapped += group.len();
        return Ok(());
    }

    let mut rc_computed = false;

    for record in group {
        if record.flags().is_unmapped() {
            *unmapped += 1;
            continue;
        }

        let chrom_name = match extract_chrom_name(record, header) {
            Some(name) => name,
            None => continue,
        };

        let chrom_kmers = match expected_kmers.get(&chrom_name) {
            Some(kmers) => kmers,
            None => continue,
        };

        let use_base_seq = record.sequence().is_empty();

        let active_seq = if use_base_seq {
            let is_rc = record.flags().is_reverse_complemented();
            if is_rc == base_is_rc {
                base_seq_buffer.as_slice()
            } else {
                if !rc_computed {
                    reverse_complement_bytes(base_seq_buffer, rc_seq_buffer);
                    rc_computed = true;
                }
                rc_seq_buffer.as_slice()
            }
        } else {
            local_seq_buffer.clear();
            local_seq_buffer.extend(record.sequence().iter().map(u8::from));
            local_seq_buffer.as_slice()
        };

        let (passes, has_zero) = evaluate_alignment(
            record,
            active_seq,
            chrom_kmers,
            ref_to_query_buffer,
            kmer_len,
            min_pct,
            min_count,
            use_base_seq,
        );

        if has_zero { *zero_kmers += 1; }
        if passes {
            writer.write_record(header, record).context("Failed to write record")?;
        }
    }

    Ok(())
}

// ==============================================================================
// CORE ENGINE: THE MATH AND BIOLOGY EVALUATION (SHARED)
// ==============================================================================

/// Evaluates a single sequence against the k-mer library.
/// Returns a tuple: (passes_threshold, has_zero_kmers)
#[allow(clippy::too_many_arguments)]
fn evaluate_alignment(
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

// ==============================================================================
// UTILITIES
// ==============================================================================

/// Helper to safely extract chromosome strings from the header map
fn extract_chrom_name(record: &bam::Record, header: &noodles::sam::Header) -> Option<String> {
    let ref_id = match record.reference_sequence_id() {
        Some(Ok(id)) => id,
        _ => return None,
    };
    match header.reference_sequences().get_index(ref_id) {
        Some((name_key, _)) => Some(String::from_utf8_lossy(name_key.as_ref()).into_owned()),
        None => None,
    }
}

/// High-speed helper to reverse-complement a DNA byte sequence
fn reverse_complement_bytes(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.reserve(seq.len());
    for &b in seq.iter().rev() {
        out.push(match b {
            b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A',
            b'a' => b't', b'c' => b'g', b'g' => b'c', b't' => b'a',
            _ => b'N',
        });
    }
}