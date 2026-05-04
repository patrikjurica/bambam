use anyhow::{Context, Result};
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::Record;
use std::fs::File;
use std::io::{BufReader, BufWriter};

use crate::distance::edit_distance;
use crate::kmer::{decode_kmer, encode_kmer};
use crate::types::KmerLibrary;

/// Filters a BAM file based on absolute positional compliance of rare k-mers.
/// ASSUMES the input BAM is name-sorted (`samtools sort -n`) to group alignments.
pub fn filter_bam(
    bam_in_path: &str,
    bam_out_path: &str,
    expected_kmers: &KmerLibrary,
    kmer_len: usize,
    min_pct: f64,
    min_count: usize,
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

    let mut zero_kmers = 0;
    let mut unmapped = 0;

    // zero-allocation buffers, hoisted out of the main loop
    let mut ref_to_query_buffer: Vec<usize> = Vec::with_capacity(100_000);
    let mut base_seq_buffer: Vec<u8> = Vec::with_capacity(100_000);
    let mut rc_seq_buffer: Vec<u8> = Vec::with_capacity(100_000);

    let mut read_group = Vec::new();
    let mut current_qname: Option<Vec<u8>> = None;

    for result in reader.records() {
        let record = result.context("Failed to read BAM record")?;

        // noodles returns the name as a BStr, we explicitly extract it as a byte slice
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
                    &read_group,
                    &mut writer,
                    &header,
                    expected_kmers,
                    kmer_len,
                    min_pct,
                    min_count,
                    &mut zero_kmers,
                    &mut unmapped,
                    &mut ref_to_query_buffer,
                    &mut base_seq_buffer,
                    &mut rc_seq_buffer,
                )?;
            }
            read_group.clear();
            current_qname = Some(qname);
        }

        read_group.push(record);
    }

    // flush the final group
    if !read_group.is_empty() {
        process_read_group(
            &read_group,
            &mut writer,
            &header,
            expected_kmers,
            kmer_len,
            min_pct,
            min_count,
            &mut zero_kmers,
            &mut unmapped,
            &mut ref_to_query_buffer,
            &mut base_seq_buffer,
            &mut rc_seq_buffer,
        )?;
    }

    println!("Number of alignments with zero rare kmers: {}", zero_kmers);
    println!("Number of unmapped/skipped alignments: {}", unmapped);

    Ok(())
}

/// Processes a group of alignments belonging to the exact same original read.
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
) -> Result<()> {
    // 1. locate the sequence inside the group (usually the primary alignment)
    base_seq_buffer.clear();
    let mut base_is_rc = false;

    for rec in group {
        if !rec.sequence().is_empty() {
            base_seq_buffer.extend(rec.sequence().iter().map(u8::from));
            base_is_rc = rec.flags().is_reverse_complemented();
            break;
        }
    }

    if base_seq_buffer.is_empty() {
        *unmapped += group.len();
        return Ok(()); // skip the whole group if no sequence is found
    }

    let mut rc_computed = false;

    // 2. process each alignment in the group independently
    for record in group {
        if record.flags().is_unmapped() {
            *unmapped += 1;
            continue;
        }

        let ref_id = match record.reference_sequence_id() {
            Some(Ok(id)) => id,
            _ => continue,
        };

        let chrom_name = match header.reference_sequences().get_index(ref_id) {
            Some((name_key, _)) => String::from_utf8_lossy(name_key.as_ref()).into_owned(),
            None => continue,
        };

        let chrom_kmers = match expected_kmers.get(&chrom_name) {
            Some(kmers) => kmers,
            None => continue,
        };

        let ref_start = match record.alignment_start() {
            Some(Ok(pos)) => usize::from(pos) - 1,
            _ => continue,
        };

        let ref_end = match record.alignment_end() {
            Some(Ok(pos)) => usize::from(pos) - 1,
            _ => continue,
        };

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
            *zero_kmers += 1;
        }

        if kmers_in_range.len() < min_count {
            continue;
        }

        let ref_len = ref_end.saturating_sub(ref_start);
        ref_to_query_buffer.clear();
        ref_to_query_buffer.resize(ref_len + 1, usize::MAX);

        let mut curr_ref = ref_start;
        let mut curr_query = 0;

        for op_result in record.cigar().iter() {
            let op = op_result?;
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
                Kind::HardClip | Kind::Pad => {}
            }
        }

        // 3. LAZY REVERSE COMPLEMENT EVALUATION
        let is_rc = record.flags().is_reverse_complemented();
        let active_seq = if is_rc == base_is_rc {
            &base_seq_buffer
        } else {
            if !rc_computed {
                reverse_complement_bytes(base_seq_buffer, rc_seq_buffer);
                rc_computed = true;
            }
            &rc_seq_buffer
        };

        let mut valid_kmer_count = 0;

        for kmer in &kmers_in_range {
            let start_idx = kmer.start.saturating_sub(ref_start);
            let end_idx = (kmer.end - 1).saturating_sub(ref_start);

            if start_idx >= ref_to_query_buffer.len() || end_idx >= ref_to_query_buffer.len() {
                continue;
            }

            let q_start = ref_to_query_buffer[start_idx];
            let q_end_inclusive = ref_to_query_buffer[end_idx];

            if q_start == usize::MAX || q_end_inclusive == usize::MAX {
                continue;
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

        if total_expected >= min_count && (survival_rate >= min_pct || min_count == 0) {
            writer.write_record(header, record).context("Failed to write record")?;
        }
    }

    Ok(())
}

/// High-speed helper to reverse-complement a DNA byte sequence
fn reverse_complement_bytes(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.reserve(seq.len());
    for &b in seq.iter().rev() {
        out.push(match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            _ => b'N',
        });
    }
}