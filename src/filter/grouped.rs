use anyhow::{Context, Result};
use noodles::bam;
use noodles::sam::alignment::Record; // Required trait for .flags() and .sequence()
use crate::types::KmerLibrary;

use super::engine::evaluate_alignment;
use super::utils::{extract_chrom_name, reverse_complement_bytes};

pub(crate) fn process_grouped_stream<R: std::io::Read, W: std::io::Write>(
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
