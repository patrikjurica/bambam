use anyhow::{Context, Result};
use noodles::bam;
use crate::types::KmerLibrary;

use super::engine::evaluate_alignment;
use super::utils::extract_chrom_name;

const PRIMARY: bool = false;

#[allow(clippy::too_many_arguments)]
pub(crate) fn process_primary_stream<R: std::io::Read, W: std::io::Write>(
    reader: &mut bam::io::Reader<R>,
    writer: &mut bam::io::Writer<W>,
    header: &noodles::sam::Header,
    expected_kmers: &KmerLibrary,
    kmer_len: usize,
    min_pct: f64,
    min_count: usize,
    ins_cost: usize,
    del_cost: usize,
    sub_cost: usize,
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
            PRIMARY, // primary alignments contains the sequence - no need to borrow
            ins_cost, // <--- Passed to engine
            del_cost, // <--- Passed to engine
            sub_cost, // <--- Passed to engine
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