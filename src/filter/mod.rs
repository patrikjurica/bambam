use anyhow::{Context, Result};
use noodles::bam;
use std::fs::File;
use std::io::{BufReader, BufWriter};

use crate::types::KmerLibrary;

// Declare the internal module files
mod primary;
mod grouped;
mod engine;
mod utils;

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
        primary::process_primary_stream(
            &mut reader, &mut writer, &header, expected_kmers, kmer_len, min_pct, min_count
        )?;
    } else {
        println!("Running in GROUPED mode. Expecting name-sorted BAM...");
        grouped::process_grouped_stream(
            &mut reader, &mut writer, &header, expected_kmers, kmer_len, min_pct, min_count
        )?;
    }

    Ok(())
}
