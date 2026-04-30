use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::kmer::decode_kmer;
use crate::types::KmerLibrary;

/// exports the rare k-mer coordinates to a standard BED file.
/// format: chrom \t start \t end \t sequence(name)
pub fn export_bed(library: &KmerLibrary, output_path: &str, kmer_len: usize) -> Result<()> {
    let file = File::create(output_path)
        .with_context(|| format!("Failed to create BED file at {}", output_path))?;

    let mut writer = BufWriter::new(file);

    for (chrom, kmers) in library.iter() {
        for kmer in kmers {
            let kmer_seq = decode_kmer(kmer.val, kmer_len);

            writeln!(
                writer,
                "{}\t{}\t{}\t{}",
                chrom, kmer.start, kmer.end, kmer_seq
            ).with_context(|| "Failed to write to BED file")?;
        }
    }

    writer.flush().context("Failed to flush BED file contents to disk")?;

    Ok(())
}