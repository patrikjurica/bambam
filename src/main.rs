use anyhow::Result;
use clap::Parser;
use std::time::Instant;

use bambam::filter::filter_bam;
use bambam::index::{apply_dynamic_threshold_tolerance, build_rare_kmers};
use bambam::io::export_bed;

/// A CLI tool to filter long-read BAM files using sequence-specific rare k-mers.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// The path to the input BAM file to be filtered.
    input_bam_file: String,

    /// The file containing the reference genome (FASTA).
    ref_file: String,

    /// The path where the filtered output BAM file will be written.
    #[arg(short, long, default_value = "./out.bam")]
    output: String,

    /// Minimal k-mer count in the reference to be considered "rare".
    #[arg(long, default_value_t = 5)]
    min: u32,

    /// Maximal k-mer count in the reference to be considered "rare".
    #[arg(long, default_value_t = 10)]
    max: u32,

    /// Minimal count of k-mers in a read for the read to be kept.
    #[arg(short = 'i', long = "inf", default_value_t = 1)]
    min_count: usize,

    /// Minimum percentage of intact expected rare k-mers required to keep a read.
    #[arg(short = 'p', long = "pct", default_value_t = 100.0)]
    pct: f64,

    /// Length of the unique k-mer.
    #[arg(short, long, default_value_t = 31)]
    len: usize,

    /// Maximum allowed edit distance (or base tolerance if --dyn-tol is used).
    #[arg(short, long, default_value_t = 0)]
    tolerance: usize,

    /// Enable dynamic, density-aware tolerance based on k-mer isolation.
    #[arg(long)]
    dyn_tol: bool,

    /// Optional: Path to output a BED file of the rare k-mer coordinates.
    #[arg(long)]
    bed: Option<String>,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let start_time = Instant::now();

    println!("Building the k-mer dictionary...");
    // build_rare_kmers applies the baseline tolerance to all k-mers initially
    let mut kmers = build_rare_kmers(&args.ref_file, args.len, args.min, args.max, args.tolerance)?;

    // ff dynamic tolerance is requested, overwrite
    if args.dyn_tol {
        println!("Applying dynamic threshold tolerance (Distance Threshold: 5000)...");
        kmers = apply_dynamic_threshold_tolerance(kmers, args.tolerance, 5000);
    } else {
        println!("Applying static tolerance of {} to all k-mers...", args.tolerance);
    }

    // calculate total unique k-mers quickly
    let total_kmers: usize = kmers.values().map(|v| v.len()).sum();
    println!("Successfully loaded {} unique rare k-mers.", total_kmers);

    // export BED if requested
    if let Some(bed_path) = &args.bed {
        println!("Exporting k-mer coordinates to {}...", bed_path);
        export_bed(&kmers, bed_path, args.len)?;
        println!("BED export complete.");
    }

    println!("Processing BAM file: {}...", args.input_bam_file);

    filter_bam(
        &args.input_bam_file,
        &args.output,
        &kmers,
        args.len,
        args.pct,
        args.min_count,
    )?;

    println!("BAM filtering complete. Saved to {}", args.output);

    println!("Took {:.2?} to run.", start_time.elapsed());

    Ok(())
}