mod parse;
mod filter;

use clap::Parser;
use std::path::PathBuf;
use rust_htslib::bam::Read;
use crate::parse::parse;
use crate::filter::filter;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    // 1. Required Positional Argument (Input File)
    /// The path to the input BAM file to be filtered.
    pub input_bam_file: PathBuf,

    // 2. Required Positional Argument (Output File)
    /// The path where the filtered output BAM file will be written.
    pub output_bam_file: PathBuf,

    // 3. Optional Flag: --min
    /// Minimal k-mer length to be considered for filtering.
    #[arg(long, default_value_t = 5)]
    pub min: u8, // Using u8 assuming k-mer length is a small positive integer

    // 4. Optional Flag: --max
    /// Maximal k-mer length to be considered for filtering.
    #[arg(long, default_value_t = 10)]
    pub max: u8,

    // 5. Optional Flag: --weight / -w
    /// Significance weight of each k-mer for filtering criteria.
    #[arg(short = 'w', long, default_value_t = 10)]
    pub weight: u32, // Using u32 for weight/significance
}

fn main() {
    let args = Args::parse();

    let mut reads = parse(&args.input_bam_file).unwrap();

    match filter(&mut reads) {
        Ok(_) => println!("BAM iteration complete."),
        Err(e) => eprintln!("An error occurred during BAM processing: {}", e),
    }
}
