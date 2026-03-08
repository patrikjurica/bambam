mod parse;
mod samtools_filter;
mod filter;
mod data_structures;
mod load;

use clap::Parser;
use std::path::PathBuf;
use rust_htslib::bam::Read;
use crate::parse::parse;
use crate::parse::load;
use crate::filter::filter;

// DNA base encoding
pub const A: u64 = 0b00; // 0
pub const C: u64 = 0b01; // 1
pub const G: u64 = 0b10; // 2
pub const T: u64 = 0b11; // 3

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
    /// Minimal k-mer count to be considered for filtering.
    #[arg(long, default_value_t = 5)]
    pub min: u8, // Using u8 assuming k-mer length is a small positive integer

    // 4. Optional Flag: --max
    /// Maximal k-mer count to be considered for filtering.
    #[arg(long, default_value_t = 10)]
    pub max: u8,

    // 5. Optional Flag: --len
    /// Length of the unique k-mer (lower yields more values)
    #[arg(long, default_value_t = 31)]
    pub len: u8,

    // 5. Optional Flag: --weight / -w
    /// Significance weight of each k-mer for filtering criteria.
    #[arg(short = 'w', long, default_value_t = 10)]
    pub weight: u32, // Using u32 for weight/significance
}

fn main() {
    let args = Args::parse();

    // let mut reads = parse(&args.input_bam_file).unwrap();

    let mut reads = load(&args.input_bam_file).unwrap();
    
    let kmers = filter(&mut reads, &args.len);

    match filter(&mut reads) {
        Ok(_) => println!("BAM iteration complete."),
        Err(e) => eprintln!("An error occurred during BAM processing: {}", e),
    }
}
