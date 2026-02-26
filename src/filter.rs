use std::fs::File;
use std::io::{self, BufRead, BufReader};

use rust_htslib::{bam, bam::Read};
use rust_htslib::errors::Error;
use std::collections::HashMap;
use std::vec;
use bio::bio_types::genome::Position;
use rust_htslib::errors::Error::BamAuxStringError;
use crate::{A, C, G, T};

pub type Kmer = u64;

#[derive(Debug)]
pub struct Location {
    chrom: String,
    coords: Coordinates,
}

impl Location {
    pub fn new(chrom: String, coords: Coordinates) -> Self {
        Self { chrom, coords }
    }

    pub fn get_chrom(&self) -> &String {
        &self.chrom
    }

    pub fn get_coords(&self) -> &Coordinates {
        &self.coords
    }
}
#[derive(Debug)]
pub struct Coordinates {
    start: u64,
    end: u64,
}

impl Coordinates {
    pub fn get_start(&self) -> usize { self.start as usize }
    pub fn get_end(&self) -> usize { self.end as usize }
}
pub type Kmers = HashMap<Location, Kmer>;

pub fn filter(reader: &mut bam::Reader, len: &u8) -> Result<Kmers, Error> {
    let mut kmers: Kmers = HashMap::new();

    for record_result in reader.records() {
        let record = record_result?;

        // if record.mapq() < 30 { continue; }

        find_kmers(&mut kmers, record, len);
    }

    Ok(kmers)
}

fn find_kmers(kmers: &mut Kmers, read: bam::Record, len: &u8) {
    let mut coords = get_kmer_coords("i").unwrap();


    let mut pos = 0;
    let seq = read.seq().encoded;

    while ((pos + len) as usize) < seq.len() {


        pos += 1;
    }

}

fn extract_kmers_from_ref(ref_file: &str, coords: &mut HashMap<String, Vec<Location>>) -> Kmers {
    let mut kmers: Kmers = HashMap::new();

    let file = File::open(ref_file).unwrap();
    let mut liner = BufReader::new(&file).lines();

    while let Some(Ok(header)) = liner.next() {
        let part_kmers : &Vec<Location> = coords.get(&header).unwrap();

        let seq = liner.next().unwrap().unwrap();
        for loc in part_kmers {
            extract_kmer(loc, seq.);

        }
    }

    kmers
}

fn extract_kmers_from_ref(ref_file: &str, coords: &mut HashMap<String, Vec<Location>>) -> Kmers {
    let mut kmers: Kmers = HashMap::new();
    let file = File::open(ref_file).expect("Could not open ref file");
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let mut current_seq = String::new();
    let mut current_header = String::new();

    // Strategy: Slurp the whole sequence for a chromosome, then extract
    while let Some(Ok(line)) = lines.next() {
        if line.starts_with('>') {
            // If we just hit a new header, process the PREVIOUS chromosome first
            if !current_seq.is_empty() {
                process_ref_region(&current_header, &current_seq, coords, &mut kmers);
            }
            // Reset for the new chromosome
            // Note: FASTA headers are ">chr1", BED usually uses "chr1". Strip the '>'
            current_header = line.trim_start_matches('>').split_whitespace().next().unwrap().to_string();
            current_seq.clear();
        } else {
            current_seq.push_str(line.trim());
        }
    }
    // Don't forget the last chromosome in the file
    process_ref_region(&current_header, &current_seq, coords, &mut kmers);

    kmers
}

// Helper to handle the actual extraction once a chromosome is fully loaded
fn process_ref_region(header: &str, seq: &str, coords: &HashMap<String, Vec<Location>>, kmers: &mut Kmers) {
    if let Some(locations) = coords.get(header) {
        for loc in locations {
            let kmer_str = extract_kmer(loc, seq);
            let packed = pack_kmer_to_u64(&kmer_str); // You'll need to implement this bit-packing

            // Note: Your Location struct doesn't implement Copy/Clone,
            // so you might need to rethink the Key of your HashMap if you want to keep 'loc'
            // For now, let's assume Kmers is HashMap<u64, String> or similar.
        }
    }
}

fn get_kmer_coords(_file: &str) -> Result<HashMap<String, Vec<Location>>, Error> {
    let file = "/Users/patrik/Desktop/uni/bam_filter/dna_data/ref_kmers.bed";
    let mut coords: HashMap<String, Vec<Location>> = HashMap::new();

    let mut file = File::open(file).unwrap();
    let mut reader = BufReader::new(&file);

    for line_result in reader.lines() {
        let line = line_result.expect("Failed to read line");

        // Skip empty lines to prevent crashes
        if line.trim().is_empty() {
            continue;
        }

        // Create an iterator over the separated words
        let mut parts = line.split_whitespace();

        if let (Some(chrom), Some(start_str), Some(end_str)) = (parts.next(), parts.next(), parts.next()) {

            // Parse the strings into integers ( u64 for genomic coordinates)
            let start: u64 = start_str.parse().expect("Start is not a valid number");
            let end: u64 = end_str.parse().expect("End is not a valid number");

            coords
                .entry(chrom.to_string())
                .or_default()
                .push(Location { chrom: chrom.to_string(), coords: Coordinates { start, end } });

        } else {

            eprintln!("Warning: Line did not have exactly 3 parts: {}", line);

        }
    }

    Ok(coords)
}

fn extract_kmer(location: &Location, refseq: &str) -> String {
    let coords = location.get_coords();
    refseq[coords.get_start()..coords.get_end()].to_owned()
}

fn encode_base(base: u8) -> Result<u64, Error> {
    match base {
        b'A' | b'a' => Ok(A),
        b'C' | b'c' => Ok(C),
        b'G' | b'g' => Ok(G),
        b'T' | b't' => Ok(T),
        _ => Err(BamAuxStringError),
    }
}