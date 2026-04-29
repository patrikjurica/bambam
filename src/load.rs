use std::collections::HashMap;
use std::error::Error;
use std::vec;
use std::path::PathBuf;

use bio::bio_types::genome::Position;
use bio::io::{bed, fasta};

use rust_htslib::{bam, bam::Read};
use rust_htslib::bam::Reader;

use crate::{A, C, G, T};
use crate::data_structures::{Kmers, Kmer, Location, Coordinates, ChromosomeTable, ChromosomeID};

pub fn load_kmers(bed_file: &PathBuf, ref_file: &PathBuf) -> Kmers {
    let (coords, table) = load_bed(bed_file).unwrap();
    let kmers = load_from_ref(ref_file, &table, &coords).unwrap();
    
    kmers
}


/// Loads the rare k-mer coordinates from a BED file
fn load_bed(file: &PathBuf) -> Result<(Vec<Location>, ChromosomeTable), Box<dyn std::error::Error>> {
    // Takes a filename
    // Uses crate bio for parsing of the BED file
    // Uses Location, Coordinates from my custom data_structures
    // Uses ChromosomeTable that assigns an u32 id for every chromosome ("ptg000001l -> 0")
    //     as using string for insertion and lookup would be wildly inefficient
    let mut reader = bed::Reader::from_file(file).expect("Mistake reading BED file");

    let mut coords: Vec<Location> = Vec::new();
    let mut lookup_table = ChromosomeTable::new();

    for rec in reader.records() {
        let record = rec?;

        coords.push(
            Location::new(lookup_table.insert_and_get(record.chrom()),
                          Coordinates::new(record.start(), record.end()
                          )
            )
        );
    }

    coords.sort();
    Ok((coords, lookup_table))
}

/// Extracts marked rare k-mers from a reference FASTA file
fn load_from_ref(file: &PathBuf, lookup_table: &ChromosomeTable, coords: &Vec<Location>) -> Result<Kmers, Box<dyn std::error::Error>> {
    // Takes a filename of the reference FASTA and a reference to a sorted vec of Locations
    // Uses crate bio for parsing of the FASTA file
    // Uses HashMap to store the k-mers and their location
    let mut reader = fasta::Reader::from_file(file).expect("Mistake reading reference file");
    let mut kmer_library: Kmers = HashMap::new();

    let mut coords_iter = coords.into_iter().copied().peekable();

    for rec in reader.records() {
        let record = rec?;
        let Some(chrom_id) = lookup_table.get_id(record.id()) else {
            // Here we skip because we dont have any rare k-mers in this part of the reference
            continue;
        };
        let seq = record.seq();

        // We peek and if the chromosome ID is still the same, we consume
        while let Some(loc) = coords_iter.peek() {

            if loc.get_chrom() != chrom_id {
                break;
            }

            let current_loc = coords_iter.next().unwrap();

            let start = current_loc.get_coords().get_start();
            let end = current_loc.get_coords().get_end();

            let current_kmer = encode_kmer(&seq[start..end]).unwrap();
            kmer_library
                    .entry(current_kmer)
                    .or_default()
                    .push(current_loc);
        }
    }

    Ok(kmer_library)
}

fn encode_kmer(slice: &[u8]) -> Option<Kmer> {
    let mut kmer: Kmer = 0;

    for &letter in slice {
        kmer = kmer << 2;
        kmer = kmer | encode_base(letter)?;
    }

    Some(kmer)
}

fn encode_base(base: u8) -> Option<Kmer> {
    match base {
        b'A' | b'a' => Some(A),
        b'C' | b'c' => Some(C),
        b'G' | b'g' => Some(G),
        b'T' | b't' => Some(T),
        _ => None,
    }
}