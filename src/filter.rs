use std::fs::File;
use std::io::{self, BufRead, BufReader};

use rust_htslib::{bam, bam::Read};
use rust_htslib::errors::Error;
use std::collections::HashMap;
use std::vec;
use bio::bio_types::genome::Position;
use rust_htslib::errors::Error::BamAuxStringError;
use crate::{A, C, G, T};
use crate::data_structures::{Kmers, Kmer, Location, Coordinates};

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
    

    let mut pos = 0;
    let seq = read.seq().encoded;

    while ((pos + len) as usize) < seq.len() {


        pos += 1;
    }

}