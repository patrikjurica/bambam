use std::collections::HashMap;

pub type Kmer = u64;

#[derive(Debug)]
pub struct Coordinates {
    start: u64,
    end: u64,
}

impl Coordinates {
    pub fn get_start(&self) -> usize { self.start as usize }
    pub fn get_end(&self) -> usize { self.end as usize }
}

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

pub type Kmers = HashMap<Location, Kmer>;