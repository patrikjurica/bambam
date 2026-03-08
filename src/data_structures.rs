use std::collections::HashMap;
use std::error::Error;

use smallvec::SmallVec;

pub type Kmer = u64;
pub type ChromosomeID = u32;
pub type LocationList = SmallVec<[Location; 10]>;

pub struct ChromosomeTable {
    table: HashMap<String, ChromosomeID>,
    reverse_table: HashMap<ChromosomeID, String>,
    lowest_free_id: ChromosomeID,
}

impl ChromosomeTable {
    pub fn new() -> Self {
        Self {
            table: HashMap::new(),
            reverse_table: HashMap::new(),
            lowest_free_id: 0,
        }
    }

    pub fn insert(&mut self, refseq: String) {
        if self.table.contains_key(&refseq) { return }

        self.table.insert(refseq.clone(), self.lowest_free_id);
        self.reverse_table.insert(self.lowest_free_id, refseq);
        self.lowest_free_id += 1;
    }

    pub fn insert_and_get(&mut self, refseq: &str) -> ChromosomeID {
        self.insert(refseq.to_string());

        self.get_id(refseq).unwrap()
    }

    pub fn get_id(&self, refseq: &str) -> Option<ChromosomeID> {
        self.table.get(refseq).cloned()
    }

    pub fn get_string(&self, id: ChromosomeID) -> Option<String> {
        self.reverse_table.get(&id).cloned()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Coordinates {
    pub start: u64,
    pub end: u64,
}

impl Coordinates {
    pub fn new(start: u64, end: u64) -> Self {
        Self {
            start,
            end
        }
    }
    pub fn get_start(&self) -> usize { self.start as usize }
    pub fn get_end(&self) -> usize { self.end as usize }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Location {
    chrom_id: ChromosomeID,
    coords: Coordinates,
}

impl Location {
    pub fn new(chrom_id: ChromosomeID, coords: Coordinates) -> Self {
        Self { chrom_id, coords }
    }

    pub fn get_chrom(&self) -> ChromosomeID {
        self.chrom_id
    }

    pub fn get_coords(&self) -> &Coordinates {
        &self.coords
    }
}

pub type Kmers = HashMap<Kmer, LocationList>;