use rustc_hash::FxHashMap;

pub type KmerVal = u128;

/// Represents a single rare k-mer anchored to the reference genome.
///
/// derive `PartialOrd` and `Ord` for sorting
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct RareKmer {
    pub start: usize,
    pub end: usize,
    pub val: KmerVal,
    pub local_tolerance: usize,
}

/// map chromosome names to a sorted list of their rare k-mers.
pub type KmerLibrary = FxHashMap<String, Vec<RareKmer>>;