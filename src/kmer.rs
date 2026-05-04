use crate::types::KmerVal;

const A: u64 = 0b00; // 0
const C: u64 = 0b01; // 1
const G: u64 = 0b10; // 2
const T: u64 = 0b11; // 3

const U: u64 = 0b11; // for RNA (T = U)

/// encodes a single nucleotide byte into its 2-bit representation
/// returns `None` if it encounters an 'N' or any invalid character
#[inline(always)]
pub fn encode_base(base: u8) -> Option<u64> {
    match base {
        b'A' | b'a' => Some(A),
        b'C' | b'c' => Some(C),
        b'G' | b'g' => Some(G),
        b'T' | b't' => Some(T),
        b'U' | b'u' => Some(U),
        _ => None,
    }
}

/// encodes a complete DNA sequence (as a byte slice) into a 2-bit integer representation
pub fn encode_kmer(seq: &[u8]) -> Option<KmerVal> {
    let mut kmer: u64 = 0;

    for &byte in seq {
        kmer <<= 2;
        match encode_base(byte) {
            Some(val) => kmer |= val,
            None => return None,
        }
    }

    Some(kmer)
}

/// decodes a 2-bit integer back into a DNA String of the specified length
pub fn decode_kmer(mut val: KmerVal, length: usize) -> String {
    // pre-allocate the vector to the exact length to avoid memory reallocations
    let mut chars = Vec::with_capacity(length);

    for _ in 0..length {
        let base_val = val & 0b11;
        let base_char = match base_val {
            A => b'A',
            C => b'C',
            G => b'G',
            T => b'T',
            _ => unreachable!(),
        };
        chars.push(base_char);
        val >>= 2;
    }

    // for cycle goes from zero to length, so we need to reverse
    chars.reverse();

    String::from_utf8(chars).unwrap()
}