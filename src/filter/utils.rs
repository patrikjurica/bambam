use noodles::bam;

/// Helper to safely extract chromosome strings from the header map
pub(crate) fn extract_chrom_name(record: &bam::Record, header: &noodles::sam::Header) -> Option<String> {
    let ref_id = match record.reference_sequence_id() {
        Some(Ok(id)) => id,
        _ => return None,
    };
    match header.reference_sequences().get_index(ref_id) {
        Some((name_key, _)) => Some(String::from_utf8_lossy(name_key.as_ref()).into_owned()),
        None => None,
    }
}

/// High-speed helper to reverse-complement a DNA byte sequence
pub(crate) fn reverse_complement_bytes(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.reserve(seq.len());
    for &b in seq.iter().rev() {
        out.push(match b {
            b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A',
            b'a' => b't', b'c' => b'g', b'g' => b'c', b't' => b'a',
            _ => b'N',
        });
    }
}