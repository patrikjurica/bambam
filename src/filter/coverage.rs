use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};

/// Takes the raw intervals of all surviving reads, merges them, and computes the
/// inverse (gaps) based on the absolute length of each chromosome.
pub fn write_coverage_gaps(
    header: &noodles::sam::Header,
    mut coverage: Vec<Vec<(usize, usize)>>,
    out_path: &str,
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(out_path).context("Failed to create coverage BED")?);

    for (ref_id, (name_bytes, ref_map)) in header.reference_sequences().iter().enumerate() {
        if ref_id >= coverage.len() { break; }

        let name = String::from_utf8_lossy(name_bytes.as_ref());
        let ref_len = usize::from(ref_map.length());

        let intervals = &mut coverage[ref_id];

        // 1. Sort intervals by start coordinate in O(N log N) time
        intervals.sort_unstable_by_key(|&(s, _)| s);

        // 2. Merge overlapping intervals in O(N) time
        let mut merged = Vec::new();
        if !intervals.is_empty() {
            let mut current = intervals[0];
            for &interval in &intervals[1..] {
                if interval.0 <= current.1 {
                    current.1 = current.1.max(interval.1);
                } else {
                    merged.push(current);
                    current = interval;
                }
            }
            merged.push(current);
        }

        // 3. Invert the merged intervals to find the gaps
        let mut current_pos = 0;
        for &(start, end) in &merged {
            if start > current_pos {
                writeln!(writer, "{}\t{}\t{}", name, current_pos, start)?;
            }
            current_pos = current_pos.max(end);
        }

        // 4. Capture the final gap extending to the end of the chromosome
        if current_pos < ref_len {
            writeln!(writer, "{}\t{}\t{}", name, current_pos, ref_len)?;
        }
    }

    Ok(())
}