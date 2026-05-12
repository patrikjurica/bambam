use anyhow::Result;
use crossbeam_channel::bounded;
use noodles::bam;
use std::thread;

use crate::types::KmerLibrary;
use super::engine::evaluate_alignment;
use super::utils::extract_chrom_name;

#[allow(clippy::too_many_arguments)]
pub(crate) fn process_grouped_stream<R: std::io::Read, W: std::io::Write + Send>(
    reader: &mut bam::io::Reader<R>,
    writer: &mut bam::io::Writer<W>,
    header: &noodles::sam::Header,
    expected_kmers: &KmerLibrary,
    kmer_len: usize,
    min_pct: f64,
    min_count: usize,
    ins_cost: usize,
    del_cost: usize,
    sub_cost: usize,
) -> Result<()> {
    // 1. Create bounded channels to prevent RAM blowouts
    // Capacity of 1000 means the main thread will pause reading if workers get backed up.
    let (tx_work, rx_work) = bounded::<Vec<bam::Record>>(1000);
    let (tx_write, rx_write) = bounded::<Vec<bam::Record>>(1000);

    // Channel to aggregate statistics from all threads before shutting down
    let (tx_stats, rx_stats) = bounded::<(usize, usize)>(50);

    let num_cores = thread::available_parallelism().map(|n| n.get()).unwrap_or(4);
    println!("Spawning {} worker threads for grouped secondary processing...", num_cores);

    // 2. Open the Thread Scope
    thread::scope(|s| {
        // --- WRITER THREAD ---
        // Because of `thread::scope`, we can safely borrow `writer` and `header`
        // without needing Arc or Mutex!
        s.spawn(|| {
            for surviving_group in rx_write {
                for record in surviving_group {
                    // Write to disk sequentially
                    let _ = writer.write_record(header, &record);
                }
            }
        });

        // --- WORKER THREADS ---
        for _ in 0..num_cores {
            let rx = rx_work.clone();
            let tx = tx_write.clone();
            let tx_stat = tx_stats.clone();

            s.spawn(move || {
                // EVERY thread gets its own isolated, zero-allocation buffers.
                // This prevents thread locking and maximizes CPU cache hits.
                let mut ref_to_query_buffer = Vec::with_capacity(100_000);
                let mut base_seq_buffer = Vec::with_capacity(100_000);
                let mut local_seq_buffer = Vec::with_capacity(100_000);
                let mut valid_group_buffer = Vec::with_capacity(10);

                let mut local_zero = 0;
                let mut local_unmapped = 0;

                for group in rx {
                    valid_group_buffer.clear();

                    evaluate_group(
                        &group, header, expected_kmers, kmer_len, min_pct, min_count,
                        &mut local_zero, &mut local_unmapped, &mut ref_to_query_buffer,
                        &mut base_seq_buffer, &mut local_seq_buffer,
                        &mut valid_group_buffer,
                        ins_cost, del_cost, sub_cost // <--- Passed to helper
                    );

                    if !valid_group_buffer.is_empty() {
                        // `mem::take` swaps the full buffer with an empty one,
                        // sending the data to the writer without cloning the Vec!
                        let _ = tx.send(std::mem::take(&mut valid_group_buffer));
                    }
                }

                // When rx channel closes, send stats back to main
                let _ = tx_stat.send((local_zero, local_unmapped));
            });
        }

        // Drop the original clones of transmitters in the main thread
        // so the channels know exactly when to close.
        drop(tx_write);
        drop(tx_stats);

        // --- PRODUCER (Main Thread) ---
        let mut read_group = Vec::new();
        let mut current_qname: Option<Vec<u8>> = None;

        for result in reader.records() {
            let record = result.expect("Failed to read BAM record");

            let qname = match record.name() {
                Some(name) => {
                    let bytes: &[u8] = name.as_ref();
                    bytes.to_vec()
                }
                None => continue,
            };

            if Some(&qname) != current_qname.as_ref() {
                if !read_group.is_empty() {
                    // Send group to workers. using mem::take reuses the allocation.
                    tx_work.send(std::mem::take(&mut read_group)).expect("Worker channel crashed");
                }
                current_qname = Some(qname);
            }
            read_group.push(record);
        }

        if !read_group.is_empty() {
            let _ = tx_work.send(read_group);
        }

        // Drop work transmitter to tell workers we are out of BAM records.
        drop(tx_work);

    }); // Scope ends here. Rust automatically waits for all threads to cleanly finish.

    // 3. Aggregate statistics
    let mut total_zero = 0;
    let mut total_unmapped = 0;
    for (zero, unmapped) in rx_stats {
        total_zero += zero;
        total_unmapped += unmapped;
    }

    println!("Number of alignments with zero rare kmers: {}", total_zero);
    println!("Number of unmapped/skipped alignments: {}", total_unmapped);
    Ok(())
}

// Extracted core loop logic to keep the thread closure clean
#[allow(clippy::too_many_arguments)]
fn evaluate_group(
    group: &[bam::Record],
    header: &noodles::sam::Header,
    expected_kmers: &KmerLibrary,
    kmer_len: usize,
    min_pct: f64,
    min_count: usize,
    local_zero: &mut usize,
    local_unmapped: &mut usize,
    ref_to_query_buffer: &mut Vec<usize>,
    base_seq_buffer: &mut Vec<u8>,
    local_seq_buffer: &mut Vec<u8>,
    valid_group_buffer: &mut Vec<bam::Record>,
    ins_cost: usize,
    del_cost: usize,
    sub_cost: usize,
) {
    base_seq_buffer.clear();

    // Find longest sequence to borrow
    for rec in group {
        let seq = rec.sequence();
        if seq.len() > base_seq_buffer.len() {
            base_seq_buffer.clear();
            base_seq_buffer.extend(seq.iter().map(u8::from));
        }
    }

    if base_seq_buffer.is_empty() {
        *local_unmapped += group.len();
        return;
    }

    for record in group {
        if record.flags().is_unmapped() {
            *local_unmapped += 1;
            continue;
        }

        let chrom_name = match extract_chrom_name(record, header) {
            Some(name) => name,
            None => continue,
        };

        let chrom_kmers = match expected_kmers.get(&chrom_name) {
            Some(kmers) => kmers,
            None => continue,
        };

        let use_base_seq = record.sequence().is_empty();

        let active_seq = if use_base_seq {
            base_seq_buffer.as_slice()
        } else {
            local_seq_buffer.clear();
            local_seq_buffer.extend(record.sequence().iter().map(u8::from));
            local_seq_buffer.as_slice()
        };

        let (passes, has_zero) = evaluate_alignment(
            record, active_seq, chrom_kmers, ref_to_query_buffer,
            kmer_len, min_pct, min_count, use_base_seq,
            ins_cost, del_cost, sub_cost // <--- Passed to engine
        );

        if has_zero { *local_zero += 1; }
        if passes {
            valid_group_buffer.push(record.clone());
        }
    }
}