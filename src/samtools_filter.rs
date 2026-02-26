use rust_htslib::{bam, bam::Read};
use rust_htslib::errors::Error;
use rust_htslib::bam::Reader;

pub fn filter(reader: &mut Reader) -> Result<(), rust_htslib::errors::Error> {
    let mut read_count = 0;
    for record_result in reader.records() {
        let mut record = record_result?;

        process_read(&mut record);
        read_count += 1;
    }

    println!("--- Processed {} reads ---", read_count);
    Ok(())
}

fn process_read(record: &mut bam::Record) -> bool {
    let qname = record.qname(); // Read name
    let flag = record.flags();  // Bitwise FLAG

    // mapping
    let pos = record.pos(); // 0-based start position
    let tid = record.tid(); // Reference sequence ID
    let mapq = record.mapq(); // Mapping quality

    // Accessing sequence and quality
    let seq_len = record.seq_len();
    let seq = record.seq().as_bytes(); // The raw sequence bytes (requires ownership)

    println!("Read: {}", String::from_utf8_lossy(qname));
    println!("  -> Chromosome ID: {}, Start Pos: {}", tid, pos);
    println!("  -> Length: {}, MAPQ: {}", seq_len, mapq);

    if record.is_unmapped() {
        println!("  -> Status: UNMAPPED");
    } else {
        println!("  -> Status: MAPPED");
    }

    true
}

fn get_sequence_bytes(record: &bam::Record) -> Vec<u8> {
    // record.seq() returns a Sequence object.
    // .as_bytes() converts the IUPAC codes from the BAM record
    // into a Vec<u8> of standard ASCII characters (A, C, G, T, N).
    record.seq().as_bytes()
}
