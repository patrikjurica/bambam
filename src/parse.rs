use rust_htslib::{bam, bam::Read};
use std::path::PathBuf;
use rust_htslib::bam::Reader;
use rust_htslib::errors::Error;

pub fn parse(file: &PathBuf) -> Result<Reader, rust_htslib::errors::Error> {
    let mut reader = bam::Reader::from_path(file)?;
    let header = bam::Header::from_template(reader.header());

    // samtools simulation
    for (key, records) in header.to_hashmap() {
        if key != "SQ" {
            continue;
        }
        for record in records {
            println!("@{}\tSN:{}\tLN:{}", key, record["SN"], record["LN"]);
        }
    }

    Result::Ok(reader)
}