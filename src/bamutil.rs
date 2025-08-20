use rust_htslib::{bam, bam::Read};
use std::str;

pub fn get_reader(input: &str) -> bam::Reader {
    match bam::Reader::from_path(input) {
        Ok(reader) => reader,
        Err(error) => {
            panic!("Error opening BAM file. {}", error);
        }
    }
}

pub fn get_header(reader: &bam::Reader) -> bam::HeaderView {
    bam::HeaderView::from_header(&bam::Header::from_template(reader.header()))
}

pub fn tid2chrom(tid: i32, header: &bam::HeaderView) -> String {
    str::from_utf8(header.tid2name(tid as u32))
        .expect("Error parsing chromosome name.")
        .to_string()
}

pub fn chrom2tid(chrom: &[u8], header: &bam::HeaderView) -> u32 {
    header.tid(chrom).unwrap()
}

pub fn is_paired_end(input: &str) -> bool {
    let mut reader = get_reader(input);
    let mut flag = false;

    if let Some(r) = reader.records().map(|r| r.unwrap()).next() {
        if r.is_paired() {
            flag = true;
        }
    }
    flag
}
