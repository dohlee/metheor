use rust_htslib::{bam, bam::Read};

pub fn get_reader(input: &str) -> bam::Reader {
    let reader = bam::Reader::from_path(&input)
                    .ok()
                    .expect("Error opening BAM file.");

    reader
}

pub fn get_header(reader: &bam::Reader) -> bam::HeaderView {
    let header = bam::HeaderView::from_header(&bam::Header::from_template(reader.header()));

    header
}