mod common;

use common::Workspace;
use rust_htslib::{bam, bam::Read};

#[test]
fn xm_pads_missing_flanks_and_keeps_records_across_chromosomes() {
    let cases = [
        ("right_one_flank", 0, "chr1", 7, "C", "U"),
        ("right_no_flanks", 0, "chr1", 8, "T", "u"),
        ("left_no_flanks", 16, "chr17", 1, "G", "U"),
        ("left_one_flank", 16, "chr17", 2, "A", "u"),
        ("interior", 0, "chr17", 4, "C", "Z"),
        ("final", 0, "chr17", 8, "C", "U"),
    ];
    for format in ["sam", "bam"] {
        let w = Workspace::new();
        let reference = w.write("reference.fa", ">chr1\nGGAAAACC\n>chr17\nGGACGTCC\n");
        rust_htslib::faidx::build(&reference).unwrap();
        let mut sam =
            String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:8\n@SQ\tSN:chr17\tLN:8\n");
        for (name, flags, chrom, pos, seq, _) in cases {
            sam.push_str(&format!(
                "{name}\t{flags}\t{chrom}\t{pos}\t60\t1M\t*\t0\t0\t{seq}\t*\n"
            ));
        }
        let mut input = w.write("input.sam", &sam);
        if format == "bam" {
            let mut reader = bam::Reader::from_path(&input).unwrap();
            input = w.path("input.bam");
            let header = bam::Header::from_template(reader.header());
            let mut writer = bam::Writer::from_path(&input, &header, bam::Format::Bam).unwrap();
            for record in reader.records() {
                writer.write(&record.unwrap()).unwrap();
            }
        }
        w.command("tag", &input, &["-g", reference.to_str().unwrap()])
            .assert()
            .success();
        let mut reader = bam::Reader::from_path(w.path("output.tsv")).unwrap();
        let header = reader.header().to_owned();
        let records: Vec<_> = reader.records().map(Result::unwrap).collect();
        assert_eq!(records.len(), cases.len(), "{format}");
        for (record, (name, flags, chrom, pos, seq, expected)) in records.iter().zip(cases) {
            assert_eq!(record.qname(), name.as_bytes(), "{format}");
            assert_eq!(record.flags(), flags);
            assert_eq!(header.tid2name(record.tid() as u32), chrom.as_bytes());
            assert_eq!(record.pos(), pos - 1);
            assert_eq!(record.seq().as_bytes(), seq.as_bytes());
            assert_eq!(
                record.aux(b"XM").unwrap(),
                bam::record::Aux::String(expected)
            );
        }
    }
}
