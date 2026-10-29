use rust_htslib::bam;
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use crate::readutil::{CpGPosition, Quartet};

pub(crate) struct ReferenceCpGs(Vec<CpGPosition>);

impl ReferenceCpGs {
    pub fn load(path: &Option<String>, header: &bam::HeaderView) -> Option<Self> {
        let input = File::open(path.as_ref()?).expect("Could not read reference CpG catalogue");
        let mut positions = Vec::new();
        for line in BufReader::new(input).lines() {
            let line = line.expect("Error reading reference CpG catalogue");
            let line = line.trim();
            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("track ")
                || line.starts_with("browser ")
            {
                continue;
            }
            let mut fields = line.split_whitespace();
            let chrom = fields.next().unwrap();
            let start: i32 = fields
                .next()
                .expect("Reference CpG BED needs three columns")
                .parse()
                .expect("Invalid reference CpG start");
            let end: i32 = fields
                .next()
                .expect("Reference CpG BED needs three columns")
                .parse()
                .expect("Invalid reference CpG end");
            assert!(
                start >= 0 && (1..=2).contains(&(i64::from(end) - i64::from(start))),
                "Reference CpG BED rows must span one or two bases"
            );
            if let Some(tid) = header.tid(chrom.as_bytes()) {
                assert!(
                    end as u64 <= header.target_len(tid).unwrap(),
                    "Reference CpG exceeds alignment reference length"
                );
                positions.push(CpGPosition::new(tid as i32, start));
            }
        }
        positions.sort_unstable();
        positions.dedup();
        Some(Self(positions))
    }

    pub fn contains(&self, quartet: &Quartet) -> bool {
        let Ok(index) = self.0.binary_search(&quartet.pos1) else {
            return false;
        };
        self.0.get(index..index + 4).is_some_and(|positions| {
            positions == [quartet.pos1, quartet.pos2, quartet.pos3, quartet.pos4]
        })
    }
}
