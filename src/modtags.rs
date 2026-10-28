use std::collections::{HashMap, HashSet};

use rust_htslib::bam::{ext::BamRecordExtensions, record::Aux, Record};

use super::{CpG, CpGPosition};

pub(super) fn get_cpgs(record: &Record) -> Vec<CpG> {
    let sequence = record.seq().as_bytes();
    let positions: Vec<_> = record.reference_positions_full().collect();
    let aligned: HashMap<_, _> = positions
        .iter()
        .zip(&sequence)
        .filter_map(|(pos, base)| pos.map(|pos| (pos, *base)))
        .collect();
    let modifications = record
        .basemods_iter()
        .expect("Invalid MM/ML modification tags");
    let cpg_position = |query_pos: usize, negative: bool| {
        let pos = positions[query_pos]?;
        let anchor = pos - i64::from(negative);
        (aligned.get(&anchor) == Some(&b'C') && aligned.get(&(anchor + 1)) == Some(&b'G'))
            .then_some(CpGPosition::new(record.tid(), anchor as i32))
    };
    let mut states = HashMap::new();
    let mm = match record.aux(b"MM").or_else(|_| record.aux(b"Mm")).unwrap() {
        Aux::String(mm) => mm,
        _ => panic!("MM must be a string"),
    };
    let mut probability_count = 0;
    for group in mm.split(';').filter(|group| !group.is_empty()) {
        let mut parts = group.split(',');
        let header = parts.next().unwrap();
        let codes = header
            .get(2..)
            .expect("Invalid MM group")
            .trim_end_matches(['.', '?']);
        assert!(!codes.is_empty(), "MM group must name a modification");
        let canonical = header.as_bytes()[0];
        let mut count: usize = 0;
        let mut canonical_positions: Vec<_> = sequence
            .iter()
            .enumerate()
            .filter_map(|(position, &base)| {
                let base = original_base(base, record.is_reverse());
                (canonical == b'N' || base == canonical || (canonical == b'U' && base == b'T'))
                    .then_some(position)
            })
            .collect();
        if record.is_reverse() {
            canonical_positions.reverse();
        }
        let mut explicit_positions = HashSet::new();
        for delta in parts {
            count = count
                .checked_add(delta.parse::<usize>().expect("Invalid MM delta"))
                .and_then(|count| count.checked_add(1))
                .expect("MM delta overflow");
            assert!(
                count <= canonical_positions.len(),
                "MM coordinates exceed read sequence"
            );
            explicit_positions.insert(canonical_positions[count - 1]);
            probability_count += if codes.bytes().all(|base| base.is_ascii_digit()) {
                1
            } else {
                codes.len()
            };
        }
        if header.ends_with('?') || !(codes == "27551" || codes.contains('m')) {
            continue;
        }
        let negative = record.is_reverse() ^ (header.as_bytes()[1] == b'-');
        for (query_pos, &base) in sequence.iter().enumerate() {
            if (canonical == b'N' || canonical == original_base(base, record.is_reverse()))
                && cpg_position(query_pos, negative).is_some()
                && !explicit_positions.contains(&query_pos)
            {
                states.insert((query_pos, negative), Some(false));
            }
        }
    }
    if let Ok(Aux::ArrayU8(probabilities)) = record.aux(b"ML").or_else(|_| record.aux(b"Ml")) {
        assert_eq!(
            probabilities.len(),
            probability_count,
            "ML length must match MM calls"
        );
    }
    let mut explicit = HashMap::new();
    for modification in modifications {
        let (query_pos, modification) =
            modification.expect("Invalid MM/ML modification coordinates");
        if !matches!(modification.modified_base, 109 | -27551) {
            continue;
        }
        let query_pos = query_pos as usize;
        let negative = record.is_reverse() ^ (modification.strand != 0);
        let value = (modification.qual >= 0).then_some(modification.qual >= 128);
        explicit
            .entry((query_pos, negative))
            .and_modify(|previous| {
                if *previous != value {
                    *previous = None;
                }
            })
            .or_insert(value);
    }
    for (position, value) in explicit {
        states
            .entry(position)
            .and_modify(|previous| {
                if *previous != value {
                    *previous = None;
                }
            })
            .or_insert(value);
    }
    let mut cpgs = Vec::new();
    for ((query_pos, negative), state) in states {
        if let (Some(position), Some(methylated)) = (cpg_position(query_pos, negative), state) {
            cpgs.push(CpG::new(
                query_pos as i32,
                position,
                if methylated { 'Z' } else { 'z' },
            ));
        }
    }
    cpgs.sort_by_key(|cpg| (cpg.abspos.pos, cpg.relpos));
    // A duplex disagreement cannot supply one binary state for this molecule.
    cpgs.chunk_by(|a, b| a.abspos == b.abspos)
        .filter(|calls| {
            calls
                .iter()
                .all(|call| call.methylated == calls[0].methylated)
        })
        .map(|calls| calls[0])
        .collect()
}

fn original_base(base: u8, reverse: bool) -> u8 {
    if !reverse {
        return base;
    }
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => base,
    }
}
