mod common;

use common::record;
use metheor::readutil::{BismarkRead, CpGPosition, ReadConcordanceState};
use rust_htslib::bam::{HeaderView, Record};

fn read(flags: u16, cigar: &str, xm: &str) -> BismarkRead {
    let header = HeaderView::from_bytes(b"@SQ\tSN:chr1\tLN:1000\n");
    let sam = record("read", flags, 101, 60, cigar, xm);
    let r = Record::from_sam(&header, sam.trim_end().as_bytes()).unwrap();
    BismarkRead::new(&r)
}

#[test]
fn cpg_strand_is_unchanged_by_duplicate_secondary_or_qc_flags() {
    for (flags, expected) in [
        (0, 100),
        (16, 99),
        (99, 100),
        (147, 100),
        (83, 99),
        (163, 99),
    ] {
        for extra in [0, 256, 512, 1024, 2048] {
            let r = read(flags | extra, "4M", "Z.z.");
            let positions: Vec<_> = r.get_cpg_positions().iter().map(|p| p.pos).collect();
            assert_eq!(
                positions,
                vec![expected, expected + 2],
                "flags={}",
                flags | extra
            );
        }
    }
}

#[test]
fn lpmd_distance_is_genomic_across_deletions_and_insertions() {
    for (cigar, xm, min, max, expected_count) in [
        ("3M10D3M", "Z..z..", 2, 4, 0),
        ("3M10D3M", "Z..z..", 13, 13, 1),
        ("3M2I3M", "Z....z..", 3, 3, 1),
        ("3M2I3M", "Z....z..", 5, 5, 0),
    ] {
        let r = read(0, cigar, xm);
        let (c, d, pairs) = r.compute_pairwise_cpg_concordance_discordance(min, max);
        assert_eq!(
            (c, d, pairs.len()),
            (0, expected_count, expected_count as usize),
            "{cigar}"
        );
    }
}

#[test]
fn cpg_positions_follow_cigar_and_ignore_unaligned_query_bases() {
    let r = read(0, "2S2M1I2M1D2M1H", "ZZZ.zZ.z.");
    let cpgs: Vec<_> = r
        .get_cpgs()
        .iter()
        .map(|c| (c.relpos, c.abspos.pos, c.methylated))
        .collect();
    assert_eq!(cpgs, vec![(2, 100, true), (5, 102, true), (7, 105, false)]);
    assert_eq!((r.get_start_pos(), r.get_end_pos()), (100, 106));
}

#[test]
fn all_sixteen_quartet_patterns_have_the_specified_encoding() {
    let patterns = [
        "z.z.z.z.", "z.z.z.Z.", "z.z.Z.z.", "z.z.Z.Z.", "z.Z.z.z.", "z.Z.z.Z.", "z.Z.Z.z.",
        "z.Z.Z.Z.", "Z.z.z.z.", "Z.z.z.Z.", "Z.z.Z.z.", "Z.z.Z.Z.", "Z.Z.z.z.", "Z.Z.z.Z.",
        "Z.Z.Z.z.", "Z.Z.Z.Z.",
    ];
    for (expected, xm) in patterns.iter().enumerate() {
        let r = read(0, "8M", xm);
        let (quartets, patterns) = r.get_cpg_quartets_and_patterns();
        assert_eq!(patterns, vec![expected]);
        let q = &quartets[0];
        assert_eq!(
            [q.pos1.pos, q.pos2.pos, q.pos3.pos, q.pos4.pos],
            [100, 102, 104, 106]
        );
    }
}

#[test]
fn quartets_slide_one_cpg_at_a_time() {
    let r = read(0, "10M", "Z.z.Z.z.Z.");
    let (q, patterns) = r.get_cpg_quartets_and_patterns();
    assert_eq!(patterns, vec![10, 5]);
    let positions: Vec<_> = q
        .iter()
        .map(|q| [q.pos1.pos, q.pos2.pos, q.pos3.pos, q.pos4.pos])
        .collect();
    assert_eq!(positions, vec![[100, 102, 104, 106], [102, 104, 106, 108]]);
    assert!(read(0, "6M", "Z.z.Z.")
        .get_cpg_quartets_and_patterns()
        .0
        .is_empty());
}

#[test]
fn reads_without_cpgs_have_no_patterns_pairs_or_stretches() {
    let r = read(0, "8M", "XxHhUu..");
    assert_eq!(r.get_num_cpgs(), 0);
    assert!(r.get_first_cpg_position().is_none());
    assert!(r.get_cpg_quartets_and_patterns().0.is_empty());
    assert!(r.get_stretch_info().is_empty());
    let (c, d, pairs) = r.compute_pairwise_cpg_concordance_discordance(2, 16);
    assert_eq!((c, d, pairs.len()), (0, 0, 0));
}

#[test]
fn target_filter_preserves_positions_states_and_chromosome_identity() {
    let mut r = read(0, "10M", "Z.z.Z.z.Z.");
    let targets = [
        CpGPosition::new(0, 102),
        CpGPosition::new(0, 106),
        CpGPosition::new(1, 100),
    ]
    .iter()
    .copied()
    .collect();
    r.filter_isin(&targets);
    assert_eq!(
        r.get_cpg_positions()
            .iter()
            .map(|p| p.pos)
            .collect::<Vec<_>>(),
        vec![102, 106]
    );
    assert!(matches!(
        r.get_concordance_state(),
        ReadConcordanceState::Concordant
    ));
    assert_eq!(r.get_first_cpg_position().unwrap().pos, 102);
}

#[test]
fn methylated_stretches_count_all_substretches_and_flush_the_last_run() {
    let r = read(0, "12M", "Z.Z.z.Z.Z.Z.");
    let expected = [(1, 5), (2, 3), (3, 1)].iter().copied().collect();
    assert_eq!(r.get_stretch_info(), expected);
    assert!(matches!(
        r.get_concordance_state(),
        ReadConcordanceState::Discordant
    ));
    assert!(read(0, "4M", "z.z.").get_stretch_info().is_empty());
}

#[test]
fn lpmd_distance_limits_are_inclusive() {
    let r = read(0, "6M", "Z.z.Z.");
    for (min, max, expected) in [
        (2, 4, (1, 2, 3)),
        (2, 2, (0, 2, 2)),
        (4, 4, (1, 0, 1)),
        (3, 3, (0, 0, 0)),
    ] {
        let (c, d, pairs) = r.compute_pairwise_cpg_concordance_discordance(min, max);
        assert_eq!((c, d, pairs.len()), expected);
    }
}

#[test]
#[should_panic(expected = "XM tag length must match read length")]
fn malformed_xm_length_is_rejected_instead_of_silently_truncating_calls() {
    let header = HeaderView::from_bytes(b"@SQ\tSN:chr1\tLN:1000\n");
    let r = Record::from_sam(
        &header,
        b"r\t0\tchr1\t101\t60\t4M\t*\t0\t0\tCGCG\t*\tXM:Z:Z.",
    )
    .unwrap();
    BismarkRead::new(&r);
}
