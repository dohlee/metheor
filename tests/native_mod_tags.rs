mod common;

use common::{assert_close, Workspace};
use metheor::readutil::BismarkRead;
use rust_htslib::bam::{HeaderView, Record};

fn native_read(flags: u16, cigar: &str, seq: &str, tags: &str) -> BismarkRead {
    let header = HeaderView::from_bytes(b"@SQ\tSN:chr1\tLN:100000\n");
    let sam = format!("read\t{flags}\tchr1\t101\t60\t{cigar}\t*\t0\t0\t{seq}\t*\t{tags}");
    let record = Record::from_sam(&header, sam.as_bytes()).unwrap();
    BismarkRead::new(&record)
}

#[test]
fn native_5mc_probabilities_produce_cpg_positions_and_quartet_states() {
    let read = native_read(
        0,
        "8M",
        "CGCGCGCG",
        "MM:Z:C+m?,0,0,0,0;\tML:B:C,255,0,128,127",
    );
    let cpgs: Vec<_> = read
        .get_cpgs()
        .iter()
        .map(|c| (c.relpos, c.abspos.pos, c.methylated))
        .collect();
    assert_eq!(
        cpgs,
        [
            (0, 100, true),
            (2, 102, false),
            (4, 104, true),
            (6, 106, false)
        ]
    );
    assert_eq!(read.get_cpg_quartets_and_patterns().1, [10]);
}

#[test]
fn skipped_cytosines_follow_the_mm_unknown_or_implicit_mode() {
    for (mode, expected) in [
        (
            "",
            vec![(100, false), (102, true), (104, false), (106, false)],
        ),
        (
            ".",
            vec![(100, false), (102, true), (104, false), (106, false)],
        ),
        ("?", vec![(102, true)]),
    ] {
        let read = native_read(
            0,
            "8M",
            "CGCGCGCG",
            &format!("MM:Z:C+m{mode},1;\tML:B:C,255"),
        );
        let calls: Vec<_> = read
            .get_cpgs()
            .iter()
            .map(|c| (c.abspos.pos, c.methylated))
            .collect();
        assert_eq!(calls, expected, "mode={mode}");
    }
}

#[test]
fn native_calls_follow_alignment_orientation_and_cigar() {
    for (flags, cigar, seq, tags, expected) in [
        (
            16,
            "8M",
            "CGCGCGCG",
            "MM:Z:C+m?,0,0,0,0;\tML:B:C,255,0,128,127",
            vec![
                (1, 100, false),
                (3, 102, true),
                (5, 104, false),
                (7, 106, true),
            ],
        ),
        (
            0,
            "8M",
            "CGCGCGCG",
            "MM:Z:G-m?,0,0,0,0;\tML:B:C,255,0,128,127",
            vec![
                (1, 100, true),
                (3, 102, false),
                (5, 104, true),
                (7, 106, false),
            ],
        ),
        (
            0,
            "1S2M1I4M1S",
            "CCGACGCGT",
            "MM:Z:C+m?,0,0,0,0;\tML:B:C,255,255,0,255",
            vec![(1, 100, true), (4, 102, false), (6, 104, true)],
        ),
    ] {
        let read = native_read(flags, cigar, seq, tags);
        let calls: Vec<_> = read
            .get_cpgs()
            .iter()
            .map(|c| (c.relpos, c.abspos.pos, c.methylated))
            .collect();
        assert_eq!(calls, expected, "flags={flags}, cigar={cigar}");
    }
}

#[test]
fn native_modification_codes_and_interleaved_probabilities_select_only_5mc() {
    for tags in [
        "MM:Z:C+27551?,0,0,0,0;\tML:B:C,255,0,128,127",
        "MM:Z:C+mh?,0,0,0,0;\tML:B:C,255,0,0,255,128,1,127,1",
        "Mm:Z:C+m?,0,0,0,0;\tMl:B:C,255,0,128,127",
    ] {
        let read = native_read(0, "8M", "CGCGCGCG", tags);
        assert_eq!(read.get_cpg_quartets_and_patterns().1, [10], "{tags}");
    }
}

#[test]
fn duplex_calls_do_not_create_repeated_or_conflicting_cpg_observations() {
    let read = native_read(
        0,
        "8M",
        "CGCGCGCG",
        "MM:Z:C+m?,0,0,0,0;G-m?,0,0,0,0;\tML:B:C,255,0,255,0,255,0,0,255",
    );
    let calls: Vec<_> = read
        .get_cpgs()
        .iter()
        .map(|c| (c.abspos.pos, c.methylated))
        .collect();
    assert_eq!(calls, [(100, true), (102, false)]);
}

#[test]
fn native_inputs_produce_the_expected_values_for_all_seven_cli_measures() {
    for metric in ["pdr", "pm", "me", "mhl", "fdrp", "qfdrp", "lpmd"] {
        let w = Workspace::new();
        let input = w.sam("methylated\t0\tchr1\t101\t60\t8M\t*\t0\t0\tCGCGCGCG\t*\tMM:Z:C+m?,0,0,0,0;\tML:B:C,255,255,255,255\nunmethylated\t0\tchr1\t101\t60\t8M\t*\t0\t0\tCGCGCGCG\t*\tMM:Z:C+m?,0,0,0,0;\tML:B:C,0,0,0,0\n");
        let args: &[&str] = match metric {
            "lpmd" => &[],
            "fdrp" | "qfdrp" => &["-d", "2", "-l", "1"],
            _ => &["-d", "2"],
        };
        w.command(metric, &input, args).assert().success();
        let rows = w.rows();
        match metric {
            "pm" | "me" => {
                assert_eq!(rows.len(), 1, "{metric}");
                assert_eq!(&rows[0][..5], &["chr1", "100", "102", "104", "106"]);
                assert_close(
                    rows[0][5].parse().unwrap(),
                    if metric == "pm" { 0.5 } else { 0.25 },
                );
            }
            "lpmd" => {
                assert_close(rows[1].last().unwrap().parse().unwrap(), 0.0);
            }
            _ => {
                assert_eq!(rows.len(), 4, "{metric}");
                let expected = match metric {
                    "pdr" => 0.0,
                    "mhl" => 0.5,
                    _ => 1.0,
                };
                for (row, position) in rows.iter().zip([100, 102, 104, 106]) {
                    assert_eq!(row[1], position.to_string());
                    assert_close(row[3].parse().unwrap(), expected);
                }
            }
        }
    }
}

#[test]
fn invalid_native_tag_coordinates_lengths_and_probabilities_preserve_output() {
    for tags in [
        "MM:Z:C+m?,4;\tML:B:C,255",
        "MM:Z:C+m?,0,0,0,0;\tML:B:C,255,0,128,127,42",
        "MM:Z:C+m?,0,0,0,0;\tML:B:C,255",
        "MM:Z:C+m?,0;\tML:B:C,255\tMN:i:7",
        "MM:Z:C+m?,0;\tML:Z:255",
    ] {
        let w = Workspace::new();
        let input = w.sam(&format!(
            "read\t0\tchr1\t101\t60\t8M\t*\t0\t0\tCGCGCGCG\t*\t{tags}\n"
        ));
        w.write("output.tsv", "previous output\n");
        w.command("pdr", &input, &["-d", "1"]).assert().failure();
        assert_eq!(
            std::fs::read_to_string(w.path("output.tsv")).unwrap(),
            "previous output\n",
            "{tags}"
        );
    }
}

#[test]
fn absent_probabilities_stay_unknown_and_existing_xm_takes_precedence() {
    let unknown = native_read(0, "8M", "CGCGCGCG", "MM:Z:C+m?,0,0,0,0;");
    assert_eq!(unknown.get_num_cpgs(), 0);
    let implicit = native_read(0, "8M", "CGCGCGCG", "MM:Z:C+m.;");
    assert_eq!(implicit.get_num_cpgs(), 4);
    assert!(implicit.get_cpgs().iter().all(|c| !c.methylated));
    let existing = native_read(0, "8M", "CGCGCGCG", "XM:Z:Z.Z.Z.Z.\tMM:Z:invalid");
    assert_eq!(existing.get_cpg_quartets_and_patterns().1, [15]);
}

#[test]
fn long_native_reads_contribute_to_fdrp_and_qfdrp() {
    let sequence = format!("CG{}CG", "A".repeat(996));
    for metric in ["fdrp", "qfdrp"] {
        let w = Workspace::new();
        let input = w.sam(&format!("methylated\t0\tchr1\t1\t60\t1000M\t*\t0\t0\t{sequence}\t*\tMM:Z:C+m?,0,0;\tML:B:C,255,255\nunmethylated\t0\tchr1\t1\t60\t1000M\t*\t0\t0\t{sequence}\t*\tMM:Z:C+m?,0,0;\tML:B:C,0,0\n"));
        w.command(metric, &input, &["-d", "2", "-l", "900"])
            .assert()
            .success();
        let rows = w.rows();
        assert_eq!(rows.len(), 2, "{metric}");
        assert_eq!(rows[0][1], "0");
        assert_eq!(rows[1][1], "998");
        for row in rows {
            assert_close(row[3].parse().unwrap(), 1.0);
        }
    }
}

#[test]
fn pairwise_measures_include_reverse_strand_cpgs_at_the_alignment_edge() {
    for (metric, expected) in [("fdrp", 1.0), ("qfdrp", 0.25)] {
        let w = Workspace::new();
        let input = w.sam("a\t16\tchr1\t101\t60\t8M\t*\t0\t0\tGAGAGAGA\t*\tXM:Z:Z.Z.Z.Z.\nb\t16\tchr1\t101\t60\t8M\t*\t0\t0\tAAGAGAGA\t*\tXM:Z:z.Z.Z.Z.\n");
        w.command(metric, &input, &["-d", "2", "-l", "1"])
            .assert()
            .success();
        let rows = w.rows();
        assert_eq!(rows.len(), 4);
        assert_eq!(rows[0][1], "99");
        for row in rows {
            assert_close(row[3].parse().unwrap(), expected);
        }
    }
}

#[test]
fn conflicting_duplicate_5mc_codes_remain_unknown() {
    let read = native_read(0, "8M", "CGCGCGCG", "MM:Z:C+m?,0;C+27551?,0;\tML:B:C,255,0");
    assert_eq!(read.get_num_cpgs(), 0);
}

#[test]
fn implicit_and_explicit_5mc_aliases_cannot_override_each_other() {
    for groups in ["C+m.;C+27551?,0;", "C+27551?,0;C+m.;"] {
        let read = native_read(0, "8M", "CGCGCGCG", &format!("MM:Z:{groups}\tML:B:C,255"));
        let calls: Vec<_> = read
            .get_cpgs()
            .iter()
            .map(|c| (c.abspos.pos, c.methylated))
            .collect();
        assert_eq!(
            calls,
            [(102, false), (104, false), (106, false)],
            "{groups}"
        );
    }
}

#[test]
fn long_native_read_lpmd_counts_all_pairs_in_a_wide_distance_window() {
    let sequence = format!("CG{}", "A".repeat(98)).repeat(100);
    let deltas = vec!["0"; 100].join(",");
    let probabilities = (0..100)
        .map(|i| if i % 2 == 0 { "255" } else { "0" })
        .collect::<Vec<_>>()
        .join(",");
    let tags = format!("MM:Z:C+m?,{deltas};\tML:B:C,{probabilities}");
    let read = native_read(0, "10000M", &sequence, &tags);
    let (concordant, discordant, pairs) =
        read.compute_pairwise_cpg_concordance_discordance(2, 10000);
    assert_eq!((concordant, discordant, pairs.len()), (2450, 2500, 4950));
    let w = Workspace::new();
    let input = w.write("input.sam", &format!("@SQ\tSN:chr1\tLN:20000\nread\t0\tchr1\t101\t60\t10000M\t*\t0\t0\t{sequence}\t*\t{tags}\n"));
    w.command("lpmd", &input, &["-m", "2", "-M", "10000"])
        .assert()
        .success();
    assert_close(w.rows()[1][1].parse().unwrap(), 50.0 / 99.0);
}
