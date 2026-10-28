mod common;

use common::{record, Workspace};
use predicates::str::contains;

#[test]
fn streaming_measures_reject_unsorted_alignments_before_replacing_output() {
    for metric in ["mhl", "pdr", "fdrp", "qfdrp"] {
        let w = Workspace::new();
        let input = w.sam(
            &(record("late", 0, 111, 60, "8M", "Z.Z.Z.Z.")
                + &record("early", 0, 101, 60, "8M", "z.z.z.z.")),
        );
        w.write("output.tsv", "previous output\n");
        w.command(metric, &input, &["-d", "1"])
            .assert()
            .failure()
            .stderr(contains("coordinate-sorted"));
        assert_eq!(
            std::fs::read_to_string(w.path("output.tsv")).unwrap(),
            "previous output\n"
        );
    }
}

#[test]
fn all_measures_skip_and_count_unmapped_and_untagged_records() {
    for metric in ["pdr", "mhl", "fdrp", "qfdrp", "pm", "me", "lpmd"] {
        let w = Workspace::new();
        let clean =
            record("m", 0, 101, 60, "8M", "Z.Z.Z.Z.") + &record("u", 0, 101, 60, "8M", "z.z.z.z.");
        let args: &[&str] = if metric == "lpmd" { &[] } else { &["-d", "2"] };
        let input = w.sam(&clean);
        w.command(metric, &input, args).assert().success();
        let mut expected = w.rows();
        expected.sort();
        let dirty = clean + "no_tag\t0\tchr1\t102\t60\t8M\t*\t0\t0\tCGCGCGCG\t*\nunmapped\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
        let input = w.sam(&dirty);
        w.command(metric, &input, args)
            .assert()
            .success()
            .stderr(contains("Skipped 1 unmapped"))
            .stderr(contains("1 mapped"));
        let mut actual = w.rows();
        actual.sort();
        assert_eq!(actual, expected, "{metric}");
    }
}
