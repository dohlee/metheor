mod common;

use common::{assert_close, record, Workspace};

#[test]
fn reference_catalogue_excludes_skipped_cpgs_but_keeps_true_neighbors_across_large_gaps() {
    for (metric, score) in [("pm", 0.5), ("me", 0.25)] {
        let w = Workspace::new();
        let records = record("first", 0, 101, 60, "1M9N1M9N1M9N1M", "Zzzz")
            + &record("second", 0, 101, 60, "1M9N1M9N1M9N1M", "zzzz")
            + &record("skipped", 0, 101, 60, "1M399N1M9N1M9N1M", "ZZZZ")
            + &record("neighbors", 0, 131, 60, "1M369N1M9N1M9N1M", "ZZZZ");
        let input = w.sam(&records);
        let bed = w.write("reference.bed", "chr1\t100\t102\nchr1\t110\t112\nchr1\t120\t122\nchr1\t130\t132\nchr1\t500\t502\nchr1\t510\t512\nchr1\t520\t522\n");
        w.command(metric, &input, &["-d", "1"]).assert().success();
        assert_eq!(w.rows().len(), 3);
        w.command(
            metric,
            &input,
            &["-d", "1", "--reference-cpgs", bed.to_str().unwrap()],
        )
        .assert()
        .success();
        let mut rows = w.rows();
        rows.sort();
        assert_eq!(rows.len(), 2);
        assert_eq!(&rows[0][..5], &["chr1", "100", "110", "120", "130"]);
        assert_close(rows[0][5].parse().unwrap(), score);
        assert_eq!(&rows[1][..5], &["chr1", "130", "500", "510", "520"]);
        assert_close(rows[1][5].parse().unwrap(), 0.0);
    }
}

#[test]
fn whole_genome_catalogues_are_independent_of_target_selection() {
    for metric in ["pm", "me"] {
        let w = Workspace::new();
        let input = w.sam(&record("read", 0, 101, 60, "1M9N1M9N1M9N1M369N1M", "ZzzZZ"));
        let reference = w.write("reference.bed", "# complete catalogue\n\nchr3\t0\t2\nchr1\t500\t502\nchr1\t100\t102\nchr1\t130\t132\nchr1\t110\t112\nchr1\t120\t122\nchr1\t110\t112\n");
        w.command(
            metric,
            &input,
            &["-d", "1", "--reference-cpgs", reference.to_str().unwrap()],
        )
        .assert()
        .success();
        assert_eq!(w.rows().len(), 2);
        let targets = w.write(
            "targets.bed",
            "chr1\t100\t102\nchr1\t120\t122\nchr1\t130\t132\nchr1\t500\t502\n",
        );
        w.command(
            metric,
            &input,
            &["-d", "1", "-c", targets.to_str().unwrap()],
        )
        .assert()
        .success();
        assert_eq!(w.rows().len(), 1);
        w.command(
            metric,
            &input,
            &[
                "-d",
                "1",
                "-c",
                targets.to_str().unwrap(),
                "--reference-cpgs",
                reference.to_str().unwrap(),
            ],
        )
        .assert()
        .success();
        assert!(w.rows().is_empty());
    }
}

#[test]
fn invalid_reference_catalogues_do_not_replace_existing_results() {
    for metric in ["pm", "me"] {
        for contents in [
            "chr1\t100\t99\n",
            "chr1\t-1\t1\n",
            "chr1\t1001\t1003\n",
            "chr1\ttext\t102\n",
            "chr1\t100\n",
        ] {
            let w = Workspace::new();
            let input = w.sam(&record("read", 0, 101, 60, "8M", "Z.Z.Z.Z."));
            let reference = w.write("reference.bed", contents);
            w.write("output.tsv", "previous output\n");
            w.command(
                metric,
                &input,
                &["-d", "1", "--reference-cpgs", reference.to_str().unwrap()],
            )
            .assert()
            .failure();
            assert_eq!(
                std::fs::read_to_string(w.path("output.tsv")).unwrap(),
                "previous output\n"
            );
        }
    }
}
