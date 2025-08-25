use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::path::PathBuf;
use tempfile::NamedTempFile;

mod data_generator;
use data_generator::*;

fn setup_output_file() -> (NamedTempFile, PathBuf) {
    let output_file = NamedTempFile::new().expect("Failed to create temp output file");
    let output_path = output_file.path().to_path_buf();
    (output_file, output_path)
}

/// Benchmark FDRP with different parallel thresholds to show performance impact
fn benchmark_fdrp_parallel_thresholds(c: &mut Criterion) {
    let mut group = c.benchmark_group("fdrp_parallel_thresholds");

    let data = generate_small_dataset().expect("Failed to generate dataset");
    let thresholds = vec![10, 50, 100, 200, 1000]; // Small dataset should have ~16 reads, so threshold > 16 forces sequential

    for threshold in thresholds {
        let label = if threshold <= 16 {
            "parallel"
        } else {
            "sequential"
        };
        group.bench_function(
            BenchmarkId::new(format!("threshold_{}", threshold), label),
            |b| {
                b.iter(|| {
                    let (_output_file, output_path) = setup_output_file();
                    let output_str = output_path.to_str().unwrap();
                    metheor::fdrp::compute_with_threshold(
                        &data.bam_path,
                        output_str,
                        10,
                        10,
                        40,
                        35,
                        &None,
                        threshold,
                    );
                });
            },
        );
    }

    group.finish();
}

/// Benchmark qFDRP with different parallel thresholds to show performance impact
fn benchmark_qfdrp_parallel_thresholds(c: &mut Criterion) {
    let mut group = c.benchmark_group("qfdrp_parallel_thresholds");

    let data = generate_small_dataset().expect("Failed to generate dataset");
    let thresholds = vec![10, 50, 100, 200, 1000]; // Small dataset should have ~16 reads

    for threshold in thresholds {
        let label = if threshold <= 16 {
            "parallel"
        } else {
            "sequential"
        };
        group.bench_function(
            BenchmarkId::new(format!("threshold_{}", threshold), label),
            |b| {
                b.iter(|| {
                    let (_output_file, output_path) = setup_output_file();
                    let output_str = output_path.to_str().unwrap();
                    metheor::qfdrp::compute_with_threshold(
                        &data.bam_path,
                        output_str,
                        10,
                        10,
                        40,
                        35,
                        &None,
                        threshold,
                    );
                });
            },
        );
    }

    group.finish();
}

/// Compare all measures to show which ones benefit from parallelization
fn benchmark_parallelization_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("parallelization_comparison");
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(8));

    let data = generate_small_dataset().expect("Failed to generate dataset");

    // PDR - no parallelization implemented
    group.bench_function(BenchmarkId::from_parameter("PDR_no_parallel"), |b| {
        b.iter(|| {
            let (_output_file, output_path) = setup_output_file();
            let output_str = output_path.to_str().unwrap();
            metheor::pdr::compute(&data.bam_path, output_str, 10, 4, 10, &None);
        });
    });

    // MHL - no parallelization implemented
    group.bench_function(BenchmarkId::from_parameter("MHL_no_parallel"), |b| {
        b.iter(|| {
            let (_output_file, output_path) = setup_output_file();
            let output_str = output_path.to_str().unwrap();
            metheor::mhl::compute(&data.bam_path, output_str, 10, 4, 10, &None);
        });
    });

    // FDRP - with smart parallel threshold
    group.bench_function(BenchmarkId::from_parameter("FDRP_smart_parallel"), |b| {
        b.iter(|| {
            let (_output_file, output_path) = setup_output_file();
            let output_str = output_path.to_str().unwrap();
            metheor::fdrp::compute_with_threshold(
                &data.bam_path,
                output_str,
                10,
                10,
                40,
                35,
                &None,
                100,
            );
        });
    });

    // qFDRP - with smart parallel threshold
    group.bench_function(BenchmarkId::from_parameter("qFDRP_smart_parallel"), |b| {
        b.iter(|| {
            let (_output_file, output_path) = setup_output_file();
            let output_str = output_path.to_str().unwrap();
            metheor::qfdrp::compute_with_threshold(
                &data.bam_path,
                output_str,
                10,
                10,
                40,
                35,
                &None,
                100,
            );
        });
    });

    group.finish();
}

/// Performance analysis with different dataset sizes
fn benchmark_dataset_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("dataset_scaling");
    group.sample_size(10);

    let datasets = vec![
        (
            "small",
            generate_small_dataset().expect("Failed to generate small dataset"),
        ),
        (
            "medium",
            generate_medium_dataset().expect("Failed to generate medium dataset"),
        ),
    ];

    for (name, data) in datasets {
        // Test FDRP scaling with sequential processing (high threshold)
        group.bench_function(BenchmarkId::new("FDRP_sequential", name), |b| {
            b.iter(|| {
                let (_output_file, output_path) = setup_output_file();
                let output_str = output_path.to_str().unwrap();
                metheor::fdrp::compute_with_threshold(
                    &data.bam_path,
                    output_str,
                    10,
                    10,
                    40,
                    35,
                    &None,
                    1000,
                );
            });
        });

        // Test FDRP scaling with parallel processing (low threshold)
        group.bench_function(BenchmarkId::new("FDRP_parallel", name), |b| {
            b.iter(|| {
                let (_output_file, output_path) = setup_output_file();
                let output_str = output_path.to_str().unwrap();
                metheor::fdrp::compute_with_threshold(
                    &data.bam_path,
                    output_str,
                    10,
                    10,
                    40,
                    35,
                    &None,
                    10,
                );
            });
        });
    }

    group.finish();
}

criterion_group!(
    benches_parallel,
    benchmark_fdrp_parallel_thresholds,
    benchmark_qfdrp_parallel_thresholds,
    benchmark_parallelization_comparison,
    benchmark_dataset_scaling
);

criterion_main!(benches_parallel);
