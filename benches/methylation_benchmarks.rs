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

fn benchmark_pdr(c: &mut Criterion) {
    let mut group = c.benchmark_group("PDR");
    group.sample_size(10);
    
    let datasets = vec![
        ("small", generate_small_dataset().expect("Failed to generate small dataset")),
        ("medium", generate_medium_dataset().expect("Failed to generate medium dataset")),
    ];
    
    for (name, data) in datasets {
        group.bench_function(BenchmarkId::from_parameter(name), |b| {
            b.iter(|| {
                let (_output_file, output_path) = setup_output_file();
                let output_str = output_path.to_str().unwrap();
                
                metheor::pdr::compute(&data.bam_path, output_str, 10, 4, 10, &None);
            });
        });
    }
    
    group.finish();
}

fn benchmark_lpmd(c: &mut Criterion) {
    let mut group = c.benchmark_group("LPMD");
    group.sample_size(10);
    
    let datasets = vec![
        ("small", generate_small_dataset().expect("Failed to generate small dataset")),
        ("medium", generate_medium_dataset().expect("Failed to generate medium dataset")),
    ];
    
    for (name, data) in datasets {
        group.bench_function(BenchmarkId::from_parameter(name), |b| {
            b.iter(|| {
                let (_output_file, output_path) = setup_output_file();
                let output_str = output_path.to_str().unwrap();
                
                metheor::lpmd::compute(&data.bam_path, output_str, 2, 16, 10, &None, &None);
            });
        });
    }
    
    group.finish();
}

fn benchmark_mhl(c: &mut Criterion) {
    let mut group = c.benchmark_group("MHL");
    group.sample_size(10);
    
    let datasets = vec![
        ("small", generate_small_dataset().expect("Failed to generate small dataset")),
        ("medium", generate_medium_dataset().expect("Failed to generate medium dataset")),
    ];
    
    for (name, data) in datasets {
        group.bench_function(BenchmarkId::from_parameter(name), |b| {
            b.iter(|| {
                let (_output_file, output_path) = setup_output_file();
                let output_str = output_path.to_str().unwrap();
                
                metheor::mhl::compute(&data.bam_path, output_str, 10, 4, 10, &None);
            });
        });
    }
    
    group.finish();
}

fn benchmark_pm(c: &mut Criterion) {
    let mut group = c.benchmark_group("PM");
    group.sample_size(10);
    
    let datasets = vec![
        ("small", generate_small_dataset().expect("Failed to generate small dataset")),
        ("medium", generate_medium_dataset().expect("Failed to generate medium dataset")),
    ];
    
    for (name, data) in datasets {
        group.bench_function(BenchmarkId::from_parameter(name), |b| {
            b.iter(|| {
                let (_output_file, output_path) = setup_output_file();
                let output_str = output_path.to_str().unwrap();
                
                metheor::pm::compute(&data.bam_path, output_str, 10, 10, &None);
            });
        });
    }
    
    group.finish();
}

fn benchmark_me(c: &mut Criterion) {
    let mut group = c.benchmark_group("ME");
    group.sample_size(10);
    
    let datasets = vec![
        ("small", generate_small_dataset().expect("Failed to generate small dataset")),
        ("medium", generate_medium_dataset().expect("Failed to generate medium dataset")),
    ];
    
    for (name, data) in datasets {
        group.bench_function(BenchmarkId::from_parameter(name), |b| {
            b.iter(|| {
                let (_output_file, output_path) = setup_output_file();
                let output_str = output_path.to_str().unwrap();
                
                metheor::me::compute(&data.bam_path, output_str, 10, 10, &None);
            });
        });
    }
    
    group.finish();
}

fn benchmark_fdrp(c: &mut Criterion) {
    let mut group = c.benchmark_group("FDRP");
    group.sample_size(10);
    
    let datasets = vec![
        ("small", generate_small_dataset().expect("Failed to generate small dataset")),
        ("medium", generate_medium_dataset().expect("Failed to generate medium dataset")),
    ];
    
    for (name, data) in datasets {
        group.bench_function(BenchmarkId::from_parameter(name), |b| {
            b.iter(|| {
                let (_output_file, output_path) = setup_output_file();
                let output_str = output_path.to_str().unwrap();
                
                metheor::fdrp::compute(&data.bam_path, output_str, 10, 10, 40, 35, &None);
            });
        });
    }
    
    group.finish();
}

fn benchmark_qfdrp(c: &mut Criterion) {
    let mut group = c.benchmark_group("qFDRP");
    group.sample_size(10);
    
    let datasets = vec![
        ("small", generate_small_dataset().expect("Failed to generate small dataset")),
        ("medium", generate_medium_dataset().expect("Failed to generate medium dataset")),
    ];
    
    for (name, data) in datasets {
        group.bench_function(BenchmarkId::from_parameter(name), |b| {
            b.iter(|| {
                let (_output_file, output_path) = setup_output_file();
                let output_str = output_path.to_str().unwrap();
                
                metheor::qfdrp::compute(&data.bam_path, output_str, 10, 10, 40, 35, &None);
            });
        });
    }
    
    group.finish();
}

fn benchmark_parameter_variations(c: &mut Criterion) {
    let mut group = c.benchmark_group("parameter_variations");
    group.sample_size(10);
    
    let data = generate_medium_dataset().expect("Failed to generate dataset");
    
    for min_depth in &[5, 10, 20, 30] {
        group.bench_function(
            BenchmarkId::new("PDR_min_depth", min_depth),
            |b| {
                b.iter(|| {
                    let (_output_file, output_path) = setup_output_file();
                    let output_str = output_path.to_str().unwrap();
                    
                    metheor::pdr::compute(&data.bam_path, output_str, *min_depth, 4, 10, &None);
                });
            }
        );
    }
    
    for min_cpgs in &[2, 4, 8, 16] {
        group.bench_function(
            BenchmarkId::new("MHL_min_cpgs", min_cpgs),
            |b| {
                b.iter(|| {
                    let (_output_file, output_path) = setup_output_file();
                    let output_str = output_path.to_str().unwrap();
                    
                    metheor::mhl::compute(&data.bam_path, output_str, 10, *min_cpgs, 10, &None);
                });
            }
        );
    }
    
    group.finish();
}

fn benchmark_methylation_patterns(c: &mut Criterion) {
    let mut group = c.benchmark_group("methylation_patterns");
    group.sample_size(10);
    
    let datasets = vec![
        ("low_methylation", generate_low_methylation_dataset().expect("Failed to generate dataset")),
        ("high_discordance", generate_high_discordance_dataset().expect("Failed to generate dataset")),
        ("high_depth", generate_high_depth_dataset().expect("Failed to generate dataset")),
    ];
    
    for (name, data) in datasets {
        group.bench_function(
            BenchmarkId::new("PDR", name),
            |b| {
                b.iter(|| {
                    let (_output_file, output_path) = setup_output_file();
                    let output_str = output_path.to_str().unwrap();
                    
                    metheor::pdr::compute(&data.bam_path, output_str, 10, 4, 10, &None);
                });
            }
        );
    }
    
    group.finish();
}

criterion_group!(
    benches,
    benchmark_pdr,
    benchmark_lpmd,
    benchmark_mhl,
    benchmark_pm,
    benchmark_me,
    benchmark_fdrp,
    benchmark_qfdrp,
    benchmark_parameter_variations,
    benchmark_methylation_patterns
);

criterion_main!(benches);