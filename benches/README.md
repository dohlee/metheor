# Metheor Benchmarks

This directory contains comprehensive benchmarks for all methylation heterogeneity measures implemented in Metheor.

## Overview

The benchmark suite tests the performance of seven methylation measures:
- **PDR** (Proportion of Discordant Reads)
- **LPMD** (Local Pairwise Methylation Disorder)
- **MHL** (Methylation Haplotype Load)
- **PM** (Epipolymorphism)
- **ME** (Methylation Entropy)
- **FDRP** (Fraction of Discordant Read Pairs)
- **qFDRP** (Quantitative FDRP)

## Running Benchmarks

### Run quiet benchmarks (minimal output, fastest)
```bash
cargo bench --bench methylation_benchmarks_quiet
```

### Run clean benchmarks (some output, comprehensive)
```bash
cargo bench --bench methylation_benchmarks_clean
```

### Run all benchmarks (comprehensive, verbose)
```bash
cargo bench --bench methylation_benchmarks
```

### Run specific benchmark group
```bash
cargo bench --bench methylation_benchmarks_quiet quiet_measures
cargo bench --bench methylation_benchmarks_quiet dataset_sizes
```

### Generate detailed HTML reports
```bash
cargo bench --bench methylation_benchmarks_quiet -- --verbose
```

Reports are saved to `target/criterion/` directory.

### Run with profiling
```bash
cargo bench -- --profile-time=10
```

## Benchmark Categories

### Quiet Benchmarks (`methylation_benchmarks_quiet`) ⚡ **RECOMMENDED**
- **Silent Measures Only**: Excludes verbose LPMD to minimize output noise
- **Quick Results**: ~5 second measurement time per benchmark  
- **Cross-Measure Comparison**: 6 clean measures (PDR, MHL, PM, ME, FDRP, qFDRP)
- **Dataset Size Tests**: Small vs medium performance comparison
- **Parameter Quick Test**: Essential parameter variations

### Clean Benchmarks (`methylation_benchmarks_clean`)
- **All 7 Measures**: Including LPMD (with verbose output)
- **Dataset Size Tests**: Small vs medium dataset performance  
- **Moderate Detail**: ~6 second measurement time per benchmark
- **Comprehensive Coverage**: All measures on identical datasets

### Comprehensive Benchmarks (`methylation_benchmarks`) 
- **Parameter Variation Tests**: Impact of `min_depth`, `min_cpgs`, etc.
- **Methylation Pattern Tests**: Different methylation characteristics  
- **Extensive Analysis**: Longer measurement times, detailed statistics
- **Full Coverage**: All edge cases and parameter combinations

## Data Generation

The benchmark suite includes synthetic data generators that create:
- Bismark-compatible BAM files
- Configurable methylation patterns
- XM tags for methylation states
- Various read depths and coverage patterns

### Generator Configuration
```rust
DataGenerator::new(num_reads)
    .with_read_length(150)
    .with_methylation_rate(0.7)
    .with_discordance_rate(0.2)
    .generate()
```

## Performance Baselines

After running benchmarks, baseline results are stored in:
- `target/criterion/*/base/` - Baseline measurements
- `target/criterion/*/report/` - HTML reports
- `target/criterion/*/estimates.json` - Statistical estimates

## Interpreting Results

### Key Metrics
- **Time**: Wall-clock time per iteration
- **Throughput**: Reads processed per second
- **Memory**: Peak memory usage (when profiling enabled)

### Statistical Analysis
Criterion provides:
- Mean execution time with confidence intervals
- Standard deviation and outlier detection
- Regression analysis for performance changes
- Statistical significance testing

## Continuous Benchmarking

### GitHub Actions Integration
Benchmarks can be run automatically on:
- Pull requests (comparison with base branch)
- Scheduled runs (nightly performance tracking)
- Release tags (performance documentation)

### Performance Regression Detection
The CI pipeline will fail if performance regresses by more than 10% (configurable).

## Advanced Usage

### Custom Benchmark Scenarios
Add new scenarios in `methylation_benchmarks.rs`:
```rust
fn benchmark_custom_scenario(c: &mut Criterion) {
    let mut group = c.benchmark_group("custom");
    // ... your benchmark code
}
```

### Memory Profiling
Use with valgrind/massif:
```bash
valgrind --tool=massif --massif-out-file=massif.out \
    target/release/metheor pdr -i input.bam -o output.tsv
```

### CPU Profiling
Use with perf:
```bash
perf record -g target/release/metheor pdr -i input.bam -o output.tsv
perf report
```

## Optimization Workflow

1. **Establish Baseline**: Run benchmarks before optimization
2. **Profile Hotspots**: Identify performance bottlenecks
3. **Implement Optimization**: Make targeted improvements
4. **Verify Improvement**: Re-run benchmarks to measure impact
5. **Check Correctness**: Ensure tests still pass

## Future Enhancements

Planned benchmark improvements:
- [ ] Real-world dataset benchmarks
- [ ] Multi-threaded performance scaling tests
- [ ] Memory usage tracking
- [ ] I/O performance benchmarks
- [ ] Cache efficiency analysis
- [ ] SIMD optimization benchmarks