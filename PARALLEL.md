# Parallel Processing in Metheor

Metheor includes intelligent parallel processing support for computationally intensive methylation heterogeneity measures, particularly FDRP and qFDRP which involve O(n²) pairwise comparisons.

## Features

### Smart Threshold-Based Parallelization
- **Automatic optimization**: Uses sequential processing for small datasets to avoid parallel overhead
- **Configurable threshold**: Default parallel threshold of 100 reads can be customized
- **Performance protection**: Prevents 35-40% performance regression on small datasets

### Supported Measures
- **FDRP** (Fraction of Discordant Read Pairs): Fully parallelized pairwise comparisons
- **qFDRP** (Quantitative FDRP): Fully parallelized with Hamming distance calculations  
- **PDR, MHL, PM, ME**: Sequential processing (complex algorithmic dependencies)

## CLI Options

### Global Options
```bash
# Control thread count (0 = auto-detect)
metheor --threads 4 fdrp --input data.bam --output result.tsv [other options]

# Customize parallel threshold (minimum reads for parallel processing)
metheor --parallel-threshold 50 fdrp --input data.bam --output result.tsv [other options]
```

### Examples
```bash
# Use 8 threads with default threshold (100 reads)
metheor --threads 8 fdrp --input large_dataset.bam --output fdrp_results.tsv --min-qual 10

# Force parallel processing on small datasets (threshold=10)  
metheor --parallel-threshold 10 qfdrp --input small_data.bam --output qfdrp_results.tsv --min-qual 10

# Auto-detect threads with high threshold (mostly sequential)
metheor --threads 0 --parallel-threshold 500 fdrp --input data.bam --output results.tsv --min-qual 10
```

## Performance Characteristics

### Benchmark Results (Small Test Dataset ~16 reads)
| Configuration | FDRP Time | qFDRP Time | Notes |
|---------------|-----------|------------|-------|
| Sequential (threshold > dataset) | ~235µs | ~254µs | Optimal for small data |
| Parallel (threshold < dataset) | ~370µs | ~390µs | 35-40% overhead |
| Smart threshold (default=100) | ~235µs | ~254µs | Auto-selects sequential |

### Scaling Expectations
- **Small datasets** (< 100 reads): Sequential processing automatically selected
- **Medium datasets** (100-1000 reads): Parallel processing begins to show benefits  
- **Large datasets** (> 1000 reads): Significant parallel speedup expected

## Technical Implementation

### Parallelization Strategy
```rust
// Smart threshold logic in FDRP/qFDRP compute functions
let fdrp = if num_reads >= parallel_threshold {
    // Parallel: Generate combinations and process with rayon
    combinations.par_iter().map(|comb| { ... }).sum()
} else {
    // Sequential: Traditional iterator for small datasets  
    (0..num_reads).combinations(2).map(|comb| { ... }).sum()
};
```

### Thread Pool Configuration
- **Auto-detection**: Uses `num_cpus::get()` when `--threads 0`
- **Manual control**: Specify exact thread count with `--threads N`
- **Global setup**: Thread pool configured once at program startup

## Benchmarking

### Available Benchmarks
```bash
# Quick performance comparison
cargo bench --bench methylation_benchmarks_quiet

# Comprehensive parallel analysis
cargo bench --bench parallel_performance

# All benchmarks
cargo bench
```

### Benchmark Categories
- **Threshold comparison**: Shows sequential vs parallel performance
- **Measure comparison**: Compares all methylation measures
- **Dataset scaling**: Performance across different data sizes

## When Parallelization Helps

### Ideal Scenarios
- **Large BAM files** with high read depth (>100 reads per CpG)
- **FDRP/qFDRP analysis** on comprehensive datasets
- **Multi-core systems** with available CPU resources
- **Batch processing** of multiple large files

### Limited Benefit Scenarios  
- **Small test datasets** (automatically handled by smart threshold)
- **PDR/MHL analysis** (sequential algorithms)
- **I/O bound workflows** where disk access is the bottleneck
- **Memory-constrained systems**

## Algorithm-Specific Notes

### FDRP & qFDRP (Parallelized)
- **Bottleneck**: O(n²) pairwise read comparisons
- **Parallelization**: Each read pair processed independently  
- **Speedup**: Scales with number of CPU cores
- **Memory**: Combinations collected then processed in parallel

### PDR & MHL (Sequential Only)
- **Bottleneck**: Sequential `retain()` operations with state dependencies
- **Complexity**: Algorithms maintain state across read processing iterations
- **Future work**: Could potentially parallelize with algorithmic restructuring

## Troubleshooting

### Performance Issues
- **Slower than expected**: Check if dataset exceeds parallel threshold
- **High memory usage**: Large datasets generate many combinations for parallel processing
- **CPU not fully utilized**: Increase `--threads` or decrease `--parallel-threshold`

### Common Solutions
```bash
# Force sequential processing if parallel is slower
metheor --parallel-threshold 999999 fdrp --input data.bam --output results.tsv

# Optimize for your specific hardware  
metheor --threads $(nproc) --parallel-threshold 200 fdrp --input data.bam --output results.tsv
```

## Development Notes

The parallel implementation prioritizes **correctness** and **performance safety** over maximum parallelization. The smart threshold approach ensures users never experience performance regressions while automatically gaining benefits when datasets are large enough to benefit from parallel processing.