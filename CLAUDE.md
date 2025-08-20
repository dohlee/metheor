# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Metheor is a Rust-based bioinformatics tool for computing DNA methylation heterogeneity levels from Bismark-aligned bisulfite sequencing data. It implements seven different methylation heterogeneity measures:

- **PDR** (Proportion of Discordant Reads): Fraction of reads with mixed methylation states
- **LPMD** (Local Pairwise Methylation Disorder): Local concordance considering CpG distances  
- **MHL** (Methylation Haplotype Load): Conservation of methylation haplotypes
- **PM** (Epipolymorphism): Heterogeneity measure similar to entropy
- **ME** (Methylation Entropy): Information theoretic measure of epiallele diversity
- **FDRP** (Fraction of Discordant Read Pairs): Single CpG resolution diversity measure
- **qFDRP** (Quantitative FDRP): Soft version using normalized hamming distance

## Build and Development Commands

```bash
# Build the project
cargo build

# Build in release mode
cargo build --release

# Run tests
cargo test

# Run tests with verbose output
cargo test --verbose

# Run a specific test
cargo test <test_name>

# Build and run
cargo run -- <subcommand> <options>
```

## Architecture

The project follows a modular Rust architecture:

- **`src/main.rs`**: Entry point that parses CLI arguments and dispatches to appropriate modules
- **`src/lib.rs`**: Library crate with CLI definition using clap derive macros
- **Core computation modules**: Each methylation measure has its own module (pdr.rs, lpmd.rs, mhl.rs, etc.)
- **Utility modules**:
  - `bamutil.rs`: BAM file handling utilities
  - `readutil.rs`: Read processing utilities  
  - `progressbar.rs`: Progress indication
  - `tag.rs`: XM tag addition functionality

### Command Structure

All commands follow the pattern:
```bash
metheor <subcommand> --input <file.bam> --output <file.tsv> [options]
```

Each subcommand is implemented as a separate module with a `compute()` function that handles the core algorithm.

### Key Dependencies

- **rust-htslib**: BAM/SAM file processing
- **clap**: Command-line argument parsing (derive API)
- **bio-types**: Bioinformatics data structures
- **indicatif**: Progress bars
- **interval-tree**: Genomic interval operations

## Testing

The project includes integration tests in the `tests/` directory that use sample genomic data files. Tests use the `assert_cmd` crate for CLI testing and include BAM files, reference sequences, and expected outputs.

Test data includes:
- Sample BAM files (test1.bam through test6.bam)
- Reference genome files (hg38.chr19.fa, tinyref.fa)
- Sample SAM files for various scenarios