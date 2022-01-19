# Metheor

[![conda](https://anaconda.org/dohlee/metheor/badges/installer/conda.svg)](https://anaconda.org/dohlee/metheor)

Compute levels of DNA methylation heterogeneity from Bismark-aligned bisulfite sequencing data.

## Installation
Install with `conda`.
```
$ conda install -c dohlee metheor
```

## Usage

### Local pairwise methylation disorder (LPMD)
```
$ metheor lpmd --input <INPUT> --output <OUTPUT> [OPTIONS]
```

*Options*

`-i, --input`: Path to input BAM file.

`-o, --output`: Path to output table file summarizing the result of LPMD calculation.

`-m, --min-distance`:

`-M, --max-distance`:

`-c, --cpg-set`:

`-q, --min-qual`:
