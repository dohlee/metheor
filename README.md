# Metheor

[![conda](https://anaconda.org/dohlee/metheor/badges/installer/conda.svg)](https://anaconda.org/dohlee/metheor)

Compute levels of DNA methylation heterogeneity from Bismark-aligned bisulfite sequencing data.

## Installation
Install with `conda`.
```
$ conda install -c dohlee metheor
```

## Usage

### Epipolymorphism (PM)
```
$ metheor pm --input <INPUT> --output <OUTPUT>
```

*Options*

`-i, --input`: Path to input BAM file.

### Methylation entropy (ME)
```
$ metheor me --input <INPUT> --output <OUTPUT>
```

*Options*

`-i, --input`: Path to input BAM file.

### Proportion of discordant reads (PDR)

```
$ metheor pdr --input <INPUT> > <OUTPUT>
```

*Options*

`-i, --input`: Path to input BAM file.

### Local pairwise methylation disorder (LPMD)
```
$ metheor lpmd --input <INPUT> --output <OUTPUT> [OPTIONS]
```

*Options*

`-i, --input`: Path to input BAM file.

`-o, --output`: Path to output table file summarizing the result of LPMD calculation.

`-p, --pairs`: (Optional) Concordance information for all CpG pairs.

`-m, --min-distance`: Minimum distance between CpG pairs to consider. [default: 2]

`-M, --max-distance`: Maximum distance between CpG pairs to consider. [default: 16]

`-q, --min-qual`: Minimum quality for a read to be considered. [default: 10]

`-c, --cpg-set`: (Optional) Specify a predefined set of CpGs (in BED file) to be analyzed.


