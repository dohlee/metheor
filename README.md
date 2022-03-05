# Metheor ☄

[![conda](https://anaconda.org/dohlee/metheor/badges/installer/conda.svg)](https://anaconda.org/dohlee/metheor)
![version](https://anaconda.org/dohlee/metheor/badges/version.svg)
![downloads](https://anaconda.org/dohlee/metheor/badges/downloads.svg)
![github actions](https://github.com/dohlee/metheor/actions/workflows/rust.yml/badge.svg)

Compute DNA methylation heterogeneity levels from Bismark-aligned bisulfite sequencing data.

## Installation
Install with `conda`.
```
conda install -c dohlee metheor
```

## Usage

### Supported methylation heterogeneity measures

`Metheor` supports seven methylation heterogeneity measures in total.

**Proportion of discordant reads (PDR)**

PDR is defined as a fraction of reads carrying CpGs in discordant methylation states (i.e. containing both methylated and unmethylated CpGs in a single read) with respect to all reads mapped to a CpG.
```
metheor pdr --input <input.bam> --output <output.tsv>
    --min-depth <min_depth> --min-cpgs <min_cpgs>
    --min-qual <min_qual> --cpg-set <cpg_set.bed>
```

*Options*

`-i, --input`: Path to input BAM file.

`-o, --output`: Path to output table file summarizing the result of PDR calculation.

`-d, --min-depth`: Minimum depth of CpG stretches to consider. [default: 10]

`-p, --min-cpgs`: Minimum number of consecutive CpGs in a CpG stretch to consider. [default: 10]

`-q, --min-qual`: Minimum quality for a read to be considered [default: 10]

`-c, --cpg-set`: (Optional) Specify a predefined set of CpGs (in BED file) to be analyzed.

**Methylation entropy (ME)**

Xie et al. proposed an information theoretic measure called methylation entropy (ME), which is calculated as the entropy of epialleles originating from a single genomic locus. 
```
metheor me --input <INPUT> --output <OUTPUT>
    --min-depth <min-depth> --min-qual <min-qual> --cpg-set <cpg-set.bed>
```

*Options*

`-i, --input`: Path to input BAM file.

`-o, --output`: Path to output table file summarizing the result of ME calculation.

`-d, --min-depth`: Minimum depth of reads covering epialleles to consider. [default: 10]

`-q, --min-qual`: Minimum quality for a read to be considered. [default: 10]

`-c, --cpg-set`: (Optional) Specify a predefined set of CpGs (in BED file) to be analyzed.


**Epipolymorphism (PM)**

Landan et al., proposed another measure named epipolymorphism (PM), which also captures the amount of heterogeneity of DNA methylation for a given genomic region. The relationship between ME and PM is analogous to that between entropy and Gini index used for decision trees, and they can be considered as similar measures of DNA methylation heterogeneity in general. Using `metheor`, PM can be calculated with the command below:
```
metheor pm --input <input.bam> --output <output.tsv>
    --min-depth <min-depth> --min-qual <min-qual> --cpg-set <cpg-set.bed>
```

*Options*

`-i, --input`: Path to input BAM file.

`-o, --output`: Path to output table file summarizing the result of epipolymorphism calculation.

`-d, --min-depth`: Minimum depth of reads covering epialleles to consider. [default: 10]

`-q, --min-qual`: Minimum quality for a read to be considered. [default: 10]

`-c, --cpg-set`: (Optional) Specify a predefined set of CpGs (in BED file) to be analyzed.

**Fraction of discordant read pairs (FDRP)**

```
metheor fdrp --input <input.bam> --output <output.tsv>
    --min-qual <min-qual> --max-depth <max-depth> --min-overlap <min-overlap>
    --cpg-set <cpg-set.bed>
```

*Options*

`-i, --input`: Path to input BAM file.

`-o, --output`: Path to output table file summarizing the result of epipolymorphism calculation.

`-q, --min-qual`: Minimum quality for a read to be considered. [default: 10]

`-d, --min-depth`: Minimum depth of reads covering epialleles to consider. [default: 10]

`-n, --max-depth`: Maximum number of reads to consider. [default: 40]

`-c, --cpg-set`: (Optional) Specify a predefined set of CpGs (in BED file) to be analyzed.

**Quantative fraction of discordant read pairs (qFDRP)**

```
metheor qfdrp --input <input.bam> --output <output.tsv>
    --min-qual <min-qual> --max-depth <max-depth> --min-overlap <min-overlap>
    --cpg-set <cpg-set.bed>
```

*Options*

`-i, --input`: Path to input BAM file.

`-o, --output`: Path to output table file summarizing the result of epipolymorphism calculation.

`-q, --min-qual`: Minimum quality for a read to be considered. [default: 10]

`-d, --min-depth`: Minimum depth of reads covering epialleles to consider. [default: 10]

`-n, --max-depth`: Maximum number of reads to consider. [default: 40]

`-c, --cpg-set`: (Optional) Specify a predefined set of CpGs (in BED file) to be analyzed.

**Local pairwise methylation disorder (LPMD)**
```
metheor lpmd --input <INPUT> --output <OUTPUT> [OPTIONS]
```

*Options*

`-i, --input`: Path to input BAM file.

`-o, --output`: Path to output table file summarizing the result of LPMD calculation.

`-p, --pairs`: (Optional) Concordance information for all CpG pairs.

`-m, --min-distance`: Minimum distance between CpG pairs to consider. [default: 2]

`-M, --max-distance`: Maximum distance between CpG pairs to consider. [default: 16]

`-q, --min-qual`: Minimum quality for a read to be considered. [default: 10]

`-c, --cpg-set`: (Optional) Specify a predefined set of CpGs (in BED file) to be analyzed.

## Miscellaneous

**Add bismark `XM` tag to BAM file created with aligners other than bismark**
```
metheor tag --input <INPUT.bam> --output <OUTPUT.bam> --genome <GENOME.fa>
```

*Options*

`-i, --input`: Path to input BAM file.

`-o, --output`: Path to output BAM file tagged with XM tag.

`-g, --genome`: Path to genome fasta file.