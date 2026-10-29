# Reference-adjacent CpG quartets

By default, PM and ME use four consecutive CpG observations within each read.
An unobserved site or an alignment gap can therefore produce a quartet that
skips CpGs seen in other reads. `--min-depth` and `--min-qual` filter coverage
and mapping quality; they do not establish adjacency in the reference genome.

To retain only reference-adjacent quartets, supply a complete CpG catalogue:

```sh
metheor me --input sample.sorted.bam --output sample.me.tsv \
    --reference-cpgs genome.cpg.bed
```

The same option is available for `metheor pm`. Each row must identify one CpG
in BED format: chromosome, zero-based cytosine start, and exclusive end.
One-base and two-base intervals are accepted. Rows may be unsorted or repeated;
blank lines, comments and BED track/browser lines are ignored. References not
present in the BAM header are ignored, allowing a whole-genome catalogue with
a subset BAM.

The catalogue must be complete for the analyzed regions. A partial list of
target sites cannot establish adjacency in the reference genome. The optional
`--cpg-set` target filter is separate: removing an intervening target does not
make its neighbors adjacent in the reference catalogue.

Quartets containing an unlisted site, or skipping a listed CpG, are excluded.
A large genomic gap remains valid when no other CpG intervenes. This is an
adjacency rule, rather than a distance cutoff. Output columns and eligible
quartet depth are unchanged. Omitting the option preserves existing results.

The catalogue stores eight bytes per CpG position, plus vector capacity. No
reference genome sequence is required or loaded.
