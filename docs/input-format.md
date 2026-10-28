# Alignment inputs and coordinates

Use reference-coordinate-sorted BAM or SAM alignments for calculations:

```sh
samtools sort -o sample.sorted.bam sample.bam
metheor mhl --input sample.sorted.bam --output sample.mhl.tsv
```

PDR, MHL, FDRP and qFDRP finalize loci as alignments advance. They reject a
decrease in reference order (the `@SQ` dictionary), or in position on the same
reference, before replacing the output file. The actual records are checked;
an `SO:coordinate` header alone does not establish that the input is sorted.
PM, ME and LPMD do not otherwise require record order for their calculations.

Unmapped records and mapped records without a methylation tag are skipped. The
seven measures report their skipped counts on stderr. XM and native MM/Mm tags
are recognized as methylation tags; a present but malformed tag remains an error.

For MHL, `--min-depth` is the minimum number of eligible reads covering a CpG.
Eligibility includes `--min-qual` and `--min-cpgs`; the latter counts observed
CpGs per read after any target-CpG filtering. Increasing the minimum depth can
only remove loci from a correctly sorted input.

CpG intervals use BED/BedGraph coordinates: zero-based start, exclusive end.
`chr17 43125376 43125378` spans two bases. In one-based inclusive coordinates,
these are positions 43125377 and 43125378. Keep this convention when comparing
with tools that report one-based positions.
