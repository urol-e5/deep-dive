`deep-dive/D-Apul/output/09-Apul-sRNAseq-miRTrace`

## Output files produced by [`09-Apul-sRNAseq-miRTrace.Rmd`](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/code/09-Apul-sRNAseq-miRTrace.Rmd).

---

- [`mirtrace.config`](https://github.com/urol-e5/deep-dive/tree/main/D-Apul/output/09-Apul-sRNAseq-miRTracemirtrace.config): Comma-separated file with one FastQ and corresponding sample name per line. Used as input for [miRTrace](https://github.com/friedlanderlab/mirtrace).

- [`mirtrace-report.html`](https://github.com/urol-e5/deep-dive/tree/main/D-Apul/output/09-Apul-sRNAseq-miRTracemirtrace-report.html): HTML-formatted report generated by [miRTrace](https://github.com/friedlanderlab/mirtrace).

- [`mirtrace-stats-contamination_basic.tsv`](https://github.com/urol-e5/deep-dive/tree/main/D-Apul/output/09-Apul-sRNAseq-miRTracemirtrace-stats-contamination_basic.tsv): Tab-delimited report with counts of sequences from each collapsed FastAs having matches to known miRNAs within each of the [miRTrace](https://github.com/friedlanderlab/mirtrace) Clades.

- [`mirtrace-stats-contamination_detailed.tsv`](https://github.com/urol-e5/deep-dive/tree/main/D-Apul/output/09-Apul-sRNAseq-miRTracemirtrace-stats-contamination_detailed.tsv): Tab-delimited report of _only_ Clades with which sequences were matched, along with the corresponding miRNA families in each clade, and the sequence counts.