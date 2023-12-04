`deep-dive/D-Apul/code`

## Directory of code used for analysis.

Each code file should have a matching directory in [`deep-dive/D-Apul/output`](https://github.com/urol-e5/deep-dive/tree/main/D-Apul/output/).

---

- `*.html`: HTML versions of corresponding `.Rmd` files.

- `*.md`: GitHub markdown versions of corresponding `.Rmd` files.

- [`08-Apul-sRNAseq-trimming.Rmd`](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/code/08-Apul-sRNAseq-trimming.Rmd): R markdown file detailing QC and trimming of _A.pulchra_ sRNAseq data.

- [`09-Apul-sRNAseq-miRTrace.Rmd`](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/code/09-Apul-sRNAseq-miRTrace.Rmd): R Markdown detailing analysis of _A.pulchra_ sRNAseq using [miRTrace](https://github.com/friedlanderlab/mirtrace) to identify taxonomic origins of miRNAs.

- [`10-Apul-sRNAseq-BLASTn.Rmd`](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/code/10-Apul-sRNAseq-BLASTn.Rmd): R Markdown detailing analysis of _A.pulchra_ sRNAseq using NCBI BLASTn against miRBase and MiRgene databases.

- [`11-Apul-sRNAseq-miRdeep2.Rmd`](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/code/11-Apul-sRNAseq-miRdeep2.Rmd): R Markdown detailing analysis of _A.pulchra_ sRNAseq using [`miRDeep2`](https://github.com/rajewsky-lab/mirdeep2) for identification of known and novel (predicted) miRNAs.

- [`12-Apul-sRNAseq-MirMachine.Rmd`](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/code/12-Apul-sRNAseq-MirMachine.Rmd): R Markdown detailing identification of miRNAs in the _A. millepora_ genome. Since _A. pulchra_ does not have a sequenced genome available, the _A. millepora_ genome was determined to produce best alignments of _A. pulchra_ RNAseq data compared to other coral genomes.

- [`13-Apul-sRNAseq-ShortStack.Rmd`](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/code/13-Apul-sRNAseq-ShortStack.Rmd): R Markdown detailing annotation of miRNAs using [`ShortStack`](https://github.com/MikeAxtell/ShortStack) with _A. pulchra_ sRNAseq data and the _A. millepora_ genome. Since _A. pulchra_ does not have a sequenced genome available, the _A. millepora_ genome was determined to produce best alignments of _A. pulchra_ RNAseq data compared to other coral genomes.