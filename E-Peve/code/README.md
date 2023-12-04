`deep-dive/E-Peve/code`

## Directory of code used for analysis.

Each code file should have a matching directory in [`deep-dive/E-Peve/output`](https://github.com/urol-e5/deep-dive/tree/main/E-Peve/output/).

---

- `*.html`: HTML versions of corresponding `.Rmd` files.

- `*.md`: GitHub markdown versions of corresponding `.Rmd` files.

- [`06-Peve-sRNAseq-trimming.Rmd`](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/code/06-Peve-sRNAseq-trimming.Rmd): R markdown file detailing QC and trimming of _P. evermanni_ sRNAseq data.

- [`07-Peve-sRNAseq-MirMachine.Rmd`](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/code/07-Peve-sRNAseq-MirMachine.Rmd): R Markdown detailing identification of miRNAs in the _P. evermanni_ genome.

- [`08-Peve-sRNAseq-ShortStack.Rmd`](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/code/08-Peve-sRNAseq-ShortStack.Rmd): R Markdown detailing annotation of miRNAs using [`ShortStack`](https://github.com/MikeAxtell/ShortStack) with _P. evermanni_ sRNAseq data and the _P. evermanni_ genome. 

- [`09-Peve-sRNAseq-miRTrace.Rmd`](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/code/09-Peve-sRNAseq-miRTrace.Rmd): R Markdown detailing analysis of _P. evermanni_ sRNAseq using [miRTrace](https://github.com/friedlanderlab/mirtrace) to identify taxonomic origins of miRNAs.

- [`10-Peve-sRNAseq-BLASTn.Rmd`](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/code/10-Peve-sRNAseq-BLASTn.Rmd): R Markdown detailing analysis of _P. evermanni_ sRNAseq using NCBI BLASTn against miRBase and MiRgene databases.

- [`11-Peve-sRNAseq-miRdeep2.Rmd`](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/code/11-Peve-sRNAseq-miRdeep2.Rmd): R Markdown detailing analysis of _P. evermanni_ sRNAseq using [`miRDeep2`](https://github.com/rajewsky-lab/mirdeep2) for identification of known and novel (predicted) miRNAs.