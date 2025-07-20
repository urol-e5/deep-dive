## deep-dive

This repository hosts analyses and datasets from [Ashey et al. 2025](https://doi.org/10.1101/2025.03.15.643469), which examines ncRNAs in three coral species. 

*Citation: Ashey J, Rodriguez Casariego J, Bengtsson Z, Huffmyer AS,  Becker DM, Durkin KM,  White S, Eirin-Lopez J, Putnam HM, Roberts SB. 2025. Non-coding RNA Repertoire in reef-building corals. bioRxiv. https://doi.org/10.1101/2025.03.15.643469.*

![](https://urol-e5.github.io/deep-dive-genome-browser/img/dive.png)

### 🔍 Explore 

- [DIVE genome browser](https://urol-e5.github.io/deep-dive-genome-browser/) - View ncRNA features across all three species in our interactive genome browser.
- [OSF Project Link](https://osf.io/aw53f/) - Access large genomic reference files in our Open Science Framework (OSF) storage directory. 
- [![SRA](https://img.shields.io/badge/SRA-PRJNA1236658-blue)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1236658/) - Access raw RNA sequences stored on NCBI. 
- [![SRA](https://img.shields.io/badge/SRA-PRJNA1236666-blue)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1236666/) - Access raw small RNA sequences stored on NCBI. 

### 📂 Repository Structure

Top level directories are organized by species.

For instance:

```
A-Pver
B-Mcap
C-Pacu
D-Apul
E-Peve
F-Ptua
```

Each species folder contains 3 directories: 

```
data # raw and processed data files 
code # Analysis code (prefix files with two digit number, eg 01-*.Rmd)
output # organized into subfolders matching two digit code names 
```

Example: Output for `01-methylation-explore.Rmd` is stored in `A-Pver/output/01-methylation-explore/`.

Please use **Relative Paths**. Commit and Push often. 

### 📑 Key Resources

The table below provides links to lncRNA, miRNA, and piRNA feature and fasta files for the three species. 


| Species  | miRNA FASTA                                                                                                                                                           | miRNA GFF3                                                                                                                                                                  | piRNA BED                                                                                                                                                           | lncRNA BED                                                                                                        | lncRNA FASTA                                                                                                          |
| -------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| **Apul** | [mir.fasta](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/mir.fasta) | [Results.gff3](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results.gff3) | [piRNA clusters](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/output/18-Apul-piRNA/0_piRNA_pipeline_proTRAC/proTRAC_results_APUL/APUL.merged.clusters.bed) | [lncRNA.bed](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/output/05.33-lncRNA-discovery/Apul_lncRNA.bed) | [lncRNA.fasta](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/output/05.33-lncRNA-discovery/Apul_lncRNA.fasta) |
| **Peve** | [mir.fasta](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/output/08.2-Peve-sRNAseq-ShortStack-31bp-fastp-merged/ShortStack_out/mir.fasta)                     | [Results.gff3](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/output/08.2-Peve-sRNAseq-ShortStack-31bp-fastp-merged/ShortStack_out/Results.gff3)                     | [piRNA clusters](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/output/18-Peve-piRNA/0_piRNA_pipeline_proTRAC/proTRAC_results_PEVE/PEVE.merged.clusters.bed) | [lncRNA.bed](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/output/Peve_lncRNA.bed)                        | [lncRNA.fasta](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/output/Peve_lncRNA.fasta)                        |
| **Ptuh** | [mir.fasta](https://github.com/urol-e5/deep-dive/blob/main/F-Pmea/output/13.2.1-Pmea-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/mir.fasta) | [Results.gff3](https://github.com/urol-e5/deep-dive/blob/main/F-Pmea/output/13.2.1-Pmea-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results.gff3) | [piRNA clusters](https://github.com/urol-e5/deep-dive/blob/main/F-Pmea/output/18-Pmea-piRNA/0_piRNA_pipeline_proTRAC/proTRAC_results_PMEA/PMEA.merged.clusters.bed) | [lncRNA.bed](https://github.com/urol-e5/deep-dive/blob/main/F-Pmea/output/02-lncRNA-discovery/Pmea_lncRNA.bed)    | [lncRNA.fasta](https://github.com/urol-e5/deep-dive/blob/main/F-Pmea/output/02-lncRNA-discovery/Pmea_lncRNA.fasta)    |

### 🔗 More Info

See our [species descriptions and genomic resources wiki page](https://github.com/urol-e5/deep-dive/wiki/Species-Characteristics-and-Genomic-Resources) for more details on reference genomes and annotation sources. 