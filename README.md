 ## deep-dive

This repository is focused on taking molecular data from three experiments where integrated epigenetic analysis has already been performed, and doing a _deeper dive_ into datasets to extract and analyze remaining points of excitement. 

### Explore ncRNA features across three species in our [DIVE genome browser](https://urol-e5.github.io/deep-dive-genome-browser/)!

### [OSF Project Link](https://osf.io/aw53f/)

Follow this link to our Open Science Framework (OSF) storage directory for more easily downloadable GitHub data and storage of larger data files (i.e., genomic data).

# Specific sub-efforts


**D**) Acropora pulchra

**E**) Porites evermanni

**F**) Pocillopora tuahiniensis

---


# How to work in this repo
### (_file structure_)

Top level directories are associated with each sub-effort categorized by species.

For instance:

```
A-Pver
B-Mcap
C-Pacu
D-Apul
E-Peve
F-Ptua
```

Within each top level directory there should be 3 directories: 

```
data
code
output
```

For any document code it should start with a 2 number prefix (eg `01-methylation-explore.Rmd`). All output from that code should be in a sub-directory of `output` named the same as the code. For example the output of `01-methylation-explore.Rmd` would be in `A-pver/output/01-methylation-explore/`

Please use **Relative Paths**. Commit and Push often. 

Links to other data types (e.g. FastQs, BAMs) can be [found in the project wiki](https://github.com/urol-e5/deep-dive/wiki).

---

## More

### Genomes of interest

All genomes of interest can be found in our [species descriptions and genomic resources wiki page](https://github.com/urol-e5/deep-dive/wiki/Species-Characteristics-and-Genomic-Resources).

### ncRNA fasta and gff files 

ncRNA fasta and gff files can be found in our [species descriptions and genomic resources wiki page](https://github.com/urol-e5/deep-dive/wiki/Species-Characteristics-and-Genomic-Resources).


