 ## deep-dive

This repository is focused on taking molecular data from three experiments where integrated epigenetic analysis has already been performed, and doing a _deeper dive_ into datasets to extract and analyze remaining points of excitement.

## [OSF Project Link](https://osf.io/aw53f/)

# Specific sub-efforts


**D**) Acropora pulcra

**E**) Porites evermanni

**F**) Pocillopora meandrina

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
F-Pmea
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

### Workflows


