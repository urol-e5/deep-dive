## deep-dive

This repository is focused on takeing molecular data from three experiments where integrated epigenetic analysis has already been performed, and doing a _deeper dive_ into datasets to extract and analyze remaining points of excitement. 

# Specific sub-efforts

Short-hand and other repos for the three experiments are

A) https://github.com/hputnam/Becker_E5 (Pocillopora verrucosa genome)

- https://gannet.fish.washington.edu/Atumefaciens/hputnam-Becker_E5/
- [OSF](https://osf.io/uayvk/)
- Nutrient exposure 


B) https://github.com/hputnam/HI_Bleaching_Timeseries/tree/main/Dec-July-2019-analysis (Montipora capitata genome)

C) https://github.com/emmastrand/Acclim_Dynamics_molecular (Pocillopora acuta genome)

- This effort is somewhat comprimised as ploidy (unintended), overwhelm signals

D) Acropora pulcra

---


# How to work in this repo
### (_file structure_)

Top level directories are associated with each sub-effort 
For instance.

```
A-Pver
B-Mcap
C-Pacu
D-Apul
```

Within each top level directory there should be 3 directories
```
data
code
output
```
For any document code it should start with a 2 number prefix (eg `01-methylation-explore.Rmd`). All output from that code should be in a sub-directory of `output` named the same as the code. For example the output of ``01-methylation-explore.Rmd` would be in `A-pver/output/01-methylation-explore/`

Please use **Relative Paths**. Commit and Push often. 


---

## More

### Genomes of interest

- *Pocillopora verrucosa* genome v1.0 : http://pver.reefgenomics.org : https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_014529365.1/
- *Montipora capitata* v3 - http://cyanophora.rutgers.edu/montipora/ : https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_006542545.1/
- *Pocillopora acuta* v2 - http://cyanophora.rutgers.edu/Pocillopora_acuta/ : 
