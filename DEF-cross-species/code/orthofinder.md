20231003 

## Orthofinder 

On Andromeda 

Make new directories

```
mkdir e5 
cd e5 
mkdir ortho
cd ortho
```

Put protein sequences of interest into a folder together 


```
mkdir protein_seqs
cd protein_seqs

# Make softlink so don't have to make copies of protein seqs and move them into new folder 
ln -s /data/putnamlab/jillashey/genome/Apul/GCF_013753865.1_Amil_v2.1.protein.faa
ln -s /data/putnamlab/jillashey/genome/Peve/Porites_evermanni_v1.annot.pep.fa
ln -s /data/putnamlab/jillashey/genome/Pmea/Pocillopora_meandrina_HIv1.genes.pep.faa
```

Orthofinder code 

```
cd ../
nano orthofinder.sh 

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=100GB
#SBATCH --export=NONE
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -D /data/putnamlab/jillashey/e5/ortho
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="OrthoFinder_out_error"
#SBATCH --output="OrthoFinder_out"

# load modules needed
module load OrthoFinder/2.5.2-intel-2019b-Python-3.7.4
#module load DIAMOND/0.9.22-foss-2018b  
#module load MCL/14.137-GCCcore-8.3.0
#module load FastME/2.1.6.1-iccifort-2019.5.281 
#module load BLAST+/2.8.1-foss-2018b

# using 10 threads, matching the SLRUM parameters above

orthofinder -f protein_seqs/ -t 10

sbatch orthofinder.sh 
```

Submitted batch job 283344. Took about an hour and a half to run. 

Results discussed below. 

OrthoFinder generates A LOT of output, which is great, but it's a lot to sift through. Here's all the folders it generated: 

```
Citation.txt                     Gene_Duplication_Events  Log.txt      Orthogroup_Sequences  Phylogenetically_Misplaced_Genes       Putative_Xenologs    Single_Copy_Orthologue_Sequences  WorkingDirectory
Comparative_Genomics_Statistics  Gene_Trees               Orthogroups  Orthologues           Phylogenetic_Hierarchical_Orthogroups  Resolved_Gene_Trees  Species_Tree
```

First, let's look at the `Comparative_Genomics_Statistics` folder. In the `Statistics_Overall.tsv` file:

```
Number of species       3
Number of genes 114089
Number of genes in orthogroups  100229
Number of unassigned genes      13860
Percentage of genes in orthogroups      87.9
Percentage of unassigned genes  12.1
Number of orthogroups   22783
Number of species-specific orthogroups  5282
Number of genes in species-specific orthogroups 24092
Percentage of genes in species-specific orthogroups     21.1
Mean orthogroup size    4.4
Median orthogroup size  3.0
G50 (assigned genes)    5
G50 (all genes) 4
O50 (assigned genes)    5789
O50 (all genes) 7507
Number of orthogroups with all species present  12489
Number of single-copy orthogroups       6401
Date    2023-10-04
Orthogroups file        Orthogroups.tsv
Unassigned genes file   Orthogroups_UnassignedGenes.tsv
Per-species statistics  Statistics_PerSpecies.tsv
Overall statistics      Statistics_Overall.tsv
Orthogroups shared between species      Orthogroups_SpeciesOverlaps.tsv

Average number of genes per-species in orthogroup       Number of orthogroups   Percentage of orthogroups       Number of genes Percentage of genes
<1      4783    21.0    9566    9.5
'1      14013   61.5    49555   49.4
'2      2260    9.9     15167   15.1
'3      786     3.4     7670    7.7
'4      383     1.7     4909    4.9
'5      189     0.8     2988    3.0
'6      112     0.5     2113    2.1
'7      70      0.3     1535    1.5
'8      43      0.2     1066    1.1
'9      33      0.1     921     0.9
'10     22      0.1     684     0.7
11-15   62      0.3     2343    2.3
16-20   19      0.1     1026    1.0
21-50   8       0.0     686     0.7
51-100  0       0.0     0       0.0
101-150 0       0.0     0       0.0
151-200 0       0.0     0       0.0
201-500 0       0.0     0       0.0
501-1000        0       0.0     0       0.0
'1001+  0       0.0     0       0.0

Number of species in orthogroup Number of orthogroups
1       5282
2       5012
3       12489
```

Almost 88% of the genes were sorted into an orthogroup! The majority of genes look like they were mostly assigned into only one orthogroup. 

In the `Statistics_PerSpecies.tsv` file:

```
        GCF_013753865.1_Amil_v2.1.protein       Pocillopora_meandrina_HIv1.genes.pep    Porites_evermanni_v1.annot.pep
Number of genes 41860   31840   40389
Number of genes in orthogroups  37278   28162   34789
Number of unassigned genes      4582    3678    5600
Percentage of genes in orthogroups      89.1    88.4    86.1
Percentage of unassigned genes  10.9    11.6    13.9
Number of orthogroups containing species        16804   17367   18602
Percentage of orthogroups containing species    73.8    76.2    81.6
Number of species-specific orthogroups  1784    1363    2135
Number of genes in species-specific orthogroups 9270    6031    8791
Percentage of genes in species-specific orthogroups     22.1    18.9    21.8
```
