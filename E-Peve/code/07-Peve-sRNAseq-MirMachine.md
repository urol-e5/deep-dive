07-Peve-sRNAseq-MirMachine
================
Sam White
2023-11-15

- <a href="#1-set-r-variables" id="toc-1-set-r-variables">1 Set R
  variables</a>
- <a href="#2-create-a-bash-variables-file"
  id="toc-2-create-a-bash-variables-file">2 Create a Bash variables
  file</a>
- <a href="#3-load-mirmachine-conda-environment"
  id="toc-3-load-mirmachine-conda-environment">3 Load MirMachine conda
  environment</a>
- <a href="#4-download-pevermanni-genome-from-our-server"
  id="toc-4-download-pevermanni-genome-from-our-server">4 Download
  <em>P.evermanni</em> genome from our server</a>
- <a href="#5-run-mirmachine" id="toc-5-run-mirmachine">5 Run
  MirMachine</a>
- <a href="#6-results" id="toc-6-results">6 Results</a>
  - <a href="#61-view-mirmachine-output-structure"
    id="toc-61-view-mirmachine-output-structure">6.1 View MirMachine output
    structure</a>
  - <a href="#62-check-fasta" id="toc-62-check-fasta">6.2 Check FastA</a>
  - <a href="#63-check-filtered-gff" id="toc-63-check-filtered-gff">6.3
    Check filtered GFF</a>
  - <a href="#64-counts-of-predicted-high-confidence-mirnas"
    id="toc-64-counts-of-predicted-high-confidence-mirnas">6.4 Counts of
    predicted high-confidence miRNAs</a>
- <a href="#7-citations" id="toc-7-citations">7 Citations</a>

Use [MirMachine](https://github.com/sinanugur/MirMachine) ([Umu et al.
2022](#ref-umu2022)) to identify potential miRNA homologs in the
*Pevermanni* genome.

Genome FastA originally downloaded from
<https://www.genoscope.cns.fr/corals/genomes.html> on 20230629.

------------------------------------------------------------------------

Inputs:

- Genome FastA.

Outputs:

- Confidence-filtered GFF of predicted miRNAs.
- Unfiltered FastA of predicted miRNAs.
- Unfiltered GFF of predicted miRNAs.

Software requirements:

- Requires a [MirMachine](https://github.com/sinanugur/MirMachine)
  Conda/Mamba environment, per the installation instructions.

------------------------------------------------------------------------

# 1 Set R variables

Replace with name of your MirMachine environment and the path to the
corresponding conda installation (find this *after* you’ve activated the
environment).

E.g.

``` r
# Activate environment
conda activate mirmachine_env

# Find conda path
which conda
```

``` r
mirmachine_conda_env_name <- c("mirmachine_env")
mirmachine_cond_path <- c("/home/sam/programs/mambaforge/condabin/conda")
```

# 2 Create a Bash variables file

This allows usage of Bash variables across R Markdown chunks.

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "# Data directories"
echo 'export deep_dive_dir=/home/shared/8TB_HDD_01/sam/gitrepos/deep-dive'
echo 'export output_dir_top=${deep_dive_dir}/E-Peve/output/07-Peve-sRNAseq-MirMachine'
echo ""

echo 'export mirmarchine_output_prefix="Peve-MirMachine"'

echo "# Input/Output files"
echo 'export genome_fasta_dir=${deep_dive_dir}/E-Peve/data'
echo 'export genome_fasta_name="Porites_evermanni_v1.fa"'
echo 'export genome_fasta="${genome_fasta_dir}/${genome_fasta_name}"'

echo ""

echo "# External data URLs"
echo 'export genome_fasta_url="https://gannet.fish.washington.edu/seashell/snaps/"'
echo ""

echo "# Set number of CPUs to use"
echo 'export threads=40'
echo ""
} > .bashvars

cat .bashvars
```

    #### Assign Variables ####

    # Data directories
    export deep_dive_dir=/home/shared/8TB_HDD_01/sam/gitrepos/deep-dive
    export output_dir_top=${deep_dive_dir}/E-Peve/output/07-Peve-sRNAseq-MirMachine

    export mirmarchine_output_prefix="Peve-MirMachine"
    # Input/Output files
    export genome_fasta_dir=${deep_dive_dir}/E-Peve/data
    export genome_fasta_name="Porites_evermanni_v1.fa"
    export genome_fasta="${genome_fasta_dir}/${genome_fasta_name}"

    # External data URLs
    export genome_fasta_url="https://gannet.fish.washington.edu/seashell/snaps/"

    # Set number of CPUs to use
    export threads=40

# 3 Load [MirMachine](https://github.com/sinanugur/MirMachine) conda environment

If this is successful, the first line of output should show that the
Python being used is the one in your
[MirMachine](https://github.com/sinanugur/MirMachine) conda environment
path.

E.g.

`python:         /home/sam/programs/mambaforge/envs/mirmachine_env/bin/python`

``` r
use_condaenv(condaenv = mirmachine_conda_env_name, conda = mirmachine_cond_path)

# Check successful env loading
py_config()
```

    python:         /home/sam/programs/mambaforge/envs/mirmachine_env/bin/python
    libpython:      /home/sam/programs/mambaforge/envs/mirmachine_env/lib/libpython3.10.so
    pythonhome:     /home/sam/programs/mambaforge/envs/mirmachine_env:/home/sam/programs/mambaforge/envs/mirmachine_env
    version:        3.10.12 | packaged by conda-forge | (main, Jun 23 2023, 22:40:32) [GCC 12.3.0]
    numpy:          /home/sam/programs/mambaforge/envs/mirmachine_env/lib/python3.10/site-packages/numpy
    numpy_version:  1.25.2

    NOTE: Python version was forced by use_python() function

# 4 Download *P.evermanni* genome from our server

``` bash
# Load bash variables into memory
source .bashvars

wget \
--directory-prefix ${genome_fasta_dir} \
--recursive \
--no-check-certificate \
--continue \
--no-host-directories \
--no-directories \
--no-parent \
--quiet \
--execute robots=off \
--accept "${genome_fasta_name}" ${genome_fasta_url}

ls -lh "${genome_fasta_dir}"
```

    total 625M
    drwxr-xr-x 3 sam sam 4.0K Nov 15 13:50 06-Peve-sRNAseq-trimming
    -rw-rw-r-- 1 sam sam  15M Oct 30 14:00 peve_bedtools_lncRNAs.fasta
    -rw-r--r-- 1 sam sam  24M Nov  6 12:54 Porites_evermanni_v1.annot.gff
    -rw-r--r-- 1 sam sam 586M Jun 30 09:21 Porites_evermanni_v1.fa
    -rw-r--r-- 1 sam sam 422K Nov 15 16:29 Porites_evermanni_v1.fa.fai
    -rw-rw-r-- 1 sam sam    0 Oct 30 14:00 README.md

# 5 Run [MirMachine](https://github.com/sinanugur/MirMachine)

The `--species` option specifies output naming structure

``` bash
# Load bash variables into memory
source .bashvars

cd "${output_dir_top}"

time \
MirMachine.py \
--node Metazoa \
--species ${mirmarchine_output_prefix} \
--genome ${genome_fasta} \
--cpu 40 \
--add-all-nodes \
&> mirmachine.log
```

# 6 Results

## 6.1 View [MirMachine](https://github.com/sinanugur/MirMachine) output structure

``` bash
# Load bash variables into memory
source .bashvars

tree -h ${output_dir_top}/results/predictions
```

    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/E-Peve/output/07-Peve-sRNAseq-MirMachine/results/predictions
    ├── [4.0K]  fasta
    │   └── [109K]  Peve-MirMachine.PRE.fasta
    ├── [4.0K]  filtered_gff
    │   └── [ 29K]  Peve-MirMachine.PRE.gff
    ├── [4.0K]  gff
    │   └── [204K]  Peve-MirMachine.PRE.gff
    └── [4.0K]  heatmap
        └── [ 30K]  Peve-MirMachine.heatmap.csv

    4 directories, 4 files

## 6.2 Check FastA

``` bash
# Load bash variables into memory
source .bashvars

echo "Total number of predicted miRNAs (high & low confidence):"
grep --count "^>" "${output_dir_top}/results/predictions/fasta/${mirmarchine_output_prefix}.PRE.fasta"
echo ""
echo "----------------------------------------------------------------------------------"
echo ""

head "${output_dir_top}/results/predictions/fasta/${mirmarchine_output_prefix}.PRE.fasta"
```

    Total number of predicted miRNAs (high & low confidence):
    794

    ----------------------------------------------------------------------------------

    >Bantam.PRE_Porites_evermani_scaffold_956_127250_127374_(+)_LOWconf
    CGGTTTTCGCATGATATCATGAATTATCAAAACCGAGGTCTGTGTTATCTGCCGAAGCCGAAGGCTGAGGCAGATAATACAGACACGAGGTTTTGATAATTCATGATATCATGCGAAAACCGAAT
    >Bantam.PRE_Porites_evermani_scaffold_4467_1811_1935_(+)_LOWconf
    CGGTTTTCGCATGATATCATGAATTATCAAAACCGAGGTCTGTGTTATCTGCCCAAGCCGAAGGCTGAGGCAGATAATACAGACACGAGGTTTTGATAATTCATGATATCATGCGAAAACCGAAT
    >Iab-4.PRE_Porites_evermani_scaffold_246_107990_108039_(+)_LOWconf
    AAGTATACTTGAAGTATACTTAATCGGAAGTATACCTCAAGTATACTTAA
    >Mir-1.PRE_Porites_evermani_scaffold_238_275734_275815_(+)_LOWconf
    GCATAATTCTTCGCATCTTAGGAGAGCCGAATTCAGAATAATTACTGAATTCAGCTCTCATAAGATGTGAAGAATTATGCAG
    >Mir-10.PRE_Porites_evermani_scaffold_49_151611_151659_(-)_HIGHconf
    CCGTAGATCCGAACTTGTGGGTATTTTCTCCACAGGTTGGGCTCTACGG

## 6.3 Check filtered GFF

Per [MirMachine](https://github.com/sinanugur/MirMachine)
recommendations:

> `filtered_gff/` High confidence miRNA family predictions after
> bitscore filtering. (This file is what you need in most cases)

So, we’ll stick with that.

``` bash
# Load bash variables into memory
source .bashvars

# Avoid 10th line - lengthy list of all miRNA families examined
head -n 9 "${output_dir_top}/results/predictions/filtered_gff/${mirmarchine_output_prefix}.PRE.gff"

echo ""
# Pull out lines that do _not_ begin with a '#'.
grep "^[^#]" "${output_dir_top}/results/predictions/filtered_gff/${mirmarchine_output_prefix}.PRE.gff" \
| head
```

    ##gff-version 3
    # MirMachine version: 0.2.12
    # CM Models: Built using MirGeneDB 2.1(2022)
    # Total families searched: 556
    # Node: Metazoa
    # Model: combined
    # Genome file: /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/E-Peve/data/Porites_evermanni_v1.fa
    # Species: Peve-MirMachine
    # Params: /home/sam/programs/mambaforge/envs/mirmachine_env/bin/MirMachine.py --node Metazoa --species Peve-MirMachine --genome /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/E-Peve/data/Porites_evermanni_v1.fa --cpu 40 --add-all-nodes

    Porites_evermani_scaffold_49    cmsearch    ncRNA   151611  151659  37.2    -   .   gene_id=Mir-10.PRE;E-value=0.0014;sequence_with_30nt=TGATTTAAAATTCTAGCAGCTCATGCGATCCCGTAGATCCGAACTTGTGGGTATTTTCTCCACAGGTTGGGCTCTACGGTCATATTTGCTGTAATATAACATGATGAGG
    Porites_evermani_scaffold_102   cmsearch    ncRNA   285408  285466  36.9    -   .   gene_id=Mir-154.PRE;E-value=0.0005;sequence_with_30nt=ATTTGGAAAAAGAGACATTCTGATGCGATAGCGGTTTTCCATATAAAGTTTACTTCCTTTATATTAAAGATATATGATTTACCTCACCTAAACACGTTCAAATAAATTTATCACTTACG
    Porites_evermani_scaffold_629   cmsearch    ncRNA   159499  159617  40.8    +   .   gene_id=Mir-216.PRE;E-value=3.7e-05;sequence_with_30nt=CTTTTAAACGCATGCGTCGCGTGAAATCGTTAATATCAGCAGATAATTCTGCCAAATTAATACGTTTAGGTGTGGTCCTCGCTCGCGGCTTAACCTTTTGCCCCACACCTAAACGTATTAATTTGGCAGAATGATCTGCTGATATTAAAACATATCTCCTTGCTCGTGGTTATAATTAT
    Porites_evermani_scaffold_522   cmsearch    ncRNA   105208  105268  33.3    +   .   gene_id=Mir-2191.PRE;E-value=0.0073;sequence_with_30nt=GAAATTGTAGGATCACACTCATAAATATGCATGTATTTATGTATGTATGCGTGTATGTATGTACATACACGCATACATACATGGTTTTTATAGTATATCAACAATTTCCAAATCCATCCCC
    Porites_evermani_scaffold_452   cmsearch    ncRNA   47661   47728   37.2    +   .   gene_id=Mir-303.PRE;E-value=0.0043;sequence_with_30nt=CCGTGTTAACTTTTGTGGTAGAAAAAATTTTTAGTTTTCCAAGGTGTCTAGTTTTTTTATAATAATTATAAAAAAACTAGACACCTTGGAAAACTAAAATTTTTTTTTTACGGTCTTTTTACATATAA
    Porites_evermani_scaffold_1642  cmsearch    ncRNA   31391   31446   37.9    -   .   gene_id=Mir-303.PRE;E-value=0.0027;sequence_with_30nt=AGGTTAAGCCCTCTCCCGATCAACAAAATGCAAGTTATCAAATGTAAGAAGTTTCATCACAAAACTTCTTACATTTGATAACTTGCATTTTGTTCCCGCCACTACGACATTCTCGC
    Porites_evermani_scaffold_79    cmsearch    ncRNA   110132  110208  42.2    +   .   gene_id=Mir-430.PRE;E-value=2.2e-05;sequence_with_30nt=AAATGCTTCTTACTTGTGTTTTTTTTAAACATCCCAACAGTGAAACACTTATTGTTGATTGTACAATTTTTGTACAATCAATAATAAGTGTTTCACTGTTAGGATGTTTAGAGCGATTTCGTATGACTTGAAATGAA
    Porites_evermani_scaffold_1000  cmsearch    ncRNA   149006  149065  36.3    -   .   gene_id=Mir-430.PRE;E-value=0.0011;sequence_with_30nt=GACCAAAGATAAAGGGAATTCGCATATCAAGTTTTGGAAAAGAAGGCGCTTTCTTTCGTGAAAGAAAGAAAGTGCTTTCTTTTGGAAAAAACAAGCACTTATCACCGTTTCGAGGCTCTC
    Porites_evermani_scaffold_132   cmsearch    ncRNA   263628  263697  47.8    -   .   gene_id=Mir-4983.PRE;E-value=1.4e-05;sequence_with_30nt=TTCATTTGCAGCATACTCTCATAACCGTCGCATTTGCTCGCTGGCATTTTTTATTGTTTTAAGCTCAAAAAACAATAAAAAATTCCAATGGGCAAATGCGACGCTTTTACGCGAATATTGCAAATGAATC
    Porites_evermani_scaffold_846   cmsearch    ncRNA   27006   27073   36.6    -   .   gene_id=Mir-4983.PRE;E-value=0.0098;sequence_with_30nt=TTCATTTGCAACATACGCGTATAGCTGTCGCAATTGCCCGTTGGTATTTTTTTGTATTTAAACGTAAAAACAATAAACAATTCTAGTGGGCAAATGCGACAGCTACACGCGTATTTTTCAAATGAATC
    grep: write error: Broken pipe

## 6.4 Counts of predicted high-confidence miRNAs

``` bash
# Load bash variables into memory
source .bashvars

# Skip lines which begin with '#' (i.e. the header components)
echo "Predicted loci:"
grep -c "^[^#]" "${output_dir_top}/results/predictions/filtered_gff/${mirmarchine_output_prefix}.PRE.gff"

echo ""
echo "Unique miRNA families:"
grep "^[^#]" "${output_dir_top}/results/predictions/filtered_gff/${mirmarchine_output_prefix}.PRE.gff" \
| awk -F"[\t=;]" '{print $10}' \
| sort -u \
| wc -l
```

    Predicted loci:
    83

    Unique miRNA families:
    15

------------------------------------------------------------------------

# 7 Citations

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-umu2022" class="csl-entry">

Umu, Sinan Uğur, Vanessa M. Paynter, Håvard Trondsen, Tilo Buschmann,
Trine B. Rounge, Kevin J. Peterson, and Bastian Fromm. 2022. “Accurate
microRNA Annotation of Animal Genomes Using Trained Covariance Models of
Curated microRNA Complements in MirMachine.”
<http://dx.doi.org/10.1101/2022.11.23.517654>.

</div>

</div>
