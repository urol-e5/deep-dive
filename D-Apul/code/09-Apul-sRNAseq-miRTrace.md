09-Apul-sRNAseq-miRTrace
================
Sam White
2023-11-01

- <a href="#1-create-a-bash-variables-file"
  id="toc-1-create-a-bash-variables-file">1 Create a Bash variables
  file</a>
- <a href="#2-create-mirtrace-config-file"
  id="toc-2-create-mirtrace-config-file">2 Create miRTrace config file</a>
- <a href="#3-run-mirtrace" id="toc-3-run-mirtrace">3 Run miRTrace</a>
- <a href="#4-results" id="toc-4-results">4 Results</a>
  - <a href="#41-read-in-table-as-data-frame"
    id="toc-41-read-in-table-as-data-frame">4.1 Read in table as data
    frame</a>
  - <a href="#42-number-of-samples-with-matches"
    id="toc-42-number-of-samples-with-matches">4.2 Number of samples with
    matches</a>
  - <a href="#43-percentage-of-samples-with-matches"
    id="toc-43-percentage-of-samples-with-matches">4.3 Percentage of samples
    with matches</a>
  - <a href="#44-number-of-clades-with-matches"
    id="toc-44-number-of-clades-with-matches">4.4 Number of clades with
    matches</a>
  - <a href="#45-mirtrace-table" id="toc-45-mirtrace-table">4.5 miRTrace
    table</a>
- <a href="#5-citations" id="toc-5-citations">5 Citations</a>

Use [miRTrace](https://github.com/friedlanderlab/mirtrace) \[@kang2018\]
to identify taxonomic origins of miRNA sequencing data.

NOTE: This requires you to have previously run
[`08-Apul-sRNAseq-trimming.Rmd`](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/code/08-Apul-sRNAseq-trimming.Rmd),
as the code relies on the trimmed reads output from that code.

------------------------------------------------------------------------

Inputs:

- Trimmed sRNAseq FastQs generated by
  [`08-Apul-sRNAseq-trimming.Rmd`](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/code/08-Apul-sRNAseq-trimming.Rmd)

  - Filenames formatted: `*flexbar_trim.25bp*.gz`

Outputs:

- `mirtrace.config`: A
  [miRTrace](https://github.com/friedlanderlab/mirtrace) config file. A
  comma-separated file with this layout (one FastQ per line):
  `/path/to/fastq,custom_sample_name`

- “Collapsed” (i.e. unique sequences only) FastA for each corresponding
  input FastQ.

- `mirtrace-report.html`: HTML-formatted report generated by
  [miRTrace](https://github.com/friedlanderlab/mirtrace).

- `mirtrace-stats-contamination_basic.tsv`: Tab-delimited report with
  counts of sequences from each collapsed FastAs having matches to known
  miRNAs within each of the
  [miRTrace](https://github.com/friedlanderlab/mirtrace) Clades.

- `mirtrace-stats-contamination_detailed.csv`: Tab-delimited report of
  *only* Clades with which sequences were matched, along with the
  corresponding miRNA families in each clade, and the sequence counts.

# 1 Create a Bash variables file

This allows usage of Bash variables across R Markdown chunks.

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "# Data directories"
echo 'export deep_dive_dir=/home/shared/8TB_HDD_01/sam/gitrepos/deep-dive'
echo 'export output_dir_top=${deep_dive_dir}/D-Apul/output/09-Apul-sRNAseq-miRTrace'
echo 'export trimmed_reads_dir=${deep_dive_dir}/D-Apul/output/08-Apul-sRNAseq-trimming/trimmed-reads'
echo ""

echo "# Paths to programs"
echo 'export mirtrace=/home/sam/programs/mambaforge/envs/miRTrace_env/bin/mirtrace'
echo ""

echo "# Set number of CPUs to use"
echo 'export threads=40'
echo ""

echo "export fastq_pattern='*flexbar_trim.25bp*.gz'"

echo "# Programs associative array"
echo "declare -A programs_array"
echo "programs_array=("
echo '[mirtrace]="${mirtrace}"'
echo ")"
} > .bashvars

cat .bashvars
```

    #### Assign Variables ####

    # Data directories
    export deep_dive_dir=/home/shared/8TB_HDD_01/sam/gitrepos/deep-dive
    export output_dir_top=${deep_dive_dir}/D-Apul/output/09-Apul-sRNAseq-miRTrace
    export trimmed_reads_dir=${deep_dive_dir}/D-Apul/output/08-Apul-sRNAseq-trimming/trimmed-reads

    # Paths to programs
    export mirtrace=/home/sam/programs/mambaforge/envs/miRTrace_env/bin/mirtrace

    # Set number of CPUs to use
    export threads=40

    export fastq_pattern='*flexbar_trim.25bp*.gz'
    # Programs associative array
    declare -A programs_array
    programs_array=(
    [mirtrace]="${mirtrace}"
    )

# 2 Create [miRTrace](https://github.com/friedlanderlab/mirtrace) config file

``` bash
# Load bash variables into memory
source .bashvars

# Declare array
fastq_array=()

# Populate array
fastq_array=(${trimmed_reads_dir}/${fastq_pattern})

# Loop through read pairs
# Increment by 2 to process next pair of FastQ files
if [ -f "${output_dir_top}/mirtrace.config" ]; then
  echo "mirtrace.config already exists. Nothing to do."
  
else

  for (( i=0; i<${#fastq_array[@]} ; i+=2 ))
  do
    # Use first three parts of filename to create short sample name
    R1_name=$(echo "${fastq_array[i]##*/}" | awk -F "-" '{print $1"-"$2"-"$3}')
    R2_name=$(echo "${fastq_array[i+1]##*/}" | awk -F "-" '{print $1"-"$2"-"$3}')
    echo "${fastq_array[i]},${R1_name}_1"
    echo "${fastq_array[i+1]},${R2_name}_2"
  done >> "${output_dir_top}/mirtrace.config"

fi

cat "${output_dir_top}/mirtrace.config"
```

    mirtrace.config already exists. Nothing to do.
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/08-Apul-sRNAseq-trimming/trimmed-reads/sRNA-ACR-140-S1-TP2.flexbar_trim.25bp_1.fastq.gz,sRNA-ACR-140_1
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/08-Apul-sRNAseq-trimming/trimmed-reads/sRNA-ACR-140-S1-TP2.flexbar_trim.25bp_2.fastq.gz,sRNA-ACR-140_2
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/08-Apul-sRNAseq-trimming/trimmed-reads/sRNA-ACR-145-S1-TP2.flexbar_trim.25bp_1.fastq.gz,sRNA-ACR-145_1
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/08-Apul-sRNAseq-trimming/trimmed-reads/sRNA-ACR-145-S1-TP2.flexbar_trim.25bp_2.fastq.gz,sRNA-ACR-145_2
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/08-Apul-sRNAseq-trimming/trimmed-reads/sRNA-ACR-150-S1-TP2.flexbar_trim.25bp_1.fastq.gz,sRNA-ACR-150_1
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/08-Apul-sRNAseq-trimming/trimmed-reads/sRNA-ACR-150-S1-TP2.flexbar_trim.25bp_2.fastq.gz,sRNA-ACR-150_2
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/08-Apul-sRNAseq-trimming/trimmed-reads/sRNA-ACR-173-S1-TP2.flexbar_trim.25bp_1.fastq.gz,sRNA-ACR-173_1
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/08-Apul-sRNAseq-trimming/trimmed-reads/sRNA-ACR-173-S1-TP2.flexbar_trim.25bp_2.fastq.gz,sRNA-ACR-173_2
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/08-Apul-sRNAseq-trimming/trimmed-reads/sRNA-ACR-178-S1-TP2.flexbar_trim.25bp_1.fastq.gz,sRNA-ACR-178_1
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/08-Apul-sRNAseq-trimming/trimmed-reads/sRNA-ACR-178-S1-TP2.flexbar_trim.25bp_2.fastq.gz,sRNA-ACR-178_2

# 3 Run [miRTrace](https://github.com/friedlanderlab/mirtrace)

``` bash
# Load bash variables into memory
source .bashvars

time \
${programs_array[mirtrace]} trace \
--config ${output_dir_top}/mirtrace.config \
--write-fasta \
--num-threads ${threads} \
--output-dir ${output_dir_top} \
--force

tree -h ${output_dir_top}
```

    miRTrace version 1.0.1 starting. Processing 10 sample(s).
    NOTE: reusing existing output directory, outdated files may be present.
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/09-Apul-sRNAseq-miRTrace/qc_passed_reads.all.collapsed/sRNA-ACR-173-S1-TP2.flexbar_trim.25bp_1.fasta (Permission denied)
    I/O Error, aborting.

    real    0m20.432s
    user    3m18.809s
    sys 0m3.585s
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/09-Apul-sRNAseq-miRTrace
    ├── [1.6K]  mirtrace.config
    ├── [302K]  mirtrace-report.html
    ├── [ 575]  mirtrace-stats-contamination_basic.tsv
    ├── [ 822]  mirtrace-stats-contamination_detailed.tsv
    ├── [4.0K]  qc_passed_reads.all.collapsed
    │   ├── [117M]  sRNA-ACR-140-S1-TP2.flexbar_trim.25bp_1.fasta
    │   ├── [163M]  sRNA-ACR-140-S1-TP2.flexbar_trim.25bp_2.fasta
    │   ├── [128M]  sRNA-ACR-145-S1-TP2.flexbar_trim.25bp_1.fasta
    │   ├── [175M]  sRNA-ACR-145-S1-TP2.flexbar_trim.25bp_2.fasta
    │   ├── [135M]  sRNA-ACR-150-S1-TP2.flexbar_trim.25bp_1.fasta
    │   ├── [179M]  sRNA-ACR-150-S1-TP2.flexbar_trim.25bp_2.fasta
    │   ├── [113M]  sRNA-ACR-173-S1-TP2.flexbar_trim.25bp_1.fasta
    │   ├── [151M]  sRNA-ACR-173-S1-TP2.flexbar_trim.25bp_2.fasta
    │   ├── [109M]  sRNA-ACR-178-S1-TP2.flexbar_trim.25bp_1.fasta
    │   └── [140M]  sRNA-ACR-178-S1-TP2.flexbar_trim.25bp_2.fasta
    └── [1.4K]  README.md

    1 directory, 15 files

# 4 Results

## 4.1 Read in table as data frame

``` r
mirtrace.detailed.df <- read.csv("../output/09-Apul-sRNAseq-miRTrace/mirtrace-stats-contamination_detailed.tsv", sep = "\t", header = TRUE)

str(mirtrace.detailed.df)
```

    'data.frame':   7 obs. of  14 variables:
     $ CLADE         : chr  "lophotrochozoa" "lophotrochozoa" "lophotrochozoa" "rodents" ...
     $ FAMILY_ID     : int  1994 1985 1984 351 618 576 576
     $ MIRBASE_IDS   : chr  "cla-miR-1994,cte-miR-1994,lgi-miR-1994a,lgi-miR-1994b" "hru-miR-1985,lgi-miR-1985" "hru-miR-1984,lgi-miR-1984" "mmu-miR-351,rno-miR-351" ...
     $ SEQ           : chr  "TGAGACAGTGTGTCCTCCCT" "TGCCATTTTTATCAGTCACT" "TGCCCTATCCGTCAGGAACT" "TCCCTGAGGAGCCCTTTGAG" ...
     $ sRNA.ACR.140_1: int  0 0 0 0 0 0 0
     $ sRNA.ACR.140_2: int  0 0 0 0 0 0 0
     $ sRNA.ACR.145_1: int  0 0 0 0 0 0 0
     $ sRNA.ACR.145_2: int  0 0 0 0 0 0 0
     $ sRNA.ACR.150_1: int  1 11 10 0 0 1 0
     $ sRNA.ACR.150_2: int  0 0 0 0 0 0 0
     $ sRNA.ACR.173_1: int  3 15 16 0 1 0 1
     $ sRNA.ACR.173_2: int  0 0 0 0 0 0 0
     $ sRNA.ACR.178_1: int  0 0 0 2 0 0 0
     $ sRNA.ACR.178_2: int  0 0 0 0 0 0 0

## 4.2 Number of samples with matches

``` r
# Select columns corresponding to sample names
sample_columns <- mirtrace.detailed.df %>%
  select(starts_with("sRNA.ACR."))

# Calculate the sum for each column
sample_sums <- colSums(sample_columns)

# Count the number of columns with a sum greater than 0
samples_with_sum_gt_0 <- sum(sample_sums > 0)

paste("Number of samples with matches: ", samples_with_sum_gt_0)
```

    [1] "Number of samples with matches:  3"

## 4.3 Percentage of samples with matches

``` r
# Total number of samples (columns)
total_samples <- ncol(sample_columns)

# Percentage of samples with sums greater than 0
percentage_samples_gt_0 <- (samples_with_sum_gt_0 / total_samples) * 100

paste("Percentage of samples with matches: ", percentage_samples_gt_0)
```

    [1] "Percentage of samples with matches:  30"

## 4.4 Number of clades with matches

``` r
unique_clade_count <- mirtrace.detailed.df %>%
  distinct(CLADE) %>%    # Get unique entries in CLADE column
  count()               # Count the number of unique entries



paste("Number of clades with matches:", unique_clade_count)
```

    [1] "Number of clades with matches: 3"

## 4.5 [miRTrace](https://github.com/friedlanderlab/mirtrace) table

To make them easier to see, counts \> 0 are highlighted in green.

``` r
mirtrace.detailed.df %>%
  mutate(
    across(
      starts_with("sRNA"),
      ~cell_spec(
        .,
        background = ifelse(
          . > 0,
          "lightgreen",
          "white"
          )
        )
      )
    ) %>%
  kable(escape = F, caption = "Clades identified as having sRNAseq matches.") %>%
  kable_styling("striped") %>% 
  scroll_box(width = "100%", height = "500px")
```

<div
style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:500px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<caption>

Clades identified as having sRNAseq matches.

</caption>
<thead>
<tr>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

CLADE

</th>
<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

FAMILY_ID

</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

MIRBASE_IDS

</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

SEQ

</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

sRNA.ACR.140_1

</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

sRNA.ACR.140_2

</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

sRNA.ACR.145_1

</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

sRNA.ACR.145_2

</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

sRNA.ACR.150_1

</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

sRNA.ACR.150_2

</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

sRNA.ACR.173_1

</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

sRNA.ACR.173_2

</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

sRNA.ACR.178_1

</th>
<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

sRNA.ACR.178_2

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

lophotrochozoa

</td>
<td style="text-align:right;">

1994

</td>
<td style="text-align:left;">

cla-miR-1994,cte-miR-1994,lgi-miR-1994a,lgi-miR-1994b

</td>
<td style="text-align:left;">

TGAGACAGTGTGTCCTCCCT

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[1\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: lightgreen !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[3\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: lightgreen !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
</tr>
<tr>
<td style="text-align:left;">

lophotrochozoa

</td>
<td style="text-align:right;">

1985

</td>
<td style="text-align:left;">

hru-miR-1985,lgi-miR-1985

</td>
<td style="text-align:left;">

TGCCATTTTTATCAGTCACT

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[11\]{style=” border-radius: 4px; padding-right: 4px; padding-left:
4px; background-color: lightgreen !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[15\]{style=” border-radius: 4px; padding-right: 4px; padding-left:
4px; background-color: lightgreen !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
</tr>
<tr>
<td style="text-align:left;">

lophotrochozoa

</td>
<td style="text-align:right;">

1984

</td>
<td style="text-align:left;">

hru-miR-1984,lgi-miR-1984

</td>
<td style="text-align:left;">

TGCCCTATCCGTCAGGAACT

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[10\]{style=” border-radius: 4px; padding-right: 4px; padding-left:
4px; background-color: lightgreen !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[16\]{style=” border-radius: 4px; padding-right: 4px; padding-left:
4px; background-color: lightgreen !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
</tr>
<tr>
<td style="text-align:left;">

rodents

</td>
<td style="text-align:right;">

351

</td>
<td style="text-align:left;">

mmu-miR-351,rno-miR-351

</td>
<td style="text-align:left;">

TCCCTGAGGAGCCCTTTGAG

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[2\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: lightgreen !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
</tr>
<tr>
<td style="text-align:left;">

primates

</td>
<td style="text-align:right;">

618

</td>
<td style="text-align:left;">

hsa-miR-618,mml-miR-618,ppy-miR-618,ptr-miR-618

</td>
<td style="text-align:left;">

AAACTCTACTTGTCCTTCTG

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[1\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: lightgreen !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
</tr>
<tr>
<td style="text-align:left;">

primates

</td>
<td style="text-align:right;">

576

</td>
<td style="text-align:left;">

hsa-miR-576,mml-miR-576,ppy-miR-576,ptr-miR-576

</td>
<td style="text-align:left;">

AAGATGTGGAAAAATTGGAA

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[1\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: lightgreen !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
</tr>
<tr>
<td style="text-align:left;">

primates

</td>
<td style="text-align:right;">

576

</td>
<td style="text-align:left;">

hsa-miR-576

</td>
<td style="text-align:left;">

ATTCTAATTTCTCCACGTCT

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[1\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: lightgreen !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
<td style="text-align:left;">

\[0\]{style=” border-radius: 4px; padding-right: 4px; padding-left: 4px;
background-color: white !important;“}

</td>
</tr>
</tbody>
</table>

</div>

# 5 Citations