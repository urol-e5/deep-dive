12-Peve-RNAseq-kallisto
================
Kathleen Durkin
2024-01-25

- <a href="#1-create-a-bash-variables-file"
  id="toc-1-create-a-bash-variables-file">1 Create a Bash variables
  file</a>
- <a href="#2-download-trimmed-rnaseq-reads"
  id="toc-2-download-trimmed-rnaseq-reads">2 Download trimmed RNAseq
  reads</a>
- <a href="#3-fastqcmultiqc-on-trimmed-reads"
  id="toc-3-fastqcmultiqc-on-trimmed-reads">3 FastQC/MultiQC on trimmed
  reads</a>
- <a href="#4-retrieve-the-reference-transcriptome-and-genome"
  id="toc-4-retrieve-the-reference-transcriptome-and-genome">4 Retrieve
  the reference transcriptome and genome</a>
  - <a href="#41-verify-transcriptomegenome-fasta-md5-checksum"
    id="toc-41-verify-transcriptomegenome-fasta-md5-checksum">4.1 Verify
    transcriptome/genome FastA MD5 checksum</a>
- <a href="#5-convert-gff-to-a-genes-fasta"
  id="toc-5-convert-gff-to-a-genes-fasta">5 Convert gff to a genes
  FASTA</a>
  - <a href="#51-extract-only-cds-from-gff"
    id="toc-51-extract-only-cds-from-gff">5.1 Extract only CDS from gff</a>
  - <a href="#52-convert-gff-to-bed-file"
    id="toc-52-convert-gff-to-bed-file">5.2 Convert gff to bed file</a>
  - <a href="#53-generate-transcriptome-fasta"
    id="toc-53-generate-transcriptome-fasta">5.3 Generate transcriptome
    fasta</a>
  - <a href="#54-check-transcriptome-fasta"
    id="toc-54-check-transcriptome-fasta">5.4 Check transcriptome fasta</a>
- <a href="#6-align-to-reference-transcriptome-kallisto-pseudoalignment"
  id="toc-6-align-to-reference-transcriptome-kallisto-pseudoalignment">6
  Align to reference transcriptome (Kallisto pseudoalignment)</a>
  - <a href="#61-building-index" id="toc-61-building-index">6.1 Building
    Index</a>
  - <a href="#62-sample-quantification"
    id="toc-62-sample-quantification">6.2 Sample Quantification</a>
  - <a href="#63-trinity-matrix-with-kallisto-output"
    id="toc-63-trinity-matrix-with-kallisto-output">6.3 Trinity Matrix with
    Kallisto Output</a>
- <a href="#7-summary" id="toc-7-summary">7 Summary</a>

**Description:**

Use kallisto to quantify *P. evermanni* RNAseq transcript abundances

**Inputs:**

Trimmed RNAseq reads (e.g. `*.fastq.gz`). *P. evermanni* reads found
[here](https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_evermanni/trimmed/).

Coding sequence gff file (e.g. `*.gff`) and scaffold genome
(e.g. `*.fa`), found [here](https://www.genoscope.cns.fr/corals/data/).

**Outputs:**

kalisto counts matrix (`kallisto.isoform.counts.matrix`)

------------------------------------------------------------------------

# 1 Create a Bash variables file

This allows usage of Bash variables (e.g. paths to common directories)
across R Markdown chunks.

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "# Data directories"
echo 'export Peve_dir=/home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve'
echo 'export output_dir_top=${Peve_dir}/output/12-Peve-RNAseq-kallisto'
echo 'export trimmed_reads_dir=${output_dir_top}/trimmed-reads'
echo 'export kallisto_output_dir=${output_dir_top}/kallisto'
echo ""

echo "# Input/Output files"
echo 'export transcriptome_dir=${Peve_dir}/data'
echo 'export transcriptome_gff_name="Porites_evermanni_v1.annot.gff"'
echo 'export transcriptome_gff=${transcriptome_dir}/${transcriptome_gff_name}'
echo 'export transcriptome_gff_filtered_name="Porites_evermanni_v1_CDS.annot.gff"'
echo 'export transcriptome_gff_filtered=${transcriptome_dir}/${transcriptome_gff_filtered_name}'
echo 'export transcriptome_bed_name="Porites_evermanni_v1_CDS.annot.bed"'
echo 'export transcriptome_bed=${transcriptome_dir}/${transcriptome_bed_name}'
echo 'export genome_fasta_name="Porites_evermanni_v1.fa"'
echo 'export genome_fasta=${transcriptome_dir}/${genome_fasta_name}'
echo 'export transcriptome_fasta_name="Porites_evermanni_CDS.fasta"'
echo 'export transcriptome_fasta=${transcriptome_dir}/${transcriptome_fasta_name}'
echo 'export kallisto_index_name="Peve_kallisto_index.idx"'

echo "# External data URLs"
echo 'export trimmed_reads_url="https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_evermanni/trimmed/"'
echo 'export transcriptome_url="https://www.genoscope.cns.fr/corals/data/Porites_evermanni_v1.annot.gff"'
echo 'export genome_url="https://www.genoscope.cns.fr/corals/data/Porites_evermanni_v1.fa"'
echo ""

echo "# Set filename patterns"
echo "export fastq_pattern='*.fastq.gz'"
echo "export R1_fastq_pattern='*_R1.fastq.gz'"
echo "export R2_fastq_pattern='*_R2.fastq.gz'"
echo ""

echo "# Paths to programs"
echo 'export kallisto=/home/shared/kallisto_linux-v0.50.1/kallisto'
echo 'export trinity_abund_to_matrix=/home/shared/trinityrnaseq-v2.12.0/util/abundance_estimates_to_matrix.pl'
echo 'export bedtools=/home/shared/bedtools2/bin/bedtools'
echo 'export bedops=/home/shared/bedops_linux_x86_64-v2.4.41/bin'
echo ""

echo "# Set number of CPUs to use"
echo 'export threads=20'
echo ""

echo "# Programs associative array"
echo "declare -A programs_array"
echo "programs_array=("
echo '[kallisto]="${kallisto}" \'
echo '[trinity_abund_to_matrix]="${trinity_abund_to_matrix}" \'
echo '[bedtools]="${bedtools}" \'
echo '[bedops]="${bedops}" \'
echo ")"

} > .bashvars

cat .bashvars
```

    #### Assign Variables ####

    # Data directories
    export Peve_dir=/home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve
    export output_dir_top=${Peve_dir}/output/12-Peve-RNAseq-kallisto
    export trimmed_reads_dir=${output_dir_top}/trimmed-reads
    export kallisto_output_dir=${output_dir_top}/kallisto

    # Input/Output files
    export transcriptome_dir=${Peve_dir}/data
    export transcriptome_gff_name="Porites_evermanni_v1.annot.gff"
    export transcriptome_gff=${transcriptome_dir}/${transcriptome_gff_name}
    export transcriptome_gff_filtered_name="Porites_evermanni_v1_CDS.annot.gff"
    export transcriptome_gff_filtered=${transcriptome_dir}/${transcriptome_gff_filtered_name}
    export transcriptome_bed_name="Porites_evermanni_v1_CDS.annot.bed"
    export transcriptome_bed=${transcriptome_dir}/${transcriptome_bed_name}
    export genome_fasta_name="Porites_evermanni_v1.fa"
    export genome_fasta=${transcriptome_dir}/${genome_fasta_name}
    export transcriptome_fasta_name="Porites_evermanni_CDS.fasta"
    export transcriptome_fasta=${transcriptome_dir}/${transcriptome_fasta_name}
    export kallisto_index_name="Peve_kallisto_index.idx"
    # External data URLs
    export trimmed_reads_url="https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_evermanni/trimmed/"
    export transcriptome_url="https://www.genoscope.cns.fr/corals/data/Porites_evermanni_v1.annot.gff"
    export genome_url="https://www.genoscope.cns.fr/corals/data/Porites_evermanni_v1.fa"

    # Set filename patterns
    export fastq_pattern='*.fastq.gz'
    export R1_fastq_pattern='*_R1.fastq.gz'
    export R2_fastq_pattern='*_R2.fastq.gz'

    # Paths to programs
    export kallisto=/home/shared/kallisto_linux-v0.50.1/kallisto
    export trinity_abund_to_matrix=/home/shared/trinityrnaseq-v2.12.0/util/abundance_estimates_to_matrix.pl
    export bedtools=/home/shared/bedtools2/bin/bedtools
    export bedops=/home/shared/bedops_linux_x86_64-v2.4.41/bin

    # Set number of CPUs to use
    export threads=20

    # Programs associative array
    declare -A programs_array
    programs_array=(
    [kallisto]="${kallisto}" \
    [trinity_abund_to_matrix]="${trinity_abund_to_matrix}" \
    [bedtools]="${bedtools}" \
    [bedops]="${bedops}" \
    )

# 2 Download trimmed RNAseq reads

Reads are downloaded from:
<https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_evermanni/trimmed/>

The `--cut-dirs 4` command cuts the preceding directory structure
(i.e. `Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_evermanni/trimmed/`)
so that we just end up with the reads.

``` bash
# Load bash variables into memory
source .bashvars

wget \
--directory-prefix ${trimmed_reads_dir} \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 4 \
--no-host-directories \
--no-parent \
--quiet \
--accept "RNA-*.fastq.gz" ${trimmed_reads_url}
```

``` bash
# Load bash variables into memory
source .bashvars

ls -lh "${trimmed_reads_dir}"
```

    total 24G
    drwxr-xr-x 2 shedurkin labmembers 4.0K Jan 25 15:20 multiqc_data
    -rw-r--r-- 1 shedurkin labmembers 2.5G May 19  2023 RNA-POR-71-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.6G May 19  2023 RNA-POR-71-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.2G May 19  2023 RNA-POR-73-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.4G May 19  2023 RNA-POR-73-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.4G May 19  2023 RNA-POR-76-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.5G May 19  2023 RNA-POR-76-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.1G May 19  2023 RNA-POR-79-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.3G May 19  2023 RNA-POR-79-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.6G May 19  2023 RNA-POR-82-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.7G May 19  2023 RNA-POR-82-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Feb  8 08:49 sample71_R1.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve/output/12-Peve-RNAseq-kallisto/trimmed-reads/RNA-POR-71-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Feb  8 08:49 sample71_R2.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve/output/12-Peve-RNAseq-kallisto/trimmed-reads/RNA-POR-71-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Feb  8 08:49 sample73_R1.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve/output/12-Peve-RNAseq-kallisto/trimmed-reads/RNA-POR-73-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Feb  8 08:49 sample73_R2.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve/output/12-Peve-RNAseq-kallisto/trimmed-reads/RNA-POR-73-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Feb  8 08:49 sample76_R1.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve/output/12-Peve-RNAseq-kallisto/trimmed-reads/RNA-POR-76-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Feb  8 08:49 sample76_R2.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve/output/12-Peve-RNAseq-kallisto/trimmed-reads/RNA-POR-76-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Feb  8 08:49 sample79_R1.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve/output/12-Peve-RNAseq-kallisto/trimmed-reads/RNA-POR-79-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Feb  8 08:49 sample79_R2.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve/output/12-Peve-RNAseq-kallisto/trimmed-reads/RNA-POR-79-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Feb  8 08:49 sample82_R1.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve/output/12-Peve-RNAseq-kallisto/trimmed-reads/RNA-POR-82-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Feb  8 08:49 sample82_R2.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve/output/12-Peve-RNAseq-kallisto/trimmed-reads/RNA-POR-82-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz

# 3 FastQC/MultiQC on trimmed reads

Already performed, can view multiqc report at
<https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_evermanni/trimmed/multiqc_report.html>

# 4 Retrieve the reference transcriptome and genome

Provided by <https://www.genoscope.cns.fr/corals/genomes.html>

Download the gene CDS (coding sequence) gff file

``` bash
# Load bash variables into memory
source .bashvars

wget \
--directory-prefix ${transcriptome_dir} \
--no-check-certificate \
--continue \
--no-host-directories \
--no-directories \
--no-parent \
--quiet \
--execute robots=off \
--accept "${transcriptome_gff_name}" ${transcriptome_url}
```

Note that this is a CDS (coding sequence) gff file, not a FASTA, so
can’t input directly into kallisto. We’ll need the reference genome as
well to convert gff to FASTA

Download the genome scaffolds FASTA file

``` bash
# Load bash variables into memory
source .bashvars

wget \
--directory-prefix ${transcriptome_dir} \
--no-check-certificate \
--continue \
--no-host-directories \
--no-directories \
--no-parent \
--quiet \
--execute robots=off \
--accept "${genome_fasta_name}" ${genome_url}
```

``` bash
# Load bash variables into memory
source .bashvars

ls -lh "${transcriptome_dir}"
```

    total 774M
    -rw------- 1 shedurkin labmembers 3.2K Jan 25 15:46 index.html.tmp
    -rw-r--r-- 1 shedurkin labmembers  15M Nov  2 08:39 peve_bedtools_lncRNAs.fasta
    -rw-r--r-- 1 shedurkin labmembers  57M Feb  7 22:12 Porites_evermanni_CDS.fasta
    -rw-r--r-- 1 shedurkin labmembers  57M Feb  6 17:32 Porites_evermanni_CDS.from_gff.fasta
    -rw-r--r-- 1 shedurkin labmembers  24M Mar 11  2022 Porites_evermanni_v1.annot.gff
    -rw-r--r-- 1 shedurkin labmembers  19M Feb  7 17:27 Porites_evermanni_v1_CDS.annot.bed
    -rw-r--r-- 1 shedurkin labmembers  18M Feb  7 17:16 Porites_evermanni_v1_CDS.annot.gff
    -rw-r--r-- 1 shedurkin labmembers 586M Mar 11  2022 Porites_evermanni_v1.fa
    -rw-r--r-- 1 shedurkin labmembers 422K Jan 29 16:55 Porites_evermanni_v1.fa.fai
    -rw-r--r-- 1 shedurkin labmembers    0 Nov  2 08:39 README.md

## 4.1 Verify transcriptome/genome FastA MD5 checksum

No checksum file(s) provided with download, so skipping this step

# 5 Convert gff to a genes FASTA

## 5.1 Extract only CDS from gff

We only want the sequences classified as “CDS”

``` bash
# Load bash variables into memory
source .bashvars

# Extract only the CDS sequence lines from the gff
grep -w 'CDS' ${transcriptome_gff} > ${transcriptome_gff_filtered}

head -n 5 ${transcriptome_gff_filtered}
```

    Porites_evermani_scaffold_1 Gmove   CDS 3107    3444    .   -   .   Parent=Peve_00000001
    Porites_evermani_scaffold_1 Gmove   CDS 4284    4488    .   -   .   Parent=Peve_00000001
    Porites_evermani_scaffold_1 Gmove   CDS 424479  425361  .   -   .   Parent=Peve_00000002
    Porites_evermani_scaffold_1 Gmove   CDS 426181  426735  .   -   .   Parent=Peve_00000002
    Porites_evermani_scaffold_1 Gmove   CDS 427013  427140  .   -   .   Parent=Peve_00000002

## 5.2 Convert gff to bed file

To extract FASTAs for each of the CDS sequences we just extracted we’ll
be using bedtools, which takes bed files as input, so we need to convert
our CDS gff file to a bed file. (Note that bedtools does accept gff
files, but since gff and bed files have slightly different coordinate
systems we’re going to convert bed just in case)

``` bash
# Load bash variables into memory
source .bashvars

# Ensure bedops can find its dependencies when running
export PATH=/home/shared/bedops_linux_x86_64-v2.4.41/bin:$PATH

${bedops}/gff2bed \
--do-not-sort \
< ${transcriptome_gff_filtered} \
> ${transcriptome_bed}

head -n 3 ${transcriptome_gff_filtered}
echo ""
head -n 3 ${transcriptome_bed}
```

    Porites_evermani_scaffold_1 Gmove   CDS 3107    3444    .   -   .   Parent=Peve_00000001
    Porites_evermani_scaffold_1 Gmove   CDS 4284    4488    .   -   .   Parent=Peve_00000001
    Porites_evermani_scaffold_1 Gmove   CDS 424479  425361  .   -   .   Parent=Peve_00000002

    Porites_evermani_scaffold_1 3106    3444    .   .   -   Gmove   CDS .   Parent=Peve_00000001
    Porites_evermani_scaffold_1 4283    4488    .   .   -   Gmove   CDS .   Parent=Peve_00000001
    Porites_evermani_scaffold_1 424478  425361  .   .   -   Gmove   CDS .   Parent=Peve_00000002

## 5.3 Generate transcriptome fasta

The below script will take as input a bed file containing information on
CDS sequences, where multiple CDS sequences may originate from the same
parent mRNA. The script will extract FASTAs for each sequence, and
concatenate and label by parent. This should output a gene FASTA that we
can use for kallisto pseudoalignment and abundance quantification!
**Warning**: This script will take a while to run – for our bed file of
231,320 CDS sequences, the script took \~4hours to output a complete
gene fasta.

``` bash
# Load bash variables into memory
source .bashvars

# Navigate to correct directory and make output file
cd ${transcriptome_dir}
echo > ${transcriptome_fasta_name}

# Helper list for processing all parent IDs
processed_ids=()

######################################################

# Helper function to concatenate and format several bedtools output sequences 
# into a single, appropriately named contig
concatenate_helper() {
    local input_bedtools_fastas="$1"
    local parent_ID="$2"
    local reference_name=""
    local positions=""
    local concatenated_sequences=""

    # Read the input line by line
    while IFS= read -r line; do
        # Check if the line starts with ">"
        if [[ "$line" == ">"* ]]; then
            # Extract reference name and position from the line
            reference_position="${line:1}"  # Remove ">"
            reference_name=$(echo "$reference_position" | cut -d: -f1)
            position=$(echo "$reference_position" | cut -d: -f2)

            # Append position to the positions variable
            positions+="$position,"
        else
            # Concatenate sequences
            concatenated_sequences+="$line"
        fi
    done <<< "$input_bedtools_fastas"

    # Remove trailing comma from positions
    positions="${positions%,}"

    # Output the reformatted result
    echo ">$parent_ID $reference_name:$positions"
    echo "$concatenated_sequences"
}

######################################################

# Process your input bed file
while IFS= read -r line; do

    # pull the parent ID number for the current line of the bed
    parentID=$(echo "$line" | grep -o 'Parent=Peve_[0-9]\+')
    
    # Only continue if you haven't already processed the CDS sequences associated with this parent ID
    if [[ ! " ${processed_ids[@]} " =~ " $parentID " ]]; then
 
        # Add the current parentID to the processed list
        processed_ids+=("$parentID")

        # Create temporary files to store intermediate results
        temp_CDS_bed_file=$(mktemp)
        temp_bedtools_fasta_file=$(mktemp)

        # Grab all of the CDS sequences with the same parent ID and write to temporary file
        grep "$parentID" ${transcriptome_bed} > "$temp_CDS_bed_file"

        # Use bedtools to extract corresponding FASTAs and write to temporary file
        ${programs_array[bedtools]} getfasta -fi ${genome_fasta} -bed "$temp_CDS_bed_file" -fo "$temp_bedtools_fasta_file"

        # Use our helper function to concatenate and format all of these CDS fastas into a single contig
        concatenated_fasta=$(concatenate_helper "$(cat "$temp_bedtools_fasta_file")" "$parentID")
 
        # Add the concatenated CDS fasta to our output file on a new line
        echo "$concatenated_fasta" >> ${transcriptome_fasta}

        # Remove the temporary files
        rm "$temp_CDS_bed_file" "$temp_bedtools_fasta_file"
    fi
done < ${transcriptome_bed}

# The output file ends up having a blank first line before the data, so delete that unwanted empty first line
sed -i '1{/^$/d}' ${transcriptome_fasta}

head -n 4 ${transcriptome_fasta}
```

## 5.4 Check transcriptome fasta

Let’s do a quick check to see whether the output file contains all the
CDS sequences we want and is grouping them appropriately

``` bash
# Load bash variables into memory
source .bashvars

# Output the first five fasta names in our transcriptome fasta (i.e. first five odd lines)
sed -nu '1~2p' ${transcriptome_fasta} | head -n 5

echo ""

# Output the first twenty CDS sequences listed in the CDS gff
head -n 20 ${transcriptome_bed}
```

    >Parent=Peve_00000001 Porites_evermani_scaffold_1:3106-3444,4283-4488
    >Parent=Peve_00000002 Porites_evermani_scaffold_1:424478-425361,426180-426735,427012-427140,427664-427724,428641-429034
    >Parent=Peve_00000003 Porites_evermani_scaffold_1:429499-429746,430884-431009,432043-432167,432627-432757,433482-433588,434246-434336,435358-435439,436216-436374,437429-437557,438130-438232,438531-438698
    >Parent=Peve_00000004 Porites_evermani_scaffold_1:441399-441851,442759-443100,446240-447172
    >Parent=Peve_00000005 Porites_evermani_scaffold_1:448044-448206,448308-448363,451515-451591,451816-451871,454946-455060,455886-456006,456781-456901,457021-457073,457639-457767,458842-458920,459554-459697,459961-460031

    Porites_evermani_scaffold_1 3106    3444    .   .   -   Gmove   CDS .   Parent=Peve_00000001
    Porites_evermani_scaffold_1 4283    4488    .   .   -   Gmove   CDS .   Parent=Peve_00000001
    Porites_evermani_scaffold_1 424478  425361  .   .   -   Gmove   CDS .   Parent=Peve_00000002
    Porites_evermani_scaffold_1 426180  426735  .   .   -   Gmove   CDS .   Parent=Peve_00000002
    Porites_evermani_scaffold_1 427012  427140  .   .   -   Gmove   CDS .   Parent=Peve_00000002
    Porites_evermani_scaffold_1 427664  427724  .   .   -   Gmove   CDS .   Parent=Peve_00000002
    Porites_evermani_scaffold_1 428641  429034  .   .   -   Gmove   CDS .   Parent=Peve_00000002
    Porites_evermani_scaffold_1 429499  429746  .   .   +   Gmove   CDS .   Parent=Peve_00000003
    Porites_evermani_scaffold_1 430884  431009  .   .   +   Gmove   CDS .   Parent=Peve_00000003
    Porites_evermani_scaffold_1 432043  432167  .   .   +   Gmove   CDS .   Parent=Peve_00000003
    Porites_evermani_scaffold_1 432627  432757  .   .   +   Gmove   CDS .   Parent=Peve_00000003
    Porites_evermani_scaffold_1 433482  433588  .   .   +   Gmove   CDS .   Parent=Peve_00000003
    Porites_evermani_scaffold_1 434246  434336  .   .   +   Gmove   CDS .   Parent=Peve_00000003
    Porites_evermani_scaffold_1 435358  435439  .   .   +   Gmove   CDS .   Parent=Peve_00000003
    Porites_evermani_scaffold_1 436216  436374  .   .   +   Gmove   CDS .   Parent=Peve_00000003
    Porites_evermani_scaffold_1 437429  437557  .   .   +   Gmove   CDS .   Parent=Peve_00000003
    Porites_evermani_scaffold_1 438130  438232  .   .   +   Gmove   CDS .   Parent=Peve_00000003
    Porites_evermani_scaffold_1 438531  438698  .   .   +   Gmove   CDS .   Parent=Peve_00000003
    Porites_evermani_scaffold_1 441399  441851  .   .   -   Gmove   CDS .   Parent=Peve_00000004
    Porites_evermani_scaffold_1 442759  443100  .   .   -   Gmove   CDS .   Parent=Peve_00000004

``` bash
# Load bash variables into memory
source .bashvars

# Output the first five fasta names in our transcriptome fasta (i.e. first five odd lines)
sed -nu '1~2p' ${transcriptome_fasta} | tail -n 5

echo ""

# Output the first twenty CDS sequences listed in the CDS gff
tail -n 20 ${transcriptome_bed}
```

    >Parent=Peve_00045353 Porites_evermani_scaffold_999:120490-120865
    >Parent=Peve_00045355 Porites_evermani_scaffold_999:124812-124870,124921-125022
    >Parent=Peve_00045356 Porites_evermani_scaffold_999:157406-158547,159689-159727
    >Parent=Peve_00045357 Porites_evermani_scaffold_999:25255-25390,26056-26534,26987-27076
    >Parent=Peve_00045358 Porites_evermani_scaffold_999:43038-43296,43486-43561,43917-44070

    Porites_evermani_scaffold_999   83896   84640   .   .   +   Gmove   CDS .   Parent=Peve_00045348
    Porites_evermani_scaffold_999   102184  102529  .   .   -   Gmove   CDS .   Parent=Peve_00045349
    Porites_evermani_scaffold_999   102652  103198  .   .   -   Gmove   CDS .   Parent=Peve_00045349
    Porites_evermani_scaffold_999   109634  109938  .   .   -   Gmove   CDS .   Parent=Peve_00045350
    Porites_evermani_scaffold_999   110551  110643  .   .   -   Gmove   CDS .   Parent=Peve_00045350
    Porites_evermani_scaffold_999   117114  117259  .   .   +   Gmove   CDS .   Parent=Peve_00045351
    Porites_evermani_scaffold_999   117363  117833  .   .   +   Gmove   CDS .   Parent=Peve_00045351
    Porites_evermani_scaffold_999   118857  119024  .   .   -   Gmove   CDS .   Parent=Peve_00045352
    Porites_evermani_scaffold_999   119372  119940  .   .   -   Gmove   CDS .   Parent=Peve_00045352
    Porites_evermani_scaffold_999   120490  120865  .   .   -   Gmove   CDS .   Parent=Peve_00045353
    Porites_evermani_scaffold_999   124812  124870  .   .   -   Gmove   CDS .   Parent=Peve_00045355
    Porites_evermani_scaffold_999   124921  125022  .   .   -   Gmove   CDS .   Parent=Peve_00045355
    Porites_evermani_scaffold_999   157406  158547  .   .   -   Gmove   CDS .   Parent=Peve_00045356
    Porites_evermani_scaffold_999   159689  159727  .   .   -   Gmove   CDS .   Parent=Peve_00045356
    Porites_evermani_scaffold_999   25255   25390   .   .   +   Gmove   CDS .   Parent=Peve_00045357
    Porites_evermani_scaffold_999   26056   26534   .   .   +   Gmove   CDS .   Parent=Peve_00045357
    Porites_evermani_scaffold_999   26987   27076   .   .   +   Gmove   CDS .   Parent=Peve_00045357
    Porites_evermani_scaffold_999   43038   43296   .   .   -   Gmove   CDS .   Parent=Peve_00045358
    Porites_evermani_scaffold_999   43486   43561   .   .   -   Gmove   CDS .   Parent=Peve_00045358
    Porites_evermani_scaffold_999   43917   44070   .   .   -   Gmove   CDS .   Parent=Peve_00045358

It looks like each grouped/concatenated FASTA in our output contains the
correct number of original gff sequences from the correct parent and
coordinates, and the output also contains sequences for all of the
parents listed in the original gff!

# 6 Align to reference transcriptome (Kallisto pseudoalignment)

## 6.1 Building Index

``` bash
# Load bash variables into memory
source .bashvars

cd "${kallisto_output_dir}"

${programs_array[kallisto]} index \
--threads=${threads} \
--index="${kallisto_index_name}" \
"${transcriptome_fasta}"

ls -lh ${kallisto_output_dir}
```

    [build] loading fasta file /home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve/data/Porites_evermanni_CDS.fasta
    [build] k-mer length: 31
    [build] warning: clipped off poly-A tail (longer than 10)
            from 2 target sequences
    [build] warning: replaced 43264 non-ACGUT characters in the input sequence
            with pseudorandom nucleotides
    KmerStream::KmerStream(): Start computing k-mer cardinality estimations (1/2)
    KmerStream::KmerStream(): Start computing k-mer cardinality estimations (1/2)
    KmerStream::KmerStream(): Finished
    CompactedDBG::build(): Estimated number of k-mers occurring at least once: 46039676
    CompactedDBG::build(): Estimated number of minimizer occurring at least once: 11275045
    CompactedDBG::filter(): Processed 52809246 k-mers in 40389 reads
    CompactedDBG::filter(): Found 46115211 unique k-mers
    CompactedDBG::filter(): Number of blocks in Bloom filter is 314726
    CompactedDBG::construct(): Extract approximate unitigs (1/2)
    CompactedDBG::construct(): Extract approximate unitigs (2/2)
    CompactedDBG::construct(): Closed all input files

    CompactedDBG::construct(): Splitting unitigs (1/2)

    CompactedDBG::construct(): Splitting unitigs (2/2)
    CompactedDBG::construct(): Before split: 489376 unitigs
    CompactedDBG::construct(): After split (1/1): 489376 unitigs
    CompactedDBG::construct(): Unitigs split: 1098
    CompactedDBG::construct(): Unitigs deleted: 0

    CompactedDBG::construct(): Joining unitigs
    CompactedDBG::construct(): After join: 467413 unitigs
    CompactedDBG::construct(): Joined 22286 unitigs
    [build] building MPHF
    [build] creating equivalence classes ... 
    [build] target de Bruijn graph has k-mer length 31 and minimizer length 23
    [build] target de Bruijn graph has 467413 contigs and contains 46155696 k-mers 

    total 83M
    -rw-r--r-- 1 shedurkin labmembers 1.5M Feb  8 09:16 kallisto.isoform.counts.matrix
    -rw-r--r-- 1 shedurkin labmembers    0 Feb  8 09:16 kallisto.isoform.TMM.EXPR.matrix
    -rw-r--r-- 1 shedurkin labmembers 2.2M Feb  8 09:16 kallisto.isoform.TPM.not_cross_norm
    -rw-r--r-- 1 shedurkin labmembers  532 Feb  8 09:16 kallisto.isoform.TPM.not_cross_norm.runTMM.R
    drwxr-xr-x 2 shedurkin labmembers 4.0K Feb  8 08:53 kallisto_quant_sample71
    drwxr-xr-x 2 shedurkin labmembers 4.0K Feb  8 08:57 kallisto_quant_sample73
    drwxr-xr-x 2 shedurkin labmembers 4.0K Feb  8 09:00 kallisto_quant_sample76
    drwxr-xr-x 2 shedurkin labmembers 4.0K Feb  8 09:04 kallisto_quant_sample79
    drwxr-xr-x 2 shedurkin labmembers 4.0K Feb  8 09:08 kallisto_quant_sample82
    -rw-r--r-- 1 shedurkin labmembers  80M Feb  8 11:12 Peve_kallisto_index.idx

Note that, when building an index, kallisto warns us that it “replaced
43264 non-ACGUT characters in the input sequence with pseudorandom
nucleotides.” This high number of identified “non-ACGUT” characters is
related to the type of reference sequences we used to build the index.
We obtained a coding sequence (CDS) gff file for P.evermanni and the
associated scaffold genome fasta, filtered the gff to retain only coding
sequences, and then used bedtools to extract the fasta sequences of
every CDS from the scaffold fasta. Notably, scaffolds are basically
fragments of known DNA sequences “stitched” together by stretches of Ns
to approximate the full sequence structure without complete sequence
data. This means some of our mRNA sequences contain long,
relatively-meaningless stretches of Ns. I’m not sure how/to what extent
this will interfere with the kallisto pseudoallignment process, since it
differs from standard alignment of a full read to reference, but we’ll
continue for now

## 6.2 Sample Quantification

Kallisto can run quantification on either single- or paired-end reads.
The default option is paired-end, which requires the input of an even
number of paired fastq files (e.g., pairA_R1.fastq, pairA_R2.fastq). To
use single-end mode, include the –single flag, as well as -l
(–fragment-length=DOUBLE, estimated avg. fragment length) and -s
(–sd=DOUBLE, estimates stand. dev. of fragment length), and a number of
fastq files. Again, gzipped files are acceptable.

Kallisto quant is rather finicky about how you input sets of paired
reads, and you can only input a single pair at a time. To circumvent,
I’ll create symlinks to each of the input files with simplified names,
create a quantification function, and apply it iteratively to each pair
using a loop.

``` bash
# Load bash variables into memory
source .bashvars

# Create sym links to each of the trimmed read files with simplified names

for file in "${trimmed_reads_dir}"/*.fastp-trim.20230519.fastq.gz; do
    # Extract sample ID and read number from the file name
    sample_id=$(echo "$file" | grep -oP 'RNA-POR-\K\d+')
    read_number=$(echo "$file" | grep -oP '_R\K\d+')

    # Create the shortened name
    shortened_name="sample${sample_id}_R${read_number}.fastq.gz"

    # Create symbolic link
    ln -s "$file" "${trimmed_reads_dir}/${shortened_name}"

done

ls -lh ${trimmed_reads_dir}
```

``` bash
# Load bash variables into memory
source .bashvars

# Function to run kallisto quant. Takes two (paired) reads as input, outputs to sample-associated directory
run_kallisto_quant() {
    source .bashvars  # Source .bashvars inside the function to make its variables accessible
    local R1_fastq=${1}
    local R2_fastq=${2}
    
    cd ${kallisto_output_dir}
    sample_num=$(basename "${R1_fastq}" "_R1.fastq.gz")
    mkdir kallisto_quant_${sample_num}

    ${programs_array[kallisto]} quant \
        --threads=${threads} \
        --index="${kallisto_output_dir}/${kallisto_index_name}" \
        --output-dir="${kallisto_output_dir}/kallisto_quant_${sample_num}" \
        --bootstrap-samples=100 \
        ${trimmed_reads_dir}/${R1_fastq} ${trimmed_reads_dir}/${R2_fastq}
}



# Iteratively apply run_kallisto_quant on each pair of input reads
for file_r1 in "${trimmed_reads_dir}"/*_R1.fastq.gz; do
    # Extract the sample name from the file name
    sample_name=$(basename "${file_r1}" "_R1.fastq.gz")

    # Form the file names (function takes input file names, not paths)
    file_r1_name="${sample_name}_R1.fastq.gz"
    file_r2_name="${sample_name}_R2.fastq.gz"

    # Check that the sample hasn't already been quantified
    if [ ! -d "${kallisto_output_dir}/kallisto_quant_${sample_name}" ]; then
    
        # Check if the corresponding R2 file exists
        if [ -e "${trimmed_reads_dir}/${file_r2}" ]; then
            # Run kallisto quant on the file pair
            run_kallisto_quant "${file_r1_name}" "${file_r2_name}" 

            echo "Processed sample: ${sample_name}"
        fi
    else
        echo "Sample already processed: ${sample_name}"
    fi
done

ls -lh ${kallisto_output_dir}
```

## 6.3 Trinity Matrix with Kallisto Output

``` bash
# Load bash variables into memory
source .bashvars

cd ${kallisto_output_dir}

${programs_array[trinity_abund_to_matrix]} \
--est_method 'kallisto' \
--gene_trans_map 'none' \
--out_prefix 'kallisto' \
--name_sample_by_basedir ${kallisto_output_dir}/kallisto_quant_*/abundance.tsv

ls -lh ${kallisto_output_dir}
```

    -reading file: /home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve/output/12-Peve-RNAseq-kallisto/kallisto/kallisto_quant_sample71/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve/output/12-Peve-RNAseq-kallisto/kallisto/kallisto_quant_sample73/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve/output/12-Peve-RNAseq-kallisto/kallisto/kallisto_quant_sample76/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve/output/12-Peve-RNAseq-kallisto/kallisto/kallisto_quant_sample79/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/deep-dive/E-Peve/output/12-Peve-RNAseq-kallisto/kallisto/kallisto_quant_sample82/abundance.tsv


    * Outputting combined matrix.

    /home/shared/trinityrnaseq-v2.12.0/util/support_scripts/run_TMM_scale_matrix.pl --matrix kallisto.isoform.TPM.not_cross_norm > kallisto.isoform.TMM.EXPR.matrixCMD: R --no-save --no-restore --no-site-file --no-init-file -q < kallisto.isoform.TPM.not_cross_norm.runTMM.R 1>&2 
    sh: 1: R: not found
    Error, cmd: R --no-save --no-restore --no-site-file --no-init-file -q < kallisto.isoform.TPM.not_cross_norm.runTMM.R 1>&2  died with ret (32512)  at /home/shared/trinityrnaseq-v2.12.0/util/support_scripts/run_TMM_scale_matrix.pl line 105.
    Error, CMD: /home/shared/trinityrnaseq-v2.12.0/util/support_scripts/run_TMM_scale_matrix.pl --matrix kallisto.isoform.TPM.not_cross_norm > kallisto.isoform.TMM.EXPR.matrix died with ret 6400 at /home/shared/trinityrnaseq-v2.12.0/util/abundance_estimates_to_matrix.pl line 385.
    total 83M
    -rw-r--r-- 1 shedurkin labmembers 1.5M Feb  8 11:12 kallisto.isoform.counts.matrix
    -rw-r--r-- 1 shedurkin labmembers    0 Feb  8 11:12 kallisto.isoform.TMM.EXPR.matrix
    -rw-r--r-- 1 shedurkin labmembers 2.2M Feb  8 11:12 kallisto.isoform.TPM.not_cross_norm
    -rw-r--r-- 1 shedurkin labmembers  532 Feb  8 11:12 kallisto.isoform.TPM.not_cross_norm.runTMM.R
    drwxr-xr-x 2 shedurkin labmembers 4.0K Feb  8 08:53 kallisto_quant_sample71
    drwxr-xr-x 2 shedurkin labmembers 4.0K Feb  8 08:57 kallisto_quant_sample73
    drwxr-xr-x 2 shedurkin labmembers 4.0K Feb  8 09:00 kallisto_quant_sample76
    drwxr-xr-x 2 shedurkin labmembers 4.0K Feb  8 09:04 kallisto_quant_sample79
    drwxr-xr-x 2 shedurkin labmembers 4.0K Feb  8 09:08 kallisto_quant_sample82
    -rw-r--r-- 1 shedurkin labmembers  80M Feb  8 11:12 Peve_kallisto_index.idx

# 7 Summary
