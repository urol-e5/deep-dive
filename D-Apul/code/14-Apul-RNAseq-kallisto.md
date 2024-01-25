14-Apul-RNAseq-kallisto
================
Kathleen Durkin
2024-01-24

Description

Inputs:

trimmed RNAseq reads transcripts FASTA (note using transcripts from
A.millepora, not A.pulchra)

Outputs:

------------------------------------------------------------------------

# Create a Bash variables file

This allows usage of Bash variables (e.g. paths to common directories)
across R Markdown chunks.

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "# Data directories"
echo 'export Apul_dir=/home/shared/8TB_HDD_02/shedurkin/deep-dive/D-Apul'
echo 'export output_dir_top=${Apul_dir}/output/14-Apul-RNAseq-kallisto'
echo 'export trimmed_reads_dir=${output_dir_top}/trimmed-reads'
echo 'export kallisto_output_dir=${output_dir_top}/kallisto'
echo ""

echo "# Input/Output files"
echo 'export Amil_ncbi_downloads_dir=${Apul_dir}/data/Amil'
echo 'export transcriptome_fasta_dir=${Amil_ncbi_downloads_dir}/ncbi_dataset/data/GCF_013753865.1'
echo 'export transcriptome_fasta_name="rna.fna"'
echo 'export transcriptome_fasta="${transcriptome_fasta_dir}/${transcriptome_fasta_name}"'
echo 'export kallisto_index_name="Amil_kallisto_index.idx"'

echo "# External data URLs"
echo 'export trimmed_reads_url="https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/A_pulchra/trimmed/"'
echo 'export transcriptome_fasta_url="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_013753865.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_013753865.1.zip"'
echo ""

echo "# Set filename patterns"
echo "export ncbi_accession='GCF_013753865.1'"
echo "export fastq_pattern='*.fastq.gz'"
echo "export R1_fastq_pattern='*_R1.fastq.gz'"
echo "export R2_fastq_pattern='*_R2.fastq.gz'"
echo ""

echo "# Paths to programs"
echo 'export kallisto=/home/shared/kallisto_linux-v0.50.1/kallisto'
echo 'export trinity_abund_to_matrix=/home/shared/trinityrnaseq-v2.12.0/util/abundance_estimates_to_matrix.pl'
echo ""

echo "# Set number of CPUs to use"
echo 'export threads=20'
echo ""

echo "# Programs associative array"
echo "declare -A programs_array"
echo "programs_array=("
echo '[kallisto]="${kallisto}" \'
echo '[trinity_abund_to_matrix]="${trinity_abund_to_matrix}" \'
echo ")"

} > .bashvars

cat .bashvars
```

    #### Assign Variables ####

    # Data directories
    export Apul_dir=/home/shared/8TB_HDD_02/shedurkin/deep-dive/D-Apul
    export output_dir_top=${Apul_dir}/output/14-Apul-RNAseq-kallisto
    export trimmed_reads_dir=${output_dir_top}/trimmed-reads
    export kallisto_output_dir=${output_dir_top}/kallisto

    # Input/Output files
    export Amil_ncbi_downloads_dir=${Apul_dir}/data/Amil
    export transcriptome_fasta_dir=${Amil_ncbi_downloads_dir}/ncbi_dataset/data/GCF_013753865.1
    export transcriptome_fasta_name="rna.fna"
    export transcriptome_fasta="${transcriptome_fasta_dir}/${transcriptome_fasta_name}"
    export kallisto_index_name="Amil_kallisto_index.idx"
    # External data URLs
    export trimmed_reads_url="https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/A_pulchra/trimmed/"
    export transcriptome_fasta_url="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_013753865.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_013753865.1.zip"

    # Set filename patterns
    export ncbi_accession='GCF_013753865.1'
    export fastq_pattern='*.fastq.gz'
    export R1_fastq_pattern='*_R1.fastq.gz'
    export R2_fastq_pattern='*_R2.fastq.gz'

    # Paths to programs
    export kallisto=/home/shared/kallisto_linux-v0.50.1/kallisto
    export trinity_abund_to_matrix=/home/shared/trinityrnaseq-v2.12.0/util/abundance_estimates_to_matrix.pl

    # Set number of CPUs to use
    export threads=20

    # Programs associative array
    declare -A programs_array
    programs_array=(
    [kallisto]="${kallisto}" \
    [trinity_abund_to_matrix]="${trinity_abund_to_matrix}" \
    )

# Download trimmed RNAseq reads

Reads are downloaded from:
<https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/A_pulchra/trimmed/>

The `--cut-dirs 4` command cuts the preceding directory structure
(i.e. `Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/A_pulchra/trimmed/`)
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

ls -lh "${trimmed_reads_dir}"
```

# FastQC/MultiQC on trimmed reads

Already performed, can view multiqc report at
<https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/A_pulchra/trimmed/multiqc_report.html>

# Retrieve the reference transcriptome

``` bash
# Load bash variables into memory
source .bashvars

curl \
-JX GET ${transcriptome_fasta_url} \
-H "Accept: application/zip" \
-o ${Amil_ncbi_downloads_dir}/${ncbi_accession}.zip

unzip ${Amil_ncbi_downloads_dir}/${ncbi_accession}.zip -d ${Amil_ncbi_downloads_dir}

rm ${Amil_ncbi_downloads_dir}/${ncbi_accession}.zip
```

## Verify transcriptome FastA MD5 checksum

No checksum file(s) provided with ncbi download, so skippping this step

# Align to reference transcriptome (Kallisto pseudoalignment)

## Building Index

``` bash
# Load bash variables into memory
source .bashvars

cd "${kallisto_output_dir}"

${programs_array[kallisto]} index \
--threads=${threads} \
--index="${kallisto_index_name}" \
"${transcriptome_fasta}"
```

## Sample Quantification

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
    sample_id=$(echo "$file" | grep -oP 'RNA-ACR-\K\d+')
    read_number=$(echo "$file" | grep -oP '_R\K\d+')

    # Create the shortened name
    shortened_name="sample${sample_id}_R${read_number}.fastq.gz"

    # Create symbolic link
    ln -s "$file" "${trimmed_reads_dir}/${shortened_name}"

    echo "Created link: ${shortened_name}"
done
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
```

## Trinity Matrix with Kallisto Output

``` bash
# Load bash variables into memory
source .bashvars

cd ${kallisto_output_dir}

${programs_array[trinity_abund_to_matrix]} \
--est_method 'kallisto' \
--gene_trans_map 'none' \
--out_prefix 'kallisto' \
--name_sample_by_basedir ${kallisto_output_dir}/kallisto_quant_*/abundance.tsv
```

# Summary
