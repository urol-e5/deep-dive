06-Peve-sRNAseq-trimming
================
Sam White
2023-11-15

- <a href="#1-create-a-bash-variables-file"
  id="toc-1-create-a-bash-variables-file">1 Create a Bash variables
  file</a>
- <a href="#2-download-raw-srnaseq-reads"
  id="toc-2-download-raw-srnaseq-reads">2 Download raw sRNAseq reads</a>
  - <a href="#21-verify-raw-read-checksums"
    id="toc-21-verify-raw-read-checksums">2.1 Verify raw read checksums</a>
- <a href="#3-create-adapters-fasta-for-use-with-flexbar-trimming"
  id="toc-3-create-adapters-fasta-for-use-with-flexbar-trimming">3 Create
  adapters FastA for use with flexbar trimming</a>
- <a href="#4-fastqcmultiqc-on-raw-reads"
  id="toc-4-fastqcmultiqc-on-raw-reads">4 FastQC/MultiQC on raw reads</a>
- <a href="#5-trimming-with-flexbar" id="toc-5-trimming-with-flexbar">5
  Trimming with flexbar</a>
- <a href="#6-fastqcmultiqc-on-trimmed-reads"
  id="toc-6-fastqcmultiqc-on-trimmed-reads">6 FastQC/MultiQC on trimmed
  reads</a>
- <a href="#7-summary" id="toc-7-summary">7 Summary</a>
- <a href="#8-citations" id="toc-8-citations">8 Citations</a>

FastQC/MultiQC ([Ewels et al. 2016](#ref-ewels2016); [Andrews,
n.d.](#ref-Andrews_undated-nf)) assessment of raw and
[flexbar](https://github.com/seqan/flexbar)-trimmed ([Dodt et al.
2012](#ref-Dodt2012-rt); [Roehr, Dieterich, and Reinert
2017](#ref-Roehr2017-dr)) sequences of E5 *P.evermanni* sRNAseq data
from
[20230515](https://robertslab.github.io/sams-notebook/posts/2023/2023-05-17-Data-Management---E5-Coral-RNA-seq-and-sRNA-seq-Reorganizing-and-Renaming/).

Inputs:

- sRNAseq gzipped FastQs (e.g. `*.fastq.gz`)

Outputs:

- [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  HTML reports for raw and trimmed reads.

- [`MultiQC`](https://multiqc.info/) HTML summaries of
  [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  for raw and trimmed reads.

- Trimmed reads with final length of 25bp: `*flexbar_trim.25bp.fastq.gz`

Libraries were prepared and sequenced by Azenta:

- Library prep: [NEB nebnext-small-rna-library-prep-set-for-illumina
  kit](https://www.neb.com/en-us/-/media/nebus/files/manuals/manuale7300_e7330_e7560_e7580.pdf?rev=d0964a2e637843b1afcb9f7d666d07b2&hash=7AC0B0EB012708EFAB0E4DBEEAF1446A)
  (PDF)

- Sequencing: Illumina HiSeq 4000, 150bp PE

------------------------------------------------------------------------

# 1 Create a Bash variables file

This allows usage of Bash variables across R Markdown chunks.

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "# Data directories"
echo 'export deep_dive_dir=/home/shared/8TB_HDD_01/sam/gitrepos/deep-dive'
echo 'export output_dir_top=${deep_dive_dir}/E-Peve/output/06-Peve-sRNAseq-trimming'
echo 'export raw_fastqc_dir=${deep_dive_dir}/E-Peve/output/06-Peve-sRNAseq-trimming/raw-fastqc'
echo 'export raw_reads_dir=${deep_dive_dir}/E-Peve/data/06-Peve-sRNAseq-trimming/raw-reads'
echo 'export raw_reads_url="https://owl.fish.washington.edu/nightingales/P_evermanni/30-852430235/"'
echo 'export trimmed_fastqc_dir=${deep_dive_dir}/E-Peve/output/06-Peve-sRNAseq-trimming/trimmed-fastqc'
echo 'export trimmed_reads_dir=${deep_dive_dir}/E-Peve/output/06-Peve-sRNAseq-trimming/trimmed-reads'
echo ""

echo "# Paths to programs"
echo 'export fastqc=/home/shared/FastQC-0.12.1/fastqc'
echo 'export flexbar=/home/shared/flexbar-3.5.0-linux/flexbar'
echo 'export multiqc=/home/sam/programs/mambaforge/bin/multiqc'
echo ""



echo "# Set FastQ filename patterns"
echo "export fastq_pattern='*.fastq.gz'"
echo "export R1_fastq_pattern='*_R1_*.fastq.gz'"
echo "export R2_fastq_pattern='*_R2_*.fastq.gz'"
echo ""

echo "# Set number of CPUs to use"
echo 'export threads=40'
echo ""

echo "# Set maximum read length"
echo 'export max_read_length=25'
echo ""

echo "# Input/output files"
echo 'export fastq_checksums=input_fastq_checksums.md5'
echo 'export trimmed_checksums=trimmed_fastq_checksums.md5'
echo 'export NEB_adapters_fasta=NEB-adapters.fasta'
echo ""

echo "## NEB nebnext-small-rna-library-prep-set-for-illumina adapters"
echo 'export first_adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"'
echo 'export second_adapter="GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT"'
echo ""

echo "## Inititalize arrays"
echo 'export fastq_array_R1=()'
echo 'export fastq_array_R2=()'
echo 'export raw_fastqs_array=()'
echo 'export R1_names_array=()'
echo 'export R2_names_array=()'
echo 'export trimmed_fastqs_array=()'
echo ""

echo "# Programs associative array"
echo "declare -A programs_array"
echo "programs_array=("
echo '[fastqc]="${fastqc}" \'
echo '[multiqc]="${multiqc}" \'
echo '[flexbar]="${flexbar}"'
echo ")"
} > .bashvars

cat .bashvars
```

    #### Assign Variables ####

    # Data directories
    export deep_dive_dir=/home/shared/8TB_HDD_01/sam/gitrepos/deep-dive
    export output_dir_top=${deep_dive_dir}/E-Peve/output/06-Peve-sRNAseq-trimming
    export raw_fastqc_dir=${deep_dive_dir}/E-Peve/output/06-Peve-sRNAseq-trimming/raw-fastqc
    export raw_reads_dir=${deep_dive_dir}/E-Peve/data/06-Peve-sRNAseq-trimming/raw-reads
    export raw_reads_url="https://owl.fish.washington.edu/nightingales/P_evermanni/30-852430235/"
    export trimmed_fastqc_dir=${deep_dive_dir}/E-Peve/output/06-Peve-sRNAseq-trimming/trimmed-fastqc
    export trimmed_reads_dir=${deep_dive_dir}/E-Peve/output/06-Peve-sRNAseq-trimming/trimmed-reads

    # Paths to programs
    export fastqc=/home/shared/FastQC-0.12.1/fastqc
    export flexbar=/home/shared/flexbar-3.5.0-linux/flexbar
    export multiqc=/home/sam/programs/mambaforge/bin/multiqc

    # Set FastQ filename patterns
    export fastq_pattern='*.fastq.gz'
    export R1_fastq_pattern='*_R1_*.fastq.gz'
    export R2_fastq_pattern='*_R2_*.fastq.gz'

    # Set number of CPUs to use
    export threads=40

    # Set maximum read length
    export max_read_length=25

    # Input/output files
    export fastq_checksums=input_fastq_checksums.md5
    export trimmed_checksums=trimmed_fastq_checksums.md5
    export NEB_adapters_fasta=NEB-adapters.fasta

    ## NEB nebnext-small-rna-library-prep-set-for-illumina adapters
    export first_adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    export second_adapter="GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT"

    ## Inititalize arrays
    export fastq_array_R1=()
    export fastq_array_R2=()
    export raw_fastqs_array=()
    export R1_names_array=()
    export R2_names_array=()
    export trimmed_fastqs_array=()

    # Programs associative array
    declare -A programs_array
    programs_array=(
    [fastqc]="${fastqc}" \
    [multiqc]="${multiqc}" \
    [flexbar]="${flexbar}"
    )

# 2 Download raw sRNAseq reads

Reads are downloaded from
<https://owl.fish.washington.edu/nightingales/P_evermanni/30-852430235/>

The `--cut-dirs 3` command cuts the preceding directory structure
(i.e. `nightingales/P_evermanni/30-852430235/`) so that we just end up
with the reads.

``` bash
# Load bash variables into memory
source .bashvars

wget \
--directory-prefix ${raw_reads_dir} \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 3 \
--no-host-directories \
--no-parent \
--quiet \
--accept "sRNA*,checksums.md5" ${raw_reads_url}

ls -lh "${raw_reads_dir}"
```

    total 5.6G
    -rw-r--r-- 1 sam sam  798 May 17 11:23 checksums.md5
    -rw-r--r-- 1 sam sam 953M May 17 10:46 sRNA-POR-73-S1-TP2_R1_001.fastq.gz
    -rw-r--r-- 1 sam sam 1.1G May 17 10:47 sRNA-POR-73-S1-TP2_R2_001.fastq.gz
    -rw-r--r-- 1 sam sam 899M May 17 10:48 sRNA-POR-79-S1-TP2_R1_001.fastq.gz
    -rw-r--r-- 1 sam sam 916M May 17 10:49 sRNA-POR-79-S1-TP2_R2_001.fastq.gz
    -rw-r--r-- 1 sam sam 934M May 17 10:51 sRNA-POR-82-S1-TP2_R1_001.fastq.gz
    -rw-r--r-- 1 sam sam 959M May 17 10:52 sRNA-POR-82-S1-TP2_R2_001.fastq.gz

## 2.1 Verify raw read checksums

``` bash
# Load bash variables into memory
source .bashvars

cd "${raw_reads_dir}"

# Checksums file contains other files, so this just looks for the sRNAseq files.
grep "sRNA" checksums.md5 | md5sum --check
```

    sRNA-POR-73-S1-TP2_R1_001.fastq.gz: OK
    sRNA-POR-73-S1-TP2_R2_001.fastq.gz: OK
    sRNA-POR-79-S1-TP2_R1_001.fastq.gz: OK
    sRNA-POR-79-S1-TP2_R2_001.fastq.gz: OK
    sRNA-POR-82-S1-TP2_R1_001.fastq.gz: OK
    sRNA-POR-82-S1-TP2_R2_001.fastq.gz: OK

# 3 Create adapters FastA for use with [flexbar](https://github.com/seqan/flexbar) trimming

``` bash
# Load bash variables into memory
source .bashvars

echo "Creating adapters FastA."
echo ""
adapter_count=0

# Check for adapters file first
# Then create adapters file if doesn't exist
if [ -f "${output_dir_top}/${NEB_adapters_fasta}" ]; then
  echo "${output_dir_top}/${NEB_adapters_fasta} already exists. Nothing to do."
else
  for adapter in "${first_adapter}" "${second_adapter}"
  do
    adapter_count=$((adapter_count + 1))
    printf ">%s\n%s\n" "adapter_${adapter_count}" "${adapter}"
  done >> "${output_dir_top}/${NEB_adapters_fasta}"
fi

echo ""
echo "Adapters FastA:"
echo ""
cat "${output_dir_top}/${NEB_adapters_fasta}"
echo ""
```

    Creating adapters FastA.

    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/E-Peve/output/06-Peve-sRNAseq-trimming/NEB-adapters.fasta already exists. Nothing to do.

    Adapters FastA:

    >adapter_1
    AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
    >adapter_2
    GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT

# 4 FastQC/MultiQC on raw reads

``` bash
# Load bash variables into memory
source .bashvars

############ RUN FASTQC ############


# Create array of trimmed FastQs
raw_fastqs_array=(${raw_reads_dir}/${fastq_pattern})

# Pass array contents to new variable as space-delimited list
raw_fastqc_list=$(echo "${raw_fastqs_array[*]}")

echo "Beginning FastQC on raw reads..."
echo ""

# Run FastQC
### NOTE: Do NOT quote raw_fastqc_list
${programs_array[fastqc]} \
--threads ${threads} \
--outdir ${raw_fastqc_dir} \
--quiet \
${raw_fastqc_list}

echo "FastQC on raw reads complete!"
echo ""

############ END FASTQC ############

############ RUN MULTIQC ############
echo "Beginning MultiQC on raw FastQC..."
echo ""

${programs_array[multiqc]} ${raw_fastqc_dir} -o ${raw_fastqc_dir}

echo ""
echo "MultiQC on raw FastQs complete."
echo ""

############ END MULTIQC ############

echo "Removing FastQC zip files."
echo ""
rm ${raw_fastqc_dir}/*.zip
echo "FastQC zip files removed."
echo ""

# View directory contents
ls -lh ${raw_fastqc_dir}
```

# 5 Trimming with [flexbar](https://github.com/seqan/flexbar)

``` bash
# Load bash variables into memory
source .bashvars

# Change to directory with raw reads
cd "${raw_reads_dir}"

# Create arrays of FastQ R1 files and sample names
# Do NOT quote R1_fastq_pattern variable
for fastq in ${R1_fastq_pattern}
do
  fastq_array_R1+=("${fastq}")

  # Use parameter substitution to remove all text up to and including last "." from
  # right side of string.
  R1_names_array+=("${fastq%%.*}")
done

# Create array of FastQ R2 files
# Do NOT quote R2_fastq_pattern variable
for fastq in ${R2_fastq_pattern}
do
  fastq_array_R2+=("${fastq}")

  # Use parameter substitution to remove all text up to and including last "." from
  # right side of string.
  R2_names_array+=("${fastq%%.*}")
done


############ RUN FLEXBAR ############
# Uses parameter substitution (e.g. ${R1_sample_name%%_*})to rm the _R[12]
# Uses NEB adapter file
# --adapter-pair-overlap ON: Recommended by NEB sRNA kit
# --qtrim-threshold 25: Minimum quality
# --qtrim-format i1.8: Sets sequencer as illumina
# --post-trim-length: Trim reads from 3' end to max length
# --target: Sets file naming patterns
# --zip-output GZ: Sets type of compression. GZ = gzip

# Run flexbar on files
echo "Beginning flexbar trimming."
echo ""

time \
for index in "${!fastq_array_R1[@]}"
do
  R1_sample_name="${R1_names_array[index]}"
  R2_sample_name="${R2_names_array[index]}"

  # Begin flexbar trimming
  ${programs_array[flexbar]} \
  --reads ${fastq_array_R1[index]} \
  --reads2 ${fastq_array_R2[index]}  \
  --adapters ${output_dir_top}/${NEB_adapters_fasta} \
  --adapter-pair-overlap ON \
  --qtrim-format i1.8 \
  --qtrim-threshold 25 \
  --post-trim-length ${max_read_length} \
  --threads ${threads} \
  --target "${trimmed_reads_dir}/${R1_sample_name%%_*}.flexbar_trim.${max_read_length}bp" \
  --zip-output GZ
        
    # Move to trimmed directory
    # This is done so checksums file doesn't include excess path
    cd ${trimmed_reads_dir}

    # Generate md5 checksums for newly trimmed files
    {
      md5sum "${R1_sample_name%%_*}.flexbar_trim.${max_read_length}bp_1.fastq.gz"
      md5sum "${R2_sample_name%%_*}.flexbar_trim.${max_read_length}bp_2.fastq.gz"
    } >> "${trimmed_checksums}"
    
    # Change back to to raw reads directory
    cd "${raw_reads_dir}"

done

echo ""
echo "flexbar trimming complete."
echo ""

echo "Trimmed FastQs MD5 checksums:"
echo ""

cat "${trimmed_reads_dir}/${trimmed_checksums}"

############ END FLEXBAR ############
```

# 6 FastQC/MultiQC on trimmed reads

``` bash
# Load bash variables into memory
source .bashvars

############ RUN FASTQC ############

### NOTE: Do NOT quote raw_fastqc_list
# Create array of trimmed FastQs
trimmed_fastqs_array=(${trimmed_reads_dir}/${fastq_pattern})

# Pass array contents to new variable as space-delimited list
trimmed_fastqc_list=$(echo "${trimmed_fastqs_array[*]}")

echo "Beginning FastQC on raw reads..."
echo ""

# Run FastQC
${programs_array[fastqc]} \
--threads ${threads} \
--outdir ${trimmed_fastqc_dir} \
--quiet \
${trimmed_fastqc_list}

echo "FastQC on trimmed reads complete!"
echo ""

############ END FASTQC ############

############ RUN MULTIQC ############
echo "Beginning MultiQC on raw FastQC..."
echo ""

${programs_array[multiqc]} ${trimmed_fastqc_dir} -o ${trimmed_fastqc_dir}

echo ""
echo "MultiQC on trimmed FastQs complete."
echo ""

############ END MULTIQC ############

echo "Removing FastQC zip files."
echo ""
rm ${trimmed_fastqc_dir}/*.zip
echo "FastQC zip files removed."
echo ""

# View directory contents
ls -lh ${trimmed_fastqc_dir}
```

# 7 Summary

A quick comparison of raw and trimmed reads to show trimming worked:

- quality is improved
- length is 25bp
- adapters removed

| RAW                                                                                                                                                                                                | TRIMMED                                                                                                                                                                                                    |
|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ![Raw MultiQC per base sequence quality plot](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/output/06-Peve-sRNAseq-trimming/raw-fastqc/fastqc_per_base_sequence_quality_plot.png?raw=true) | ![Trimmed MultiQC per base sequence quality plot](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/output/06-Peve-sRNAseq-trimming/trimmed-fastqc/fastqc_per_base_sequence_quality_plot.png?raw=true) |
| ![Raw MultiQC adapter content plot](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/output/06-Peve-sRNAseq-trimming/raw-fastqc/fastqc_adapter_content_plot.png?raw=true)                     | ![Trimmed MultiQC adapter content plot](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/output/06-Peve-sRNAseq-trimming/trimmed-fastqc/fastqc_adapter_content_plot.png?raw=true)                     |

# 8 Citations

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Andrews_undated-nf" class="csl-entry">

Andrews, Simon. n.d. “FastQC.” Github.

</div>

<div id="ref-Dodt2012-rt" class="csl-entry">

Dodt, Matthias, Johannes T Roehr, Rina Ahmed, and Christoph Dieterich.
2012. “FLEXBAR-Flexible Barcode and Adapter Processing for
Next-Generation Sequencing Platforms.” *Biology* 1 (3): 895–905.

</div>

<div id="ref-ewels2016" class="csl-entry">

Ewels, Philip, Måns Magnusson, Sverker Lundin, and Max Käller. 2016.
“MultiQC: Summarize Analysis Results for Multiple Tools and Samples in a
Single Report.” *Bioinformatics* 32 (19): 3047–48.
<https://doi.org/10.1093/bioinformatics/btw354>.

</div>

<div id="ref-Roehr2017-dr" class="csl-entry">

Roehr, Johannes T, Christoph Dieterich, and Knut Reinert. 2017. “Flexbar
3.0 - SIMD and Multicore Parallelization.” *Bioinformatics* 33 (18):
2941–42.

</div>

</div>
