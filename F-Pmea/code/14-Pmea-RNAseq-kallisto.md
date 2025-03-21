14-Pmea-RNAseq-kallisto
================
Kathleen Durkin
2024-01-30

Description

Inputs:

Trimmed *P. meandrina RNAseq* reads (found
[here](https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_meandrina/trimmed/))

Transcripts FASTA (found
[here](http://cyanophora.rutgers.edu/Pocillopora_meandrina/))

Outputs:

Counts matrix

------------------------------------------------------------------------

# Create a Bash variables file

This allows usage of Bash variables (e.g. paths to common directories)
across R Markdown chunks.

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "# Data directories"
echo 'export Pmea_dir=/home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea'
echo 'export output_dir_top=${Pmea_dir}/output/14-Pmea-RNAseq-kallisto'
echo 'export transcriptome_fasta_dir=${Pmea_dir}/data/Pmea'
echo 'export trimmed_reads_dir=${output_dir_top}/trimmed-reads'
echo 'export kallisto_output_dir=${output_dir_top}/kallisto'
echo ""

echo "# Input/Output files"
echo 'export transcriptome_fasta_name="Pocillopora_meandrina_HIv1.genes.cds.fna.gz"'
echo 'export transcriptome_fasta="${transcriptome_fasta_dir}/${transcriptome_fasta_name}"'
echo 'export kallisto_index_name="Pmea_kallisto_index.idx"'

echo "# External data URLs"
echo 'export trimmed_reads_url="https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_meandrina/trimmed/"'
echo 'export transcriptome_fasta_url="http://cyanophora.rutgers.edu/Pocillopora_meandrina/Pocillopora_meandrina_HIv1.genes.cds.fna.gz"'
echo ""

echo "# Set filename patterns"
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
    export Pmea_dir=/home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea
    export output_dir_top=${Pmea_dir}/output/14-Pmea-RNAseq-kallisto
    export transcriptome_fasta_dir=${Pmea_dir}/data/Pmea
    export trimmed_reads_dir=${output_dir_top}/trimmed-reads
    export kallisto_output_dir=${output_dir_top}/kallisto

    # Input/Output files
    export transcriptome_fasta_name="Pocillopora_meandrina_HIv1.genes.cds.fna.gz"
    export transcriptome_fasta="${transcriptome_fasta_dir}/${transcriptome_fasta_name}"
    export kallisto_index_name="Pmea_kallisto_index.idx"
    # External data URLs
    export trimmed_reads_url="https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_meandrina/trimmed/"
    export transcriptome_fasta_url="http://cyanophora.rutgers.edu/Pocillopora_meandrina/Pocillopora_meandrina_HIv1.genes.cds.fna.gz"

    # Set filename patterns
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
(i.e. `Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_meandrina/trimmed/`)
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

    total 31G
    drwxr-xr-x 2 shedurkin labmembers 4.0K Jan 30 18:57 multiqc_data
    -rw-r--r-- 1 shedurkin labmembers 3.2G May 19  2023 RNA-POC-47-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.3G May 19  2023 RNA-POC-47-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.0G May 19  2023 RNA-POC-48-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.1G May 19  2023 RNA-POC-48-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.3G May 19  2023 RNA-POC-50-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.4G May 19  2023 RNA-POC-50-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.2G May 19  2023 RNA-POC-53-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.3G May 19  2023 RNA-POC-53-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.5G May 19  2023 RNA-POC-57-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.6G May 19  2023 RNA-POC-57-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample47_R1.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-47-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample47_R2.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-47-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample48_R1.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-48-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample48_R2.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-48-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample50_R1.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-50-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample50_R2.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-50-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample53_R1.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-53-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample53_R2.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-53-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample57_R1.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-57-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample57_R2.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-57-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz

# FastQC/MultiQC on trimmed reads

Already performed, can view multiqc report at
<https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_meandrina/trimmed/multiqc_report.html>

# Retrieve the reference transcriptome

``` bash
# Load bash variables into memory
source .bashvars

wget \
--directory-prefix ${transcriptome_fasta_dir} \
--no-check-certificate \
--continue \
--no-host-directories \
--no-directories \
--no-parent \
--quiet \
--execute robots=off \
--accept "${transcriptome_fasta_name}" ${transcriptome_fasta_url}

ls -lh "${transcriptome_fasta_dir}"
```

    total 2.0G
    -rw-r--r-- 1 shedurkin labmembers 360M Dec  1 08:51 Pocillopora_meandrina_HIv1.assembly.fa
    -rw-r--r-- 1 shedurkin labmembers 107M Nov 30 15:52 Pocillopora_meandrina_HIv1.assembly.fa.1.ebwt
    -rw-r--r-- 1 shedurkin labmembers  45M Nov 30 15:52 Pocillopora_meandrina_HIv1.assembly.fa.2.ebwt
    -rw-r--r-- 1 shedurkin labmembers 1.9K Nov 30 15:50 Pocillopora_meandrina_HIv1.assembly.fa.3.ebwt
    -rw-r--r-- 1 shedurkin labmembers  90M Nov 30 15:50 Pocillopora_meandrina_HIv1.assembly.fa.4.ebwt
    -rw-r--r-- 1 shedurkin labmembers  15K Nov 30 15:50 Pocillopora_meandrina_HIv1.assembly.fa.fai
    -rw-r--r-- 1 shedurkin labmembers 107M Nov 30 15:54 Pocillopora_meandrina_HIv1.assembly.fa.rev.1.ebwt
    -rw-r--r-- 1 shedurkin labmembers  45M Nov 30 15:54 Pocillopora_meandrina_HIv1.assembly.fa.rev.2.ebwt
    -rw-r--r-- 1 shedurkin labmembers 360M May 23  2022 Pocillopora_meandrina_HIv1.assembly.fasta
    -rw-r--r-- 1 shedurkin labmembers 101M May 23  2022 Pocillopora_meandrina_HIv1.assembly.fasta.gz
    -rw-r--r-- 1 shedurkin labmembers  13M May 23  2022 Pocillopora_meandrina_HIv1.genes.cds.fna.gz
    -rw-r--r-- 1 shedurkin labmembers 107M Nov 30 08:55 Pocillopora_meandrina_HIv1_nospaces.assembly.1.ebwt
    -rw-r--r-- 1 shedurkin labmembers  45M Nov 30 08:55 Pocillopora_meandrina_HIv1_nospaces.assembly.2.ebwt
    -rw-r--r-- 1 shedurkin labmembers 1.9K Nov 30 08:53 Pocillopora_meandrina_HIv1_nospaces.assembly.3.ebwt
    -rw-r--r-- 1 shedurkin labmembers  90M Nov 30 08:53 Pocillopora_meandrina_HIv1_nospaces.assembly.4.ebwt
    -rw-r--r-- 1 shedurkin labmembers 360M Dec  5 11:39 Pocillopora_meandrina_HIv1_nospaces.assembly.fasta
    -rw-r--r-- 1 shedurkin labmembers 107M Nov 30 08:57 Pocillopora_meandrina_HIv1_nospaces.assembly.rev.1.ebwt
    -rw-r--r-- 1 shedurkin labmembers  45M Nov 30 08:57 Pocillopora_meandrina_HIv1_nospaces.assembly.rev.2.ebwt

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

ls -lh ${kallisto_output_dir}
```


    [build] loading fasta file /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/data/Pmea/Pocillopora_meandrina_HIv1.genes.cds.fna.gz
    [build] k-mer length: 31
    KmerStream::KmerStream(): Start computing k-mer cardinality estimations (1/2)
    KmerStream::KmerStream(): Start computing k-mer cardinality estimations (1/2)
    KmerStream::KmerStream(): Finished
    CompactedDBG::build(): Estimated number of k-mers occurring at least once: 38035051
    CompactedDBG::build(): Estimated number of minimizer occurring at least once: 9356951
    CompactedDBG::filter(): Processed 42887685 k-mers in 31840 reads
    CompactedDBG::filter(): Found 37946767 unique k-mers
    CompactedDBG::filter(): Number of blocks in Bloom filter is 260007
    CompactedDBG::construct(): Extract approximate unitigs (1/2)
    CompactedDBG::construct(): Extract approximate unitigs (2/2)
    CompactedDBG::construct(): Closed all input files

    CompactedDBG::construct(): Splitting unitigs (1/2)

    CompactedDBG::construct(): Splitting unitigs (2/2)
    CompactedDBG::construct(): Before split: 268901 unitigs
    CompactedDBG::construct(): After split (1/1): 268901 unitigs
    CompactedDBG::construct(): Unitigs split: 854
    CompactedDBG::construct(): Unitigs deleted: 0

    CompactedDBG::construct(): Joining unitigs
    CompactedDBG::construct(): After join: 251521 unitigs
    CompactedDBG::construct(): Joined 17529 unitigs
    [build] building MPHF
    [build] creating equivalence classes ... 
    [build] target de Bruijn graph has k-mer length 31 and minimizer length 23
    [build] target de Bruijn graph has 251521 contigs and contains 37979144 k-mers 

    total 56M
    -rw-r--r-- 1 shedurkin labmembers 2.0M Jan 30 18:56 kallisto.isoform.counts.matrix
    -rw-r--r-- 1 shedurkin labmembers    0 Jan 30 18:56 kallisto.isoform.TMM.EXPR.matrix
    -rw-r--r-- 1 shedurkin labmembers 2.4M Jan 30 18:56 kallisto.isoform.TPM.not_cross_norm
    -rw-r--r-- 1 shedurkin labmembers  532 Jan 30 18:56 kallisto.isoform.TPM.not_cross_norm.runTMM.R
    drwxr-xr-x 2 shedurkin labmembers 4.0K Jan 30 18:38 kallisto_quant_sample47
    drwxr-xr-x 2 shedurkin labmembers 4.0K Jan 30 18:42 kallisto_quant_sample48
    drwxr-xr-x 2 shedurkin labmembers 4.0K Jan 30 18:46 kallisto_quant_sample50
    drwxr-xr-x 2 shedurkin labmembers 4.0K Jan 30 18:50 kallisto_quant_sample53
    drwxr-xr-x 2 shedurkin labmembers 4.0K Jan 30 18:53 kallisto_quant_sample57
    -rw-r--r-- 1 shedurkin labmembers  51M Jan 30 18:57 Pmea_kallisto_index.idx

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
    sample_id=$(echo "$file" | grep -oP 'RNA-POC-\K\d+')
    read_number=$(echo "$file" | grep -oP '_R\K\d+')

    # Create the shortened name
    shortened_name="sample${sample_id}_R${read_number}.fastq.gz"

    # Create symbolic link
    ln -s "$file" "${trimmed_reads_dir}/${shortened_name}"

    echo "Created link: ${shortened_name}"
done

ls -lh ${trimmed_reads_dir}
```

    ln: failed to create symbolic link '/home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/sample47_R1.fastq.gz': File exists
    Created link: sample47_R1.fastq.gz
    ln: failed to create symbolic link '/home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/sample47_R2.fastq.gz': File exists
    Created link: sample47_R2.fastq.gz
    ln: failed to create symbolic link '/home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/sample48_R1.fastq.gz': File exists
    Created link: sample48_R1.fastq.gz
    ln: failed to create symbolic link '/home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/sample48_R2.fastq.gz': File exists
    Created link: sample48_R2.fastq.gz
    ln: failed to create symbolic link '/home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/sample50_R1.fastq.gz': File exists
    Created link: sample50_R1.fastq.gz
    ln: failed to create symbolic link '/home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/sample50_R2.fastq.gz': File exists
    Created link: sample50_R2.fastq.gz
    ln: failed to create symbolic link '/home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/sample53_R1.fastq.gz': File exists
    Created link: sample53_R1.fastq.gz
    ln: failed to create symbolic link '/home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/sample53_R2.fastq.gz': File exists
    Created link: sample53_R2.fastq.gz
    ln: failed to create symbolic link '/home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/sample57_R1.fastq.gz': File exists
    Created link: sample57_R1.fastq.gz
    ln: failed to create symbolic link '/home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/sample57_R2.fastq.gz': File exists
    Created link: sample57_R2.fastq.gz
    total 31G
    drwxr-xr-x 2 shedurkin labmembers 4.0K Jan 30 18:57 multiqc_data
    -rw-r--r-- 1 shedurkin labmembers 3.2G May 19  2023 RNA-POC-47-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.3G May 19  2023 RNA-POC-47-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.0G May 19  2023 RNA-POC-48-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.1G May 19  2023 RNA-POC-48-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.3G May 19  2023 RNA-POC-50-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.4G May 19  2023 RNA-POC-50-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.2G May 19  2023 RNA-POC-53-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.3G May 19  2023 RNA-POC-53-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.5G May 19  2023 RNA-POC-57-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.6G May 19  2023 RNA-POC-57-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample47_R1.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-47-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample47_R2.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-47-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample48_R1.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-48-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample48_R2.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-48-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample50_R1.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-50-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample50_R2.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-50-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample53_R1.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-53-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample53_R2.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-53-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample57_R1.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-57-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    lrwxrwxrwx 1 shedurkin labmembers  149 Jan 30 18:33 sample57_R2.fastq.gz -> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/trimmed-reads/RNA-POC-57-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz

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
# This will take a while!(20-30 min for 5 samples)
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

ls -lh ${kallisto_output_dir}
```

    -reading file: /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/kallisto/kallisto_quant_sample47/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/kallisto/kallisto_quant_sample48/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/kallisto/kallisto_quant_sample50/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/kallisto/kallisto_quant_sample53/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/14-Pmea-RNAseq-kallisto/kallisto/kallisto_quant_sample57/abundance.tsv


    * Outputting combined matrix.

    /home/shared/trinityrnaseq-v2.12.0/util/support_scripts/run_TMM_scale_matrix.pl --matrix kallisto.isoform.TPM.not_cross_norm > kallisto.isoform.TMM.EXPR.matrixCMD: R --no-save --no-restore --no-site-file --no-init-file -q < kallisto.isoform.TPM.not_cross_norm.runTMM.R 1>&2 
    sh: 1: R: not found
    Error, cmd: R --no-save --no-restore --no-site-file --no-init-file -q < kallisto.isoform.TPM.not_cross_norm.runTMM.R 1>&2  died with ret (32512)  at /home/shared/trinityrnaseq-v2.12.0/util/support_scripts/run_TMM_scale_matrix.pl line 105.
    Error, CMD: /home/shared/trinityrnaseq-v2.12.0/util/support_scripts/run_TMM_scale_matrix.pl --matrix kallisto.isoform.TPM.not_cross_norm > kallisto.isoform.TMM.EXPR.matrix died with ret 6400 at /home/shared/trinityrnaseq-v2.12.0/util/abundance_estimates_to_matrix.pl line 385.
    total 56M
    -rw-r--r-- 1 shedurkin labmembers 2.0M Jan 30 18:57 kallisto.isoform.counts.matrix
    -rw-r--r-- 1 shedurkin labmembers    0 Jan 30 18:57 kallisto.isoform.TMM.EXPR.matrix
    -rw-r--r-- 1 shedurkin labmembers 2.4M Jan 30 18:57 kallisto.isoform.TPM.not_cross_norm
    -rw-r--r-- 1 shedurkin labmembers  532 Jan 30 18:57 kallisto.isoform.TPM.not_cross_norm.runTMM.R
    drwxr-xr-x 2 shedurkin labmembers 4.0K Jan 30 18:38 kallisto_quant_sample47
    drwxr-xr-x 2 shedurkin labmembers 4.0K Jan 30 18:42 kallisto_quant_sample48
    drwxr-xr-x 2 shedurkin labmembers 4.0K Jan 30 18:46 kallisto_quant_sample50
    drwxr-xr-x 2 shedurkin labmembers 4.0K Jan 30 18:50 kallisto_quant_sample53
    drwxr-xr-x 2 shedurkin labmembers 4.0K Jan 30 18:53 kallisto_quant_sample57
    -rw-r--r-- 1 shedurkin labmembers  51M Jan 30 18:57 Pmea_kallisto_index.idx

# Summary
