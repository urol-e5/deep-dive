10-Peve-sRNAseq-BLASTn
================
Sam White
2023-11-16

- <a href="#1-create-a-bash-variables-file"
  id="toc-1-create-a-bash-variables-file">1 Create a Bash variables
  file</a>
- <a href="#2-download-mirgenedb-fasta"
  id="toc-2-download-mirgenedb-fasta">2 Download MirGeneDB Fasta</a>
  - <a href="#21-inspect-mirna-fastas" id="toc-21-inspect-mirna-fastas">2.1
    Inspect miRNA FastAs</a>
- <a href="#3-convert-u-to-t-in-mirna-fastas"
  id="toc-3-convert-u-to-t-in-mirna-fastas">3 Convert <code>U</code> to
  <code>T</code> in miRNA FastAs</a>
- <a href="#4-create-blast-databases" id="toc-4-create-blast-databases">4
  Create BLAST Databases</a>
- <a href="#5-prepare-reads-for-blasting"
  id="toc-5-prepare-reads-for-blasting">5 Prepare reads for BLASTing</a>
  - <a href="#51-concatenate-all-trimmed-reads"
    id="toc-51-concatenate-all-trimmed-reads">5.1 Concatenate all trimmed
    reads</a>
  - <a href="#52-collapse-reads-to-fasta"
    id="toc-52-collapse-reads-to-fasta">5.2 Collapse reads to FastA</a>
- <a href="#6-run-blastn-default-e-value"
  id="toc-6-run-blastn-default-e-value">6 Run BLASTn Default E-value</a>
  - <a href="#61-mirbase-blastn-default-e-value"
    id="toc-61-mirbase-blastn-default-e-value">6.1 miRBase BLASTn Default
    e-value</a>
  - <a href="#62-mirgene-blastn-default-e-value"
    id="toc-62-mirgene-blastn-default-e-value">6.2 MirGene BLASTn Default
    e-value</a>
- <a href="#7-blastn-e-value--10" id="toc-7-blastn-e-value--10">7 BLASTn
  E-value = 10</a>
  - <a href="#71-mirbase-blastn-e-value--10"
    id="toc-71-mirbase-blastn-e-value--10">7.1 miRBase BLASTn e-value =
    10</a>
  - <a href="#72-mirgene-blastn-e-value--10"
    id="toc-72-mirgene-blastn-e-value--10">7.2 MirGene BLASTn e-value =
    10</a>
- <a href="#8-results" id="toc-8-results">8 Results</a>
  - <a href="#81-check-blastn-default-e-value-results"
    id="toc-81-check-blastn-default-e-value-results">8.1 Check BLASTn
    Default e-value results</a>
  - <a href="#82-check-blastn-e-value--10-results"
    id="toc-82-check-blastn-e-value--10-results">8.2 Check BLASTn e-value =
    10 results</a>
- <a href="#9-citations" id="toc-9-citations">9 Citations</a>

This notebook performs a simple [NCBI
BLASTn](https://www.ncbi.nlm.nih.gov/books/NBK279690/) \[@altschul1990\]
against two miRNA databases to attempt to identify miRNA in *A.pulchra*
sRNAseq:

- [miRBase](https://mirbase.org/download/)

- [MirGeneDB](https://www.mirgenedb.org/download)

Relies on the following software:

- [fastx_toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)

  - `fastx_collapser`: Collapses duplicate sequences in FastA/Q into
    single sequence.

# 1 Create a Bash variables file

This allows usage of Bash variables across R Markdown chunks.

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "# Trimmed FastQ naming pattern"
echo "export trimmed_fastqs_pattern='*flexbar_trim.25bp*.fastq.gz'"

echo "# Data directories"
echo 'export deep_dive_dir=/home/shared/8TB_HDD_01/sam/gitrepos/deep-dive'
echo 'export deep_dive_data_dir="${deep_dive_dir}/data"'
echo 'export output_dir_top=${deep_dive_dir}/E-Peve/output/10-Peve-sRNAseq-BLASTn'
echo 'export trimmed_fastqs_dir="${deep_dive_dir}/E-Peve/output/06-Peve-sRNAseq-trimming/trimmed-reads"'
echo 'export blast_dbs_dir="${deep_dive_dir}/data/blast_dbs"'
echo ""

echo "# Input/Output files"
echo 'export collapsed_reads_fasta="collapsed-reads-all.fasta"'
echo 'export concatenated_trimmed_reads_fastq="concatenated-trimmed-reads-all.fastq.gz"'
echo 'export mirbase_mature_fasta_name="mirbase-mature-v22.1.fa"'
echo 'export mirbase_mature_fasta_no_U="mirbase-mature-v22.1-no_U.fa"'
echo 'export mirgene_mature_fasta_name="mirgene-mature-all-v2.1.fa"'
echo 'export mirgene_mature_fasta_no_U="mirgene-mature-all-v2.1-no_U.fa"'
echo ""

echo "# External data URLs"
echo 'export mirgenedb_fasta_url="https://www.mirgenedb.org/fasta/ALL?mat=1"'
echo ""

echo "# Paths to programs"
echo 'export ncbi_blast_dir="/home/shared/ncbi-blast-2.15.0+/bin/"'
echo 'export ncbi_blastn="${ncbi_blast_dir}/blastn"'
echo 'export ncbi_makeblast_db="${ncbi_blast_dir}/makeblastdb"'
echo 'export fastx_collapser="/home/shared/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/bin/fastx_collapser"'

echo "# Set number of CPUs to use"
echo 'export threads=46'
echo ""

echo "# Initialize arrays"
echo 'export trimmed_fastqs_array=()'


} > .bashvars

cat .bashvars
```

    #### Assign Variables ####

    # Trimmed FastQ naming pattern
    export trimmed_fastqs_pattern='*flexbar_trim.25bp*.fastq.gz'
    # Data directories
    export deep_dive_dir=/home/shared/8TB_HDD_01/sam/gitrepos/deep-dive
    export deep_dive_data_dir="${deep_dive_dir}/data"
    export output_dir_top=${deep_dive_dir}/E-Peve/output/10-Peve-sRNAseq-BLASTn
    export trimmed_fastqs_dir="${deep_dive_dir}/E-Peve/output/06-Peve-sRNAseq-trimming/trimmed-reads"
    export blast_dbs_dir="${deep_dive_dir}/data/blast_dbs"

    # Input/Output files
    export collapsed_reads_fasta="collapsed-reads-all.fasta"
    export concatenated_trimmed_reads_fastq="concatenated-trimmed-reads-all.fastq.gz"
    export mirbase_mature_fasta_name="mirbase-mature-v22.1.fa"
    export mirbase_mature_fasta_no_U="mirbase-mature-v22.1-no_U.fa"
    export mirgene_mature_fasta_name="mirgene-mature-all-v2.1.fa"
    export mirgene_mature_fasta_no_U="mirgene-mature-all-v2.1-no_U.fa"

    # External data URLs
    export mirgenedb_fasta_url="https://www.mirgenedb.org/fasta/ALL?mat=1"

    # Paths to programs
    export ncbi_blast_dir="/home/shared/ncbi-blast-2.15.0+/bin/"
    export ncbi_blastn="${ncbi_blast_dir}/blastn"
    export ncbi_makeblast_db="${ncbi_blast_dir}/makeblastdb"
    export fastx_collapser="/home/shared/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/bin/fastx_collapser"
    # Set number of CPUs to use
    export threads=46

    # Initialize arrays
    export trimmed_fastqs_array=()

# 2 Download MirGeneDB Fasta

``` bash
# Load bash variables into memory
source .bashvars

wget \
--no-check-certificate \
--continue \
--no-host-directories \
--no-directories \
--no-parent \
--quiet \
--execute robots=off \
--output-document ${deep_dive_data_dir}/${mirgene_mature_fasta_name} \
 ${mirgenedb_fasta_url}
 
 ls -lh ${deep_dive_data_dir}
```

    total 13M
    drwxr-xr-x 2 sam sam 4.0K Nov  7 14:36 blast_dbs
    -rw-r--r-- 1 sam sam 3.7M Nov 15 21:30 mirbase-mature-v22.1.fa
    -rw-r--r-- 1 sam sam 3.7M Nov 14 09:39 mirbase-mature-v22.1-no_spaces.fa
    -rw-r--r-- 1 sam sam 3.7M Nov 16 09:42 mirbase-mature-v22.1-no_U.fa
    -rw-r--r-- 1 sam sam 726K Nov  7 14:36 mirgene-mature-all-v2.1.fa
    -rw-r--r-- 1 sam sam 726K Nov 16 09:42 mirgene-mature-all-v2.1-no_U.fa

## 2.1 Inspect miRNA FastAs

``` bash
# Load bash variables into memory
source .bashvars

head "${deep_dive_data_dir}"/mir*.fa
```

    ==> /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/data/mirbase-mature-v22.1.fa <==
    >cel-let-7-5p MIMAT0000001 Caenorhabditis elegans let-7-5p
    UGAGGUAGUAGGUUGUAUAGUU
    >cel-let-7-3p MIMAT0015091 Caenorhabditis elegans let-7-3p
    CUAUGCAAUUUUCUACCUUACC
    >cel-lin-4-5p MIMAT0000002 Caenorhabditis elegans lin-4-5p
    UCCCUGAGACCUCAAGUGUGA
    >cel-lin-4-3p MIMAT0015092 Caenorhabditis elegans lin-4-3p
    ACACCUGGGCUCUCCGGGUACC
    >cel-miR-1-5p MIMAT0020301 Caenorhabditis elegans miR-1-5p
    CAUACUUCCUUACAUGCCCAUA

    ==> /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/data/mirbase-mature-v22.1-no_spaces.fa <==
    >cel-let-7-5p_MIMAT0000001_Caenorhabditis_elegans_let-7-5p
    UGAGGUAGUAGGUUGUAUAGUU
    >cel-let-7-3p_MIMAT0015091_Caenorhabditis_elegans_let-7-3p
    CUAUGCAAUUUUCUACCUUACC
    >cel-lin-4-5p_MIMAT0000002_Caenorhabditis_elegans_lin-4-5p
    UCCCUGAGACCUCAAGUGUGA
    >cel-lin-4-3p_MIMAT0015092_Caenorhabditis_elegans_lin-4-3p
    ACACCUGGGCUCUCCGGGUACC
    >cel-miR-1-5p_MIMAT0020301_Caenorhabditis_elegans_miR-1-5p
    CAUACUUCCUUACAUGCCCAUA

    ==> /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/data/mirbase-mature-v22.1-no_U.fa <==
    >cel-let-7-5p MIMAT0000001 Caenorhabditis elegans let-7-5p
    TGAGGTAGTAGGTTGTATAGTT
    >cel-let-7-3p MIMAT0015091 Caenorhabditis elegans let-7-3p
    CTATGCAATTTTCTACCTTACC
    >cel-lin-4-5p MIMAT0000002 Caenorhabditis elegans lin-4-5p
    TCCCTGAGACCTCAAGTGTGA
    >cel-lin-4-3p MIMAT0015092 Caenorhabditis elegans lin-4-3p
    ACACCTGGGCTCTCCGGGTACC
    >cel-miR-1-5p MIMAT0020301 Caenorhabditis elegans miR-1-5p
    CATACTTCCTTACATGCCCATA

    ==> /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/data/mirgene-mature-all-v2.1.fa <==
    >Aae-Bantam_3p
    UGAGAUCAUUUUGAAAGCUGAUU
    >Bge-Bantam_3p
    UGAGAUCAUUGUGAAAGCUGAUU
    >Bpl-Bantam_3p
    UGAGAUCAUUGUGAAAACUGAU
    >Cgi-Bantam_3p
    UGAGAUCAUUGUGAAAACUGAUU
    >Cte-Bantam_3p
    UGAGAUCAUUGUGAAAACUAAUC

    ==> /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/data/mirgene-mature-all-v2.1-no_U.fa <==
    >Aae-Bantam_3p
    TGAGATCATTTTGAAAGCTGATT
    >Bge-Bantam_3p
    TGAGATCATTGTGAAAGCTGATT
    >Bpl-Bantam_3p
    TGAGATCATTGTGAAAACTGAT
    >Cgi-Bantam_3p
    TGAGATCATTGTGAAAACTGATT
    >Cte-Bantam_3p
    TGAGATCATTGTGAAAACTAATC

# 3 Convert `U` to `T` in miRNA FastAs

This is needed because the sRNAseq sequences do *not* have uracils
(`U`) - they have thymines (`T`).

``` bash
# Load bash variables into memory
source .bashvars

# Convert miRBase FastA
sed '/^[^>]/s/U/T/g' "${deep_dive_data_dir}/${mirbase_mature_fasta_name}" \
> "${deep_dive_data_dir}/${mirbase_mature_fasta_no_U}"

# Convert MirGene FastA
sed '/^[^>]/s/U/T/g' "${deep_dive_data_dir}/${mirgene_mature_fasta_name}" \
> "${deep_dive_data_dir}/${mirgene_mature_fasta_no_U}"

head ${deep_dive_data_dir}/*.fa
  
```

    ==> /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/data/mirbase-mature-v22.1.fa <==
    >cel-let-7-5p MIMAT0000001 Caenorhabditis elegans let-7-5p
    UGAGGUAGUAGGUUGUAUAGUU
    >cel-let-7-3p MIMAT0015091 Caenorhabditis elegans let-7-3p
    CUAUGCAAUUUUCUACCUUACC
    >cel-lin-4-5p MIMAT0000002 Caenorhabditis elegans lin-4-5p
    UCCCUGAGACCUCAAGUGUGA
    >cel-lin-4-3p MIMAT0015092 Caenorhabditis elegans lin-4-3p
    ACACCUGGGCUCUCCGGGUACC
    >cel-miR-1-5p MIMAT0020301 Caenorhabditis elegans miR-1-5p
    CAUACUUCCUUACAUGCCCAUA

    ==> /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/data/mirbase-mature-v22.1-no_spaces.fa <==
    >cel-let-7-5p_MIMAT0000001_Caenorhabditis_elegans_let-7-5p
    UGAGGUAGUAGGUUGUAUAGUU
    >cel-let-7-3p_MIMAT0015091_Caenorhabditis_elegans_let-7-3p
    CUAUGCAAUUUUCUACCUUACC
    >cel-lin-4-5p_MIMAT0000002_Caenorhabditis_elegans_lin-4-5p
    UCCCUGAGACCUCAAGUGUGA
    >cel-lin-4-3p_MIMAT0015092_Caenorhabditis_elegans_lin-4-3p
    ACACCUGGGCUCUCCGGGUACC
    >cel-miR-1-5p_MIMAT0020301_Caenorhabditis_elegans_miR-1-5p
    CAUACUUCCUUACAUGCCCAUA

    ==> /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/data/mirbase-mature-v22.1-no_U.fa <==
    >cel-let-7-5p MIMAT0000001 Caenorhabditis elegans let-7-5p
    TGAGGTAGTAGGTTGTATAGTT
    >cel-let-7-3p MIMAT0015091 Caenorhabditis elegans let-7-3p
    CTATGCAATTTTCTACCTTACC
    >cel-lin-4-5p MIMAT0000002 Caenorhabditis elegans lin-4-5p
    TCCCTGAGACCTCAAGTGTGA
    >cel-lin-4-3p MIMAT0015092 Caenorhabditis elegans lin-4-3p
    ACACCTGGGCTCTCCGGGTACC
    >cel-miR-1-5p MIMAT0020301 Caenorhabditis elegans miR-1-5p
    CATACTTCCTTACATGCCCATA

    ==> /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/data/mirgene-mature-all-v2.1.fa <==
    >Aae-Bantam_3p
    UGAGAUCAUUUUGAAAGCUGAUU
    >Bge-Bantam_3p
    UGAGAUCAUUGUGAAAGCUGAUU
    >Bpl-Bantam_3p
    UGAGAUCAUUGUGAAAACUGAU
    >Cgi-Bantam_3p
    UGAGAUCAUUGUGAAAACUGAUU
    >Cte-Bantam_3p
    UGAGAUCAUUGUGAAAACUAAUC

    ==> /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/data/mirgene-mature-all-v2.1-no_U.fa <==
    >Aae-Bantam_3p
    TGAGATCATTTTGAAAGCTGATT
    >Bge-Bantam_3p
    TGAGATCATTGTGAAAGCTGATT
    >Bpl-Bantam_3p
    TGAGATCATTGTGAAAACTGAT
    >Cgi-Bantam_3p
    TGAGATCATTGTGAAAACTGATT
    >Cte-Bantam_3p
    TGAGATCATTGTGAAAACTAATC

# 4 Create BLAST Databases

``` bash
# Load bash variables into memory
source .bashvars



# miRBase BLAST DB
## Make sure output directory exists
if [ ! -d "${blast_dbs_dir}" ]; then
  mkdir --parents "${blast_dbs_dir}"
fi

## Check for pre-exising database
if [ ! -f "${blast_dbs_dir}/${mirgene_mature_fasta_no_U%.*}.blastdb.log" ]; then
  ${ncbi_makeblast_db} \
  -in ${deep_dive_data_dir}/${mirgene_mature_fasta_no_U} \
  -title ${mirgene_mature_fasta_no_U%.*} \
  -dbtype nucl \
  -out ${blast_dbs_dir}/${mirgene_mature_fasta_no_U%.*} \
  -logfile ${blast_dbs_dir}/${mirgene_mature_fasta_no_U%.*}.blastdb.log
fi

# miRBase BLAST DB
## Make sure output directory exists
if [ ! -d "${blast_dbs_dir}" ]; then
  mkdir --parents "${blast_dbs_dir}"
fi

## Check for pre-exising database
if [ ! -f "${blast_dbs_dir}/${mirbase_mature_fasta_no_U%.*}.blastdb.log" ]; then
  ${ncbi_makeblast_db} \
  -in ${deep_dive_data_dir}/${mirbase_mature_fasta_no_U} \
  -title ${mirbase_mature_fasta_no_U%.*} \
  -dbtype nucl \
  -out ${blast_dbs_dir}/${mirbase_mature_fasta_no_U%.*} \
  -logfile ${blast_dbs_dir}/${mirbase_mature_fasta_no_U%.*}.blastdb.log
fi
```

# 5 Prepare reads for BLASTing

## 5.1 Concatenate all trimmed reads

``` bash
# Load bash variables into memory
source .bashvars

# Check for existence of concatenated FastA before running
if [ ! -f "${output_dir_top}/${concatenated_trimmed_reads_fastq}" ]; then
  cat ${trimmed_fastqs_dir}/*.fastq.gz \
  > "${output_dir_top}/${concatenated_trimmed_reads_fastq}"
fi

ls -lh "${output_dir_top}/${concatenated_trimmed_reads_fastq}"
```

    -rw-r--r-- 1 sam sam 1.3G Nov 16 08:32 /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/E-Peve/output/10-Peve-sRNAseq-BLASTn/concatenated-trimmed-reads-all.fastq.gz

## 5.2 Collapse reads to FastA

Uses `fastx_collapser` to collapse to unique reads.

Requires undocumented quality setting. Have selected `30` as cuttoff:
`-Q30`.

``` bash
# Load bash variables into memory
source .bashvars

# Check for existence of collapsed FastA before running
time \
if [ ! -f "${output_dir_top}/${collapsed_reads_fasta}" ]; then
  zcat ${output_dir_top}/${concatenated_trimmed_reads_fastq} \
  | ${fastx_collapser} \
  -Q30 \
  -o "${output_dir_top}/${collapsed_reads_fasta}"
fi

head "${output_dir_top}/${collapsed_reads_fasta}"

echo ""
echo ""

total_reads=$(grep -c "^>" ${output_dir_top}/${collapsed_reads_fasta})
echo "Total number of reads after collapse: ${total_reads}"
```

    real    0m0.000s
    user    0m0.000s
    sys 0m0.000s
    >1-2463526
    GCACTGGTGGTTCAGTGGTAGAATT
    >2-2045325
    TGAAAATCTTTGCTCTGAAGTGGAA
    >3-2026343
    TTCCACTTCAGAGCAAAGATTTTCA
    >4-1419489
    GGCGAGAATTCTACCACTGAACCAC
    >5-1161968
    ACAAATCTTAGAACAAAGGCTTAAT


    Total number of reads after collapse: 8870343

# 6 Run BLASTn Default E-value

- 1000 for blastn-short

Runs BLASTn using the `blastn-short` task for sequences \< 30bp.

Look for top match (`-max_hsps 1` & `-max_target_seqs 1`) for each
query.

- Suppress subsequent warning
  `Examining 5 or more matches is recommended` by redirecting stdout:
  `2> /dev/null`

## 6.1 miRBase BLASTn Default e-value

``` bash
# Load bash variables into memory
source .bashvars

time \
${ncbi_blastn} \
-db ${blast_dbs_dir}/${mirbase_mature_fasta_no_U%.*} \
-query ${output_dir_top}/${collapsed_reads_fasta} \
-out ${output_dir_top}/miRBase-BLASTn-eval_1000.outfmt6 \
-task blastn-short \
-max_hsps 1 \
-max_target_seqs 1 \
-outfmt 6 \
-num_threads ${threads} \
2> /dev/null
```

## 6.2 MirGene BLASTn Default e-value

``` bash
# Load bash variables into memory
source .bashvars

time \
${ncbi_blastn} \
-db ${blast_dbs_dir}/${mirgene_mature_fasta_no_U%.*} \
-query ${output_dir_top}/${collapsed_reads_fasta} \
-out ${output_dir_top}/MirGene-BLASTn-eval_1000.outfmt6 \
-task blastn-short \
-max_hsps 1 \
-max_target_seqs 1 \
-outfmt 6 \
-num_threads ${threads} \
2> /dev/null
```

# 7 BLASTn E-value = 10

Running this for simple comparison to the default `blastn-short` value
of 1000.

## 7.1 miRBase BLASTn e-value = 10

``` bash
# Load bash variables into memory
source .bashvars

time \
${ncbi_blastn} \
-db ${blast_dbs_dir}/${mirbase_mature_fasta_no_U%.*} \
-query ${output_dir_top}/${collapsed_reads_fasta} \
-out ${output_dir_top}/miRBase-BLASTn-eval_10.outfmt6 \
-task blastn-short \
-evalue 10 \
-max_hsps 1 \
-max_target_seqs 1 \
-outfmt 6 \
-num_threads ${threads} \
2> /dev/null
```

## 7.2 MirGene BLASTn e-value = 10

``` bash
# Load bash variables into memory
source .bashvars

time \
${ncbi_blastn} \
-db ${blast_dbs_dir}/${mirgene_mature_fasta_no_U%.*} \
-query ${output_dir_top}/${collapsed_reads_fasta} \
-out ${output_dir_top}/MirGene-BLASTn-eval_10.outfmt6 \
-task blastn-short \
-evalue 10 \
-max_hsps 1 \
-max_target_seqs 1 \
-outfmt 6 \
-num_threads ${threads} \
2> /dev/null
```

# 8 Results

## 8.1 Check BLASTn Default e-value results

``` bash
# Load bash variables into memory
source .bashvars

head ${output_dir_top}/*eval_1000.outfmt6

echo ""
echo "-----------------------------------------"
echo ""

wc -l ${output_dir_top}/*eval_1000.outfmt6
```

    ==> /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/E-Peve/output/10-Peve-sRNAseq-BLASTn/miRBase-BLASTn-eval_1000.outfmt6 <==
    1-2463526   ppc-miR-8214-5p 100.000 13  0   0   12  24  18  6   0.092   26.3
    2-2045325   ssc-miR-7143-5p 100.000 11  0   0   12  22  8   18  1.4 22.3
    3-2026343   ssc-miR-7143-5p 100.000 11  0   0   4   14  18  8   1.4 22.3
    4-1419489   ppc-miR-8214-5p 100.000 13  0   0   8   20  6   18  0.092   26.3
    5-1161968   cin-miR-4068-3p 100.000 12  0   0   11  22  10  21  0.36    24.3
    6-1144435   lja-miR11150-3p 100.000 10  0   0   13  22  10  19  5.7 20.3
    7-1040735   hsa-miR-3168    100.000 13  0   0   4   16  13  1   0.092   26.3
    8-942821    efu-miR-9230a   100.000 12  0   0   14  25  2   13  0.36    24.3
    9-658343    ppc-miR-8214-5p 100.000 13  0   0   8   20  6   18  0.092   26.3
    10-653354   gma-miR9757 100.000 10  0   0   13  22  11  2   5.7 20.3

    ==> /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/E-Peve/output/10-Peve-sRNAseq-BLASTn/MirGene-BLASTn-eval_1000.outfmt6 <==
    1-2463526   Pma-Mir-96-P3o2_5p  100.000 10  0   0   16  25  9   18  2.2 20.3
    2-2045325   Tca-Mir-87-P2_3p    100.000 11  0   0   5   15  13  3   0.54    22.3
    3-2026343   Tca-Mir-87-P2_3p    100.000 11  0   0   11  21  3   13  0.54    22.3
    4-1419489   Pma-Mir-96-P3o2_5p  100.000 11  0   0   6   16  19  9   0.54    22.3
    5-1161968   Cin-Mir-4062-P1_3p  100.000 12  0   0   11  22  10  21  0.14    24.3
    6-1144435   Dno-Mir-1271_5p 100.000 9   0   0   1   9   21  13  8.5 18.3
    7-1040735   Hme-Mir-2-P7_3p 100.000 11  0   0   9   19  19  9   0.54    22.3
    8-942821    Cbr-Mir-279-o10_5p  100.000 10  0   0   13  22  11  2   2.2 20.3
    9-658343    Pma-Mir-96-P3o2_5p  100.000 11  0   0   6   16  19  9   0.54    22.3
    10-653354   Gga-Mir-7441_5p 100.000 10  0   0   13  22  3   12  2.2 20.3

    -----------------------------------------

       8870343 /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/E-Peve/output/10-Peve-sRNAseq-BLASTn/miRBase-BLASTn-eval_1000.outfmt6
       8870343 /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/E-Peve/output/10-Peve-sRNAseq-BLASTn/MirGene-BLASTn-eval_1000.outfmt6
      17740686 total

## 8.2 Check BLASTn e-value = 10 results

``` bash
# Load bash variables into memory
source .bashvars

head ${output_dir_top}/*eval_10.outfmt6

echo ""
echo "-----------------------------------------"
echo ""

wc -l ${output_dir_top}/*eval_10.outfmt6
```

    ==> /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/E-Peve/output/10-Peve-sRNAseq-BLASTn/miRBase-BLASTn-eval_10.outfmt6 <==
    1-2463526   ppc-miR-8214-5p 100.000 13  0   0   12  24  18  6   0.092   26.3
    2-2045325   ssc-miR-7143-5p 100.000 11  0   0   12  22  8   18  1.4 22.3
    3-2026343   ssc-miR-7143-5p 100.000 11  0   0   4   14  18  8   1.4 22.3
    4-1419489   ppc-miR-8214-5p 100.000 13  0   0   8   20  6   18  0.092   26.3
    5-1161968   cin-miR-4068-3p 100.000 12  0   0   11  22  10  21  0.36    24.3
    6-1144435   lja-miR11150-3p 100.000 10  0   0   13  22  10  19  5.7 20.3
    7-1040735   hsa-miR-3168    100.000 13  0   0   4   16  13  1   0.092   26.3
    8-942821    efu-miR-9230a   100.000 12  0   0   14  25  2   13  0.36    24.3
    9-658343    ppc-miR-8214-5p 100.000 13  0   0   8   20  6   18  0.092   26.3
    10-653354   gma-miR9757 100.000 10  0   0   13  22  11  2   5.7 20.3

    ==> /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/E-Peve/output/10-Peve-sRNAseq-BLASTn/MirGene-BLASTn-eval_10.outfmt6 <==
    1-2463526   Pma-Mir-96-P3o2_5p  100.000 10  0   0   16  25  9   18  2.2 20.3
    2-2045325   Tca-Mir-87-P2_3p    100.000 11  0   0   5   15  13  3   0.54    22.3
    3-2026343   Tca-Mir-87-P2_3p    100.000 11  0   0   11  21  3   13  0.54    22.3
    4-1419489   Pma-Mir-96-P3o2_5p  100.000 11  0   0   6   16  19  9   0.54    22.3
    5-1161968   Cin-Mir-4062-P1_3p  100.000 12  0   0   11  22  10  21  0.14    24.3
    6-1144435   Dno-Mir-1271_5p 100.000 9   0   0   1   9   21  13  8.5 18.3
    7-1040735   Hme-Mir-2-P7_3p 100.000 11  0   0   9   19  19  9   0.54    22.3
    8-942821    Cbr-Mir-279-o10_5p  100.000 10  0   0   13  22  11  2   2.2 20.3
    9-658343    Pma-Mir-96-P3o2_5p  100.000 11  0   0   6   16  19  9   0.54    22.3
    10-653354   Gga-Mir-7441_5p 100.000 10  0   0   13  22  3   12  2.2 20.3

    -----------------------------------------

       8824359 /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/E-Peve/output/10-Peve-sRNAseq-BLASTn/miRBase-BLASTn-eval_10.outfmt6
       8783659 /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/E-Peve/output/10-Peve-sRNAseq-BLASTn/MirGene-BLASTn-eval_10.outfmt6
      17608018 total

The default e-value results in alignments for all query sequences, which
is likely not what we’d expect.

Decreasing the e-value resulted in fewer query alignments, as we’d
expect.

For further analysis, we should probably discuss a reasonable e-value to
use for filtering.

------------------------------------------------------------------------

# 9 Citations
