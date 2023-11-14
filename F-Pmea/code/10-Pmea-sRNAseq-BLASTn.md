10-Pmea-sRNAseq-BLASTn
================
Sam White (modified by K Durkin for P. meandrina analysis)
2023-11-14

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
BLASTn](https://www.ncbi.nlm.nih.gov/books/NBK279690/) ([Altschul et al.
1990](#ref-altschul1990)) against two miRNA databases to attempt to
identify miRNA in *P.meandrina* sRNAseq:

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
echo 'export deep_dive_dir=/home/shared/8TB_HDD_02/shedurkin/deep-dive'
echo 'export deep_dive_data_dir="${deep_dive_dir}/data"'
echo 'export output_dir_top=${deep_dive_dir}/F-Pmea/output/10-Pmea-sRNAseq-BLASTn'
echo 'export trimmed_fastqs_dir="${deep_dive_dir}/F-Pmea/output/08-Pmea-sRNAseq-trimming/trimmed-reads"'
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
    export deep_dive_dir=/home/shared/8TB_HDD_02/shedurkin/deep-dive
    export deep_dive_data_dir="${deep_dive_dir}/data"
    export output_dir_top=${deep_dive_dir}/F-Pmea/output/10-Pmea-sRNAseq-BLASTn
    export trimmed_fastqs_dir="${deep_dive_dir}/F-Pmea/output/08-Pmea-sRNAseq-trimming/trimmed-reads"
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
    drwxr-xr-x 2 shedurkin labmembers 4.0K Nov 14 09:39 blast_dbs
    -rw-r--r-- 1 shedurkin labmembers 3.7M Nov 14 09:39 mirbase-mature-v22.1.fa
    -rw-r--r-- 1 shedurkin labmembers 3.7M Nov 14 09:39 mirbase-mature-v22.1-no_spaces.fa
    -rw-r--r-- 1 shedurkin labmembers 3.7M Nov 14 15:33 mirbase-mature-v22.1-no_U.fa
    -rw-r--r-- 1 shedurkin labmembers 726K Nov 14 09:39 mirgene-mature-all-v2.1.fa
    -rw-r--r-- 1 shedurkin labmembers 726K Nov 14 15:33 mirgene-mature-all-v2.1-no_U.fa

## 2.1 Inspect miRNA FastAs

``` bash
# Load bash variables into memory
source .bashvars

head "${deep_dive_data_dir}"/mir*.fa
```

    ==> /home/shared/8TB_HDD_02/shedurkin/deep-dive/data/mirbase-mature-v22.1.fa <==
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

    ==> /home/shared/8TB_HDD_02/shedurkin/deep-dive/data/mirbase-mature-v22.1-no_spaces.fa <==
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

    ==> /home/shared/8TB_HDD_02/shedurkin/deep-dive/data/mirbase-mature-v22.1-no_U.fa <==
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

    ==> /home/shared/8TB_HDD_02/shedurkin/deep-dive/data/mirgene-mature-all-v2.1.fa <==
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

    ==> /home/shared/8TB_HDD_02/shedurkin/deep-dive/data/mirgene-mature-all-v2.1-no_U.fa <==
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

    ==> /home/shared/8TB_HDD_02/shedurkin/deep-dive/data/mirbase-mature-v22.1.fa <==
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

    ==> /home/shared/8TB_HDD_02/shedurkin/deep-dive/data/mirbase-mature-v22.1-no_spaces.fa <==
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

    ==> /home/shared/8TB_HDD_02/shedurkin/deep-dive/data/mirbase-mature-v22.1-no_U.fa <==
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

    ==> /home/shared/8TB_HDD_02/shedurkin/deep-dive/data/mirgene-mature-all-v2.1.fa <==
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

    ==> /home/shared/8TB_HDD_02/shedurkin/deep-dive/data/mirgene-mature-all-v2.1-no_U.fa <==
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

    -rw-r--r-- 1 shedurkin labmembers 2.0G Nov 14 11:56 /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/10-Pmea-sRNAseq-BLASTn/concatenated-trimmed-reads-all.fastq.gz

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
```

    real    0m0.000s
    user    0m0.000s
    sys 0m0.000s
    >1-10280803
    TGAAAATCTTTGCTCTGAAGTGGAA
    >2-10220688
    TTCCACTTCAGAGCAAAGATTTTCA
    >3-4878221
    GCACTGGTGGTTCAGTGGTAGAATT
    >4-2964293
    CTCGGGCTGAGACTTGAAGCG
    >5-2941994
    CGCTTCAAGTCTCAGCCCGAG

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

Running this for simple comparison to the defaul `blastn-short` value of
1000.

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

    ==> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/10-Pmea-sRNAseq-BLASTn/miRBase-BLASTn-eval_1000.outfmt6 <==
    1-10280803  ssc-miR-7143-5p 100.000 11  0   0   12  22  8   18  1.4 22.3
    2-10220688  ssc-miR-7143-5p 100.000 11  0   0   4   14  18  8   1.4 22.3
    3-4878221   ppc-miR-8214-5p 100.000 13  0   0   12  24  18  6   0.092   26.3
    4-2964293   cpo-miR-1379-5p 100.000 12  0   0   1   12  18  7   0.26    24.3
    5-2941994   cpo-miR-1379-5p 100.000 12  0   0   10  21  7   18  0.26    24.3
    6-2745200   pab-miR11514    100.000 10  0   0   4   13  16  7   5.7 20.3
    7-2697913   rno-miR-3561-3p 100.000 12  0   0   2   13  10  21  0.36    24.3
    8-2084846   hsa-miR-3168    100.000 13  0   0   4   16  13  1   0.092   26.3
    9-2015242   efu-miR-9230a   100.000 12  0   0   14  25  2   13  0.36    24.3
    10-1645687  ppc-miR-8214-5p 100.000 13  0   0   8   20  6   18  0.092   26.3

    ==> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/10-Pmea-sRNAseq-BLASTn/MirGene-BLASTn-eval_1000.outfmt6 <==
    1-10280803  Tca-Mir-87-P2_3p    100.000 11  0   0   5   15  13  3   0.54    22.3
    2-10220688  Tca-Mir-87-P2_3p    100.000 11  0   0   11  21  3   13  0.54    22.3
    3-4878221   Pma-Mir-96-P3o2_5p  100.000 10  0   0   16  25  9   18  2.2 20.3
    4-2964293   Ovu-Novel-55_5p 100.000 9   0   0   10  18  14  6   7.3 18.3
    5-2941994   Ovu-Novel-55_5p 100.000 9   0   0   4   12  6   14  7.3 18.3
    6-2745200   Ebu-Mir-1329_5p 100.000 10  0   0   4   13  14  23  2.2 20.3
    7-2697913   Cin-Mir-4047_5p 100.000 11  0   0   3   13  3   13  0.54    22.3
    8-2084846   Hme-Mir-2-P7_3p 100.000 11  0   0   9   19  19  9   0.54    22.3
    9-2015242   Cbr-Mir-279-o10_5p  100.000 10  0   0   13  22  11  2   2.2 20.3
    10-1645687  Pma-Mir-96-P3o2_5p  100.000 11  0   0   6   16  19  9   0.54    22.3

    -----------------------------------------

      13777166 /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/10-Pmea-sRNAseq-BLASTn/miRBase-BLASTn-eval_1000.outfmt6
      13777166 /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/10-Pmea-sRNAseq-BLASTn/MirGene-BLASTn-eval_1000.outfmt6
      27554332 total

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

    ==> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/10-Pmea-sRNAseq-BLASTn/miRBase-BLASTn-eval_10.outfmt6 <==
    1-10280803  ssc-miR-7143-5p 100.000 11  0   0   12  22  8   18  1.4 22.3
    2-10220688  ssc-miR-7143-5p 100.000 11  0   0   4   14  18  8   1.4 22.3
    3-4878221   ppc-miR-8214-5p 100.000 13  0   0   12  24  18  6   0.092   26.3
    4-2964293   cpo-miR-1379-5p 100.000 12  0   0   1   12  18  7   0.26    24.3
    5-2941994   cpo-miR-1379-5p 100.000 12  0   0   10  21  7   18  0.26    24.3
    6-2745200   pab-miR11514    100.000 10  0   0   4   13  16  7   5.7 20.3
    7-2697913   rno-miR-3561-3p 100.000 12  0   0   2   13  10  21  0.36    24.3
    8-2084846   hsa-miR-3168    100.000 13  0   0   4   16  13  1   0.092   26.3
    9-2015242   efu-miR-9230a   100.000 12  0   0   14  25  2   13  0.36    24.3
    10-1645687  ppc-miR-8214-5p 100.000 13  0   0   8   20  6   18  0.092   26.3

    ==> /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/10-Pmea-sRNAseq-BLASTn/MirGene-BLASTn-eval_10.outfmt6 <==
    1-10280803  Tca-Mir-87-P2_3p    100.000 11  0   0   5   15  13  3   0.54    22.3
    2-10220688  Tca-Mir-87-P2_3p    100.000 11  0   0   11  21  3   13  0.54    22.3
    3-4878221   Pma-Mir-96-P3o2_5p  100.000 10  0   0   16  25  9   18  2.2 20.3
    4-2964293   Ovu-Novel-55_5p 100.000 9   0   0   10  18  14  6   7.3 18.3
    5-2941994   Ovu-Novel-55_5p 100.000 9   0   0   4   12  6   14  7.3 18.3
    6-2745200   Ebu-Mir-1329_5p 100.000 10  0   0   4   13  14  23  2.2 20.3
    7-2697913   Cin-Mir-4047_5p 100.000 11  0   0   3   13  3   13  0.54    22.3
    8-2084846   Hme-Mir-2-P7_3p 100.000 11  0   0   9   19  19  9   0.54    22.3
    9-2015242   Cbr-Mir-279-o10_5p  100.000 10  0   0   13  22  11  2   2.2 20.3
    10-1645687  Pma-Mir-96-P3o2_5p  100.000 11  0   0   6   16  19  9   0.54    22.3

    -----------------------------------------

      13708946 /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/10-Pmea-sRNAseq-BLASTn/miRBase-BLASTn-eval_10.outfmt6
      13614786 /home/shared/8TB_HDD_02/shedurkin/deep-dive/F-Pmea/output/10-Pmea-sRNAseq-BLASTn/MirGene-BLASTn-eval_10.outfmt6
      27323732 total

The default e-value results in alignments for all query sequences, which
is likely not what we’d expect.

Decreasing the e-value resulted in fewer queary alignments, as we’d
expect.

For further analysis, we should probably discuss a reasonable e-value to
use for filtering.

------------------------------------------------------------------------

# 9 Citations

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-altschul1990" class="csl-entry">

Altschul, Stephen F., Warren Gish, Webb Miller, Eugene W. Myers, and
David J. Lipman. 1990. “Basic Local Alignment Search Tool.” *Journal of
Molecular Biology* 215 (3): 403–10.
<https://doi.org/10.1016/s0022-2836(05)80360-2>.

</div>

</div>
