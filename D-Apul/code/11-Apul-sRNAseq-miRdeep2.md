11-Apul-sRNAseq-miRdeep2
================
Sam White
2023-11-08

- <a href="#1-create-a-bash-variables-file"
  id="toc-1-create-a-bash-variables-file">1 Create a Bash variables
  file</a>
- <a href="#2-prepare-reads-for-mirdeep2"
  id="toc-2-prepare-reads-for-mirdeep2">2 Prepare reads for miRDeep2</a>
- <a href="#3-reformat-genome-fasta-description-lines"
  id="toc-3-reformat-genome-fasta-description-lines">3 Reformat genome
  FastA description lines</a>
- <a href="#4-reformat-mirbase-fasta-description-lines"
  id="toc-4-reformat-mirbase-fasta-description-lines">4 Reformat miRBase
  FastA description lines</a>
- <a href="#5-bowtie-v1-genome-index" id="toc-5-bowtie-v1-genome-index">5
  Bowtie v1 genome index</a>
- <a href="#6-map-reads-to-genome" id="toc-6-map-reads-to-genome">6 Map
  reads to genome</a>
- <a href="#7-run-mirdeep2" id="toc-7-run-mirdeep2">7 Run miRDeep2</a>
  - <a href="#71-check-runtime" id="toc-71-check-runtime">7.1 Check
    runtime</a>
- <a href="#8-move-output-files-to-output-directory"
  id="toc-8-move-output-files-to-output-directory">8 Move output files to
  output directory</a>
  - <a href="#81-move-output-files" id="toc-81-move-output-files">8.1 Move
    output files</a>
  - <a href="#82-identify-directories" id="toc-82-identify-directories">8.2
    Identify directories</a>
  - <a href="#83-rsync-directories" id="toc-83-rsync-directories">8.3 Rsync
    directories</a>
  - <a href="#84-confirm-deletion-patterns-work-before-deletion"
    id="toc-84-confirm-deletion-patterns-work-before-deletion">8.4 Confirm
    deletion patterns work <em>before</em> deletion!</a>
  - <a href="#85-remove-directores-from-code-directory"
    id="toc-85-remove-directores-from-code-directory">8.5 Remove directores
    from code directory</a>
  - <a href="#86-check-to-make-sure-theyre-gone"
    id="toc-86-check-to-make-sure-theyre-gone">8.6 Check to make sure
    they’re gone</a>
- <a href="#9-results" id="toc-9-results">9 Results</a>
  - <a href="#91-peep-output-file-format"
    id="toc-91-peep-output-file-format">9.1 Peep output file format</a>
  - <a href="#92-create-more-easily-parasable-results-file"
    id="toc-92-create-more-easily-parasable-results-file">9.2 Create more
    easily parasable results file</a>
  - <a href="#93-read-in-output-csv" id="toc-93-read-in-output-csv">9.3 Read
    in output CSV</a>
  - <a href="#94-mirnas-count-data" id="toc-94-mirnas-count-data">9.4 miRNAs
    count data</a>
- <a href="#10-citations" id="toc-10-citations">10 Citations</a>

Use [miRDeep2](https://github.com/rajewsky-lab/mirdeep2) ([Friedländer
et al. 2011](#ref-friedländer2011)) to identify potential miRNAs using
*A.pulchra* sRNAseq reads. The *A.millepora* genome will be used as the
reference genome for *A.pulchra*, as *A.pulchra* does not currently have
a sequenced genome and *A.millepora* had highest alignment rates for
standard RNAseq data compared to other published genomes tested. The
mature miRBase miRNA database will also be used.

------------------------------------------------------------------------

Inputs:

- Requires collapsed reads (i.e. concatenated, unique reads) in FastA
  format. See
  [10-Apul-sRNAseq-BLASTn.Rmd](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/code/10-Apul-sRNAseq-BLASTn.Rmd)
  for code.

- *A.millepora* genome FastA. See
  [12-Apul-sRNAseq-MirMachine.Rmd](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/code/12-Apul-sRNAseq-MirMachine.Rmd)
  for download info if needed.

- [miRBase](https://mirbase.org/download/) mature miRNAs FastA.

Outputs:

- Primary outputs are a result table in BED, CSV (tab-delimited), and
  HTML formats.

- Due to the nature of mirDeep2’s naming, trying to use variable names
  is challenging. As such, the chunks processing those files will
  require manual intervention to identify and provide the output
  filename(s); they are *not* handled at in the `.bashvars` file at the
  top of this script.

- Other output files are too large for GitHub (and some (all?) are not
  needed for this analysis). Please find a full backup here:

<https://gannet.fish.washington.edu/Atumefaciens/gitrepos/deep-dive/D-Apul/output/11-Apul-sRNAseq-miRdeep2/>

------------------------------------------------------------------------

# 1 Create a Bash variables file

This allows usage of Bash variables across R Markdown chunks.

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "# Trimmed FastQ naming pattern"
echo "export trimmed_fastqs_pattern='*flexbar_trim.25bp*.fastq.gz'"
echo ""

echo "# miRTrace FastA naming pattern"
echo "export mirtrace_fasta_pattern='*flexbar_trim.25bp*.fasta'"

echo "# Data directories"
echo 'export deep_dive_dir=/home/shared/8TB_HDD_01/sam/gitrepos/deep-dive'
echo 'export deep_dive_data_dir="${deep_dive_dir}/data"'
echo 'export output_dir_top=${deep_dive_dir}/D-Apul/output/11-Apul-sRNAseq-miRdeep2'
echo 'export genome_fasta_dir=${deep_dive_dir}/D-Apul/data/Amil/ncbi_dataset/data/GCF_013753865.1'
echo 'export trimmed_fastqs_dir="${deep_dive_dir}/D-Apul/output/08-Apul-sRNAseq-trimming/trimmed-reads"'
echo 'export collapsed_reads_dir="${deep_dive_dir}/D-Apul/output/10-Apul-sRNAseq-BLASTn"'
echo ""

echo "# Input/Output files"
echo 'export collapsed_reads_fasta="collapsed-reads-all.fasta"'
echo 'export collapsed_reads_mirdeep2="collapsed-reads-all-mirdeep2.fasta"'
echo 'export concatenated_trimmed_reads_fastq="concatenated-trimmed-reads-all.fastq.gz"'
echo 'export genome_fasta_name="GCF_013753865.1_Amil_v2.1_genomic.fna"'
echo 'export genome_fasta_no_spaces="GCF_013753865.1_Amil_v2.1_genomic-no_spaces.fna"'
echo 'export mirdeep2_mapping_file="Apul-mirdeep2-mapping.arf"'
echo 'export mirbase_mature_fasta_name="mirbase-mature-v22.1.fa"'
echo 'export mirbase_mature_fasta_no_spaces="mirbase-mature-v22.1-no_spaces.fa"'
echo 'export concatenated_mirtrace_reads_fasta="concatenated-mirtrace-reads.fasta"'
echo ""


echo "# Paths to programs"
echo 'export mirdeep2_mapper="mapper.pl"'
echo 'export mirdeep2="miRDeep2.pl"'
echo 'export bowtie_build="/home/shared/bowtie-1.3.1-linux-x86_64/bowtie-build"'
echo ""

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

    # miRTrace FastA naming pattern
    export mirtrace_fasta_pattern='*flexbar_trim.25bp*.fasta'
    # Data directories
    export deep_dive_dir=/home/shared/8TB_HDD_01/sam/gitrepos/deep-dive
    export deep_dive_data_dir="${deep_dive_dir}/data"
    export output_dir_top=${deep_dive_dir}/D-Apul/output/11-Apul-sRNAseq-miRdeep2
    export genome_fasta_dir=${deep_dive_dir}/D-Apul/data/Amil/ncbi_dataset/data/GCF_013753865.1
    export trimmed_fastqs_dir="${deep_dive_dir}/D-Apul/output/08-Apul-sRNAseq-trimming/trimmed-reads"
    export collapsed_reads_dir="${deep_dive_dir}/D-Apul/output/10-Apul-sRNAseq-BLASTn"

    # Input/Output files
    export collapsed_reads_fasta="collapsed-reads-all.fasta"
    export collapsed_reads_mirdeep2="collapsed-reads-all-mirdeep2.fasta"
    export concatenated_trimmed_reads_fastq="concatenated-trimmed-reads-all.fastq.gz"
    export genome_fasta_name="GCF_013753865.1_Amil_v2.1_genomic.fna"
    export genome_fasta_no_spaces="GCF_013753865.1_Amil_v2.1_genomic-no_spaces.fna"
    export mirdeep2_mapping_file="Apul-mirdeep2-mapping.arf"
    export mirbase_mature_fasta_name="mirbase-mature-v22.1.fa"
    export mirbase_mature_fasta_no_spaces="mirbase-mature-v22.1-no_spaces.fa"
    export concatenated_mirtrace_reads_fasta="concatenated-mirtrace-reads.fasta"

    # Paths to programs
    export mirdeep2_mapper="mapper.pl"
    export mirdeep2="miRDeep2.pl"
    export bowtie_build="/home/shared/bowtie-1.3.1-linux-x86_64/bowtie-build"

    # Set number of CPUs to use
    export threads=46

    # Initialize arrays
    export trimmed_fastqs_array=()

# 2 Prepare reads for miRDeep2

Per miRDeep2 documentation:

> The readID must end with \_xNumber and is not allowed to contain
> whitespaces. has to have the format name_uniqueNumber_xnumber

``` bash
# Load bash variables into memory
source .bashvars


sed '/^>/ s/-/_x/g' "${collapsed_reads_dir}/${collapsed_reads_fasta}" \
| sed '/^>/ s/>/>seq_/' \
> "${output_dir_top}/${collapsed_reads_mirdeep2}"

grep "^>" ${collapsed_reads_dir}/${collapsed_reads_fasta} \
| head

echo ""
echo "--------------------------------------------------"
echo ""

grep "^>" ${output_dir_top}/${collapsed_reads_mirdeep2} \
| head
```

    >1-5404223
    >2-2739211
    >3-2687482
    >4-2632527
    >5-2254991
    >6-2220700
    >7-2187568
    >8-1457045
    >9-1086004
    >10-1079942

    --------------------------------------------------

    >seq_1_x5404223
    >seq_2_x2739211
    >seq_3_x2687482
    >seq_4_x2632527
    >seq_5_x2254991
    >seq_6_x2220700
    >seq_7_x2187568
    >seq_8_x1457045
    >seq_9_x1086004
    >seq_10_x1079942

# 3 Reformat genome FastA description lines

miRDeep2 can’t process genome FastAs with spaces in the description
lines.

So, I’m replacing spaces with underscores.

And, for aesthetics, I’m also removing commas.

``` bash
# Load bash variables into memory
source .bashvars


sed '/^>/ s/ /_/g' "${genome_fasta_dir}/${genome_fasta_name}" \
| sed '/^>/ s/,//g' \
> "${genome_fasta_dir}/${genome_fasta_no_spaces}"

grep "^>" ${genome_fasta_dir}/${genome_fasta_name} \
| head

echo ""
echo "--------------------------------------------------"
echo ""

grep "^>" ${genome_fasta_dir}/${genome_fasta_no_spaces} \
| head
```

    >NC_058066.1 Acropora millepora isolate JS-1 chromosome 1, Amil_v2.1, whole genome shotgun sequence
    >NC_058067.1 Acropora millepora isolate JS-1 chromosome 2, Amil_v2.1, whole genome shotgun sequence
    >NC_058068.1 Acropora millepora isolate JS-1 chromosome 3, Amil_v2.1, whole genome shotgun sequence
    >NC_058069.1 Acropora millepora isolate JS-1 chromosome 4, Amil_v2.1, whole genome shotgun sequence
    >NC_058070.1 Acropora millepora isolate JS-1 chromosome 5, Amil_v2.1, whole genome shotgun sequence
    >NC_058071.1 Acropora millepora isolate JS-1 chromosome 6, Amil_v2.1, whole genome shotgun sequence
    >NC_058072.1 Acropora millepora isolate JS-1 chromosome 7, Amil_v2.1, whole genome shotgun sequence
    >NC_058073.1 Acropora millepora isolate JS-1 chromosome 8, Amil_v2.1, whole genome shotgun sequence
    >NC_058074.1 Acropora millepora isolate JS-1 chromosome 9, Amil_v2.1, whole genome shotgun sequence
    >NC_058075.1 Acropora millepora isolate JS-1 chromosome 10, Amil_v2.1, whole genome shotgun sequence

    --------------------------------------------------

    >NC_058066.1_Acropora_millepora_isolate_JS-1_chromosome_1_Amil_v2.1_whole_genome_shotgun_sequence
    >NC_058067.1_Acropora_millepora_isolate_JS-1_chromosome_2_Amil_v2.1_whole_genome_shotgun_sequence
    >NC_058068.1_Acropora_millepora_isolate_JS-1_chromosome_3_Amil_v2.1_whole_genome_shotgun_sequence
    >NC_058069.1_Acropora_millepora_isolate_JS-1_chromosome_4_Amil_v2.1_whole_genome_shotgun_sequence
    >NC_058070.1_Acropora_millepora_isolate_JS-1_chromosome_5_Amil_v2.1_whole_genome_shotgun_sequence
    >NC_058071.1_Acropora_millepora_isolate_JS-1_chromosome_6_Amil_v2.1_whole_genome_shotgun_sequence
    >NC_058072.1_Acropora_millepora_isolate_JS-1_chromosome_7_Amil_v2.1_whole_genome_shotgun_sequence
    >NC_058073.1_Acropora_millepora_isolate_JS-1_chromosome_8_Amil_v2.1_whole_genome_shotgun_sequence
    >NC_058074.1_Acropora_millepora_isolate_JS-1_chromosome_9_Amil_v2.1_whole_genome_shotgun_sequence
    >NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence

# 4 Reformat miRBase FastA description lines

``` bash
# Load bash variables into memory
source .bashvars


sed '/^>/ s/ /_/g' "${deep_dive_data_dir}/${mirbase_mature_fasta_name}" \
| sed '/^>/ s/,//g' \
> "${deep_dive_data_dir}/${mirbase_mature_fasta_no_spaces}"

grep "^>" ${deep_dive_data_dir}/${mirbase_mature_fasta_name} \
| head

echo ""
echo "--------------------------------------------------"
echo ""

grep "^>" ${deep_dive_data_dir}/${mirbase_mature_fasta_no_spaces} \
| head
```

    >cel-let-7-5p MIMAT0000001 Caenorhabditis elegans let-7-5p
    >cel-let-7-3p MIMAT0015091 Caenorhabditis elegans let-7-3p
    >cel-lin-4-5p MIMAT0000002 Caenorhabditis elegans lin-4-5p
    >cel-lin-4-3p MIMAT0015092 Caenorhabditis elegans lin-4-3p
    >cel-miR-1-5p MIMAT0020301 Caenorhabditis elegans miR-1-5p
    >cel-miR-1-3p MIMAT0000003 Caenorhabditis elegans miR-1-3p
    >cel-miR-2-5p MIMAT0020302 Caenorhabditis elegans miR-2-5p
    >cel-miR-2-3p MIMAT0000004 Caenorhabditis elegans miR-2-3p
    >cel-miR-34-5p MIMAT0000005 Caenorhabditis elegans miR-34-5p
    >cel-miR-34-3p MIMAT0015093 Caenorhabditis elegans miR-34-3p

    --------------------------------------------------

    >cel-let-7-5p_MIMAT0000001_Caenorhabditis_elegans_let-7-5p
    >cel-let-7-3p_MIMAT0015091_Caenorhabditis_elegans_let-7-3p
    >cel-lin-4-5p_MIMAT0000002_Caenorhabditis_elegans_lin-4-5p
    >cel-lin-4-3p_MIMAT0015092_Caenorhabditis_elegans_lin-4-3p
    >cel-miR-1-5p_MIMAT0020301_Caenorhabditis_elegans_miR-1-5p
    >cel-miR-1-3p_MIMAT0000003_Caenorhabditis_elegans_miR-1-3p
    >cel-miR-2-5p_MIMAT0020302_Caenorhabditis_elegans_miR-2-5p
    >cel-miR-2-3p_MIMAT0000004_Caenorhabditis_elegans_miR-2-3p
    >cel-miR-34-5p_MIMAT0000005_Caenorhabditis_elegans_miR-34-5p
    >cel-miR-34-3p_MIMAT0015093_Caenorhabditis_elegans_miR-34-3p

# 5 Bowtie v1 genome index

miRDeep2 requires a Bowtie v1 genome index - cannot use Bowtie2 index

``` bash
# Load bash variables into memory
source .bashvars

# Check for existence of genome index first
if [ ! -f "${genome_fasta_no_spaces%.*}.ebwt" ]; then
  ${bowtie_build} \
  ${genome_fasta_dir}/${genome_fasta_no_spaces} \
  ${genome_fasta_dir}/${genome_fasta_no_spaces%.*} \
  --threads ${threads} \
  --quiet
fi
```

# 6 Map reads to genome

Requires genome to be previously indexed with Bowtie.

Additionally, requires user to enter path to their `mirdeep2` directory
as well as their `perl5` installation.

``` bash
# Load bash variables into memory
source .bashvars

# Append miRDeep2 to system PATH and set PERL5LIB
export PATH=$PATH:/home/shared/mirdeep2/bin
export PERL5LIB=$PERL5LIB:/home/shared/mirdeep2/lib/perl5

# Run miRDeep2 mapping
time \
${mirdeep2_mapper} \
${output_dir_top}/${collapsed_reads_mirdeep2} \
-c \
-p ${genome_fasta_dir}/${genome_fasta_no_spaces%.*} \
-t ${output_dir_top}/${mirdeep2_mapping_file} \
-o ${threads}
```

# 7 Run miRDeep2

Recommendation is to use the closest related species in the miRDeep2
options, even if the species isn’t very closely related. The
documentation indicates that miRDeep2 is always more accurate when at
least a species is provided.

The options provided to the command are as follows:

- `none`: Known miRNAs of the species being analyzed.
- `-t S.pupuratus`: Related species.
- `none`: Known miRNA precursors in this species.
- `-P`: Specifies miRBase version \> 18.
- `-v`: Remove temporary files after completion.
- `-g -1`: Number of precursors to anlayze. A setting of `-1` will
  analyze all. Default is 50,000 I set this to `-1` after multiple
  attempts to run using the default kept failing.

NOTE: This will take an *extremely* long time to run (days). Could
possible by shortened by excluding `randfold` analysis.

``` bash
# Load bash variables into memory
source .bashvars

# Append miRDeep2 to system PATH and set PERL5LIB
export PATH=$PATH:/home/shared/mirdeep2/bin
export PERL5LIB=$PERL5LIB:/home/shared/mirdeep2/lib/perl5

time \
${mirdeep2} \
${output_dir_top}/${collapsed_reads_mirdeep2} \
${genome_fasta_dir}/${genome_fasta_no_spaces} \
${output_dir_top}/${mirdeep2_mapping_file} \
none \
${deep_dive_data_dir}/${mirbase_mature_fasta_no_spaces} \
none \
-t S.purpuratus \
-P \
-v \
-g -1 \
2>${output_dir_top}/miRDeep2-S.purpuratus-report.log
```

## 7.1 Check runtime

``` bash
# Load bash variables into memory
source .bashvars

tail -n 6 ${output_dir_top}/miRDeep2-S.purpuratus-report.log
```

    miRDeep runtime: 

    started: 14:47:36
    ended: 21:33:35
    total:78h:45m:59s

# 8 Move output files to output directory

MiRDeep2 outputs all files to the current working directly with no way
to redirect so want to move to intended output directory.

## 8.1 Move output files

Output files will be in the format of `result_*` and `error_*`

``` bash
# Load bash variables into memory
source .bashvars

for file in result_* error_*
do
  mv "${file}" "${output_dir_top}/"
done
```

## 8.2 Identify directories

``` bash
# Load bash variables into memory
source .bashvars

ls -l | grep "^d"
```

    drwxrwxr-x 3 sam sam    4096 Oct 30 14:00 07-Apul-lncRNA-dist_files
    drwxr-xr-x 4 sam sam    4096 Nov  7 14:29 10-Apul-sRNAseq-BLASTn_cache
    drwxr-xr-x 3 sam sam    4096 Nov  2 13:10 12-Apul-sRNAseq-MirMachine_cache
    drwxr-xr-x 4 sam sam    4096 Nov  5 17:39 13-Apul-sRNAseq-ShortStack_cache
    drwxrwxr-x 3 sam sam    4096 May  5  2023 rsconnect

## 8.3 Rsync directories

``` bash
# Load bash variables into memory
source .bashvars

rsync -aP . --include='dir***' --exclude='*' --quiet "${output_dir_top}/"

rsync -aP . --include='map***' --exclude='*' --quiet "${output_dir_top}/"

rsync -aP . --exclude='mirgene*' --include='mir***' --exclude='*' --quiet "${output_dir_top}/"

rsync -aP . --include='pdfs***' --exclude='*' --quiet "${output_dir_top}/"

echo ""
echo "Check new location:"
echo ""

ls -l "${output_dir_top}/" | grep "^d"
```

## 8.4 Confirm deletion patterns work *before* deletion!

FYI - `eval=FALSE` is set because the following command will only work
once…

``` bash
ls --directory dir_* mapper_logs mirdeep_runs mirna_results* pdfs_*
```

## 8.5 Remove directores from code directory

``` bash
# Load bash variables into memory
source .bashvars

rm -rf dir_* mapper_logs mirdeep_runs mirna_results* pdfs_*
```

## 8.6 Check to make sure they’re gone

``` bash
ls --directory dir_* mapper_logs mirdeep_runs mirna_results* pdfs_*
```

    ls: cannot access 'dir_*': No such file or directory
    ls: cannot access 'mapper_logs': No such file or directory
    ls: cannot access 'mirdeep_runs': No such file or directory
    ls: cannot access 'mirna_results*': No such file or directory
    ls: cannot access 'pdfs_*': No such file or directory

# 9 Results

## 9.1 Peep output file format

The formatting of this CSV is terrible. It has a 3-column table on top
of a 17-column table. This makes parsing a bit of a pain in its raw
format.

``` bash
# Load bash variables into memory
source .bashvars


head -n 30 "${output_dir_top}/result_08_11_2023_t_14_47_36.csv"
```

    miRDeep2 score  estimated signal-to-noise   excision gearing
    10  1.4 1
    9   1.4 1
    8   1.4 1
    7   1.4 1
    6   1.4 1
    5   1.4 1
    4   1.4 1
    3   1.5 1
    2   1.3 1
    1   1.2 1
    0   1.2 1
    -1  1.1 1
    -2  1   1
    -3  1   1
    -4  0.9 1
    -5  0.9 1
    -6  1   1
    -7  1   1
    -8  1   1
    -9  1   1
    -10 1   1



    novel miRNAs predicted by miRDeep2
    provisional id  miRDeep2 score  estimated probability that the miRNA candidate is a true positive   rfam alert  total read count    mature read count   loop read count star read count significant randfold p-value    miRBase miRNA   example miRBase miRNA with the same seed    UCSC browser    NCBI blastn consensus mature sequence   consensus star sequence consensus precursor sequence    precursor coordinate
    NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence_1178449   1732928.9       -   3399069 2519474 7036    872559  no  -   -   -   -   gagguccggacuuggggaggguuau   caccccuccccgcacauugaccucu   gagguccggacuuggggaggguuauccacucuuaaauucaguugauuaccaccccuccccgcacauugaccucu  NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence:2422830..2422904:+
    NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence_1234230   1723153     -   3379895 2515531 105 864259  no  -   -   -   -   gagguccggacuuggggaggguuau   caccccuccccgcacauugaccucu   gagguccggacuuggggaggguuaucacccuugaauucaguugauuaccaccccuccccgcacauugaccucu   NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence:4009032..4009105:-
    NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence_1231100   1722860.2       -   3379320 2515059 2   864259  no  -   -   -   -   gagguccggacuuggggaggguuau   caccccuccccgcacauugaccucu   gagguccggacuuggggaggguuaucacauuugaauuugauuaccaccccuccccgcacauugaccucu   NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence:2420994..2421063:-

## 9.2 Create more easily parasable results file

This will match the line beginning with `provisional id` and print to
the end of the file (represented by the `$p`. `$` = end, `p` = print)

``` bash
# Load bash variables into memory
source .bashvars


sed --quiet '/provisional id/,$p' "${output_dir_top}/result_08_11_2023_t_14_47_36.csv" \
> "${output_dir_top}/parsable-result_08_11_2023_t_14_47_36.csv"

head "${output_dir_top}/parsable-result_08_11_2023_t_14_47_36.csv"
```

    provisional id  miRDeep2 score  estimated probability that the miRNA candidate is a true positive   rfam alert  total read count    mature read count   loop read count star read count significant randfold p-value    miRBase miRNA   example miRBase miRNA with the same seed    UCSC browser    NCBI blastn consensus mature sequence   consensus star sequence consensus precursor sequence    precursor coordinate
    NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence_1178449   1732928.9       -   3399069 2519474 7036    872559  no  -   -   -   -   gagguccggacuuggggaggguuau   caccccuccccgcacauugaccucu   gagguccggacuuggggaggguuauccacucuuaaauucaguugauuaccaccccuccccgcacauugaccucu  NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence:2422830..2422904:+
    NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence_1234230   1723153     -   3379895 2515531 105 864259  no  -   -   -   -   gagguccggacuuggggaggguuau   caccccuccccgcacauugaccucu   gagguccggacuuggggaggguuaucacccuugaauucaguugauuaccaccccuccccgcacauugaccucu   NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence:4009032..4009105:-
    NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence_1231100   1722860.2       -   3379320 2515059 2   864259  no  -   -   -   -   gagguccggacuuggggaggguuau   caccccuccccgcacauugaccucu   gagguccggacuuggggaggguuaucacauuugaauuugauuaccaccccuccccgcacauugaccucu   NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence:2420994..2421063:-
    NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence_1178451   1413566.8       -   2772654 2515531 4332    252791  no  -   -   -   -   gagguccggacuuggggaggguuau   ccaccccuccccgcacauugaccuc   gagguccggacuuggggaggguuaucacccuuaaauucaguugauuaccaccccuccccgcacauugaccuc    NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence:2423277..2423349:+
    NC_058066.1_Acropora_millepora_isolate_JS-1_chromosome_1_Amil_v2.1_whole_genome_shotgun_sequence_157850 127338.1        -   249763  249278  0   485 no  -   nve-miR-9437_MIMAT0035398_Nematostella_vectensis_miR-9437   -   -   uuaacgaguagauaaaugaagaga    cuucguuuauucacucguucaua cuucguuuauucacucguucauauuuauuauuaacgaguagauaaaugaagaga  NC_058066.1_Acropora_millepora_isolate_JS-1_chromosome_1_Amil_v2.1_whole_genome_shotgun_sequence:20346247..20346301:-
    NC_058077.1_Acropora_millepora_isolate_JS-1_chromosome_12_Amil_v2.1_whole_genome_shotgun_sequence_1444410   74034.8     -   145223  145219  1   3   no  -   bta-miR-2424_MIMAT0011993_Bos_taurus_miR-2424   -   -   ugaucuuucgucgcgugguugcucg   uggccaacaaauuugaugaaauaug   uggccaacaaauuugaugaaauaugaauugaugcaaaaccuaccgugaaggcuuuuggugauuucuuucggaugaucuuucgucgcgugguugcucg   NC_058077.1_Acropora_millepora_isolate_JS-1_chromosome_12_Amil_v2.1_whole_genome_shotgun_sequence:7740108..7740205:-
    NC_058079.1_Acropora_millepora_isolate_JS-1_chromosome_14_Amil_v2.1_whole_genome_shotgun_sequence_1622729   59166       -   116075  115529  351 195 no  -   mmu-miR-380-5p_MIMAT0000744_Mus_musculus_miR-380-5p -   -   uugguugaauauaauaacagcgaga   uuguugacaaaacuccacauaauga   uugguugaauauaauaacagcgagaacaguuuuugaaaaugaguccaaauuuaagaugcauuugaaauagucuuguugacaaaacuccacauaauga   NC_058079.1_Acropora_millepora_isolate_JS-1_chromosome_14_Amil_v2.1_whole_genome_shotgun_sequence:11870080..11870177:+
    NC_058079.1_Acropora_millepora_isolate_JS-1_chromosome_14_Amil_v2.1_whole_genome_shotgun_sequence_1622673   59166       -   116075  115529  351 195 no  -   mmu-miR-380-5p_MIMAT0000744_Mus_musculus_miR-380-5p -   -   uugguugaauauaauaacagcgaga   uuguugacaaaacuccacauaauga   uugguugaauauaauaacagcgagaacaguuuuugaaaaugaguccaaauuuaagaugcauuugaaauagucuuguugacaaaacuccacauaauga   NC_058079.1_Acropora_millepora_isolate_JS-1_chromosome_14_Amil_v2.1_whole_genome_shotgun_sequence:11865521..11865618:+
    NC_058079.1_Acropora_millepora_isolate_JS-1_chromosome_14_Amil_v2.1_whole_genome_shotgun_sequence_1622645   59166       -   116075  115529  351 195 no  -   mmu-miR-380-5p_MIMAT0000744_Mus_musculus_miR-380-5p -   -   uugguugaauauaauaacagcgaga   uuguugacaaaacuccacauaauga   uugguugaauauaauaacagcgagaacaguuuuugaaaaugaguccaaauuuaagaugcauuugaaauagucuuguugacaaaacuccacauaauga   NC_058079.1_Acropora_millepora_isolate_JS-1_chromosome_14_Amil_v2.1_whole_genome_shotgun_sequence:11863245..11863342:+

## 9.3 Read in output CSV

This chunk provides a more convise overview of the data and it’s
columns.

``` r
mirdeep_result.df <- read.csv("../output/11-Apul-sRNAseq-miRdeep2/parsable-result_08_11_2023_t_14_47_36.csv",
                              header = TRUE,
                              sep = "\t")

str(mirdeep_result.df)
```

    'data.frame':   4553 obs. of  17 variables:
     $ provisional.id                                                   : chr  "NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence_1178449" "NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence_1234230" "NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence_1231100" "NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence_1178451" ...
     $ miRDeep2.score                                                   : num  1732929 1723153 1722860 1413567 127338 ...
     $ estimated.probability.that.the.miRNA.candidate.is.a.true.positive: logi  NA NA NA NA NA NA ...
     $ rfam.alert                                                       : chr  "-" "-" "-" "-" ...
     $ total.read.count                                                 : int  3399069 3379895 3379320 2772654 249763 145223 116075 116075 116075 116075 ...
     $ mature.read.count                                                : int  2519474 2515531 2515059 2515531 249278 145219 115529 115529 115529 115529 ...
     $ loop.read.count                                                  : int  7036 105 2 4332 0 1 351 351 351 351 ...
     $ star.read.count                                                  : int  872559 864259 864259 252791 485 3 195 195 195 195 ...
     $ significant.randfold.p.value                                     : chr  "no" "no" "no" "no" ...
     $ miRBase.miRNA                                                    : chr  "-" "-" "-" "-" ...
     $ example.miRBase.miRNA.with.the.same.seed                         : chr  "-" "-" "-" "-" ...
     $ UCSC.browser                                                     : chr  "-" "-" "-" "-" ...
     $ NCBI.blastn                                                      : chr  "-" "-" "-" "-" ...
     $ consensus.mature.sequence                                        : chr  "gagguccggacuuggggaggguuau" "gagguccggacuuggggaggguuau" "gagguccggacuuggggaggguuau" "gagguccggacuuggggaggguuau" ...
     $ consensus.star.sequence                                          : chr  "caccccuccccgcacauugaccucu" "caccccuccccgcacauugaccucu" "caccccuccccgcacauugaccucu" "ccaccccuccccgcacauugaccuc" ...
     $ consensus.precursor.sequence                                     : chr  "gagguccggacuuggggaggguuauccacucuuaaauucaguugauuaccaccccuccccgcacauugaccucu" "gagguccggacuuggggaggguuaucacccuugaauucaguugauuaccaccccuccccgcacauugaccucu" "gagguccggacuuggggaggguuaucacauuugaauuugauuaccaccccuccccgcacauugaccucu" "gagguccggacuuggggaggguuaucacccuuaaauucaguugauuaccaccccuccccgcacauugaccuc" ...
     $ precursor.coordinate                                             : chr  "NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence:2422830..2422904:+" "NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence:4009032..4009105:-" "NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence:2420994..2421063:-" "NC_058075.1_Acropora_millepora_isolate_JS-1_chromosome_10_Amil_v2.1_whole_genome_shotgun_sequence:2423277..2423349:+" ...

## 9.4 miRNAs count data

This provides some rudimentary numbers for the miRDeep2 output.

Further analysis is possibly desired to evaluate score thresholds, miRNA
families, etc.

``` bash
# Load bash variables into memory
source .bashvars

# Total predicted miRNAS
total_miRNAs=$(awk 'NR > 1' ${output_dir_top}/parsable-result_08_11_2023_t_14_47_36.csv \
| wc -l
)

echo "Total of predicted miRNAs: ${total_miRNAs}"
echo ""

# Matches to known mature miRNAs
mature_miRNAs=$(awk -F'\t' '$11 != "-" && $11 != "" {print $11}' ${output_dir_top}/parsable-result_08_11_2023_t_14_47_36.csv \
| wc -l
)

echo "Number of seed matches to known miRNAS: ${mature_miRNAs}"
echo ""

# Novel miRNAs
novel_miRNAs=$(awk -F "\t" '$11 == "-" || $11 == "" {print $11}' ${output_dir_top}/parsable-result_08_11_2023_t_14_47_36.csv \
| awk 'NR > 1' \
| wc -l
)

echo "Number of novel miRNAs: ${novel_miRNAs}"
```

    Total of predicted miRNAs: 4553

    Number of seed matches to known miRNAS: 4137

    Number of novel miRNAs: 416

------------------------------------------------------------------------

# 10 Citations

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-friedländer2011" class="csl-entry">

Friedländer, Marc R., Sebastian D. Mackowiak, Na Li, Wei Chen, and
Nikolaus Rajewsky. 2011. “miRDeep2 Accurately Identifies Known and
Hundreds of Novel microRNA Genes in Seven Animal Clades.” *Nucleic Acids
Research* 40 (1): 37–52. <https://doi.org/10.1093/nar/gkr688>.

</div>

</div>
