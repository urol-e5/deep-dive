11-Peve-sRNAseq-miRdeep2
================
Sam White
2023-12-01

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
*P.evermanni* sRNAseq reads.

------------------------------------------------------------------------

Inputs:

- Requires collapsed reads (i.e. concatenated, unique reads) in FastA
  format. See
  [10-Peve-sRNAseq-BLASTn.Rmd](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/code/10-Peve-sRNAseq-BLASTn.Rmd)
  for code.

- Genome FastA. See
  [07-Peve-sRNAseq-MirMachine.Rmd](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/code/07-Peve-sRNAseq-MirMachine.Rmd)
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

<https://gannet.fish.washington.edu/Atumefaciens/gitrepos/deep-dive/E-Peve/output/11-Peve-sRNAseq-miRdeep2/>

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
echo 'export output_dir_top=${deep_dive_dir}/E-Peve/output/11-Peve-sRNAseq-miRdeep2'
echo 'export genome_fasta_dir=${deep_dive_dir}/E-Peve/data'
echo 'export trimmed_fastqs_dir="${deep_dive_dir}/E-Peve/output/06-Peve-sRNAseq-trimming/trimmed-reads"'
echo 'export collapsed_reads_dir="${deep_dive_dir}/E-Peve/output/10-Peve-sRNAseq-BLASTn"'
echo ""

echo "# Input/Output files"
echo 'export collapsed_reads_fasta="collapsed-reads-all.fasta"'
echo 'export collapsed_reads_mirdeep2="collapsed-reads-all-mirdeep2.fasta"'
echo 'export concatenated_trimmed_reads_fastq="concatenated-trimmed-reads-all.fastq.gz"'
echo 'export genome_fasta_name="Porites_evermanni_v1.fa"'
echo 'export genome_fasta_no_spaces="Porites_evermanni_v1_genomic-no_spaces.fna"'
echo 'export mirdeep2_mapping_file="Peve-mirdeep2-mapping.arf"'
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
    export output_dir_top=${deep_dive_dir}/E-Peve/output/11-Peve-sRNAseq-miRdeep2
    export genome_fasta_dir=${deep_dive_dir}/E-Peve/data
    export trimmed_fastqs_dir="${deep_dive_dir}/E-Peve/output/06-Peve-sRNAseq-trimming/trimmed-reads"
    export collapsed_reads_dir="${deep_dive_dir}/E-Peve/output/10-Peve-sRNAseq-BLASTn"

    # Input/Output files
    export collapsed_reads_fasta="collapsed-reads-all.fasta"
    export collapsed_reads_mirdeep2="collapsed-reads-all-mirdeep2.fasta"
    export concatenated_trimmed_reads_fastq="concatenated-trimmed-reads-all.fastq.gz"
    export genome_fasta_name="Porites_evermanni_v1.fa"
    export genome_fasta_no_spaces="Porites_evermanni_v1_genomic-no_spaces.fna"
    export mirdeep2_mapping_file="Peve-mirdeep2-mapping.arf"
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

    >1-2463526
    >2-2045325
    >3-2026343
    >4-1419489
    >5-1161968
    >6-1144435
    >7-1040735
    >8-942821
    >9-658343
    >10-653354

    --------------------------------------------------

    >seq_1_x2463526
    >seq_2_x2045325
    >seq_3_x2026343
    >seq_4_x1419489
    >seq_5_x1161968
    >seq_6_x1144435
    >seq_7_x1040735
    >seq_8_x942821
    >seq_9_x658343
    >seq_10_x653354

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

    >Porites_evermani_scaffold_1
    >Porites_evermani_scaffold_2
    >Porites_evermani_scaffold_3
    >Porites_evermani_scaffold_4
    >Porites_evermani_scaffold_5
    >Porites_evermani_scaffold_6
    >Porites_evermani_scaffold_7
    >Porites_evermani_scaffold_8
    >Porites_evermani_scaffold_9
    >Porites_evermani_scaffold_10

    --------------------------------------------------

    >Porites_evermani_scaffold_1
    >Porites_evermani_scaffold_2
    >Porites_evermani_scaffold_3
    >Porites_evermani_scaffold_4
    >Porites_evermani_scaffold_5
    >Porites_evermani_scaffold_6
    >Porites_evermani_scaffold_7
    >Porites_evermani_scaffold_8
    >Porites_evermani_scaffold_9
    >Porites_evermani_scaffold_10

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

    started: 10:13:39
    ended: 6:24:29
    total:68h:10m:50s

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

    drwxrwxr-x 3 sam sam    4096 Oct 30 14:00 03-Peve-lncRNA-dist_files
    drwxr-xr-x 4 sam sam    4096 Nov 15 14:49 06-Peve-sRNAseq-trimming_cache
    drwxr-xr-x 4 sam sam    4096 Nov 15 17:51 07-Peve-sRNAseq-MirMachine_cache
    drwxr-xr-x 4 sam sam    4096 Nov 15 21:27 08-Peve-sRNAseq-ShortStack_cache
    drwxr-xr-x 4 sam sam    4096 Nov 16 10:46 10-Peve-sRNAseq-BLASTn_cache
    drwxr-xr-x 4 sam sam    4096 Nov 15 15:16 analyses
    drwxr-xr-x 4 sam sam    4096 Nov 15 15:16 data

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


head -n 30 "${output_dir_top}/result_01_12_2023_t_10_13_39.csv"
```

    miRDeep2 score  estimated signal-to-noise   excision gearing
    10  1.9 1
    9   1.9 1
    8   1.9 1
    7   1.8 1
    6   2   1
    5   3.1 1
    4   3   1
    3   2.2 1
    2   1.6 1
    1   1.3 1
    0   1.3 1
    -1  1.2 1
    -2  0.9 1
    -3  0.7 1
    -4  0.7 1
    -5  0.7 1
    -6  0.7 1
    -7  0.7 1
    -8  0.8 1
    -9  0.8 1
    -10 0.8 1



    novel miRNAs predicted by miRDeep2
    provisional id  miRDeep2 score  estimated probability that the miRNA candidate is a true positive   rfam alert  total read count    mature read count   loop read count star read count significant randfold p-value    miRBase miRNA   example miRBase miRNA with the same seed    UCSC browser    NCBI blastn consensus mature sequence   consensus star sequence consensus precursor sequence    precursor coordinate
    Porites_evermani_scaffold_875_1053054   136688.2        -   268112  268092  0   20  no  -   gga-miR-1629_MIMAT0007500_Gallus_gallus_miR-1629    -   -   agcugucgacuccugcaccaagaag   uauguuguaagcgcagcucagaugu   uauguuguaagcgcagcucagauguugaauuaagcugucgacuccugcaccaagaag   Porites_evermani_scaffold_875:120669..120726:+
    Porites_evermani_scaffold_1503_1364289  110310.8        -   216372  216235  0   137 yes -   prd-miR-7961-5p_MIMAT0030820_Panagrellus_redivivus_miR-7961-5p  -   -   aaacuaaccagccuagaccuga  aaggacuacacugguaua  aaacuaaccagccuagaccugacaagaaacaugcaaggacuacacugguaua    Porites_evermani_scaffold_1503:47771..47823:-
    Porites_evermani_scaffold_1503_1364283  110255.5        -   216257  216235  0   22  yes -   prd-miR-7961-5p_MIMAT0030820_Panagrellus_redivivus_miR-7961-5p  -   -   aaacuaaccagccuagaccuga  aggacuagacugguauag  aaacuaaccagccuagaccugacaagaaacaugccaggacuagacugguauag   Porites_evermani_scaffold_1503:46868..46921:-

## 9.2 Create more easily parasable results file

This will match the line beginning with `provisional id` and print to
the end of the file (represented by the `$p`. `$` = end, `p` = print)

``` bash
# Load bash variables into memory
source .bashvars


sed --quiet '/provisional id/,$p' "${output_dir_top}/result_01_12_2023_t_10_13_39.csv" \
> "${output_dir_top}/parsable-result_01_12_2023_t_10_13_39.csv"

head "${output_dir_top}/parsable-result_01_12_2023_t_10_13_39.csv"
```

    provisional id  miRDeep2 score  estimated probability that the miRNA candidate is a true positive   rfam alert  total read count    mature read count   loop read count star read count significant randfold p-value    miRBase miRNA   example miRBase miRNA with the same seed    UCSC browser    NCBI blastn consensus mature sequence   consensus star sequence consensus precursor sequence    precursor coordinate
    Porites_evermani_scaffold_875_1053054   136688.2        -   268112  268092  0   20  no  -   gga-miR-1629_MIMAT0007500_Gallus_gallus_miR-1629    -   -   agcugucgacuccugcaccaagaag   uauguuguaagcgcagcucagaugu   uauguuguaagcgcagcucagauguugaauuaagcugucgacuccugcaccaagaag   Porites_evermani_scaffold_875:120669..120726:+
    Porites_evermani_scaffold_1503_1364289  110310.8        -   216372  216235  0   137 yes -   prd-miR-7961-5p_MIMAT0030820_Panagrellus_redivivus_miR-7961-5p  -   -   aaacuaaccagccuagaccuga  aaggacuacacugguaua  aaacuaaccagccuagaccugacaagaaacaugcaaggacuacacugguaua    Porites_evermani_scaffold_1503:47771..47823:-
    Porites_evermani_scaffold_1503_1364283  110255.5        -   216257  216235  0   22  yes -   prd-miR-7961-5p_MIMAT0030820_Panagrellus_redivivus_miR-7961-5p  -   -   aaacuaaccagccuagaccuga  aggacuagacugguauag  aaacuaaccagccuagaccugacaagaaacaugccaggacuagacugguauag   Porites_evermani_scaffold_1503:46868..46921:-
    Porites_evermani_scaffold_72_223878 61997       -   121607  116699  1763    3145    yes -   -   -   -   gagguccggacgguugaggguuauc   caccccucauccaccaacuugaccu   gagguccggacgguugaggguuaucaauuuauacuagucugcucaacuggaauuucugaaccaccccucauccaccaacuugaccu  Porites_evermani_scaffold_72:198220..198306:+
    Porites_evermani_scaffold_910_1073823   56249.2     -   110321  84745   0   25576   yes -   hsa-miR-33a-3p_MIMAT0004506_Homo_sapiens_miR-33a-3p -   -   caauguuucggcuuguucccg   gggaacaagccgaaacauuu    caauguuucggcuuguucccguuuucgggaacaagccgaaacauuu  Porites_evermani_scaffold_910:118741..118787:+
    Porites_evermani_scaffold_26_108334 46988.5     -   92157   89552   0   2605    yes -   gga-miR-1467-5p_MIMAT0007345_Gallus_gallus_miR-1467-5p  -   -   ucucagcucaccaaucucugcu  cagggacuggugagcugauguc  cagggacuggugagcugaugucauuuacugaucucagcucaccaaucucugcu   Porites_evermani_scaffold_26:382571..382624:-
    Porites_evermani_scaffold_26_106464 46167.2     -   90546   90411   0   135 yes -   gma-miR4340_MIMAT0018228_Glycine_max_miR4340    -   -   agcagagauuggugagcugaga  caucagcucaccagucccug    agcagagauuggugagcugagaucaguaaaugacaucagcucaccagucccug   Porites_evermani_scaffold_26:382571..382624:+
    Porites_evermani_scaffold_910_1073772   43684       -   85677   84745   0   932 yes -   hsa-miR-33a-3p_MIMAT0004506_Homo_sapiens_miR-33a-3p -   -   caauguuucggcuuguucccg   aaacaaaccgaaacauuu  caauguuucggcuuguucccguuuucggaaacaaaccgaaacauuu  Porites_evermani_scaffold_910:99762..99808:+
    Porites_evermani_scaffold_1503_1364155  27757.3     -   54437   54313   0   124 yes -   mmu-miR-710_MIMAT0003500_Mus_musculus_miR-710   -   -   ucaagucuaggcugguuaguuu  cuacaccaguguagucuuggca  cuacaccaguguagucuuggcaugcuucuugucaagucuaggcugguuaguuu   Porites_evermani_scaffold_1503:47579..47632:+

## 9.3 Read in output CSV

This chunk provides a more convise overview of the data and it’s
columns.

``` r
mirdeep_result.df <- read.csv("../output/11-Peve-sRNAseq-miRdeep2/parsable-result_01_12_2023_t_10_13_39.csv",
                              header = TRUE,
                              sep = "\t")

str(mirdeep_result.df)
```

    'data.frame':   5812 obs. of  17 variables:
     $ provisional.id                                                   : chr  "Porites_evermani_scaffold_875_1053054" "Porites_evermani_scaffold_1503_1364289" "Porites_evermani_scaffold_1503_1364283" "Porites_evermani_scaffold_72_223878" ...
     $ miRDeep2.score                                                   : num  136688 110311 110256 61997 56249 ...
     $ estimated.probability.that.the.miRNA.candidate.is.a.true.positive: logi  NA NA NA NA NA NA ...
     $ rfam.alert                                                       : chr  "-" "-" "-" "-" ...
     $ total.read.count                                                 : int  268112 216372 216257 121607 110321 92157 90546 85677 54437 54249 ...
     $ mature.read.count                                                : int  268092 216235 216235 116699 84745 89552 90411 84745 54313 54225 ...
     $ loop.read.count                                                  : int  0 0 0 1763 0 0 0 0 0 0 ...
     $ star.read.count                                                  : int  20 137 22 3145 25576 2605 135 932 124 24 ...
     $ significant.randfold.p.value                                     : chr  "no" "yes" "yes" "yes" ...
     $ miRBase.miRNA                                                    : chr  "-" "-" "-" "-" ...
     $ example.miRBase.miRNA.with.the.same.seed                         : chr  "gga-miR-1629_MIMAT0007500_Gallus_gallus_miR-1629" "prd-miR-7961-5p_MIMAT0030820_Panagrellus_redivivus_miR-7961-5p" "prd-miR-7961-5p_MIMAT0030820_Panagrellus_redivivus_miR-7961-5p" "-" ...
     $ UCSC.browser                                                     : chr  "-" "-" "-" "-" ...
     $ NCBI.blastn                                                      : chr  "-" "-" "-" "-" ...
     $ consensus.mature.sequence                                        : chr  "agcugucgacuccugcaccaagaag" "aaacuaaccagccuagaccuga" "aaacuaaccagccuagaccuga" "gagguccggacgguugaggguuauc" ...
     $ consensus.star.sequence                                          : chr  "uauguuguaagcgcagcucagaugu" "aaggacuacacugguaua" "aggacuagacugguauag" "caccccucauccaccaacuugaccu" ...
     $ consensus.precursor.sequence                                     : chr  "uauguuguaagcgcagcucagauguugaauuaagcugucgacuccugcaccaagaag" "aaacuaaccagccuagaccugacaagaaacaugcaaggacuacacugguaua" "aaacuaaccagccuagaccugacaagaaacaugccaggacuagacugguauag" "gagguccggacgguugaggguuaucaauuuauacuagucugcucaacuggaauuucugaaccaccccucauccaccaacuugaccu" ...
     $ precursor.coordinate                                             : chr  "Porites_evermani_scaffold_875:120669..120726:+" "Porites_evermani_scaffold_1503:47771..47823:-" "Porites_evermani_scaffold_1503:46868..46921:-" "Porites_evermani_scaffold_72:198220..198306:+" ...

## 9.4 miRNAs count data

This provides some rudimentary numbers for the miRDeep2 output.

Further analysis is possibly desired to evaluate score thresholds, miRNA
families, etc.

``` bash
# Load bash variables into memory
source .bashvars

# Total predicted miRNAS
total_miRNAs=$(awk 'NR > 1' ${output_dir_top}/parsable-result_01_12_2023_t_10_13_39.csv \
| wc -l
)

echo "Total of predicted miRNAs: ${total_miRNAs}"
echo ""

# Matches to known mature miRNAs
mature_miRNAs=$(awk -F'\t' '$11 != "-" && $11 != "" {print $11}' ${output_dir_top}/parsable-result_01_12_2023_t_10_13_39.csv \
| wc -l
)

echo "Number of seed matches to known miRNAS: ${mature_miRNAs}"
echo ""

# Novel miRNAs
novel_miRNAs=$(awk -F "\t" '$11 == "-" || $11 == "" {print $11}' ${output_dir_top}/parsable-result_01_12_2023_t_10_13_39.csv \
| awk 'NR > 1' \
| wc -l
)

echo "Number of novel miRNAs: ${novel_miRNAs}"
```

    Total of predicted miRNAs: 5812

    Number of seed matches to known miRNAS: 5096

    Number of novel miRNAs: 716

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
