11-Pmea-sRNAseq-miRdeep2
================
Sam White (modified by K Durkin for P. meandrina analysis)
2023-11-14

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
  - <a href="#85-remove-directories-from-code-directory"
    id="toc-85-remove-directories-from-code-directory">8.5 Remove
    directories from code directory</a>
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
*P.meandrina* sRNAseq reads. The *P.meandrina* genome will be used as
the reference genome. The mature miRBase miRNA database will also be
used.

------------------------------------------------------------------------

Inputs:

- Requires collapsed reads (i.e. concatenated, unique reads) in FastA
  format. See
  [10-Pmea-sRNAseq-BLASTn.Rmd](https://github.com/urol-e5/deep-dive/blob/main/F-Pmea/code/10-Pmea-sRNAseq-BLASTn.Rmd)
  for code.

- *P.meandrina* genome FastA. See
  [12-Pmea-sRNAseq-MirMachine.Rmd](https://github.com/urol-e5/deep-dive/blob/main/F-Pmea/code/12-Pmea-sRNAseq-MirMachine.Rmd)
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

<https://gannet.fish.washington.edu/Atumefaciens/gitrepos/deep-dive/F-Pmea/output/11-Pmea-sRNAseq-miRdeep2/>

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
echo 'export deep_dive_dir=/home/shared/8TB_HDD_02/shedurkin/deep-dive'
echo 'export deep_dive_data_dir="${deep_dive_dir}/data"'
echo 'export output_dir_top=${deep_dive_dir}/F-Pmea/output/11-Pmea-sRNAseq-miRdeep2'
echo 'export genome_fasta_dir=${deep_dive_dir}/F-Pmea/data/Pmea'
echo 'export trimmed_fastqs_dir="${deep_dive_dir}/F-Pmea/output/08-Pmea-sRNAseq-trimming/trimmed-reads"'
echo 'export collapsed_reads_dir="${deep_dive_dir}/F-Pmea/output/10-Pmea-sRNAseq-BLASTn"'
echo ""

echo "# Input/Output files"
echo 'export collapsed_reads_fasta="collapsed-reads-all.fasta"'
echo 'export collapsed_reads_mirdeep2="collapsed-reads-all-mirdeep2.fasta"'
echo 'export concatenated_trimmed_reads_fastq="concatenated-trimmed-reads-all.fastq.gz"'
echo 'export genome_fasta_name="Pocillopora_meandrina_HIv1.assembly.fasta"'
echo 'export genome_fasta_no_spaces="Pocillopora_meandrina_HIv1_nospaces.assembly.fasta"'
echo 'export mirdeep2_mapping_file="Pmea-mirdeep2-mapping.arf"'
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
    export deep_dive_dir=/home/shared/8TB_HDD_02/shedurkin/deep-dive
    export deep_dive_data_dir="${deep_dive_dir}/data"
    export output_dir_top=${deep_dive_dir}/F-Pmea/output/11-Pmea-sRNAseq-miRdeep2
    export genome_fasta_dir=${deep_dive_dir}/F-Pmea/data/Pmea
    export trimmed_fastqs_dir="${deep_dive_dir}/F-Pmea/output/08-Pmea-sRNAseq-trimming/trimmed-reads"
    export collapsed_reads_dir="${deep_dive_dir}/F-Pmea/output/10-Pmea-sRNAseq-BLASTn"

    # Input/Output files
    export collapsed_reads_fasta="collapsed-reads-all.fasta"
    export collapsed_reads_mirdeep2="collapsed-reads-all-mirdeep2.fasta"
    export concatenated_trimmed_reads_fastq="concatenated-trimmed-reads-all.fastq.gz"
    export genome_fasta_name="Pocillopora_meandrina_HIv1.assembly.fasta"
    export genome_fasta_no_spaces="Pocillopora_meandrina_HIv1_nospaces.assembly.fasta"
    export mirdeep2_mapping_file="Pmea-mirdeep2-mapping.arf"
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

    >1-10280803
    >2-10220688
    >3-4878221
    >4-2964293
    >5-2941994
    >6-2745200
    >7-2697913
    >8-2084846
    >9-2015242
    >10-1645687

    --------------------------------------------------

    >seq_1_x10280803
    >seq_2_x10220688
    >seq_3_x4878221
    >seq_4_x2964293
    >seq_5_x2941994
    >seq_6_x2745200
    >seq_7_x2697913
    >seq_8_x2084846
    >seq_9_x2015242
    >seq_10_x1645687

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

    >Pocillopora_meandrina_HIv1___Sc0000000
    >Pocillopora_meandrina_HIv1___Sc0000001
    >Pocillopora_meandrina_HIv1___Sc0000002
    >Pocillopora_meandrina_HIv1___Sc0000003
    >Pocillopora_meandrina_HIv1___Sc0000004
    >Pocillopora_meandrina_HIv1___Sc0000005
    >Pocillopora_meandrina_HIv1___Sc0000006
    >Pocillopora_meandrina_HIv1___Sc0000007
    >Pocillopora_meandrina_HIv1___Sc0000008
    >Pocillopora_meandrina_HIv1___Sc0000009

    --------------------------------------------------

    >Pocillopora_meandrina_HIv1___Sc0000000
    >Pocillopora_meandrina_HIv1___Sc0000001
    >Pocillopora_meandrina_HIv1___Sc0000002
    >Pocillopora_meandrina_HIv1___Sc0000003
    >Pocillopora_meandrina_HIv1___Sc0000004
    >Pocillopora_meandrina_HIv1___Sc0000005
    >Pocillopora_meandrina_HIv1___Sc0000006
    >Pocillopora_meandrina_HIv1___Sc0000007
    >Pocillopora_meandrina_HIv1___Sc0000008
    >Pocillopora_meandrina_HIv1___Sc0000009

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

nohup time \
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
2>${output_dir_top}/miRDeep2-S.purpuratus-report.log \
2>${output_dir_top}/nohup.out \
&
```

## 7.1 Check runtime

``` bash
# Load bash variables into memory
source .bashvars

tail -n 6 ${output_dir_top}/miRDeep2-S.purpuratus-report.log
```

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

    drwxr-xr-x 4 shedurkin labmembers    4096 Nov 14 15:33 10-Pmea-sRNAseq-BLASTn_cache
    drwxr-xr-x 4 shedurkin labmembers    4096 Dec  1 08:51 13-Pmea-sRNAseq-ShortStack_cache
    drwxr-xr-x 4 shedurkin labmembers    4096 Nov 17 08:33 analyses
    drwxr-xr-x 3 shedurkin labmembers    4096 Nov 17 08:58 data

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

## 8.5 Remove directories from code directory

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


head -n 30 "${output_dir_top}/result_30_11_2023_t_09_23_06.csv"
```

    miRDeep2 score  estimated signal-to-noise   excision gearing
    10  1.8 1
    9   1.8 1
    8   1.8 1
    7   1.8 1
    6   1.8 1
    5   2.5 1
    4   2.6 1
    3   2.1 1
    2   1.6 1
    1   1.3 1
    0   1.4 1
    -1  1.5 1
    -2  1.2 1
    -3  0.8 1
    -4  0.6 1
    -5  0.6 1
    -6  0.6 1
    -7  0.6 1
    -8  0.6 1
    -9  0.7 1
    -10 0.7 1



    novel miRNAs predicted by miRDeep2
    provisional id  miRDeep2 score  estimated probability that the miRNA candidate is a true positive   rfam alert  total read count    mature read count   loop read count star read count significant randfold p-value    miRBase miRNA   example miRBase miRNA with the same seed    UCSC browser    NCBI blastn consensus mature sequence   consensus star sequence consensus precursor sequence    precursor coordinate
    Pocillopora_meandrina_HIv1___Sc0000026_2028756  2808293.9       rRNA/tRNA   5508340 5508187 118 35  yes -   cel-miR-359_MIMAT0000701_Caenorhabditis_elegans_miR-359 -   -   gcacuggugguucagugguagaauu   ucgauucccggccagugca gcacuggugguucagugguagaauucucgccugccacgcgggaggcccggguucgauucccggccagugca Pocillopora_meandrina_HIv1___Sc0000026:5391695..5391766:+
    Pocillopora_meandrina_HIv1___Sc0000000_151108   266720.9        -   523157  521464  0   1693    yes -   mmu-miR-710_MIMAT0003500_Mus_musculus_miR-710   -   -   ucaagucuaggcugguuaguuu  cuaaaccagacuaggcuucagc  cuaaaccagacuaggcuucagcauauuuauuuugucaagucuaggcugguuaguuu    Pocillopora_meandrina_HIv1___Sc0000000:20372434..20372490:-
    Pocillopora_meandrina_HIv1___Sc0000000_72708    265125      -   520035  519878  0   157 no  -   prd-miR-7961-5p_MIMAT0030820_Panagrellus_redivivus_miR-7961-5p  -   -   aaacuaaccagccuagacuuga  gaagccuagucugguuuag aaacuaaccagccuagacuugacaaaauaaauaugcugaagccuagucugguuuag    Pocillopora_meandrina_HIv1___Sc0000000:20372434..20372490:+

## 9.2 Create more easily parasable results file

This will match the line beginning with `provisional id` and print to
the end of the file (represented by the `$p`. `$` = end, `p` = print)

``` bash
# Load bash variables into memory
source .bashvars


sed --quiet '/provisional id/,$p' "${output_dir_top}/result_30_11_2023_t_09_23_06.csv" \
> "${output_dir_top}/parsable-result_30_11_2023_t_09_23_06.csv"

head "${output_dir_top}/parsable-result_30_11_2023_t_09_23_06.csv"
```

    provisional id  miRDeep2 score  estimated probability that the miRNA candidate is a true positive   rfam alert  total read count    mature read count   loop read count star read count significant randfold p-value    miRBase miRNA   example miRBase miRNA with the same seed    UCSC browser    NCBI blastn consensus mature sequence   consensus star sequence consensus precursor sequence    precursor coordinate
    Pocillopora_meandrina_HIv1___Sc0000026_2028756  2808293.9       rRNA/tRNA   5508340 5508187 118 35  yes -   cel-miR-359_MIMAT0000701_Caenorhabditis_elegans_miR-359 -   -   gcacuggugguucagugguagaauu   ucgauucccggccagugca gcacuggugguucagugguagaauucucgccugccacgcgggaggcccggguucgauucccggccagugca Pocillopora_meandrina_HIv1___Sc0000026:5391695..5391766:+
    Pocillopora_meandrina_HIv1___Sc0000000_151108   266720.9        -   523157  521464  0   1693    yes -   mmu-miR-710_MIMAT0003500_Mus_musculus_miR-710   -   -   ucaagucuaggcugguuaguuu  cuaaaccagacuaggcuucagc  cuaaaccagacuaggcuucagcauauuuauuuugucaagucuaggcugguuaguuu    Pocillopora_meandrina_HIv1___Sc0000000:20372434..20372490:-
    Pocillopora_meandrina_HIv1___Sc0000000_72708    265125      -   520035  519878  0   157 no  -   prd-miR-7961-5p_MIMAT0030820_Panagrellus_redivivus_miR-7961-5p  -   -   aaacuaaccagccuagacuuga  gaagccuagucugguuuag aaacuaaccagccuagacuugacaaaauaaauaugcugaagccuagucugguuuag    Pocillopora_meandrina_HIv1___Sc0000000:20372434..20372490:+
    Pocillopora_meandrina_HIv1___Sc0000000_72707    265036.4        -   519879  519878  0   1   no  -   prd-miR-7961-5p_MIMAT0030820_Panagrellus_redivivus_miR-7961-5p  -   -   aaacuaaccagccuagacuuga  augucuaguggauugaagugugaag   augucuaguggauugaagugugaagugagaaaacuaaccagccuagacuuga    Pocillopora_meandrina_HIv1___Sc0000000:20372404..20372456:+
    Pocillopora_meandrina_HIv1___Sc0000017_1663179  67005       -   131419  81165   0   50254   yes -   tca-miR-6008-5p_MIMAT0023557_Tribolium_castaneum_miR-6008-5p    -   -   aacugcugagauucuauggauuu auccauagaacuucugcauuuga aacugcugagauucuauggauuuaauuuaaauccauagaacuucugcauuuga   Pocillopora_meandrina_HIv1___Sc0000017:5050927..5050980:-
    Pocillopora_meandrina_HIv1___Sc0000005_690479   66008.5     -   129464  128045  0   1419    yes -   nve-miR-2023-3p_MIMAT0009756_Nematostella_vectensis_miR-2023-3p -   -   aaagaaguacaagugguaggg   cugccacucguauuuucuuuca  cugccacucguauuuucuuucacguuuaucgaugaaagaaguacaagugguaggg Pocillopora_meandrina_HIv1___Sc0000005:601591..601646:+
    Pocillopora_meandrina_HIv1___Sc0000002_415752   54244       -   106389  106109  0   280 yes -   tca-miR-3860-3p_MIMAT0018793_Tribolium_castaneum_miR-3860-3p    -   -   uuguguaacucccuaaggaagg  cucuuucgggugucacacaacg  cucuuucgggugucacacaacgucgucaaggagcguuguguaacucccuaaggaagg   Pocillopora_meandrina_HIv1___Sc0000002:3841965..3842022:-
    Pocillopora_meandrina_HIv1___Sc0000002_342336   53884       -   105682  105663  0   19  yes -   zma-miR482-3p_MIMAT0014027_Zea_mays_miR482-3p   -   -   ccuuccuuagggaguuacacaa  ugugugacacccgaaagag ccuuccuuagggaguuacacaacgcuccuugacgacguugugugacacccgaaagag   Pocillopora_meandrina_HIv1___Sc0000002:3841965..3842022:+
    Pocillopora_meandrina_HIv1___Sc0000016_1576340  45696.9     -   89625   84067   0   5558    yes -   rlcv-miR-rL1-29-5p_MIMAT0019194_Rhesus_lymphocryptovirus_miR-rL1-29-5p  -   -   ucagucccaccaucucaccaau  ggugagcuguuuggacuuaua   ggugagcuguuuggacuuauauuauugguaucagucccaccaucucaccaau    Pocillopora_meandrina_HIv1___Sc0000016:7550614..7550666:+

## 9.3 Read in output CSV

This chunk provides a more convise overview of the data and it’s
columns.

``` r
mirdeep_result.df <- read.csv("../output/11-Pmea-sRNAseq-miRdeep2/parsable-result_30_11_2023_t_09_23_06.csv",
                              header = TRUE,
                              sep = "\t")

str(mirdeep_result.df)
```

    'data.frame':   2429 obs. of  17 variables:
     $ provisional.id                                                   : chr  "Pocillopora_meandrina_HIv1___Sc0000026_2028756" "Pocillopora_meandrina_HIv1___Sc0000000_151108" "Pocillopora_meandrina_HIv1___Sc0000000_72708" "Pocillopora_meandrina_HIv1___Sc0000000_72707" ...
     $ miRDeep2.score                                                   : num  2808294 266721 265125 265036 67005 ...
     $ estimated.probability.that.the.miRNA.candidate.is.a.true.positive: logi  NA NA NA NA NA NA ...
     $ rfam.alert                                                       : chr  "rRNA/tRNA" "-" "-" "-" ...
     $ total.read.count                                                 : int  5508340 523157 520035 519879 131419 129464 106389 105682 89625 89424 ...
     $ mature.read.count                                                : int  5508187 521464 519878 519878 81165 128045 106109 105663 84067 84067 ...
     $ loop.read.count                                                  : int  118 0 0 0 0 0 0 0 0 0 ...
     $ star.read.count                                                  : int  35 1693 157 1 50254 1419 280 19 5558 5357 ...
     $ significant.randfold.p.value                                     : chr  "yes" "yes" "no" "no" ...
     $ miRBase.miRNA                                                    : chr  "-" "-" "-" "-" ...
     $ example.miRBase.miRNA.with.the.same.seed                         : chr  "cel-miR-359_MIMAT0000701_Caenorhabditis_elegans_miR-359" "mmu-miR-710_MIMAT0003500_Mus_musculus_miR-710" "prd-miR-7961-5p_MIMAT0030820_Panagrellus_redivivus_miR-7961-5p" "prd-miR-7961-5p_MIMAT0030820_Panagrellus_redivivus_miR-7961-5p" ...
     $ UCSC.browser                                                     : chr  "-" "-" "-" "-" ...
     $ NCBI.blastn                                                      : chr  "-" "-" "-" "-" ...
     $ consensus.mature.sequence                                        : chr  "gcacuggugguucagugguagaauu" "ucaagucuaggcugguuaguuu" "aaacuaaccagccuagacuuga" "aaacuaaccagccuagacuuga" ...
     $ consensus.star.sequence                                          : chr  "ucgauucccggccagugca" "cuaaaccagacuaggcuucagc" "gaagccuagucugguuuag" "augucuaguggauugaagugugaag" ...
     $ consensus.precursor.sequence                                     : chr  "gcacuggugguucagugguagaauucucgccugccacgcgggaggcccggguucgauucccggccagugca" "cuaaaccagacuaggcuucagcauauuuauuuugucaagucuaggcugguuaguuu" "aaacuaaccagccuagacuugacaaaauaaauaugcugaagccuagucugguuuag" "augucuaguggauugaagugugaagugagaaaacuaaccagccuagacuuga" ...
     $ precursor.coordinate                                             : chr  "Pocillopora_meandrina_HIv1___Sc0000026:5391695..5391766:+" "Pocillopora_meandrina_HIv1___Sc0000000:20372434..20372490:-" "Pocillopora_meandrina_HIv1___Sc0000000:20372434..20372490:+" "Pocillopora_meandrina_HIv1___Sc0000000:20372404..20372456:+" ...

## 9.4 miRNAs count data

This provides some rudimentary numbers for the miRDeep2 output.

Further analysis is possibly desired to evaluate score thresholds, miRNA
families, etc.

``` bash
# Load bash variables into memory
source .bashvars

# Total predicted miRNAS
total_miRNAs=$(awk 'NR > 1' ${output_dir_top}/parsable-result_30_11_2023_t_09_23_06.csv \
| wc -l
)

echo "Total of predicted miRNAs: ${total_miRNAs}"
echo ""

# Matches to known mature miRNAs
mature_miRNAs=$(awk -F'\t' '$11 != "-" && $11 != "" {print $11}' ${output_dir_top}/parsable-result_30_11_2023_t_09_23_06.csv \
| wc -l
)

echo "Number of seed matches to known miRNAS: ${mature_miRNAs}"
echo ""

# Novel miRNAs
novel_miRNAs=$(awk -F "\t" '$11 == "-" || $11 == "" {print $11}' ${output_dir_top}/parsable-result_30_11_2023_t_09_23_06.csv \
| awk 'NR > 1' \
| wc -l
)

echo "Number of novel miRNAs: ${novel_miRNAs}"
```

    Total of predicted miRNAs: 2429

    Number of seed matches to known miRNAS: 2143

    Number of novel miRNAs: 286

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
