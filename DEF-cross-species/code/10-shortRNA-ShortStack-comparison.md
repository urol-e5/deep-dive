10-shortRNA-ShortStack-comparison
================
Kathleen Durkin
2024-05-28

- <a href="#1-prep-data" id="toc-1-prep-data">1 Prep data</a>
  - <a href="#11-isolate-mature-mirna-sequences"
    id="toc-11-isolate-mature-mirna-sequences">1.1 Isolate mature miRNA
    sequences</a>
  - <a href="#12-merge-the-three-mature-mirna-fastas"
    id="toc-12-merge-the-three-mature-mirna-fastas">1.2 Merge the three
    mature miRNA FASTAs</a>
  - <a href="#13-blasts" id="toc-13-blasts">1.3 BLASTs</a>
    - <a href="#131-make-database-for-each-species"
      id="toc-131-make-database-for-each-species">1.3.1 Make database for each
      species:</a>
    - <a href="#132-run-blastn" id="toc-132-run-blastn">1.3.2 Run Blastn</a>
- <a href="#2-join-blast-tables" id="toc-2-join-blast-tables">2 Join BLAST
  tables</a>
- <a href="#3-identify-conserved-mirnas"
  id="toc-3-identify-conserved-mirnas">3 Identify conserved miRNAs</a>
  - <a href="#31-conserved-across-all-three-species-apul-peve-and-pmea"
    id="toc-31-conserved-across-all-three-species-apul-peve-and-pmea">3.1
    Conserved across all three species (Apul, Peve, and Pmea)</a>
  - <a href="#32-conserved-among-subsets-of-the-three-species"
    id="toc-32-conserved-among-subsets-of-the-three-species">3.2 Conserved
    among subsets of the three species</a>
    - <a href="#321-apul-and-peve" id="toc-321-apul-and-peve">3.2.1 Apul and
      Peve</a>
    - <a href="#322-apul-and-pmea" id="toc-322-apul-and-pmea">3.2.2 Apul and
      Pmea</a>
    - <a href="#323-peve-and-pmea" id="toc-323-peve-and-pmea">3.2.3 Peve and
      Pmea</a>

I want to find miRNAs that are conserved among either a subset of or all
three species of interest (*A.pulchra*, *P.evermanni*, and
*P.meandrina*) using Blastn.

# 1 Prep data

## 1.1 Isolate mature miRNA sequences

Our ShortStack output contains sequences for the mature, star, and
precursor sequences for each identified miRNA. We just want to look at
the mature miRNA sequences right now, so let’s isolate those.

``` bash
cd ../data/10-shortRNA-ShortStack-comparison

# Copy all sequences whose headers contain "mature"
awk '/^>/ {p = /mature/} p' ../../../D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/mir.fasta > Apul_ShortStack_mature.fasta

awk '/^>/ {p = /mature/} p' ../../../E-Peve/output/08.2-Peve-sRNAseq-ShortStack-31bp-fastp-merged/ShortStack_out/mir.fasta > Peve_ShortStack_mature.fasta

awk '/^>/ {p = /mature/} p' ../../../F-Pmea/output/13.2.1-Pmea-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/mir.fasta > Pmea_ShortStack_mature.fasta

grep "^>" Apul_ShortStack_mature.fasta | wc -l 
grep "^>" Peve_ShortStack_mature.fasta | wc -l
grep "^>" Pmea_ShortStack_mature.fasta | wc -l
```

    38
    46
    36

## 1.2 Merge the three mature miRNA FASTAs

``` bash
cd ../data/10-shortRNA-ShortStack-comparison

cat Apul_ShortStack_mature.fasta Peve_ShortStack_mature.fasta Pmea_ShortStack_mature.fasta > merged_all_ShortStack_mature.fasta

head merged_all_ShortStack_mature.fasta
tail merged_all_ShortStack_mature.fasta
```

    >Cluster_316.mature::NC_058066.1:12757146-12757168(-)
    TGATCTCTGCAATAGCCTGCCT
    >Cluster_514.mature::NC_058066.1:20088678-20088700(+)
    ACGCTAGGAAGGGATGCCGGGA
    >Cluster_548.mature::NC_058066.1:20346248-20346271(-)
    TTAACGAGTAGATAAATGAAGAG
    >Cluster_1506.mature::NC_058067.1:5656213-5656236(-)
    TTTGCTAGTTGCTTTTGTCCCGT
    >Cluster_1900.mature::NC_058067.1:16118269-16118291(-)
    aaaaatgtcggttgcttaagct
    >Cluster_4838.mature::Pocillopora_meandrina_HIv1___Sc0000018:6855520-6855542(+)
    TCACCCAACAGTTTTAATCTGA
    >Cluster_5273.mature::Pocillopora_meandrina_HIv1___Sc0000021:4351838-4351860(+)
    ACTGATATTCACCAAGTGATTA
    >Cluster_5633.mature::Pocillopora_meandrina_HIv1___Sc0000024:4808687-4808708(+)
    AGAACCCAAGAATCTCGAAGG
    >Cluster_5761.mature::Pocillopora_meandrina_HIv1___Sc0000026:1154771-1154793(-)
    TGTACTATGTTCATGATCTTGC
    >Cluster_6425.mature::Pocillopora_meandrina_HIv1___Sc0000035:1989841-1989863(+)
    TATTTACAACTCTCAAAACAAC

## 1.3 BLASTs

### 1.3.1 Make database for each species:

Apul

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/makeblastdb \
-in ../data/10-shortRNA-ShortStack-comparison/Apul_ShortStack_mature.fasta \
-dbtype nucl \
-out ../output/10-shortRNA-ShortStack-comparison/blasts/Apul-db/Apul_ShortStack_mature
```

    Building a new DB, current time: 05/29/2024 12:34:56
    New DB name:   /home/shared/8TB_HDD_02/shedurkin/deep-dive/DEF-cross-species/output/10-shortRNA-ShortStack-comparison/blasts/Apul-db/Apul_ShortStack_mature
    New DB title:  ../data/10-shortRNA-ShortStack-comparison/Apul_ShortStack_mature.fasta
    Sequence type: Nucleotide
    Deleted existing Nucleotide BLAST database named /home/shared/8TB_HDD_02/shedurkin/deep-dive/DEF-cross-species/output/10-shortRNA-ShortStack-comparison/blasts/Apul-db/Apul_ShortStack_mature
    Keep MBits: T
    Maximum file size: 1000000000B
    Adding sequences from FASTA; added 38 sequences in 0.00743389 seconds.

Peve

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/makeblastdb \
-in ../data/10-shortRNA-ShortStack-comparison/Peve_ShortStack_mature.fasta \
-dbtype nucl \
-out ../output/10-shortRNA-ShortStack-comparison/blasts/Peve-db/Peve_ShortStack_mature
```

    Building a new DB, current time: 05/29/2024 12:34:57
    New DB name:   /home/shared/8TB_HDD_02/shedurkin/deep-dive/DEF-cross-species/output/10-shortRNA-ShortStack-comparison/blasts/Peve-db/Peve_ShortStack_mature
    New DB title:  ../data/10-shortRNA-ShortStack-comparison/Peve_ShortStack_mature.fasta
    Sequence type: Nucleotide
    Deleted existing Nucleotide BLAST database named /home/shared/8TB_HDD_02/shedurkin/deep-dive/DEF-cross-species/output/10-shortRNA-ShortStack-comparison/blasts/Peve-db/Peve_ShortStack_mature
    Keep MBits: T
    Maximum file size: 1000000000B
    Adding sequences from FASTA; added 46 sequences in 0.00186086 seconds.

Pmea

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/makeblastdb \
-in ../data/10-shortRNA-ShortStack-comparison/Pmea_ShortStack_mature.fasta \
-dbtype nucl \
-out ../output/10-shortRNA-ShortStack-comparison/blasts/Pmea-db/Pmea_ShortStack_mature
```

    Building a new DB, current time: 05/29/2024 12:34:58
    New DB name:   /home/shared/8TB_HDD_02/shedurkin/deep-dive/DEF-cross-species/output/10-shortRNA-ShortStack-comparison/blasts/Pmea-db/Pmea_ShortStack_mature
    New DB title:  ../data/10-shortRNA-ShortStack-comparison/Pmea_ShortStack_mature.fasta
    Sequence type: Nucleotide
    Deleted existing Nucleotide BLAST database named /home/shared/8TB_HDD_02/shedurkin/deep-dive/DEF-cross-species/output/10-shortRNA-ShortStack-comparison/blasts/Pmea-db/Pmea_ShortStack_mature
    Keep MBits: T
    Maximum file size: 1000000000B
    Adding sequences from FASTA; added 36 sequences in 0.00192285 seconds.

### 1.3.2 Run Blastn

Apul to all

``` bash

/home/shared/ncbi-blast-2.11.0+/bin/blastn \
-task blastn \
-query ../data/10-shortRNA-ShortStack-comparison/merged_all_ShortStack_mature.fasta \
-db ../output/10-shortRNA-ShortStack-comparison/blasts/Apul-db/Apul_ShortStack_mature \
-out ../output/10-shortRNA-ShortStack-comparison/Apul_to_all_blastn.tab \
-evalue 1E-2 \
-num_threads 40 \
-max_target_seqs 1 \
-max_hsps 1 \
-outfmt 6

wc -l ../output/10-shortRNA-ShortStack-comparison/Apul_to_all_blastn.tab
```

    Warning: [blastn] Examining 5 or more matches is recommended
    48 ../output/10-shortRNA-ShortStack-comparison/Apul_to_all_blastn.tab

Peve to all

``` bash

/home/shared/ncbi-blast-2.11.0+/bin/blastn \
-task blastn \
-query ../data/10-shortRNA-ShortStack-comparison/merged_all_ShortStack_mature.fasta \
-db ../output/10-shortRNA-ShortStack-comparison/blasts/Peve-db/Peve_ShortStack_mature \
-out ../output/10-shortRNA-ShortStack-comparison/Peve_to_all_blastn.tab \
-evalue 1E-2 \
-num_threads 40 \
-max_target_seqs 1 \
-max_hsps 1 \
-outfmt 6

wc -l ../output/10-shortRNA-ShortStack-comparison/Peve_to_all_blastn.tab
```

    Warning: [blastn] Examining 5 or more matches is recommended
    57 ../output/10-shortRNA-ShortStack-comparison/Peve_to_all_blastn.tab

Pmea to all

``` bash

/home/shared/ncbi-blast-2.11.0+/bin/blastn \
-task blastn \
-query ../data/10-shortRNA-ShortStack-comparison/merged_all_ShortStack_mature.fasta \
-db ../output/10-shortRNA-ShortStack-comparison/blasts/Pmea-db/Pmea_ShortStack_mature \
-out ../output/10-shortRNA-ShortStack-comparison/Pmea_to_all_blastn.tab \
-evalue 1E-2 \
-num_threads 40 \
-max_target_seqs 1 \
-max_hsps 1 \
-outfmt 6

wc -l ../output/10-shortRNA-ShortStack-comparison/Pmea_to_all_blastn.tab
```

    Warning: [blastn] Examining 5 or more matches is recommended
    49 ../output/10-shortRNA-ShortStack-comparison/Pmea_to_all_blastn.tab

# 2 Join BLAST tables

``` r
apul_to_all_blastn <- read.table("../output/10-shortRNA-ShortStack-comparison/Apul_to_all_blastn.tab", sep="\t", header=FALSE)
peve_to_all_blastn <- read.table("../output/10-shortRNA-ShortStack-comparison/Peve_to_all_blastn.tab", sep="\t", header=FALSE)
pmea_to_all_blastn <- read.table("../output/10-shortRNA-ShortStack-comparison/Pmea_to_all_blastn.tab", sep="\t", header=FALSE)
```

Column labels: qseqid: Query sequence ID sseqid: Subject (database)
sequence ID pident: Percentage of identical matches length: Alignment
length (number of base pairs or amino acids) mismatch: Number of
mismatches gapopen: Number of gap openings qstart: Start of alignment in
the query qend: End of alignment in the query sstart: Start of alignment
in the subject send: End of alignment in the subject evalue: Expect
value (number of hits expected by chance) bitscore: Bit score

``` r
# Combine the three blast tables
combined_blastn <- rbind(apul_to_all_blastn, peve_to_all_blastn, pmea_to_all_blastn)

# Assign informative column labels
colnames(combined_blastn) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Remove instances where a sequence matched to itself (e.g. from querying an Apul sequence against our combined database which contained all Apul sequences)
corrected_combined_blastn <- filter(combined_blastn, qseqid != sseqid)

# Remove any hits that have mismatches
filtered_combined_blastn <- filter(corrected_combined_blastn, mismatch==0)

# View
nrow(filtered_combined_blastn)
```

    [1] 34

``` r
head(filtered_combined_blastn)
```

                                                                   qseqid
    1                   Cluster_2521.mature::NC_058068.1:597877-597899(+)
    2                Cluster_10726.mature::NC_058079.1:5232178-5232199(+)
    3               Cluster_16040.mature::NW_025322765.1:722296-722318(+)
    4  Cluster_1153.mature::Porites_evermani_scaffold_49:151639-151661(-)
    5 Cluster_5540.mature::Porites_evermani_scaffold_430:205886-205909(-)
    6   Cluster_6211.mature::Porites_evermani_scaffold_502:58996-59018(-)
                                                     sseqid pident length mismatch
    1     Cluster_2522.mature::NC_058068.1:598173-598195(+)    100     22        0
    2  Cluster_10729.mature::NC_058079.1:5261442-5261463(+)    100     21        0
    3 Cluster_16041.mature::NW_025322765.1:799181-799203(-)    100     22        0
    4 Cluster_6977.mature::NC_058073.1:12437102-12437124(+)    100     22        0
    5 Cluster_7025.mature::NC_058073.1:15996878-15996901(-)    100     23        0
    6   Cluster_7077.mature::NC_058074.1:2966522-2966545(+)    100     22        0
      gapopen qstart qend sstart send   evalue bitscore
    1       0      1   22      1   22 3.52e-09     41.0
    2       0      1   21      1   21 1.14e-08     39.2
    3       0      1   22      1   22 3.52e-09     41.0
    4       0      1   22      1   22 3.52e-09     41.0
    5       0      1   23      1   23 1.08e-09     42.8
    6       0      1   22      1   22 3.52e-09     41.0

``` r
write.table(filtered_combined_blastn, "../output/10-shortRNA-ShortStack-comparison/filtered_combined_blast.tab", sep="\t", row.names=FALSE, quote=FALSE)
```

# 3 Identify conserved miRNAs

Ok now we can start identifying conserved miRNAs. Keep in mind that this
list of filtered, combined blastn hits contains duplicates because, for
example, querying Apul sequences against a database containing Peve
sequences is functionally the same as querying those Peve sequences
against a databse which contains Apul. So, for example, this list would
contain a hit matching Apul.seq1 to Peve.seq2, *and* a hit matching
Peve.seq2 to Apul.seq1.

## 3.1 Conserved across all three species (Apul, Peve, and Pmea)

First, lets find miRNAs conserved among all three species. These would
show up as an miRNA from one species that has hits from both other
species (e.g., Apul.seq1 has a hit from Peve *and* a hit from Pmea).

``` r
# Find Apul miRNAs that have matches from both Peve and Pmea
present_in_all <- filtered_combined_blastn %>%
  # isolate Apul miRNAs with hits
  filter(!grepl("Porites_evermani|Pocillopora_meandrina", sseqid)) %>%
  group_by(sseqid) %>%
  filter(any(grepl("Porites_evermani", qseqid)) & any(grepl("Pocillopora_meandrina", qseqid)))

# View the miRNAs that match across all three species
# (recall this will include two entries for each conserved miRNA, it's Apul match in Peve, and its Apul match to Pmea)
head(present_in_all, nrow(present_in_all))
```

    # A tibble: 8 × 12
    # Groups:   sseqid [4]
      qseqid sseqid pident length mismatch gapopen qstart  qend sstart  send  evalue
      <chr>  <chr>   <dbl>  <int>    <int>   <int>  <int> <int>  <int> <int>   <dbl>
    1 Clust… Clust…    100     22        0       0      1    22      1    22 3.52e-9
    2 Clust… Clust…    100     23        0       0      1    23      1    23 1.08e-9
    3 Clust… Clust…    100     21        0       0      1    21      1    21 1.14e-8
    4 Clust… Clust…    100     22        0       0      1    22      1    22 3.52e-9
    5 Clust… Clust…    100     22        0       0      1    22      1    22 3.52e-9
    6 Clust… Clust…    100     21        0       0      2    22      2    22 1.23e-8
    7 Clust… Clust…    100     21        0       0      1    21      1    21 1.14e-8
    8 Clust… Clust…    100     22        0       0      1    22      1    22 3.52e-9
    # ℹ 1 more variable: bitscore <dbl>

``` r
# Count the number of miRNAs conserved across all three species
paste("Number of miRNAs conserved across all three species:", nrow(distinct(present_in_all, sseqid)))
```

    [1] "Number of miRNAs conserved across all three species: 4"

## 3.2 Conserved among subsets of the three species

Now we want to find miRNAs that are conserved withing subsets of the
three species

### 3.2.1 Apul and Peve

Find Apul miRNAs that have hits to Peve miRNAs but *not* hits to Pmea
miRNAs (that would make them conserved among all three species, which
we’ve already identified)

``` r
# Find Apul miRNAs that have matches from only Peve
present_in_apul_peve <- filtered_combined_blastn %>%
  # isolate Apul miRNAs with hits
  filter(!grepl("Porites_evermani|Pocillopora_meandrina", sseqid)) %>%
  group_by(sseqid) %>%
  # filter for hits to Peve only
  filter(any(grepl("Porites_evermani", qseqid)) & !any(grepl("Pocillopora_meandrina", qseqid)))

# View the miRNAs that match between Apul and Peve
head(present_in_apul_peve, nrow(present_in_apul_peve))
```

    # A tibble: 1 × 12
    # Groups:   sseqid [1]
      qseqid sseqid pident length mismatch gapopen qstart  qend sstart  send  evalue
      <chr>  <chr>   <dbl>  <int>    <int>   <int>  <int> <int>  <int> <int>   <dbl>
    1 Clust… Clust…    100     22        0       0      1    22      1    22 3.52e-9
    # ℹ 1 more variable: bitscore <dbl>

``` r
# Count the number of miRNAs conserved across the two species
paste("Number of miRNAs conserved in Apul and Peve:", nrow(distinct(present_in_apul_peve, sseqid)))
```

    [1] "Number of miRNAs conserved in Apul and Peve: 1"

### 3.2.2 Apul and Pmea

Find Apul miRNAs that have hits to Pmea miRNAs but *not* hits to Peve
miRNAs

``` r
# Find Apul miRNAs that have matches from only Pmea
present_in_apul_pmea <- filtered_combined_blastn %>%
  # isolate Apul miRNAs with hits
  filter(!grepl("Porites_evermani|Pocillopora_meandrina", sseqid)) %>%
  group_by(sseqid) %>%
  # filter for hits to Pmea only
  filter(!any(grepl("Porites_evermani", qseqid)) & any(grepl("Pocillopora_meandrina", qseqid)))

# View the miRNAs that match between Apul and Pmea
head(present_in_apul_pmea, nrow(present_in_apul_pmea))
```

    # A tibble: 1 × 12
    # Groups:   sseqid [1]
      qseqid sseqid pident length mismatch gapopen qstart  qend sstart  send  evalue
      <chr>  <chr>   <dbl>  <int>    <int>   <int>  <int> <int>  <int> <int>   <dbl>
    1 Clust… Clust…    100     22        0       0      1    22      1    22 3.52e-9
    # ℹ 1 more variable: bitscore <dbl>

``` r
# Count the number of miRNAs conserved across the two species
paste("Number of miRNAs conserved in Apul and Pmea:", nrow(distinct(present_in_apul_pmea, sseqid)))
```

    [1] "Number of miRNAs conserved in Apul and Pmea: 1"

### 3.2.3 Peve and Pmea

Find Peve miRNAs that have hits to Pmea miRNAs but *not* hits to Apul
miRNAs

``` r
# Find Peve miRNAs that have matches from only Pmea
present_in_peve_pmea <- filtered_combined_blastn %>%
  # isolate Peve miRNAs with hits
  filter(grepl("Porites_evermani", sseqid)) %>%
  group_by(sseqid) %>%
  # filter for hits to Pmea only (note the Apul sequence IDs don't contain the species name, so we have to use a non-descriptive unique identifier for filtering)
  filter(!any(grepl("mature::N", qseqid)) & any(grepl("Pocillopora_meandrina", qseqid)))

# View the miRNAs that match between Peve and Pmea
head(present_in_peve_pmea, nrow(present_in_peve_pmea))
```

    # A tibble: 1 × 12
    # Groups:   sseqid [1]
      qseqid  sseqid pident length mismatch gapopen qstart  qend sstart  send evalue
      <chr>   <chr>   <dbl>  <int>    <int>   <int>  <int> <int>  <int> <int>  <dbl>
    1 Cluste… Clust…    100     12        0       0     11    22      4    15  0.001
    # ℹ 1 more variable: bitscore <dbl>

``` r
# Count the number of miRNAs conserved across the two species
paste("Number of miRNAs conserved in Peve and Pmea:", nrow(distinct(present_in_peve_pmea, sseqid)))
```

    [1] "Number of miRNAs conserved in Peve and Pmea: 1"
