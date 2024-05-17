17-Apul-ShortStack-miRdeep2-comparison
================
Kathleen Durkin
2024-05-17

- <a href="#01-intersectbed" id="toc-01-intersectbed">0.1 intersectBed</a>
- <a href="#02-results" id="toc-02-results">0.2 Results</a>

## 0.1 intersectBed

Examine our input files (intersectBed accepts .bed and .gff files)

``` bash
head -5 ../output/11.1-Apul-sRNAseq-miRdeep2-31bp-fastp-merged/mirna_results_03_04_2024_t_13_00_39/novel_mature_03_04_2024_t_13_00_39_score-50_to_na.bed
echo ""
head -5 ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results.gff3
```

    browser position 
    browser hide all
    track name="notTrackname.novel_miRNAs" description="novel miRNAs detected by miRDeep2 for notTrackname" visibility=2
    itemRgb="On";
    NC_058066.1_Acropora_millepora_isolate_JS-1_chromosome_1_Amil_v2.1_whole_genome_shotgun_sequence    20346247    20346271    NC_058066.1_Acropora_millepora_isolate_JS-1_chromosome_1_Amil_v2.1_whole_genome_shotgun_sequence_49320  127017.5    -   20346247    20346271    0,0,255

    NC_058066.1 ShortStack  Unknown_sRNA_locus  152483  152910  140 -   .   ID=Cluster_1;DicerCall=N;MIRNA=N
    NC_058066.1 ShortStack  Unknown_sRNA_locus  161064  161674  549 .   .   ID=Cluster_2;DicerCall=N;MIRNA=N
    NC_058066.1 ShortStack  Unknown_sRNA_locus  172073  172496  105 -   .   ID=Cluster_3;DicerCall=N;MIRNA=N
    NC_058066.1 ShortStack  Unknown_sRNA_locus  203242  203651  100 .   .   ID=Cluster_4;DicerCall=N;MIRNA=N
    NC_058066.1 ShortStack  Unknown_sRNA_locus  204535  205150  313 .   .   ID=Cluster_5;DicerCall=N;MIRNA=N

We need to get two input files that contain only mature miRNAs and have
the same sequence naming conventions for matching purposes. That means
the miRdeep2 mature miRNAs file need to be reformatted, and the
ShortStack full results file needs to be filtered.

``` bash
# remove header lines, shorten chromosome ref names to match those used in ShortStack output, and ensure final file is tab-delimited instead of space-delimited
tail -n +5 ../output/11.1-Apul-sRNAseq-miRdeep2-31bp-fastp-merged/mirna_results_03_04_2024_t_13_00_39/novel_mature_03_04_2024_t_13_00_39_score-50_to_na.bed | \
awk 'BEGIN { FS = "\t" } { sub(/\.1.*/, ".1", $1); print }' | \
tr ' ' '\t' > ../output/11.1-Apul-sRNAseq-miRdeep2-31bp-fastp-merged/mirna_results_03_04_2024_t_13_00_39/novel_mature_03_04_2024_t_13_00_39_score-50_to_na_formatted.bed

# filter full results to obtain a gff file of only the mature miRNAs
awk -F'\t' '$3 == "mature_miRNA"' ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results.gff3 > ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results_mature.gff3
```

Check the files

``` bash
head -5 ../output/11.1-Apul-sRNAseq-miRdeep2-31bp-fastp-merged/mirna_results_03_04_2024_t_13_00_39/novel_mature_03_04_2024_t_13_00_39_score-50_to_na_formatted.bed
echo ""
head -5 ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results_mature.gff3
```

    NC_058066.1 20346247    20346271    NC_058066.1_Acropora_millepora_isolate_JS-1_chromosome_1_Amil_v2.1_whole_genome_shotgun_sequence_49320  127017.5    -   20346247    20346271    0,0,255
    NC_058070.1 11598966    11598988    NC_058070.1_Acropora_millepora_isolate_JS-1_chromosome_5_Amil_v2.1_whole_genome_shotgun_sequence_208678 40616.1 +   11598966    11598988    255,0,0
    NC_058068.1 597877  597899  NC_058068.1_Acropora_millepora_isolate_JS-1_chromosome_3_Amil_v2.1_whole_genome_shotgun_sequence_123927 40339.4 +   597877  597899  255,0,0
    NC_058071.1 8840802 8840823 NC_058071.1_Acropora_millepora_isolate_JS-1_chromosome_6_Amil_v2.1_whole_genome_shotgun_sequence_256047 39769.6 +   8840802 8840823 255,0,0
    NC_058068.1 598173  598195  NC_058068.1_Acropora_millepora_isolate_JS-1_chromosome_3_Amil_v2.1_whole_genome_shotgun_sequence_123929 38710.3 +   598173  598195  255,0,0

    NC_058066.1 ShortStack  mature_miRNA    12757147    12757168    1413    -   .   ID=Cluster_316.mature;Parent=Cluster_316
    NC_058066.1 ShortStack  mature_miRNA    20088679    20088700    102 +   .   ID=Cluster_514.mature;Parent=Cluster_514
    NC_058066.1 ShortStack  mature_miRNA    20346249    20346271    27205   -   .   ID=Cluster_548.mature;Parent=Cluster_548
    NC_058067.1 ShortStack  mature_miRNA    5656214 5656236 1338    -   .   ID=Cluster_1506.mature;Parent=Cluster_1506
    NC_058067.1 ShortStack  mature_miRNA    16118270    16118291    7486    -   .   ID=Cluster_1900.mature;Parent=Cluster_1900

Looks good! Now we can input into intersectBed.

``` bash
# intersectBed to ID sequences in the miRdeep2 mature miRNA output that match mature miRNAs ID'd by ShortStack
# -wa and -wb ensure we recieve full annotations from both input files in the output
/home/shared/bedtools2/bin/intersectBed \
-wa \
-wb \
-a ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results_mature.gff3 \
-b ../output/11.1-Apul-sRNAseq-miRdeep2-31bp-fastp-merged/mirna_results_03_04_2024_t_13_00_39/novel_mature_03_04_2024_t_13_00_39_score-50_to_na_formatted.bed \
&> ../output/17-Apul-ShortStack-miRdeep2-comparison/ShortStack_miRdeep2_mature_intersect.bed
```

## 0.2 Results

``` bash
echo "Number of ShortStack mature miRNAs:"
wc -l < ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results_mature.gff3
echo ""
echo "Number of ShortStack mature miRNAs also identified by miRdeep2:"
wc -l < ../output/17-Apul-ShortStack-miRdeep2-comparison/ShortStack_miRdeep2_mature_intersect.bed
```

    Number of ShortStack mature miRNAs:
    38

    Number of ShortStack mature miRNAs also identified by miRdeep2:
    35

While looking through the output file I noticed that two of the
intersects originate from the same cluster… not really sure what that’s
about…

``` bash
head -21 ../output/17-Apul-ShortStack-miRdeep2-comparison/ShortStack_miRdeep2_mature_intersect.bed | tail -2
```

    NC_058072.1 ShortStack  mature_miRNA    19030618    19030639    21560   +   .   ID=Cluster_6376.mature;Parent=Cluster_6376  NC_058072.1 19030617    19030639    NC_058072.1_Acropora_millepora_isolate_JS-1_chromosome_7_Amil_v2.1_whole_genome_shotgun_sequence_295614 12623.6 +   19030617    19030639    255,0,0
    NC_058072.1 ShortStack  mature_miRNA    19030618    19030639    21560   +   .   ID=Cluster_6376.mature;Parent=Cluster_6376  NC_058072.1 19030617    19030639    NC_058072.1_Acropora_millepora_isolate_JS-1_chromosome_7_Amil_v2.1_whole_genome_shotgun_sequence_295613 1.1 +   19030617    19030639    255,0,0

It looks like they have pretty much identical entries (cluster,
coordinates, etc.), except the end of the miRdeep locus name
(NC_058072.1_Acropora_millepora_isolate_JS-1_chromosome_7\_Amil_v2.1_whole_genome_shotgun_sequence\_***295614***
vs
NC_058072.1_Acropora_millepora_isolate_JS-1_chromosome_7\_Amil_v2.1_whole_genome_shotgun_sequence\_***295613***)
and the miRdeep2 “score” value assigned to them (this is the column
following the miRdeep2 locus name, 12623.6 vs 1.1). We can also check
these two loci in the full miRdeep2 output.

``` bash

# View full mirdeep2 output for these two loci
awk -F'\t' '$1 == "NC_058072.1_Acropora_millepora_isolate_JS-1_chromosome_7_Amil_v2.1_whole_genome_shotgun_sequence_295614"' ../output/11.1-Apul-sRNAseq-miRdeep2-31bp-fastp-merged/parsable-result_03_04_2024_t_13_00_39.csv
echo""
awk -F'\t' '$1 == "NC_058072.1_Acropora_millepora_isolate_JS-1_chromosome_7_Amil_v2.1_whole_genome_shotgun_sequence_295613"' ../output/11.1-Apul-sRNAseq-miRdeep2-31bp-fastp-merged/parsable-result_03_04_2024_t_13_00_39.csv

echo""
echo""

echo "mature miRNA sequences for these two loci:"
awk -F'\t' '$1 == "NC_058072.1_Acropora_millepora_isolate_JS-1_chromosome_7_Amil_v2.1_whole_genome_shotgun_sequence_295614"' ../output/11.1-Apul-sRNAseq-miRdeep2-31bp-fastp-merged/parsable-result_03_04_2024_t_13_00_39.csv | awk '{print $13}'

awk -F'\t' '$1 == "NC_058072.1_Acropora_millepora_isolate_JS-1_chromosome_7_Amil_v2.1_whole_genome_shotgun_sequence_295613"' ../output/11.1-Apul-sRNAseq-miRdeep2-31bp-fastp-merged/parsable-result_03_04_2024_t_13_00_39.csv | awk '{print $13}'

echo ""
echo ""

echo "precursor miRNA sequences for these two loci:"
awk -F'\t' '$1 == "NC_058072.1_Acropora_millepora_isolate_JS-1_chromosome_7_Amil_v2.1_whole_genome_shotgun_sequence_295614"' ../output/11.1-Apul-sRNAseq-miRdeep2-31bp-fastp-merged/parsable-result_03_04_2024_t_13_00_39.csv | awk '{print $15}'

awk -F'\t' '$1 == "NC_058072.1_Acropora_millepora_isolate_JS-1_chromosome_7_Amil_v2.1_whole_genome_shotgun_sequence_295613"' ../output/11.1-Apul-sRNAseq-miRdeep2-31bp-fastp-merged/parsable-result_03_04_2024_t_13_00_39.csv | awk '{print $15}'
```

    NC_058072.1_Acropora_millepora_isolate_JS-1_chromosome_7_Amil_v2.1_whole_genome_shotgun_sequence_295614 12623.6     -   24755   24435   0   320 no  -   tca-miR-11646-3p_MIMAT0045620_Tribolium_castaneum_miR-11646-3p  -   -   ugggugucaucuauuauguuuu  aacauaaaagauggcacc  ugggugucaucuauuauguuuuugcuuguuaaaacauaaaagauggcacc  NC_058072.1_Acropora_millepora_isolate_JS-1_chromosome_7_Amil_v2.1_whole_genome_shotgun_sequence:19030617..19030667:+

    NC_058072.1_Acropora_millepora_isolate_JS-1_chromosome_7_Amil_v2.1_whole_genome_shotgun_sequence_295613 1.1     -   24449   24435   14  0   no  -   tca-miR-11646-3p_MIMAT0045620_Tribolium_castaneum_miR-11646-3p  -   -   ugggugucaucuauuauguuuu  aauguaacaaaauugacggccaga    aauguaacaaaauugacggccagaagccguacguauguagaaaauguggggugagugccugggugucaucuauuauguuuu   NC_058072.1_Acropora_millepora_isolate_JS-1_chromosome_7_Amil_v2.1_whole_genome_shotgun_sequence:19030558..19030639:+


    mature miRNA sequences for these two loci:
    ugggugucaucuauuauguuuu
    ugggugucaucuauuauguuuu


    precursor miRNA sequences for these two loci:
    ugggugucaucuauuauguuuuugcuuguuaaaacauaaaagauggcacc
    aauguaacaaaauugacggccagaagccguacguauguagaaaauguggggugagugccugggugucaucuauuauguuuu

Interesting… The two loci have identical mature miRNAs, but different
precursors!

Let’s set that aside for now and do some quick investigation of the
miRdeep2 evaluation criteria for all of these ShortStack/miRdeep2 shared
miRNAs. This could give us an idea of what thresholds may be appropriate
for filtering the miRdeep2 output.

``` bash
mirdeepIDs=$(awk '{print $13}' ../output/17-Apul-ShortStack-miRdeep2-comparison/ShortStack_miRdeep2_mature_intersect.bed)

head -1 ../output/11.1-Apul-sRNAseq-miRdeep2-31bp-fastp-merged/parsable-result_03_04_2024_t_13_00_39.csv > ../output/17-Apul-ShortStack-miRdeep2-comparison/intersect_miRdeep2_stats.txt
while IFS= read -r id; do
  # Use awk to process fileA and match column 1 with the current ID
  awk -F'\t' -v id="$id" '$1 == id {print}' ../output/11.1-Apul-sRNAseq-miRdeep2-31bp-fastp-merged/parsable-result_03_04_2024_t_13_00_39.csv
done <<< "$mirdeepIDs" >> ../output/17-Apul-ShortStack-miRdeep2-comparison/intersect_miRdeep2_stats.txt
```

``` r
intersect_miRdeep2_stats <- read.csv("../output/17-Apul-ShortStack-miRdeep2-comparison/intersect_miRdeep2_stats.txt", sep="\t")

intersect_miRdeep2_stats %>% 
  select(miRDeep2.score, significant.randfold.p.value) %>% 
  arrange(desc(miRDeep2.score)) %>%
  kable()
```

| miRDeep2.score | significant.randfold.p.value |
|---------------:|:-----------------------------|
|       127017.5 | no                           |
|        40616.1 | no                           |
|        40339.4 | no                           |
|        39769.6 | no                           |
|        38710.3 | no                           |
|        14705.5 | no                           |
|        12860.5 | no                           |
|        12623.6 | no                           |
|         7816.8 | no                           |
|         6367.0 | no                           |
|         5869.1 | no                           |
|         5325.7 | no                           |
|         4346.5 | no                           |
|         4203.4 | no                           |
|         2769.9 | no                           |
|         2323.0 | no                           |
|         1422.1 | no                           |
|          957.0 | no                           |
|          928.6 | no                           |
|          637.6 | no                           |
|          494.9 | no                           |
|          296.0 | no                           |
|          214.5 | no                           |
|          202.9 | no                           |
|          178.3 | no                           |
|          171.5 | no                           |
|          117.8 | no                           |
|           97.5 | no                           |
|           97.3 | no                           |
|           78.4 | no                           |
|           23.2 | no                           |
|           23.1 | no                           |
|            3.7 | no                           |
|            3.4 | no                           |
|            1.1 | no                           |

The miRdeep2 score is “the log-odds score assigned to the hairpin” and
is essentially a probability that the locus is a miRNA (presumably based
on the hairpin structure) on, with higher values indicating higher
probability. MiRdeep2’s default score threshold for miRNA classification
is 0, but coral miRNA papers we’ve seen use myriad thresholds (e.g., 4,
10).

The randfold is also an evaluation of ncrna secondary structure, and we
want a significant randfold value (this would indicate high likelihood
of miRNA precursory structure)

For some reason there are no miRdeep2 miRNAs with significant randfold
p-values for A. pulchra, even though the other species have a bunch of
loci with significant p-values. The miRdeep2 score is interesting
though! All but 3 of the shared miRNAs have scores \>10 – maybe that
could be a good cutoff?
