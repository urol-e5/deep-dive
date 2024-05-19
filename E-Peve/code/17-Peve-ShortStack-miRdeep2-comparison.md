17-Peve-ShortStack-miRdeep2-comparison
================
Kathleen Durkin
2024-05-17

- <a href="#01-intersectbed" id="toc-01-intersectbed">0.1 intersectBed</a>
- <a href="#02-results" id="toc-02-results">0.2 Results</a>

## 0.1 intersectBed

Examine our input files (intersectBed accepts .bed and .gff files)

``` bash
head -5 ../output/11.1-Peve-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/mirna_results_22_04_2024_t_15_10_16/novel_mature_22_04_2024_t_15_10_16_score-50_to_na.bed
echo ""
head -5 ../output/08.2-Peve-sRNAseq-ShortStack-31bp-fastp-merged/ShortStack_out/Results.gff3
```

    browser position 
    browser hide all
    track name="notTrackname.novel_miRNAs" description="novel miRNAs detected by miRDeep2 for notTrackname" visibility=2
    itemRgb="On";
    Porites_evermani_scaffold_1503  46899   46921   Porites_evermani_scaffold_1503_530575   110319.9    +   46899   46921   255,0,0

    Porites_evermani_scaffold_1 ShortStack  Unknown_sRNA_locus  45711   46131   88  +   .   ID=Cluster_1;DicerCall=N;MIRNA=N
    Porites_evermani_scaffold_1 ShortStack  Unknown_sRNA_locus  201507  201931  58  -   .   ID=Cluster_2;DicerCall=N;MIRNA=N
    Porites_evermani_scaffold_1 ShortStack  Unknown_sRNA_locus  313446  313846  50  -   .   ID=Cluster_3;DicerCall=N;MIRNA=N
    Porites_evermani_scaffold_1 ShortStack  Unknown_sRNA_locus  406146  406734  175 -   .   ID=Cluster_4;DicerCall=N;MIRNA=N
    Porites_evermani_scaffold_1 ShortStack  Unknown_sRNA_locus  409839  410269  169 -   .   ID=Cluster_5;DicerCall=N;MIRNA=N

We need to get two input files that contain only mature miRNAs and are
correctly formatted. That means we need to remove the header lines of
the miRdeep2 mature miRNAs file, and the ShortStack full results file
needs to be filtered.

``` bash
# remove header lines
tail -n +5 ../output/11.1-Peve-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/mirna_results_22_04_2024_t_15_10_16/novel_mature_22_04_2024_t_15_10_16_score-50_to_na.bed > ../output/11.1-Peve-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/mirna_results_22_04_2024_t_15_10_16/novel_mature_22_04_2024_t_15_10_16_score-50_to_na._formatted.bed

# filter full results to obtain a gff file of only the mature miRNAs
awk -F'\t' '$3 == "mature_miRNA"' ../output/08.2-Peve-sRNAseq-ShortStack-31bp-fastp-merged/ShortStack_out/Results.gff3 > ../output/08.2-Peve-sRNAseq-ShortStack-31bp-fastp-merged/ShortStack_out/Results_mature.gff3
```

Check the files

``` bash
head -5 ../output/11.1-Peve-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/mirna_results_22_04_2024_t_15_10_16/novel_mature_22_04_2024_t_15_10_16_score-50_to_na._formatted.bed
echo ""
head -5 ../output/08.2-Peve-sRNAseq-ShortStack-31bp-fastp-merged/ShortStack_out/Results_mature.gff3
```

    Porites_evermani_scaffold_1503  46899   46921   Porites_evermani_scaffold_1503_530575   110319.9    +   46899   46921   255,0,0
    Porites_evermani_scaffold_26    382571  382593  Porites_evermani_scaffold_26_42156  46827.4 -   382571  382593  0,0,255
    Porites_evermani_scaffold_910   118741  118762  Porites_evermani_scaffold_910_418426    43145.7 +   118741  118762  255,0,0
    Porites_evermani_scaffold_72    198220  198245  Porites_evermani_scaffold_72_87796  31859.5 +   198220  198245  255,0,0
    Porites_evermani_scaffold_1503  47610   47632   Porites_evermani_scaffold_1503_530579   27800.7 +   47610   47632   255,0,0

    Porites_evermani_scaffold_1 ShortStack  mature_miRNA    1404272 1404293 3403    -   .   ID=Cluster_29.mature;Parent=Cluster_29
    Porites_evermani_scaffold_16    ShortStack  mature_miRNA    383437  383458  1456    -   .   ID=Cluster_578.mature;Parent=Cluster_578
    Porites_evermani_scaffold_26    ShortStack  mature_miRNA    382572  382593  69226   -   .   ID=Cluster_786.mature;Parent=Cluster_786
    Porites_evermani_scaffold_47    ShortStack  mature_miRNA    475994  476015  240 -   .   ID=Cluster_1125.mature;Parent=Cluster_1125
    Porites_evermani_scaffold_49    ShortStack  mature_miRNA    151640  151661  35368   -   .   ID=Cluster_1153.mature;Parent=Cluster_1153

Looks good! Now we can input into intersectBed.

``` bash
# intersectBed to ID sequences in the miRdeep2 mature miRNA output that match mature miRNAs ID'd by ShortStack
# -wa and -wb ensure we recieve full annotations from both input files in the output
/home/shared/bedtools2/bin/intersectBed \
-wa \
-wb \
-a ../output/08.2-Peve-sRNAseq-ShortStack-31bp-fastp-merged/ShortStack_out/Results_mature.gff3 \
-b ../output/11.1-Peve-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/mirna_results_22_04_2024_t_15_10_16/novel_mature_22_04_2024_t_15_10_16_score-50_to_na._formatted.bed \
&> ../output/17-Peve-ShortStack-miRdeep2-comparison/ShortStack_miRdeep2_mature_intersect.bed
```

## 0.2 Results

``` bash
echo "Number of ShortStack mature miRNAs:"
wc -l < ../output/08.2-Peve-sRNAseq-ShortStack-31bp-fastp-merged/ShortStack_out/Results_mature.gff3
echo ""
echo "Number of ShortStack mature miRNAs also identified by miRdeep2:"
wc -l < ../output/17-Peve-ShortStack-miRdeep2-comparison/ShortStack_miRdeep2_mature_intersect.bed
```

    Number of ShortStack mature miRNAs:
    46

    Number of ShortStack mature miRNAs also identified by miRdeep2:
    30

While looking through the output file I noticed that two of the
intersects originate from the same cluster… not really sure what that’s
about…

``` bash
head -9 ../output/17-Peve-ShortStack-miRdeep2-comparison/ShortStack_miRdeep2_mature_intersect.bed | tail -2
```

    Porites_evermani_scaffold_334   ShortStack  mature_miRNA    153605  153626  142 -   .   ID=Cluster_4722.mature;Parent=Cluster_4722  Porites_evermani_scaffold_334   153606  153626  Porites_evermani_scaffold_334_234019    5.6 -   153606  153626  0,0,255
    Porites_evermani_scaffold_334   ShortStack  mature_miRNA    153605  153626  142 -   .   ID=Cluster_4722.mature;Parent=Cluster_4722  Porites_evermani_scaffold_334   153605  153625  Porites_evermani_scaffold_334_233889    5.5 +   153605  153625  255,0,0

It looks like they have pretty much identical entries (cluster,
coordinates, etc.), except the end of the miRdeep locus name
(Porites_evermani_scaffold_334_234019 vs
Porites_evermani_scaffold_334_233889) and the miRdeep2 “score” value
assigned to them (this is the column following the miRdeep2 locus name,
5.6 vs 5.1). We can also check these two loci in the full miRdeep2
output.

``` bash

# View full mirdeep2 output for these two loci
awk -F'\t' '$1 == "Porites_evermani_scaffold_334_234019"' ../output/11.1-Peve-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_10_16.csv
echo""
awk -F'\t' '$1 == "Porites_evermani_scaffold_334_233889"' ../output/11.1-Peve-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_10_16.csv

echo""
echo""

echo "mature miRNA sequences for these two loci:"
awk -F'\t' '$1 == "Porites_evermani_scaffold_334_234019"' ../output/11.1-Peve-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_10_16.csv | awk '{print $13}'

awk -F'\t' '$1 == "Porites_evermani_scaffold_334_233889"' ../output/11.1-Peve-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_10_16.csv | awk '{print $13}'

echo ""
echo ""

echo "miRNA* sequences for these two loci:"
awk -F'\t' '$1 == "Porites_evermani_scaffold_334_234019"' ../output/11.1-Peve-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_10_16.csv | awk '{print $14}'

awk -F'\t' '$1 == "Porites_evermani_scaffold_334_233889"' ../output/11.1-Peve-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_10_16.csv | awk '{print $14}'

echo ""
echo ""

echo "precursor miRNA sequences for these two loci:"
awk -F'\t' '$1 == "Porites_evermani_scaffold_334_234019"' ../output/11.1-Peve-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_10_16.csv | awk '{print $15}'

awk -F'\t' '$1 == "Porites_evermani_scaffold_334_233889"' ../output/11.1-Peve-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_10_16.csv | awk '{print $15}'
```

    Porites_evermani_scaffold_334_234019    5.6     -   969 686 0   283 yes -   gga-miR-12259-5p_MIMAT0050009_Gallus_gallus_miR-12259-5p    -   -   ugcagguacaguuauaaggu    accuuauaacuguaccugccaa  ugcagguacaguuauaagguccccuugguggaccuuauaacuguaccugccaa   Porites_evermani_scaffold_334:153573..153626:-

    Porites_evermani_scaffold_334_233889    5.5     -   111 96  0   15  yes -   cpi-miR-9592-5p_MIMAT0037980_Chrysemys_picta_miR-9592-5p    -   -   gaccuuauaacuguaccugc    gcagguacaguuauaaggucc   gcagguacaguuauaagguccaccaaggggaccuuauaacuguaccugc   Porites_evermani_scaffold_334:153576..153625:+


    mature miRNA sequences for these two loci:
    ugcagguacaguuauaaggu
    gaccuuauaacuguaccugc


    miRNA* sequences for these two loci:
    accuuauaacuguaccugccaa
    gcagguacaguuauaaggucc


    precursor miRNA sequences for these two loci:
    ugcagguacaguuauaagguccccuugguggaccuuauaacuguaccugccaa
    gcagguacaguuauaagguccaccaaggggaccuuauaacuguaccugc

Interesting… The two loci have *very* similar precursor sequences and
reversed mature and star sequences! In other words, the mature miRNA
sequence for Porites_evermani_scaffold_334_234019 is almost identical to
the miRNA\* sequence of Porites_evermani_scaffold_334_233889, and vice
versa! I’m not exactly sure what this means though… is miRdeep2 just
incorrectly distinguishing the mature and star sequences for one of
these loci?

Let’s set that aside for now and do some quick investigation of the
miRdeep2 evaluation criteria for all of these ShortStack/miRdeep2 shared
miRNAs. This could give us an idea of what thresholds may be appropriate
for filtering the miRdeep2 output.

``` bash
mirdeepIDs=$(awk '{print $13}' ../output/17-Peve-ShortStack-miRdeep2-comparison/ShortStack_miRdeep2_mature_intersect.bed)

head -1 ../output/11.1-Peve-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_10_16.csv > ../output/17-Peve-ShortStack-miRdeep2-comparison/intersect_miRdeep2_stats.txt
while IFS= read -r id; do
  # Use awk to process fileA and match column 1 with the current ID
  awk -F'\t' -v id="$id" '$1 == id {print}' ../output/11.1-Peve-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_10_16.csv
done <<< "$mirdeepIDs" >> ../output/17-Peve-ShortStack-miRdeep2-comparison/intersect_miRdeep2_stats.txt
```

``` r
intersect_miRdeep2_stats <- read.csv("../output/17-Peve-ShortStack-miRdeep2-comparison/intersect_miRdeep2_stats.txt", sep="\t")

intersect_miRdeep2_stats %>% 
  select(miRDeep2.score, significant.randfold.p.value) %>% 
  arrange(desc(miRDeep2.score)) %>%
  kable()
```

| miRDeep2.score | significant.randfold.p.value |
|---------------:|:-----------------------------|
|        46827.4 | yes                          |
|        43145.7 | yes                          |
|        22726.6 | yes                          |
|        17914.4 | yes                          |
|        16354.3 | yes                          |
|        13472.1 | yes                          |
|         6732.5 | yes                          |
|         4846.3 | yes                          |
|         3020.5 | yes                          |
|         2121.3 | yes                          |
|         2121.2 | yes                          |
|         1893.3 | yes                          |
|         1686.5 | yes                          |
|         1635.5 | yes                          |
|         1553.3 | yes                          |
|         1455.3 | yes                          |
|         1371.1 | yes                          |
|         1205.2 | yes                          |
|          998.8 | yes                          |
|          861.7 | yes                          |
|          468.6 | yes                          |
|          204.7 | yes                          |
|          184.6 | yes                          |
|          182.4 | yes                          |
|            5.6 | yes                          |
|            5.5 | yes                          |
|            5.5 | yes                          |
|            5.1 | yes                          |
|            4.4 | yes                          |
|            0.4 | yes                          |

The miRdeep2 score is “the log-odds score assigned to the hairpin” and
is essentially a probability that the locus is a miRNA (presumably based
on the hairpin structure) on, with higher values indicating higher
probability. MiRdeep2’s default score threshold for miRNA classification
is 0, but coral miRNA papers we’ve seen use myriad thresholds (e.g., 4,
10).

The randfold is also an evaluation of ncrna secondary structure, and we
want a significant randfold value (this would indicate high likelihood
of miRNA precursory structure)

25 out of 30 shared miRNAs have miRdeep2 scores \>10, and all 30 have
significant randfold p-values, which is good support for using these
thresholds.
