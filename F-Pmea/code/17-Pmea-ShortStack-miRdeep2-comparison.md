17-Pmea-ShortStack-miRdeep2-comparison
================
Kathleen Durkin
2024-05-17

- <a href="#01-intersectbed" id="toc-01-intersectbed">0.1 intersectBed</a>
- <a href="#02-results" id="toc-02-results">0.2 Results</a>

## 0.1 intersectBed

Examine our input files (intersectBed accepts .bed and .gff files)

``` bash
head -5 ../output/11.1-Pmea-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/mirna_results_22_04_2024_t_15_27_22/novel_mature_22_04_2024_t_15_27_22_score-50_to_na.bed
echo ""
head -5 ../output/13.2.1-Pmea-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results.gff3
```

    browser position 
    browser hide all
    track name="notTrackname.novel_miRNAs" description="novel miRNAs detected by miRDeep2 for notTrackname" visibility=2
    itemRgb="On";
    Pocillopora_meandrina_HIv1___Sc0000000  20372434    20372456    Pocillopora_meandrina_HIv1___Sc0000000_60202    266239.8    -   20372434    20372456    0,0,255

    Pocillopora_meandrina_HIv1___Sc0000000  ShortStack  Unknown_sRNA_locus  9092    9521    10746   +   .   ID=Cluster_1;DicerCall=N;MIRNA=N
    Pocillopora_meandrina_HIv1___Sc0000000  ShortStack  Unknown_sRNA_locus  53578   53997   286 +   .   ID=Cluster_2;DicerCall=N;MIRNA=N
    Pocillopora_meandrina_HIv1___Sc0000000  ShortStack  Unknown_sRNA_locus  150243  150718  2549    -   .   ID=Cluster_3;DicerCall=N;MIRNA=N
    Pocillopora_meandrina_HIv1___Sc0000000  ShortStack  siRNA22_locus   173728  174150  1257    +   .   ID=Cluster_4;DicerCall=22;MIRNA=N
    Pocillopora_meandrina_HIv1___Sc0000000  ShortStack  Unknown_sRNA_locus  187562  188076  182 .   .   ID=Cluster_5;DicerCall=N;MIRNA=N

We need to get two input files that contain only mature miRNAs and are
correctly formatted. That means we need to remove the header lines of
the miRdeep2 mature miRNAs file, and the ShortStack full results file
needs to be filtered.

``` bash
# remove header lines
tail -n +5 ../output/11.1-Pmea-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/mirna_results_22_04_2024_t_15_27_22/novel_mature_22_04_2024_t_15_27_22_score-50_to_na.bed > ../output/11.1-Pmea-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/mirna_results_22_04_2024_t_15_27_22/novel_mature_22_04_2024_t_15_27_22_score-50_to_na._formatted.bed

# filter full results to obtain a gff file of only the mature miRNAs
awk -F'\t' '$3 == "mature_miRNA"' ../output/13.2.1-Pmea-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results.gff3 > ../output/13.2.1-Pmea-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results_mature.gff3
```

Check the files

``` bash
head -5 ../output/11.1-Pmea-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/mirna_results_22_04_2024_t_15_27_22/novel_mature_22_04_2024_t_15_27_22_score-50_to_na._formatted.bed
echo ""
head -5 ../output/13.2.1-Pmea-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results_mature.gff3
```

    Pocillopora_meandrina_HIv1___Sc0000000  20372434    20372456    Pocillopora_meandrina_HIv1___Sc0000000_60202    266239.8    -   20372434    20372456    0,0,255
    Pocillopora_meandrina_HIv1___Sc0000017  5050957 5050980 Pocillopora_meandrina_HIv1___Sc0000017_655645   67675.4 -   5050957 5050980 0,0,255
    Pocillopora_meandrina_HIv1___Sc0000005  601625  601646  Pocillopora_meandrina_HIv1___Sc0000005_273871   65863.2 +   601625  601646  255,0,0
    Pocillopora_meandrina_HIv1___xfSc0000017    39556   39581   Pocillopora_meandrina_HIv1___xfSc0000017_985353 65601.9 -   39556   39581   0,0,255
    Pocillopora_meandrina_HIv1___Sc0000009  11978171    11978196    Pocillopora_meandrina_HIv1___Sc0000009_442963   64606.4 -   11978171    11978196    0,0,255

    Pocillopora_meandrina_HIv1___Sc0000000  ShortStack  mature_miRNA    818049  818070  3240    +   .   ID=Cluster_21.mature;Parent=Cluster_21
    Pocillopora_meandrina_HIv1___Sc0000000  ShortStack  mature_miRNA    2872041 2872061 110 +   .   ID=Cluster_37.mature;Parent=Cluster_37
    Pocillopora_meandrina_HIv1___Sc0000000  ShortStack  mature_miRNA    20372469    20372490    769 -   .   ID=Cluster_361.mature;Parent=Cluster_361
    Pocillopora_meandrina_HIv1___Sc0000001  ShortStack  mature_miRNA    19145788    19145809    1102    +   .   ID=Cluster_759.mature;Parent=Cluster_759
    Pocillopora_meandrina_HIv1___Sc0000002  ShortStack  mature_miRNA    3841966 3841987 98704   -   .   ID=Cluster_918.mature;Parent=Cluster_918

Looks good! Now we can input into intersectBed.

``` bash
# intersectBed to ID sequences in the miRdeep2 mature miRNA output that match mature miRNAs ID'd by ShortStack
# -wa and -wb ensure we recieve full annotations from both input files in the output
/home/shared/bedtools2/bin/intersectBed \
-wa \
-wb \
-a ../output/13.2.1-Pmea-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results_mature.gff3 \
-b ../output/11.1-Pmea-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/mirna_results_22_04_2024_t_15_27_22/novel_mature_22_04_2024_t_15_27_22_score-50_to_na._formatted.bed \
&> ../output/17-Pmea-ShortStack-miRdeep2-comparison/ShortStack_miRdeep2_mature_intersect.bed
```

## 0.2 Results

``` bash
echo "Number of ShortStack mature miRNAs:"
wc -l < ../output/13.2.1-Pmea-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results_mature.gff3
echo ""
echo "Number of ShortStack mature miRNAs also identified by miRdeep2:"
wc -l < ../output/17-Pmea-ShortStack-miRdeep2-comparison/ShortStack_miRdeep2_mature_intersect.bed
```

    Number of ShortStack mature miRNAs:
    36

    Number of ShortStack mature miRNAs also identified by miRdeep2:
    27

While looking through the output file I noticed that two of the
intersects originate from the same cluster… not really sure what that’s
about…

``` bash
head -2 ../output/17-Pmea-ShortStack-miRdeep2-comparison/ShortStack_miRdeep2_mature_intersect.bed
```

    Pocillopora_meandrina_HIv1___Sc0000000  ShortStack  mature_miRNA    818049  818070  3240    +   .   ID=Cluster_21.mature;Parent=Cluster_21  Pocillopora_meandrina_HIv1___Sc0000000  818048  818071  Pocillopora_meandrina_HIv1___Sc0000000_1750 6153    +   818048  818071  255,0,0
    Pocillopora_meandrina_HIv1___Sc0000000  ShortStack  mature_miRNA    818049  818070  3240    +   .   ID=Cluster_21.mature;Parent=Cluster_21  Pocillopora_meandrina_HIv1___Sc0000000  818046  818069  Pocillopora_meandrina_HIv1___Sc0000000_34562    10  -   818046  818069  0,0,255

It looks like they have pretty much identical entries (cluster,
coordinates, etc.), except the end of the miRdeep locus name
(Pocillopora_meandrina_HIv1\_**Sc0000000_1750 vs
Pocillopora_meandrina_HIv1**\_Sc0000000_34562) and the miRdeep2 “score”
value assigned to them (this is the column following the miRdeep2 locus
name, 6153 vs 10). We can also check these two loci in the full miRdeep2
output.

``` bash

# View full mirdeep2 output for these two loci
awk -F'\t' '$1 == "Pocillopora_meandrina_HIv1___Sc0000000_1750"' ../output/11.1-Pmea-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_27_22.csv
echo""
awk -F'\t' '$1 == "Pocillopora_meandrina_HIv1___Sc0000000_34562"' ../output/11.1-Pmea-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_27_22.csv

echo""
echo""

echo "mature miRNA sequences for these two loci:"
awk -F'\t' '$1 == "Pocillopora_meandrina_HIv1___Sc0000000_1750"' ../output/11.1-Pmea-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_27_22.csv | awk '{print $13}'

awk -F'\t' '$1 == "Pocillopora_meandrina_HIv1___Sc0000000_34562"' ../output/11.1-Pmea-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_27_22.csv | awk '{print $13}'

echo ""
echo ""

echo "miRNA* sequences for these two loci:"
awk -F'\t' '$1 == "Pocillopora_meandrina_HIv1___Sc0000000_1750"' ../output/11.1-Pmea-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_27_22.csv | awk '{print $14}'

awk -F'\t' '$1 == "Pocillopora_meandrina_HIv1___Sc0000000_34562"' ../output/11.1-Pmea-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_27_22.csv | awk '{print $14}'

echo ""
echo ""

echo "precursor miRNA sequences for these two loci:"
awk -F'\t' '$1 == "Pocillopora_meandrina_HIv1___Sc0000000_1750"' ../output/11.1-Pmea-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_27_22.csv | awk '{print $15}'

awk -F'\t' '$1 == "Pocillopora_meandrina_HIv1___Sc0000000_34562"' ../output/11.1-Pmea-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_27_22.csv | awk '{print $15}'
```

    Pocillopora_meandrina_HIv1___Sc0000000_1750 6153        -   12060   11275   0   785 yes -   hsa-miR-2117_MIMAT0011162_Homo_sapiens_miR-2117 -   -   uguucucucugcaguaagcaugu augcuugcuguaaagagaacuug uguucucucugcaguaagcauguuuugacaugcuugcuguaaagagaacuug    Pocillopora_meandrina_HIv1___Sc0000000:818048..818100:+

    Pocillopora_meandrina_HIv1___Sc0000000_34562    10      -   11  9   0   2   yes -   egr-miR-153-5p_MIMAT0037428_Echinococcus_granulosus_miR-153-5p  -   -   augcuuacugcagagagaacaug aaguucucuuuacagcaagcaugucaaa    aaguucucuuuacagcaagcaugucaaaacaugcuuacugcagagagaacaug   Pocillopora_meandrina_HIv1___Sc0000000:818046..818099:-


    mature miRNA sequences for these two loci:
    uguucucucugcaguaagcaugu
    augcuuacugcagagagaacaug


    miRNA* sequences for these two loci:
    augcuugcuguaaagagaacuug
    aaguucucuuuacagcaagcaugucaaa


    precursor miRNA sequences for these two loci:
    uguucucucugcaguaagcauguuuugacaugcuugcuguaaagagaacuug
    aaguucucuuuacagcaagcaugucaaaacaugcuuacugcagagagaacaug

Interesting… The two loci have *very* similar precursor sequences and
reversed mature and star sequences! In other words, the mature miRNA
sequence for Pocillopora_meandrina_HIv1\_**Sc0000000_1750 is almost
identical to the miRNA\* sequence of
Pocillopora_meandrina_HIv1**\_Sc0000000_34562, and vice versa! I’m not
exactly sure what this means though… is miRdeep2 just incorrectly
distinguishing the mature and star sequences for one of these loci? Does
the much higher miRdeep2 score of
Pocillopora_meandrina_HIv1\_\_\_Sc0000000_1750 mean we should be more
confident in it being correctly distinguished?

Let’s set that aside for now and do some quick investigation of the
miRdeep2 evaluation criteria for all of these ShortStack/miRdeep2 shared
miRNAs. This could give us an idea of what thresholds may be appropriate
for filtering the miRdeep2 output.

``` bash
mirdeepIDs=$(awk '{print $13}' ../output/17-Pmea-ShortStack-miRdeep2-comparison/ShortStack_miRdeep2_mature_intersect.bed)

head -1 ../output/11.1-Pmea-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_27_22.csv > ../output/17-Pmea-ShortStack-miRdeep2-comparison/intersect_miRdeep2_stats.txt
while IFS= read -r id; do
  # Use awk to process fileA and match column 1 with the current ID
  awk -F'\t' -v id="$id" '$1 == id {print}' ../output/11.1-Pmea-sRNAseq-miRdeep2-31bp-fastp-merged-cnidarian_miRBase/parsable-result_22_04_2024_t_15_27_22.csv
done <<< "$mirdeepIDs" >> ../output/17-Pmea-ShortStack-miRdeep2-comparison/intersect_miRdeep2_stats.txt
```

``` r
intersect_miRdeep2_stats <- read.csv("../output/17-Pmea-ShortStack-miRdeep2-comparison/intersect_miRdeep2_stats.txt", sep="\t")

intersect_miRdeep2_stats %>% 
  select(miRDeep2.score, significant.randfold.p.value) %>% 
  arrange(desc(miRDeep2.score)) %>%
  kable()
```

| miRDeep2.score | significant.randfold.p.value |
|---------------:|:-----------------------------|
|        67675.4 | yes                          |
|        65863.2 | yes                          |
|        54159.5 | yes                          |
|        29656.8 | yes                          |
|        15853.5 | yes                          |
|        13573.1 | yes                          |
|         8138.8 | yes                          |
|         6153.0 | yes                          |
|         5845.2 | yes                          |
|         5647.4 | yes                          |
|         2344.8 | yes                          |
|         1515.4 | yes                          |
|         1188.9 | yes                          |
|          775.3 | yes                          |
|          712.0 | yes                          |
|          526.9 | yes                          |
|          509.1 | yes                          |
|          184.8 | yes                          |
|          170.7 | yes                          |
|          152.7 | yes                          |
|          136.5 | yes                          |
|           91.4 | yes                          |
|           59.1 | yes                          |
|           10.0 | yes                          |
|            5.3 | yes                          |
|            5.2 | yes                          |
|            5.1 | yes                          |

The miRdeep2 score is “the log-odds score assigned to the hairpin” and
is essentially a probability that the locus is a miRNA (presumably based
on the hairpin structure) on, with higher values indicating higher
probability. MiRdeep2’s default score threshold for miRNA classification
is 0, but coral miRNA papers we’ve seen use myriad thresholds (e.g., 4,
10).

The randfold is also an evaluation of ncrna secondary structure, and we
want a significant randfold value (this would indicate high likelihood
of miRNA precursory structure)

23 out of 27 shared miRNAs have miRdeep2 scores \>10, and all 30 have
significant randfold p-values, which is good support for using these
thresholds.
