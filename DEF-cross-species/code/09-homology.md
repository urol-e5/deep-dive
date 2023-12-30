08-homology
================
Steven Roberts
2023-12-29

``` bash
grep ">" -c ../../D-Apul/output/05.33-lncRNA-discovery/Apul_lncRNA.fasta

head -2 ../../D-Apul/output/05.33-lncRNA-discovery/Apul_lncRNA.fasta 
```

    16206
    >transcript::NC_058066.1:468618-469943
    taactgatcaaaacgtatcttcctacaacattaatttgacagtggcgtttctcaactgac

``` bash
grep ">" -c ../../E-Peve/output/05-lncRNA-discovery/Peve_lncRNA.fasta

head -2 ../../E-Peve/output/05-lncRNA-discovery/Peve_lncRNA.fasta 
```

    7378
    >transcript::Porites_evermani_scaffold_1:422643-423512
    GGCAAAGCCACAATCCATGATAAATGAGGGCATAAGCCCGAGGAAGAAGAACTCATAGAT

``` bash
grep ">" -c ../../E-Peve/output/Peve_lncRNA.fasta

head -2 ../../E-Peve/output/Peve_lncRNA.fasta 
```

    7378
    >transcript::Porites_evermani_scaffold_1:422643-423512
    GGCAAAGCCACAATCCATGATAAATGAGGGCATAAGCCCGAGGAAGAAGAACTCATAGAT

``` bash
grep ">" -c ../../F-Pmea/output/02-lncRNA-discovery/Pmea_lncRNA.fasta

head -2 ../../F-Pmea/output/02-lncRNA-discovery/Pmea_lncRNA.fasta 
```

    12395
    >transcript::Pocillopora_meandrina_HIv1___Sc0000000:164391-165433
    TGTCACGTTTATCTTCATGTAAAATGTTTTCGATTTCTTGTGAGCGACGAAAACCATCTG

``` bash
cat \
../../D-Apul/output/05.33-lncRNA-discovery/Apul_lncRNA.fasta \
../../E-Peve/output/05-lncRNA-discovery/Peve_lncRNA.fasta \
../../F-Pmea/output/02-lncRNA-discovery/Pmea_lncRNA.fasta \
> ../output/09-homology/merged.fasta
```

``` bash
grep ">" -c ../output/09-homology/merged.fasta

head -4 ../output/09-homology/merged.fasta
```

    35979
    >transcript::NC_058066.1:468618-469943
    taactgatcaaaacgtatcttcctacaacattaatttgacagtggcgtttctcaactgac
    caatcaaaacttacatttgaaaatttggtgATGGTgcgtttacaactcgtgtatctttac
    gtcacacaaccatgtttgcATACTCTCTTGCaaccacgcctctcggccaatcagagcgcg

``` bash
wc -l ../output/08-comparative-BLASTs/*tab

head -2 ../output/08-comparative-BLASTs/*tab
```

       18513 ../output/08-comparative-BLASTs/Apul.tab
       42080 ../output/08-comparative-BLASTs/combined_blast.tab
        9430 ../output/08-comparative-BLASTs/Peve.tab
       14136 ../output/08-comparative-BLASTs/Pmea.tab
       84159 total
    ==> ../output/08-comparative-BLASTs/Apul.tab <==
    transcript::NC_058066.1:468618-469943   transcript::NC_058066.1:468618-469943   100.000 1325    0   0   1   1325    1   1325    0.0 2390
    transcript::NC_058066.1:1135315-1144814 transcript::NC_058066.1:1135315-1144814 100.000 9499    0   0   1   9499    1   9499    0.0 17131

    ==> ../output/08-comparative-BLASTs/combined_blast.tab <==
    V1  V2  V3  V4  V5  V6  V7  V8  V9  V10 V11 V12
    transcript::NC_058066.1:468618-469943   transcript::NC_058066.1:468618-469943   100 1325    0   0   1   1325    1   1325    0   2390

    ==> ../output/08-comparative-BLASTs/Peve.tab <==
    transcript::NC_058066.1:1135315-1144814 transcript::Porites_evermani_scaffold_50:639161-645099  93.939  198 12  0   8637    8834    3079    2882    4.36e-81    304
    transcript::NC_058066.1:1153398-1165634 transcript::Porites_evermani_scaffold_3148:28244-29656  83.389  301 31  6   6297    6597    724 443 4.61e-82    308

    ==> ../output/08-comparative-BLASTs/Pmea.tab <==
    transcript::NC_058066.1:1135315-1144814 transcript::Pocillopora_meandrina_HIv1___Sc0000040:1878921-1903778  91.549  355 26  4   8482    8835    16484   16133   2.87e-137   491
    transcript::NC_058066.1:2836056-2837849 transcript::Pocillopora_meandrina_HIv1___Sc0000001:4735350-4735856  83.128  243 41  0   902 1144    264 22  2.44e-66    254

Joining table in R.

``` bash
perl -e '$count=0; $len=0; while(<>) {s/\r?\n//; s/\t/ /g; if (s/^>//) { if ($. != 1) {print "\n"} s/ |$/\t/; $count++; $_ .= "\t";} else {s/ //g; $len += length($_)} print $_;} print "\n"; warn "\nConverted $count FASTA records in $. lines to tabular format\nTotal sequence length: $len\n\n";' \
../output/09-homology/merged.fasta > ../output/09-homology/merged.tab
```

``` r
query <- read.csv("../output/09-homology/merged.tab", sep = '\t', header = FALSE, row.names=NULL)
```

``` r
apul <- read.csv("../output/08-comparative-BLASTs/Apul.tab", sep = '\t', header = FALSE, row.names=NULL)
```

``` r
peve <- read.csv("../output/08-comparative-BLASTs/Peve.tab", sep = '\t', header = FALSE, row.names=NULL)
```

``` r
pmea <- read.csv("../output/08-comparative-BLASTs/Pmea.tab", sep = '\t', header = FALSE, row.names=NULL)
```

``` r
comp <- left_join(query, apul, by = "V1") %>%
  left_join(peve, by = "V1") %>%
  left_join(pmea, by = "V1") %>%
  select(V1, apul_hit = V2.y, apul_evalue = V11.x, peve_hit = V2.x.x, peve_evalue = V11.y, pmea_hit = V2.y.y, pmea_evalue = V11) %>%
   rowwise() %>%
  mutate(Hits = sum(!is.na(c_across(c(apul_hit, peve_hit, pmea_hit)))))
```

``` r
datatable(comp)
```

``` r
count_table <- table(comp$Hits)
print(count_table)
```


        0     1     2     3 
       16 32266  5216  1967 

Based on this

1967 sequences are found in all species and 5216 sequences are present
in 2 species

It is worth doing some taxonomy comparisons to see relatedness of
species,etc
