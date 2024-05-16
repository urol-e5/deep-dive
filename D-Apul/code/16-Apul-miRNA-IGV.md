16-miRNA vis
================
Steven Roberts
16 May, 2024

- <a href="#1-genome" id="toc-1-genome">1 Genome</a>
- <a href="#2-short-stack" id="toc-2-short-stack">2 short stack</a>

Apul

# 1 Genome

<https://gannet.fish.washington.edu/seashell/bu-github/deep-dive/D-Apul/data/GCF_013753865.1_Amil_v2.1_genomic.fna>

Get stuff from NCBI

<https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_013753865.1/>

``` bash
/home/shared/datasets download genome accession GCF_013753865.1 --include gff3,gtf,gbff
```

``` bash

unzip ../data/ncbi_dataset.zip
```

# 2 short stack

/seashell/bu-github/deep-dive/D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out

See:
<https://gannet.fish.washington.edu/seashell/bu-github/deep-dive/D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/>

``` bash
head ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/*
```

    ==> ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/alignment_details.tsv <==
    readfile    mapping_type    read_length count
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/sRNA-ACR-140-S1-TP2-fastp-adapters-polyG-31bp-merged.bam U   <21 120657
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/sRNA-ACR-140-S1-TP2-fastp-adapters-polyG-31bp-merged.bam U   21  43442
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/sRNA-ACR-140-S1-TP2-fastp-adapters-polyG-31bp-merged.bam U   22  106124
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/sRNA-ACR-140-S1-TP2-fastp-adapters-polyG-31bp-merged.bam U   23  65290
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/sRNA-ACR-140-S1-TP2-fastp-adapters-polyG-31bp-merged.bam U   24  84916
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/sRNA-ACR-140-S1-TP2-fastp-adapters-polyG-31bp-merged.bam U   >24 3225093
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/sRNA-ACR-140-S1-TP2-fastp-adapters-polyG-31bp-merged.bam P   <21 197276
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/sRNA-ACR-140-S1-TP2-fastp-adapters-polyG-31bp-merged.bam P   21  60763
    /home/shared/8TB_HDD_01/sam/gitrepos/deep-dive/D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/sRNA-ACR-140-S1-TP2-fastp-adapters-polyG-31bp-merged.bam P   22  73569

    ==> ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Counts.txt <==
    Coords  Name    MIRNA   sRNA-ACR-140-S1-TP2-fastp-adapters-polyG-31bp-merged    sRNA-ACR-145-S1-TP2-fastp-adapters-polyG-31bp-merged    sRNA-ACR-150-S1-TP2-fastp-adapters-polyG-31bp-merged    sRNA-ACR-173-S1-TP2-fastp-adapters-polyG-31bp-merged    sRNA-ACR-178-S1-TP2-fastp-adapters-polyG-31bp-merged
    NC_058066.1:152483-152910   Cluster_1   N   2   131 2   1   4
    NC_058066.1:161064-161674   Cluster_2   N   57  48  219 32  193
    NC_058066.1:172073-172496   Cluster_3   N   36  31  0   36  2
    NC_058066.1:203242-203651   Cluster_4   N   14  28  3   17  38
    NC_058066.1:204535-205150   Cluster_5   N   3   234 17  13  46
    NC_058066.1:205745-206966   Cluster_6   N   914 432 78  247 259
    NC_058066.1:210841-211344   Cluster_7   N   315 200 78  2   652
    NC_058066.1:349655-351297   Cluster_8   N   497 1236    757 767 22
    NC_058066.1:351491-353439   Cluster_9   N   3181    1995    1498    2030    185

    ==> ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/known_miRNAs.gff3 <==
    NC_058066.1 ShortStack  sRNA    172680  172697  0   +   .   ID=cgr-miR-1260
    NC_058066.1 ShortStack  sRNA    924864  924883  0   -   .   ID=mmu-miR-466i-5p
    NC_058066.1 ShortStack  sRNA    2130443 2130462 0   +   .   ID=mmu-miR-466i-5p
    NC_058066.1 ShortStack  sRNA    2130445 2130464 0   +   .   ID=mmu-miR-466i-5p
    NC_058066.1 ShortStack  sRNA    2130447 2130466 0   +   .   ID=mmu-miR-466i-5p
    NC_058066.1 ShortStack  sRNA    2130449 2130468 0   +   .   ID=mmu-miR-466i-5p
    NC_058066.1 ShortStack  sRNA    2130451 2130470 0   +   .   ID=mmu-miR-466i-5p
    NC_058066.1 ShortStack  sRNA    2130453 2130472 0   +   .   ID=mmu-miR-466i-5p
    NC_058066.1 ShortStack  sRNA    2424062 2424080 0   +   .   ID=gma-miR1533
    NC_058066.1 ShortStack  sRNA    2582838 2582857 0   +   .   ID=mmu-miR-466i-5p

    ==> ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/mir.fasta <==
    >Cluster_316::NC_058066.1:12757124-12757218(-)
    ATGCTTTACTCCTTTGGGAGGGAGGTTAGTGCAGAGGTCATCGTTATTGATGATCTCTGCAATAGCCTGCCTCCCAAAGGAGTTCTACTAGTCC
    >Cluster_316.mature::NC_058066.1:12757146-12757168(-)
    TGATCTCTGCAATAGCCTGCCT
    >Cluster_316.star::NC_058066.1:12757176-12757198(-)
    GGAGGTTAGTGCAGAGGTCATC
    >Cluster_514::NC_058066.1:20088629-20088720(+)
    TTGAATTGTCTGTGGCCTTCACGGCACCCCTTGCTGGAGTTATTTAATAACGCTAGGAAGGGATGCCGGGAAGGAGATGGTACAATGCAAA
    >Cluster_514.mature::NC_058066.1:20088678-20088700(+)
    ACGCTAGGAAGGGATGCCGGGA

    ==> ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results.gff3 <==
    NC_058066.1 ShortStack  Unknown_sRNA_locus  152483  152910  140 -   .   ID=Cluster_1;DicerCall=N;MIRNA=N
    NC_058066.1 ShortStack  Unknown_sRNA_locus  161064  161674  549 .   .   ID=Cluster_2;DicerCall=N;MIRNA=N
    NC_058066.1 ShortStack  Unknown_sRNA_locus  172073  172496  105 -   .   ID=Cluster_3;DicerCall=N;MIRNA=N
    NC_058066.1 ShortStack  Unknown_sRNA_locus  203242  203651  100 .   .   ID=Cluster_4;DicerCall=N;MIRNA=N
    NC_058066.1 ShortStack  Unknown_sRNA_locus  204535  205150  313 .   .   ID=Cluster_5;DicerCall=N;MIRNA=N
    NC_058066.1 ShortStack  Unknown_sRNA_locus  205745  206966  1930    .   .   ID=Cluster_6;DicerCall=N;MIRNA=N
    NC_058066.1 ShortStack  Unknown_sRNA_locus  210841  211344  1247    .   .   ID=Cluster_7;DicerCall=N;MIRNA=N
    NC_058066.1 ShortStack  Unknown_sRNA_locus  349655  351297  3279    +   .   ID=Cluster_8;DicerCall=N;MIRNA=N
    NC_058066.1 ShortStack  Unknown_sRNA_locus  351491  353439  8889    .   .   ID=Cluster_9;DicerCall=N;MIRNA=N
    NC_058066.1 ShortStack  Unknown_sRNA_locus  598651  599068  114 +   .   ID=Cluster_10;DicerCall=N;MIRNA=N

    ==> ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results.txt <==
    Locus   Name    Chrom   Start   End Length  Reads   DistinctSequences   FracTop Strand  MajorRNA    MajorRNAReads   Short   Long    21  22  23  24  DicerCall   MIRNA   known_miRNAs
    NC_058066.1:152483-152910   Cluster_1   NC_058066.1 152483  152910  428 140 32  0.05    -   UAAGUACUUUAUCAACUAACUCUAGGCA    75  1   130 0   2   0   7   N   N   NA
    NC_058066.1:161064-161674   Cluster_2   NC_058066.1 161064  161674  611 549 247 0.2987249544626594  .   UUUUAGCCUAGUGCGGGUUUCCAGACGU    43  25  479 16  4   4   21  N   N   NA
    NC_058066.1:172073-172496   Cluster_3   NC_058066.1 172073  172496  424 105 40  0.12380952380952381 -   GCGAUUAUUAACGGCUGGAACGACAGGCGA  16  1   88  1   1   0   14  N   N   NA
    NC_058066.1:203242-203651   Cluster_4   NC_058066.1 203242  203651  410 100 45  0.56    .   UUCUGACUCUAUUAGCAACGAAGACUUU    26  1   96  0   1   0   2   N   N   NA
    NC_058066.1:204535-205150   Cluster_5   NC_058066.1 204535  205150  616 313 157 0.7763578274760383  .   UCCCAACACGUCUAGACUGUACAAUUUCU   32  3   304 1   1   2   2   N   N   NA
    NC_058066.1:205745-206966   Cluster_6   NC_058066.1 205745  206966  1222    1930    416 0.35544041450777203 .   CAAAAGAGCGGACAAAAUAGUCGACAGAUU  716 3   1882    5   10  7   23  N   N   NA
    NC_058066.1:210841-211344   Cluster_7   NC_058066.1 210841  211344  504 1247    333 0.7457898957497995  .   UAAUACUUGUAGUGAAGGUUCAAUCUCGA   95  10  1133    7   7   20  70  N   N   NA
    NC_058066.1:349655-351297   Cluster_8   NC_058066.1 349655  351297  1643    3279    1165    0.8127477889600488  +   UCAGCUUGGAAAUGACAGCUUUUGACGU    255 27  3141    10  22  17  62  N   N   NA
    NC_058066.1:351491-353439   Cluster_9   NC_058066.1 351491  353439  1949    8889    1615    0.4114073574080324  .   UUUCAAAUCAAAGAUCUUCGCAACGAUGA   780 82  8503    34  34  114 122 N   N   NA

``` bash
wc -l ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/*
```

        151 ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/alignment_details.tsv
      18896 ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Counts.txt
       1446 ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/known_miRNAs.gff3
        228 ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/mir.fasta
      18971 ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results.gff3
      18896 ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results.txt
      58588 total

``` bash
grep -c ">" ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/*
```

    ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/alignment_details.tsv:25
    ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Counts.txt:0
    ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/known_miRNAs.gff3:0
    ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/mir.fasta:114
    ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results.gff3:0
    ../output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results.txt:0
