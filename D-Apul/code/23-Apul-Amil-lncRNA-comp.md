# 23-Apul Amil lncRNA compare
Steven Roberts
2025-07-19

``` bash
head ../output/05.33-lncRNA-discovery/Apul_lncRNA.bed
wc -l ../output/05.33-lncRNA-discovery/Apul_lncRNA.bed
```

    NC_058066.1 468617  469943
    NC_058066.1 1135314 1144814
    NC_058066.1 1144882 1148491
    NC_058066.1 1153397 1165634
    NC_058066.1 1153397 1165634
    NC_058066.1 1153402 1165634
    NC_058066.1 1153408 1165634
    NC_058066.1 1154205 1155609
    NC_058066.1 1155785 1165634
    NC_058066.1 1222538 1225166
    16206 ../output/05.33-lncRNA-discovery/Apul_lncRNA.bed

``` bash
grep -v -- '-[0-9]' ../output/05.33-lncRNA-discovery/Apul_lncRNA.bed > ../output/23-Apul-Amil-lncRNA-comp/Apul_lncRNA.no-negatives.bed
```

``` bash
cat ../data/Amil/GCF_013753865.1_Amil_v2.1_genomic.gff | \
awk '$0 !~ /^#/ {count[$3]++} END {for (f in count) print f, count[f]}' | sort
```

    cDNA_match 22387
    CDS 317969
    exon 390533
    gene 36904
    guide_RNA 1
    lnc_RNA 6128
    mRNA 41860
    pseudogene 5871
    region 854
    rRNA 283
    snoRNA 62
    snRNA 170
    transcript 2066
    tRNA 1413

``` bash
grep -i "lnc_RNA" ../data/Amil/GCF_013753865.1_Amil_v2.1_genomic.gff > ../output/23-Apul-Amil-lncRNAlncRNA.gff
```

``` bash
head ../output/23-Apul-Amil-lncRNA-comp/lncRNA.gff
wc -l ../output/23-Apul-Amil-lncRNA-comp/lncRNA.gff
```

    NC_058066.1 Gnomon  lnc_RNA 1962    23221   .   -   .   ID=rna-XR_003825913.2;Parent=gene-LOC114963522;Dbxref=GeneID:114963522,Genbank:XR_003825913.2;Name=XR_003825913.2;gbkey=ncRNA;gene=LOC114963522;model_evidence=Supporting evidence includes similarity to: 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 2 samples with support for all annotated introns;product=uncharacterized LOC114963522;transcript_id=XR_003825913.2
    NC_058066.1 Gnomon  lnc_RNA 97370   99364   .   +   .   ID=rna-XR_003827594.2;Parent=gene-LOC114972404;Dbxref=GeneID:114972404,Genbank:XR_003827594.2;Name=XR_003827594.2;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=ncRNA;gene=LOC114972404;model_evidence=Supporting evidence includes similarity to: 3 mRNAs%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 56 samples with support for all annotated introns;product=uncharacterized LOC114972404;transcript_id=XR_003827594.2
    NC_058066.1 Gnomon  lnc_RNA 196438  197763  .   +   .   ID=rna-XR_006394789.1;Parent=gene-LOC122957587;Dbxref=GeneID:122957587,Genbank:XR_006394789.1;Name=XR_006394789.1;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=ncRNA;gene=LOC122957587;model_evidence=Supporting evidence includes similarity to: 1 mRNA%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 77 samples with support for all annotated introns;product=uncharacterized LOC122957587;transcript_id=XR_006394789.1
    NC_058066.1 Gnomon  lnc_RNA 407057  411847  .   +   .   ID=rna-XR_003823912.2;Parent=gene-LOC114952940;Dbxref=GeneID:114952940,Genbank:XR_003823912.2;Name=XR_003823912.2;gbkey=ncRNA;gene=LOC114952940;model_evidence=Supporting evidence includes similarity to: 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 1 sample with support for all annotated introns;product=uncharacterized LOC114952940;transcript_id=XR_003823912.2
    NC_058066.1 Gnomon  lnc_RNA 475161  486418  .   +   .   ID=rna-XR_003823898.2;Parent=gene-LOC114952897;Dbxref=GeneID:114952897,Genbank:XR_003823898.2;Name=XR_003823898.2;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=ncRNA;gene=LOC114952897;model_evidence=Supporting evidence includes similarity to: 11 mRNAs%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 81 samples with support for all annotated introns;product=uncharacterized LOC114952897;transcript_id=XR_003823898.2
    NC_058066.1 Gnomon  lnc_RNA 486421  493748  .   -   .   ID=rna-XR_003823897.2;Parent=gene-LOC114952896;Dbxref=GeneID:114952896,Genbank:XR_003823897.2;Name=XR_003823897.2;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=ncRNA;gene=LOC114952896;model_evidence=Supporting evidence includes similarity to: 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 14 samples with support for all annotated introns;product=uncharacterized LOC114952896%2C transcript variant X2;transcript_id=XR_003823897.2
    NC_058066.1 Gnomon  lnc_RNA 486421  493745  .   -   .   ID=rna-XR_003823896.2;Parent=gene-LOC114952896;Dbxref=GeneID:114952896,Genbank:XR_003823896.2;Name=XR_003823896.2;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=ncRNA;gene=LOC114952896;model_evidence=Supporting evidence includes similarity to: 2 mRNAs%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 90 samples with support for all annotated introns;product=uncharacterized LOC114952896%2C transcript variant X1;transcript_id=XR_003823896.2
    NC_058066.1 Gnomon  lnc_RNA 516911  518397  .   +   .   ID=rna-XR_006394821.1;Parent=gene-LOC114952899;Dbxref=GeneID:114952899,Genbank:XR_006394821.1;Name=XR_006394821.1;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=ncRNA;gene=LOC114952899;model_evidence=Supporting evidence includes similarity to: 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 6 samples with support for all annotated introns;product=uncharacterized LOC114952899;transcript_id=XR_006394821.1
    NC_058066.1 Gnomon  lnc_RNA 536049  540472  .   -   .   ID=rna-XR_006394820.1;Parent=gene-LOC122957618;Dbxref=GeneID:122957618,Genbank:XR_006394820.1;Name=XR_006394820.1;experiment=COORDINATES: polyA evidence [ECO:0006239];gbkey=ncRNA;gene=LOC122957618;model_evidence=Supporting evidence includes similarity to: 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 9 samples with support for all annotated introns;product=uncharacterized LOC122957618;transcript_id=XR_006394820.1
    NC_058066.1 Gnomon  lnc_RNA 540506  556085  .   -   .   ID=rna-XR_003823899.2;Parent=gene-LOC114952898;Dbxref=GeneID:114952898,Genbank:XR_003823899.2;Name=XR_003823899.2;gbkey=ncRNA;gene=LOC114952898;model_evidence=Supporting evidence includes similarity to: 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 15 samples with support for all annotated introns;product=uncharacterized LOC114952898;transcript_id=XR_003823899.2
    6128 ../output/23-Apul-Amil-lncRNA-comp/lncRNA.gff

``` bash
wc -l ../output/23-Apul-Amil-lncRNA-comp/Apul_lncRNA.no-negatives.bed
```

    16201 ../output/23-Apul-Amil-lncRNA-comp/Apul_lncRNA.no-negatives.bed

``` bash
bedtools intersect \
-a ../output/23-Apul-Amil-lncRNA-comp/Apul_lncRNA.no-negatives.bed \
-b ../output/23-Apul-Amil-lncRNA-comp/lncRNA.gff
```

    NC_058066.1 1135314 1136572
    NC_058066.1 17109219    17112344
    NC_058066.1 17109219    17112344
    NC_058066.1 17112380    17117165
    NC_058066.1 17112384    17117165
    NC_058066.1 17112384    17117165
    NC_058066.1 17112384    17117165
    NC_058066.1 17112384    17117165
    NC_058066.1 17112384    17116095
    NC_058066.1 17113007    17120760
    NC_058066.1 17113007    17120760
    NC_058066.1 17113007    17118453
    NC_058066.1 17113007    17118452
    NC_058066.1 17113007    17117845
    NC_058066.1 17113007    17116095
    NC_058066.1 17123556    17125188
    NC_058066.1 17110302    17112344
    NC_058066.1 17110302    17112344
    NC_058066.1 17112380    17120760
    NC_058066.1 17112384    17120760
    NC_058066.1 17112384    17118453
    NC_058066.1 17112384    17118452
    NC_058066.1 17112384    17117845
    NC_058066.1 17112384    17116095
    NC_058066.1 17123556    17125188
    NC_058066.1 17123556    17124972
    NC_058066.1 7062475 7065195
    NC_058066.1 8321013 8323206
    NC_058066.1 28992922    28993631
    NC_058067.1 13320337    13322742
    NC_058067.1 31534387    31536144
    NC_058067.1 31540314    31542285
    NC_058068.1 319076  323165
    NC_058068.1 319076  321021
    NC_058068.1 322238  323165
    NC_058068.1 322844  323165
    NC_058068.1 3026486 3029014
    NC_058068.1 3026486 3028994
    NC_058068.1 3026486 3028478
    NC_058068.1 3032927 3036783
    NC_058068.1 3036812 3040669
    NC_058068.1 3029041 3032898
    NC_058068.1 3029041 3032453
    NC_058068.1 3886678 3889941
    NC_058068.1 3886678 3889100
    NC_058068.1 3893757 3894856
    NC_058068.1 3893757 3894856
    NC_058068.1 3893757 3894856
    NC_058068.1 7652181 7653209
    NC_058068.1 8241987 8245462
    NC_058068.1 9544021 9545810
    NC_058068.1 8542539 8542540
    NC_058068.1 8542539 8542540
    NC_058068.1 8542539 8542540
    NC_058069.1 3399469 3400117
    NC_058069.1 3412085 3414170
    NC_058069.1 3399469 3400117
    NC_058069.1 3412085 3414170
    NC_058069.1 3412085 3414170
    NC_058069.1 3412085 3414170
    NC_058069.1 13244095    13244870
    NC_058069.1 22939993    22940885
    NC_058069.1 22940902    22952293
    NC_058069.1 22941485    22952293
    NC_058069.1 22942099    22952293
    NC_058069.1 22940902    22952293
    NC_058069.1 25495577    25497139
    NC_058070.1 4935615 4938195
    NC_058070.1 4935615 4938195
    NC_058070.1 4263104 4271536
    NC_058071.1 4641288 4644865
    NC_058071.1 4641290 4644865
    NC_058071.1 4641292 4644865
    NC_058071.1 19823711    19825047
    NC_058072.1 15531281    15532457
    NC_058072.1 15531281    15532457
    NC_058073.1 10767410    10770132
    NC_058073.1 10767410    10770132
    NC_058073.1 10767410    10770132
    NC_058074.1 4791534 4798294
    NC_058074.1 4791534 4798294
    NC_058074.1 4791534 4798294
    NC_058074.1 4791534 4798294
    NC_058074.1 8569313 8574051
    NC_058075.1 13051322    13052645
    NC_058075.1 13051322    13052645
    NC_058075.1 13051322    13052645
    NC_058075.1 17724946    17725407
    NC_058076.1 13591680    13593520
    NC_058077.1 16247679    16249210
    NC_058077.1 1285371 1287536
    NC_058078.1 1552006 1552801
    NC_058079.1 557889  558200
    NC_058079.1 546204  558200
    NW_025322618.1  947160  949000
    NW_025322629.1  662140  662447
    NW_025322644.1  131508  132600
    NW_025322719.1  29162   29891
    NW_025322771.1  492953  494210
    NW_025322771.1  488550  492928
    NW_025322771.1  509000  509643
    NW_025322813.1  46492   47401
    NW_025322847.1  85586   87924
    NW_025322879.1  24056   25418
    NW_025322879.1  15446   17268
    NW_025322914.1  49224   51766
    NW_025322914.1  49224   51766
    NW_025322914.1  49224   51766
    NW_025323135.1  129760  130906

``` bash
bedtools intersect \
-a ../output/23-Apul-Amil-lncRNA-comp/Apul_lncRNA.no-negatives.bed \
-b ../output/23-Apul-Amil-lncRNA-comp/lncRNA.gff \
-u | wc -l
```

    59

``` bash
bedtools intersect \
-a ../output/23-Apul-Amil-lncRNA-comp/lncRNA.gff \
-b ../output/23-Apul-Amil-lncRNA-comp/Apul_lncRNA.no-negatives.bed \
-u | wc -l
```

    79
