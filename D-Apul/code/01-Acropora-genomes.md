# Looking at genomes

<https://www.ncbi.nlm.nih.gov/data-hub/genome/?taxon=6127>

![ncbi](https://gannet.fish.washington.edu/seashell/snaps/Genome_-_NCBI_-_NLM_2023-03-02_11-52-28.png)

## millipora

    head ../data/Amil/ncbi_dataset/data/GCF_013753865.1/GCF_013753865.1_Amil_v2.1_genomic.fna

    ## >NC_058066.1 Acropora millepora isolate JS-1 chromosome 1, Amil_v2.1, whole genome shotgun sequence
    ## CCATTTGACTGACTGAGCTGTAATTCGGGCAGTTAGCCCGCTGACGTCAACGTTAGAAAATGCAGTGATAACGCACGCAG
    ## AGTCGTAGACATTAAGAGGGAGGTAAAACTGAAGACGCTGGTACTCTAAAGCAAGCAGCAGCTACCAAACAGGAAACATA
    ## AATACAGTGCCGTTTTCCAGCTTTTATCTCTTTCATTTCCTCATTACTTCTAACAGTAATCGTAACCTTATTTCTCAGTG
    ## AGAGTAAACAGGATTTTAAACTTTTATCGTTATCTTAGGGTTGGGACCTTACTCCATGTCAGACTCAGTTAACTGCCTGT
    ## GAGTGTTTGCGCGCGCTACAGCTACTGGACGAGTTAATTTTGCATACTCCATTTTGCTGcgtgctttttttttatgtttt
    ## atttttatagtAATGCAATGAATCATTGGTCTGTTGCAAGATCATGAATTAACAGTCAAAAcaataacatgtttttcatt
    ## CAGTCAGGCGTTTCGAACGACCActctaccaactgagccaaTTCATTAGAAAAAATGTATGAAcacctaaaccctaaccc
    ## ttaccCTAACACCATTAACTACCAGACTGCTATCAAATTGAATGATATATAGCCCTTACGGAACTACGCCTGAAGCTACC
    ## TTTGAATTTCACAATTGCCGCTCCAACGTTTTATTACTTAGCATTCAATCTCATTCAGGCTGCAGCTTTAGTCCTGCTTC

    grep '>' ../data/Amil/ncbi_dataset/data/GCF_013753865.1/GCF_013753865.1_Amil_v2.1_genomic.fna | wc -l

    ## 854

    grep '>' ../data/Adig/ncbi_dataset/data/GCF_000222465.1/GCF_000222465.1_Adig_1.1_genomic.fna | head -40

    ## >NW_015441057.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970692.1, whole genome shotgun sequence
    ## >NW_015441058.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970693.1, whole genome shotgun sequence
    ## >NW_015441059.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970694.1, whole genome shotgun sequence
    ## >NW_015441060.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970695.1, whole genome shotgun sequence
    ## >NW_015441061.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970696.1, whole genome shotgun sequence
    ## >NW_015441062.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970697.1, whole genome shotgun sequence
    ## >NW_015441063.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970698.1, whole genome shotgun sequence
    ## >NW_015441064.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970699.1, whole genome shotgun sequence
    ## >NW_015441065.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970700.1, whole genome shotgun sequence
    ## >NW_015441066.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970701.1, whole genome shotgun sequence
    ## >NW_015441067.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970702.1, whole genome shotgun sequence
    ## >NW_015441068.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970703.1, whole genome shotgun sequence
    ## >NW_015441069.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970704.1, whole genome shotgun sequence
    ## >NW_015441070.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970705.1, whole genome shotgun sequence
    ## >NW_015441071.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970706.1, whole genome shotgun sequence
    ## >NW_015441072.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970707.1, whole genome shotgun sequence
    ## >NW_015441073.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970708.1, whole genome shotgun sequence
    ## >NW_015441074.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970709.1, whole genome shotgun sequence
    ## >NW_015441075.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970710.1, whole genome shotgun sequence
    ## >NW_015441076.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970711.1, whole genome shotgun sequence
    ## >NW_015441077.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970712.1, whole genome shotgun sequence
    ## >NW_015441078.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970713.1, whole genome shotgun sequence
    ## >NW_015441079.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970714.1, whole genome shotgun sequence
    ## >NW_015441080.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970715.1, whole genome shotgun sequence
    ## >NW_015441081.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970716.1, whole genome shotgun sequence
    ## >NW_015441082.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970717.1, whole genome shotgun sequence
    ## >NW_015441083.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970718.1, whole genome shotgun sequence
    ## >NW_015441084.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970719.1, whole genome shotgun sequence
    ## >NW_015441085.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970720.1, whole genome shotgun sequence
    ## >NW_015441086.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970721.1, whole genome shotgun sequence
    ## >NW_015441087.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970722.1, whole genome shotgun sequence
    ## >NW_015441088.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970723.1, whole genome shotgun sequence
    ## >NW_015441089.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970724.1, whole genome shotgun sequence
    ## >NW_015441090.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970725.1, whole genome shotgun sequence
    ## >NW_015441091.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970726.1, whole genome shotgun sequence
    ## >NW_015441092.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970727.1, whole genome shotgun sequence
    ## >NW_015441093.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970728.1, whole genome shotgun sequence
    ## >NW_015441094.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970729.1, whole genome shotgun sequence
    ## >NW_015441095.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970730.1, whole genome shotgun sequence
    ## >NW_015441096.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970731.1, whole genome shotgun sequence

## A digitera

    head ../data/Adig/ncbi_dataset/data/GCF_000222465.1/GCF_000222465.1_Adig_1.1_genomic.fna

    ## >NW_015441057.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970692.1, whole genome shotgun sequence
    ## tttgaaaactgctctaattgctttggttttggtttttcgacagtcatttgaaaaccgctctattataACAATTCTTACGA
    ## TTACGTTTTTCGATTCACATTGGATTTTAACGCTCGTTGTTTTTGTCGTTTAGTTGTCAATGACTTTTAAGGTAACTTTC
    ## CTAGTTTTGAAAATGTAATCTtcaggcctcagttgttcaaagggtgggtaAGTGCTGTCcgctggataaatcactatcca
    ## gtggaCGAGTGCTATCGAAATCATTTGCGTTacccagtggatagtgatttatccaatggatattgAACGACTGGGGCCAG
    ## AAGTGTTCTGTAGCATCTTTTAATTGGCTATTTCTCCTTGTAGACTTACGCACCTTTACGTTAGTTAATCGGTgcaaagt
    ## ttgaatcgagtttcattttcttcaggTGTGATGGTGGGGGATGAAATCATTGCAGTTAATGACATTGATGTAACAGAAGT
    ## AGAGAATGCTGTGGAAGATCTGAGGGAAGCTTTGAAAGGTGAGGATGCCATCTAGTGAATCAATAATAGCACTGGATACA
    ## TTTCATTTCCAGATTTATGAAGGCTCGGAGAGAGGTTATTATAAATATAGTTGATGTTCTGTAAGACCTCAGTAATTTTG
    ## GCTATATTATGTAAGGGTAGTTTTCTTACAAAATCTGCAGGCCTTGATGACTATATAAGGTTTTGAGATGCATTGAGTTT

    grep '>' ../data/Adig/ncbi_dataset/data/GCF_000222465.1/GCF_000222465.1_Adig_1.1_genomic.fna | wc -l

    ## 2421

    grep '>' ../data/Adig/ncbi_dataset/data/GCF_000222465.1/GCF_000222465.1_Adig_1.1_genomic.fna | head -40

    ## >NW_015441057.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970692.1, whole genome shotgun sequence
    ## >NW_015441058.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970693.1, whole genome shotgun sequence
    ## >NW_015441059.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970694.1, whole genome shotgun sequence
    ## >NW_015441060.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970695.1, whole genome shotgun sequence
    ## >NW_015441061.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970696.1, whole genome shotgun sequence
    ## >NW_015441062.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970697.1, whole genome shotgun sequence
    ## >NW_015441063.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970698.1, whole genome shotgun sequence
    ## >NW_015441064.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970699.1, whole genome shotgun sequence
    ## >NW_015441065.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970700.1, whole genome shotgun sequence
    ## >NW_015441066.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970701.1, whole genome shotgun sequence
    ## >NW_015441067.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970702.1, whole genome shotgun sequence
    ## >NW_015441068.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970703.1, whole genome shotgun sequence
    ## >NW_015441069.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970704.1, whole genome shotgun sequence
    ## >NW_015441070.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970705.1, whole genome shotgun sequence
    ## >NW_015441071.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970706.1, whole genome shotgun sequence
    ## >NW_015441072.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970707.1, whole genome shotgun sequence
    ## >NW_015441073.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970708.1, whole genome shotgun sequence
    ## >NW_015441074.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970709.1, whole genome shotgun sequence
    ## >NW_015441075.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970710.1, whole genome shotgun sequence
    ## >NW_015441076.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970711.1, whole genome shotgun sequence
    ## >NW_015441077.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970712.1, whole genome shotgun sequence
    ## >NW_015441078.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970713.1, whole genome shotgun sequence
    ## >NW_015441079.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970714.1, whole genome shotgun sequence
    ## >NW_015441080.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970715.1, whole genome shotgun sequence
    ## >NW_015441081.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970716.1, whole genome shotgun sequence
    ## >NW_015441082.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970717.1, whole genome shotgun sequence
    ## >NW_015441083.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970718.1, whole genome shotgun sequence
    ## >NW_015441084.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970719.1, whole genome shotgun sequence
    ## >NW_015441085.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970720.1, whole genome shotgun sequence
    ## >NW_015441086.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970721.1, whole genome shotgun sequence
    ## >NW_015441087.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970722.1, whole genome shotgun sequence
    ## >NW_015441088.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970723.1, whole genome shotgun sequence
    ## >NW_015441089.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970724.1, whole genome shotgun sequence
    ## >NW_015441090.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970725.1, whole genome shotgun sequence
    ## >NW_015441091.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970726.1, whole genome shotgun sequence
    ## >NW_015441092.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970727.1, whole genome shotgun sequence
    ## >NW_015441093.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970728.1, whole genome shotgun sequence
    ## >NW_015441094.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970729.1, whole genome shotgun sequence
    ## >NW_015441095.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970730.1, whole genome shotgun sequence
    ## >NW_015441096.1 Acropora digitifera unplaced genomic scaffold, Adig_1.1 DF970731.1, whole genome shotgun sequence

## A hyacinthus

    head ../data/Ahya/ncbi_dataset/data/GCA_020536085.1/GCA_020536085.1_Ahyacinthus.chrsV1_genomic.fna

    ## >CM035865.1 Acropora hyacinthus isolate CA74 chromosome 1, whole genome shotgun sequence
    ## GATTAGCATCTGCTATCTGTGGGCATGTCCAATTTGTTGACAATACTTGCTAGAATTAGAAAAATAAAGACACAGATAAA
    ## GACACATAAGTGAGGGAATGTTCTGAAATTAGTCATAATTGTTGTGAAAGCACTGTCAAAACTAAATGAATGTGATCTTT
    ## CTATAGTGCACCAGTGTGATATGAAAGGCTGCAAAGAGAACCTTGTTATTAATCGAAGTTTTGACAATGTCACATGGAGT
    ## GCTGCTTGGCAAAGGAAAGATTAGTGGGTAATATTGAATATGATTCCCTTCCTGGACTTCCTGGGCAGATCAGAACGGGC
    ## TGCATAAGAACGGCAAAACTTGGGAGCCAATCATGTGAAGACCATCACCAAGATCACAAAGAATTAGACAAGGTTGGTAT
    ## AAACTGAAAGCAGAAACGTAGGTGATATAtgtattatattatattatattatattaataatatcaCTCAAGCATGCATTT
    ## AGGAGCTGTTAGAGGATAGTTCAGGACTTGTGTAACTGTGCAATGAATGCTAAAATCGGAGCATCAGATATGAAATGGGA
    ## gttcaaatttttgttatcgttactaaattaaatgtttttgtCAGTGTAACAAAAACAATGGCATTATTTGTACTCATCAT
    ## GGTCGTAATGGTGATACTAATAATGTTCCATGACAGGTGCATGAACTCGAtgtaaaagaaattgaaaaggacTTAGGTGA

    grep '>' ../data/Ahya/ncbi_dataset/data/GCA_020536085.1/GCA_020536085.1_Ahyacinthus.chrsV1_genomic.fna | wc -l

    ## 907

    grep '>' ../data/Ahya/ncbi_dataset/data/GCA_020536085.1/GCA_020536085.1_Ahyacinthus.chrsV1_genomic.fna | head -40

    ## >CM035865.1 Acropora hyacinthus isolate CA74 chromosome 1, whole genome shotgun sequence
    ## >CM035866.1 Acropora hyacinthus isolate CA74 chromosome 2, whole genome shotgun sequence
    ## >CM035867.1 Acropora hyacinthus isolate CA74 chromosome 3, whole genome shotgun sequence
    ## >CM035868.1 Acropora hyacinthus isolate CA74 chromosome 4, whole genome shotgun sequence
    ## >CM035869.1 Acropora hyacinthus isolate CA74 chromosome 5, whole genome shotgun sequence
    ## >CM035870.1 Acropora hyacinthus isolate CA74 chromosome 6, whole genome shotgun sequence
    ## >CM035871.1 Acropora hyacinthus isolate CA74 chromosome 7, whole genome shotgun sequence
    ## >CM035872.1 Acropora hyacinthus isolate CA74 chromosome 8, whole genome shotgun sequence
    ## >CM035873.1 Acropora hyacinthus isolate CA74 chromosome 9, whole genome shotgun sequence
    ## >CM035874.1 Acropora hyacinthus isolate CA74 chromosome 10, whole genome shotgun sequence
    ## >CM035875.1 Acropora hyacinthus isolate CA74 chromosome 11, whole genome shotgun sequence
    ## >CM035876.1 Acropora hyacinthus isolate CA74 chromosome 12, whole genome shotgun sequence
    ## >CM035877.1 Acropora hyacinthus isolate CA74 chromosome 13, whole genome shotgun sequence
    ## >CM035878.1 Acropora hyacinthus isolate CA74 chromosome 14, whole genome shotgun sequence
    ## >JAIFHZ010000015.1 Acropora hyacinthus isolate CA74 sc001, whole genome shotgun sequence
    ## >JAIFHZ010000016.1 Acropora hyacinthus isolate CA74 sc002, whole genome shotgun sequence
    ## >JAIFHZ010000017.1 Acropora hyacinthus isolate CA74 sc003, whole genome shotgun sequence
    ## >JAIFHZ010000018.1 Acropora hyacinthus isolate CA74 sc004, whole genome shotgun sequence
    ## >JAIFHZ010000019.1 Acropora hyacinthus isolate CA74 sc005, whole genome shotgun sequence
    ## >JAIFHZ010000020.1 Acropora hyacinthus isolate CA74 sc006, whole genome shotgun sequence
    ## >JAIFHZ010000021.1 Acropora hyacinthus isolate CA74 sc007, whole genome shotgun sequence
    ## >JAIFHZ010000022.1 Acropora hyacinthus isolate CA74 sc008, whole genome shotgun sequence
    ## >JAIFHZ010000023.1 Acropora hyacinthus isolate CA74 sc009, whole genome shotgun sequence
    ## >JAIFHZ010000024.1 Acropora hyacinthus isolate CA74 sc010, whole genome shotgun sequence
    ## >JAIFHZ010000025.1 Acropora hyacinthus isolate CA74 sc011, whole genome shotgun sequence
    ## >JAIFHZ010000026.1 Acropora hyacinthus isolate CA74 sc012, whole genome shotgun sequence
    ## >JAIFHZ010000027.1 Acropora hyacinthus isolate CA74 sc013, whole genome shotgun sequence
    ## >JAIFHZ010000028.1 Acropora hyacinthus isolate CA74 sc014, whole genome shotgun sequence
    ## >JAIFHZ010000029.1 Acropora hyacinthus isolate CA74 sc015, whole genome shotgun sequence
    ## >JAIFHZ010000030.1 Acropora hyacinthus isolate CA74 sc016, whole genome shotgun sequence
    ## >JAIFHZ010000031.1 Acropora hyacinthus isolate CA74 sc017, whole genome shotgun sequence
    ## >JAIFHZ010000032.1 Acropora hyacinthus isolate CA74 sc018, whole genome shotgun sequence
    ## >JAIFHZ010000033.1 Acropora hyacinthus isolate CA74 sc019, whole genome shotgun sequence
    ## >JAIFHZ010000034.1 Acropora hyacinthus isolate CA74 sc020, whole genome shotgun sequence
    ## >JAIFHZ010000035.1 Acropora hyacinthus isolate CA74 sc021, whole genome shotgun sequence
    ## >JAIFHZ010000036.1 Acropora hyacinthus isolate CA74 sc022, whole genome shotgun sequence
    ## >JAIFHZ010000037.1 Acropora hyacinthus isolate CA74 sc023, whole genome shotgun sequence
    ## >JAIFHZ010000038.1 Acropora hyacinthus isolate CA74 sc024, whole genome shotgun sequence
    ## >JAIFHZ010000039.1 Acropora hyacinthus isolate CA74 sc025, whole genome shotgun sequence
    ## >JAIFHZ010000040.1 Acropora hyacinthus isolate CA74 sc026, whole genome shotgun sequence
