# Looking at genomes

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

    grep '>' ../data/Amil/ncbi_dataset/data/GCF_013753865.1/GCF_013753865.1_Amil_v2.1_genomic.fna | head -40

    ## >NC_058066.1 Acropora millepora isolate JS-1 chromosome 1, Amil_v2.1, whole genome shotgun sequence
    ## >NC_058067.1 Acropora millepora isolate JS-1 chromosome 2, Amil_v2.1, whole genome shotgun sequence
    ## >NC_058068.1 Acropora millepora isolate JS-1 chromosome 3, Amil_v2.1, whole genome shotgun sequence
    ## >NC_058069.1 Acropora millepora isolate JS-1 chromosome 4, Amil_v2.1, whole genome shotgun sequence
    ## >NC_058070.1 Acropora millepora isolate JS-1 chromosome 5, Amil_v2.1, whole genome shotgun sequence
    ## >NC_058071.1 Acropora millepora isolate JS-1 chromosome 6, Amil_v2.1, whole genome shotgun sequence
    ## >NC_058072.1 Acropora millepora isolate JS-1 chromosome 7, Amil_v2.1, whole genome shotgun sequence
    ## >NC_058073.1 Acropora millepora isolate JS-1 chromosome 8, Amil_v2.1, whole genome shotgun sequence
    ## >NC_058074.1 Acropora millepora isolate JS-1 chromosome 9, Amil_v2.1, whole genome shotgun sequence
    ## >NC_058075.1 Acropora millepora isolate JS-1 chromosome 10, Amil_v2.1, whole genome shotgun sequence
    ## >NC_058076.1 Acropora millepora isolate JS-1 chromosome 11, Amil_v2.1, whole genome shotgun sequence
    ## >NC_058077.1 Acropora millepora isolate JS-1 chromosome 12, Amil_v2.1, whole genome shotgun sequence
    ## >NC_058078.1 Acropora millepora isolate JS-1 chromosome 13, Amil_v2.1, whole genome shotgun sequence
    ## >NC_058079.1 Acropora millepora isolate JS-1 chromosome 14, Amil_v2.1, whole genome shotgun sequence
    ## >NW_025322615.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000015, whole genome shotgun sequence
    ## >NW_025322616.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000028, whole genome shotgun sequence
    ## >NW_025322617.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000030, whole genome shotgun sequence
    ## >NW_025322618.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000038, whole genome shotgun sequence
    ## >NW_025322619.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000050, whole genome shotgun sequence
    ## >NW_025322620.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000055, whole genome shotgun sequence
    ## >NW_025322621.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000059, whole genome shotgun sequence
    ## >NW_025322622.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000061, whole genome shotgun sequence
    ## >NW_025322623.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000066, whole genome shotgun sequence
    ## >NW_025322624.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000070, whole genome shotgun sequence
    ## >NW_025322625.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000072, whole genome shotgun sequence
    ## >NW_025322626.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000077, whole genome shotgun sequence
    ## >NW_025322627.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000080, whole genome shotgun sequence
    ## >NW_025322628.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000086, whole genome shotgun sequence
    ## >NW_025322629.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000089, whole genome shotgun sequence
    ## >NW_025322630.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000091, whole genome shotgun sequence
    ## >NW_025322631.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000094, whole genome shotgun sequence
    ## >NW_025322632.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000095, whole genome shotgun sequence
    ## >NW_025322633.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000098, whole genome shotgun sequence
    ## >NW_025322634.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000102, whole genome shotgun sequence
    ## >NW_025322635.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000104, whole genome shotgun sequence
    ## >NW_025322636.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000107, whole genome shotgun sequence
    ## >NW_025322637.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000109, whole genome shotgun sequence
    ## >NW_025322638.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000115, whole genome shotgun sequence
    ## >NW_025322639.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000118, whole genome shotgun sequence
    ## >NW_025322640.1 Acropora millepora isolate JS-1 unplaced genomic scaffold, Amil_v2.1 Sc0000120, whole genome shotgun sequence
