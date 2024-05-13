15-Hisat
================
Steven Roberts
13 May, 2024

- <a href="#1-run-hisat-on-rna-seq" id="toc-1-run-hisat-on-rna-seq">1 Run
  HiSat on RNA-seq</a>
  - <a href="#11-grab-trimmed-rna-seq-reads"
    id="toc-11-grab-trimmed-rna-seq-reads">1.1 Grab Trimmed RNA-seq
    Reads</a>
  - <a href="#12-genome" id="toc-12-genome">1.2 Genome</a>
  - <a href="#13-hisat" id="toc-13-hisat">1.3 HiSat</a>
  - <a href="#14-convert-to-bams" id="toc-14-convert-to-bams">1.4 convert to
    bams</a>
  - <a href="#15-looking-at-bams" id="toc-15-looking-at-bams">1.5 Looking at
    Bams</a>
- <a href="#2-stringtie" id="toc-2-stringtie">2 StringTie</a>
- <a href="#3-gffcompare" id="toc-3-gffcompare">3 GFFcompare</a>
- <a href="#4-filter" id="toc-4-filter">4 Filter</a>
- <a href="#5-bedtools" id="toc-5-bedtools">5 Bedtools</a>
- <a href="#6-cpc2" id="toc-6-cpc2">6 CPC2</a>
- <a href="#7-subsetting-fasta" id="toc-7-subsetting-fasta">7 subsetting
  fasta</a>
- <a href="#8-getting-genome-feature-track"
  id="toc-8-getting-genome-feature-track">8 Getting genome feature
  track</a>

# 1 Run HiSat on RNA-seq

Will end up with 5 sorted bam files.

## 1.1 Grab Trimmed RNA-seq Reads

``` bash

ls ../data/fastq/ 
```

    RNA-ACR-140-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    RNA-ACR-140-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    RNA-ACR-145-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    RNA-ACR-145-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    RNA-ACR-150-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    RNA-ACR-150-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    RNA-ACR-173-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    RNA-ACR-173-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    RNA-ACR-178-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    RNA-ACR-178-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz

## 1.2 Genome

``` bash
ls ../data/GCF_013753865.1_Amil_v2.1_genomic.fna
```

    ../data/GCF_013753865.1_Amil_v2.1_genomic.fna

``` bash
head ../data/Amil/ncbi_dataset/data/GCF_013753865.1/genomic.gtf

wc -l ../data/Amil/ncbi_dataset/data/GCF_013753865.1/genomic.gtf

grep -v '^#' ../data/Amil/ncbi_dataset/data/GCF_013753865.1/genomic.gtf | cut -f3 | sort | uniq

grep -c "transcript" ../data/Amil/ncbi_dataset/data/GCF_013753865.1/genomic.gtf

grep -c "gene" ../data/Amil/ncbi_dataset/data/GCF_013753865.1/genomic.gtf
```

    #gtf-version 2.2
    #!genome-build Amil_v2.1
    #!genome-build-accession NCBI_Assembly:GCF_013753865.1
    #!annotation-source NCBI Acropora millepora Annotation Release 101
    NC_058066.1 Gnomon  gene    1962    23221   .   -   .   gene_id "LOC114963522"; transcript_id ""; db_xref "GeneID:114963522"; gbkey "Gene"; gene "LOC114963522"; gene_biotype "lncRNA"; 
    NC_058066.1 Gnomon  transcript  1962    23221   .   -   .   gene_id "LOC114963522"; transcript_id "XR_003825913.2"; db_xref "GeneID:114963522"; gbkey "ncRNA"; gene "LOC114963522"; model_evidence "Supporting evidence includes similarity to: 100% coverage of the annotated genomic feature by RNAseq alignments, including 2 samples with support for all annotated introns"; product "uncharacterized LOC114963522"; transcript_biotype "lnc_RNA"; 
    NC_058066.1 Gnomon  exon    23085   23221   .   -   .   gene_id "LOC114963522"; transcript_id "XR_003825913.2"; db_xref "GeneID:114963522"; gene "LOC114963522"; model_evidence "Supporting evidence includes similarity to: 100% coverage of the annotated genomic feature by RNAseq alignments, including 2 samples with support for all annotated introns"; product "uncharacterized LOC114963522"; transcript_biotype "lnc_RNA"; exon_number "1"; 
    NC_058066.1 Gnomon  exon    21001   21093   .   -   .   gene_id "LOC114963522"; transcript_id "XR_003825913.2"; db_xref "GeneID:114963522"; gene "LOC114963522"; model_evidence "Supporting evidence includes similarity to: 100% coverage of the annotated genomic feature by RNAseq alignments, including 2 samples with support for all annotated introns"; product "uncharacterized LOC114963522"; transcript_biotype "lnc_RNA"; exon_number "2"; 
    NC_058066.1 Gnomon  exon    18711   18775   .   -   .   gene_id "LOC114963522"; transcript_id "XR_003825913.2"; db_xref "GeneID:114963522"; gene "LOC114963522"; model_evidence "Supporting evidence includes similarity to: 100% coverage of the annotated genomic feature by RNAseq alignments, including 2 samples with support for all annotated introns"; product "uncharacterized LOC114963522"; transcript_biotype "lnc_RNA"; exon_number "3"; 
    NC_058066.1 Gnomon  exon    1962    2119    .   -   .   gene_id "LOC114963522"; transcript_id "XR_003825913.2"; db_xref "GeneID:114963522"; gene "LOC114963522"; model_evidence "Supporting evidence includes similarity to: 100% coverage of the annotated genomic feature by RNAseq alignments, including 2 samples with support for all annotated introns"; product "uncharacterized LOC114963522"; transcript_biotype "lnc_RNA"; exon_number "4"; 
    874259 ../data/Amil/ncbi_dataset/data/GCF_013753865.1/genomic.gtf
    CDS
    exon
    gene
    start_codon
    stop_codon
    transcript
    874254
    874254

``` bash
head ../data/Amil/ncbi_dataset/data/GCF_013753865.1/genomic.gff

wc -l ../data/Amil/ncbi_dataset/data/GCF_013753865.1/genomic.gff

grep -v '^#' ../data/Amil/ncbi_dataset/data/GCF_013753865.1/genomic.gff | cut -f3 | sort | uniq

grep -c "transcript" ../data/Amil/ncbi_dataset/data/GCF_013753865.1/genomic.gff

grep -c "gene" ../data/Amil/ncbi_dataset/data/GCF_013753865.1/genomic.gff
```

    ##gff-version 3
    #!gff-spec-version 1.21
    #!processor NCBI annotwriter
    #!genome-build Amil_v2.1
    #!genome-build-accession NCBI_Assembly:GCF_013753865.1
    #!annotation-source NCBI Acropora millepora Annotation Release 101
    ##sequence-region NC_058066.1 1 39361238
    ##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=45264
    NC_058066.1 RefSeq  region  1   39361238    .   +   .   ID=NC_058066.1:1..39361238;Dbxref=taxon:45264;Name=1;chromosome=1;collection-date=2017;country=Indonesia;gbkey=Src;genome=chromosome;isolate=JS-1;isolation-source=Whole tissue;mol_type=genomic DNA;tissue-type=Adult tissue
    NC_058066.1 Gnomon  gene    1962    23221   .   -   .   ID=gene-LOC114963522;Dbxref=GeneID:114963522;Name=LOC114963522;gbkey=Gene;gene=LOC114963522;gene_biotype=lncRNA
    828216 ../data/Amil/ncbi_dataset/data/GCF_013753865.1/genomic.gff
    cDNA_match
    CDS
    exon
    gene
    guide_RNA
    lnc_RNA
    mRNA
    pseudogene
    region
    rRNA
    snoRNA
    snRNA
    transcript
    tRNA
    431975
    803260

## 1.3 HiSat

``` bash
/home/shared/hisat2-2.2.1/hisat2-build \
../data/GCF_013753865.1_Amil_v2.1_genomic.fna \
../output/05.33-lncRNA-discovery/GCF_013753865.1_Amil_v2.1.index \
-p 24 \
../data/Amil/ncbi_dataset/data/GCF_013753865.1/genomic.gtf \
2> ../output/05.33-lncRNA-discovery/hisat2-build_stats.txt
```

``` bash
cat ../output/05.33-lncRNA-discovery/hisat2-build_stats.txt
```

    Settings:
      Output files: "../output/05.33-lncRNA-discovery/GCF_013753865.1_Amil_v2.1.index.*.ht2"
      Line rate: 6 (line is 64 bytes)
      Lines per side: 1 (side is 64 bytes)
      Offset rate: 4 (one in 16)
      FTable chars: 10
      Strings: unpacked
      Local offset rate: 3 (one in 8)
      Local fTable chars: 6
      Local sequence length: 57344
      Local sequence overlap between two consecutive indexes: 1024
      Endianness: little
      Actual local endianness: little
      Sanity checking: disabled
      Assertions: disabled
      Random seed: 0
      Sizeofs: void*:8, int:4, long:8, size_t:8
    Input files DNA, FASTA:
      ../data/GCF_013753865.1_Amil_v2.1_genomic.fna
    Reading reference sizes
      Time reading reference sizes: 00:00:02
    Calculating joined length
    Writing header
    Reserving space for joined string
    Joining reference sequences
      Time to join reference sequences: 00:00:02
      Time to read SNPs and splice sites: 00:00:00
    Using parameters --bmax 3713613 --dcv 1024
      Doing ahead-of-time memory usage test
      Passed!  Constructing with these parameters: --bmax 3713613 --dcv 1024
    Constructing suffix-array element generator
    Converting suffix-array elements to index image
    Allocating ftab, absorbFtab
    Entering GFM loop
    Exited GFM loop
    fchr[A]: 0
    fchr[C]: 144834195
    fchr[G]: 237637860
    fchr[T]: 330514140
    fchr[$]: 475342477
    Exiting GFM::buildToDisk()
    Returning from initFromVector
    Wrote 162768226 bytes to primary GFM file: ../output/05.33-lncRNA-discovery/GCF_013753865.1_Amil_v2.1.index.1.ht2
    Wrote 118835624 bytes to secondary GFM file: ../output/05.33-lncRNA-discovery/GCF_013753865.1_Amil_v2.1.index.2.ht2
    Re-opening _in1 and _in2 as input streams
    Returning from GFM constructor
    Returning from initFromVector
    Wrote 212701387 bytes to primary GFM file: ../output/05.33-lncRNA-discovery/GCF_013753865.1_Amil_v2.1.index.5.ht2
    Wrote 120924662 bytes to secondary GFM file: ../output/05.33-lncRNA-discovery/GCF_013753865.1_Amil_v2.1.index.6.ht2
    Re-opening _in5 and _in5 as input streams
    Returning from HGFM constructor
    Headers:
        len: 475342477
        gbwtLen: 475342478
        nodes: 475342478
        sz: 118835620
        gbwtSz: 118835620
        lineRate: 6
        offRate: 4
        offMask: 0xfffffff0
        ftabChars: 10
        eftabLen: 0
        eftabSz: 0
        ftabLen: 1048577
        ftabSz: 4194308
        offsLen: 29708905
        offsSz: 118835620
        lineSz: 64
        sideSz: 64
        sideGbwtSz: 48
        sideGbwtLen: 192
        numSides: 2475743
        numLines: 2475743
        gbwtTotLen: 158447552
        gbwtTotSz: 158447552
        reverse: 0
        linearFM: Yes
    Total time for call to driver() for forward index: 00:02:37

``` bash
find ../data/fastq/*gz
```

    ../data/fastq/RNA-ACR-140-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-ACR-140-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-ACR-145-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-ACR-145-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-ACR-150-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-ACR-150-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-ACR-173-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-ACR-173-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-ACR-178-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-ACR-178-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz

``` bash
find ../data/fastq/*R2_001.fastp-trim.20230519.fastq.gz \
| xargs basename -s -S1-TP2_R2_001.fastp-trim.20230519.fastq.gz | xargs -I{} \
echo {}
```

    RNA-ACR-140
    RNA-ACR-145
    RNA-ACR-150
    RNA-ACR-173
    RNA-ACR-178

``` bash
find ../data/fastq/*R2_001.fastp-trim.20230519.fastq.gz \
| xargs basename -s -S1-TP2_R2_001.fastp-trim.20230519.fastq.gz | xargs -I{} \
/home/shared/hisat2-2.2.1/hisat2 \
-x ../output/05.33-lncRNA-discovery/GCF_013753865.1_Amil_v2.1.index\ \
-p 48 \
-1 ../data/fastq/{}-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz \
-2 ../data/fastq/{}-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz \
-S ../output/05.33-lncRNA-discovery/{}.sam \
2> ../output/05.33-lncRNA-discovery/hisat.out
```

``` bash
cat ../output/05.33-lncRNA-discovery/hisat.out
```

    47710408 reads; of these:
      47710408 (100.00%) were paired; of these:
        27060558 (56.72%) aligned concordantly 0 times
        19176285 (40.19%) aligned concordantly exactly 1 time
        1473565 (3.09%) aligned concordantly >1 times
        ----
        27060558 pairs aligned concordantly 0 times; of these:
          1274633 (4.71%) aligned discordantly 1 time
        ----
        25785925 pairs aligned 0 times concordantly or discordantly; of these:
          51571850 mates make up the pairs; of these:
            41230080 (79.95%) aligned 0 times
            9296317 (18.03%) aligned exactly 1 time
            1045453 (2.03%) aligned >1 times
    56.79% overall alignment rate
    42864294 reads; of these:
      42864294 (100.00%) were paired; of these:
        23640036 (55.15%) aligned concordantly 0 times
        17633629 (41.14%) aligned concordantly exactly 1 time
        1590629 (3.71%) aligned concordantly >1 times
        ----
        23640036 pairs aligned concordantly 0 times; of these:
          1207857 (5.11%) aligned discordantly 1 time
        ----
        22432179 pairs aligned 0 times concordantly or discordantly; of these:
          44864358 mates make up the pairs; of these:
            35417854 (78.94%) aligned 0 times
            8221441 (18.33%) aligned exactly 1 time
            1225063 (2.73%) aligned >1 times
    58.69% overall alignment rate
    43712298 reads; of these:
      43712298 (100.00%) were paired; of these:
        30353070 (69.44%) aligned concordantly 0 times
        12098419 (27.68%) aligned concordantly exactly 1 time
        1260809 (2.88%) aligned concordantly >1 times
        ----
        30353070 pairs aligned concordantly 0 times; of these:
          826489 (2.72%) aligned discordantly 1 time
        ----
        29526581 pairs aligned 0 times concordantly or discordantly; of these:
          59053162 mates make up the pairs; of these:
            51066335 (86.48%) aligned 0 times
            6685083 (11.32%) aligned exactly 1 time
            1301744 (2.20%) aligned >1 times
    41.59% overall alignment rate
    47501524 reads; of these:
      47501524 (100.00%) were paired; of these:
        27841628 (58.61%) aligned concordantly 0 times
        18249998 (38.42%) aligned concordantly exactly 1 time
        1409898 (2.97%) aligned concordantly >1 times
        ----
        27841628 pairs aligned concordantly 0 times; of these:
          1281752 (4.60%) aligned discordantly 1 time
        ----
        26559876 pairs aligned 0 times concordantly or discordantly; of these:
          53119752 mates make up the pairs; of these:
            42721011 (80.42%) aligned 0 times
            9197577 (17.31%) aligned exactly 1 time
            1201164 (2.26%) aligned >1 times
    55.03% overall alignment rate
    42677752 reads; of these:
      42677752 (100.00%) were paired; of these:
        25633048 (60.06%) aligned concordantly 0 times
        15651560 (36.67%) aligned concordantly exactly 1 time
        1393144 (3.26%) aligned concordantly >1 times
        ----
        25633048 pairs aligned concordantly 0 times; of these:
          1075688 (4.20%) aligned discordantly 1 time
        ----
        24557360 pairs aligned 0 times concordantly or discordantly; of these:
          49114720 mates make up the pairs; of these:
            38244722 (77.87%) aligned 0 times
            9400848 (19.14%) aligned exactly 1 time
            1469150 (2.99%) aligned >1 times
    55.19% overall alignment rate

## 1.4 convert to bams

``` bash
for samfile in ../output/05.33-lncRNA-discovery/*.sam; do
  bamfile="${samfile%.sam}.bam"
  sorted_bamfile="${samfile%.sam}.sorted.bam"
  
  # Convert SAM to BAM
  /home/shared/samtools-1.12/samtools view -bS -@ 20 "$samfile" > "$bamfile"
  
  # Sort BAM
  /home/shared/samtools-1.12/samtools sort -@ 20 "$bamfile" -o "$sorted_bamfile"
  
  # Index sorted BAM
  /home/shared/samtools-1.12/samtools index -@ 20 "$sorted_bamfile"
done
```

## 1.5 Looking at Bams

![igv](http//gannet.fish.washington.edu/seashell/snaps/Monosnap_IGV_2023-11-05_09-57-16.png)

# 2 StringTie

StringTie uses the sorted BAM files to assemble transcripts for each
sample, outputting them as GTF (Gene Transfer Format) files. And then
merges all individual GTF assemblies into a single merged GTF file. This
step extracts transcript information and merges GTFs from all samples–an
important step in creating a canonical list of lncRNAs across all
samples included in the pipeline.

``` bash
find ../output/05.33-lncRNA-discovery/*sorted.bam \
| xargs basename -s .sorted.bam | xargs -I{} \
/home/shared/stringtie-2.2.1.Linux_x86_64/stringtie \
-p 16 \
-G ../data/Amil/ncbi_dataset/data/GCF_013753865.1/genomic.gff \
-o ../output/05.33-lncRNA-discovery/{}.gtf \
../output/05.33-lncRNA-discovery/{}.sorted.bam
```

``` bash
wc -l ../output/05.33-lncRNA-discovery/RNA*.gtf
ls ../output/05.33-lncRNA-discovery/RNA*.gtf
#head ../output/05.33-lncRNA-discovery/RNA*.gtf
```

       359447 ../output/05.33-lncRNA-discovery/RNA-ACR-140.gtf
       308507 ../output/05.33-lncRNA-discovery/RNA-ACR-145.gtf
       326448 ../output/05.33-lncRNA-discovery/RNA-ACR-150.gtf
       354558 ../output/05.33-lncRNA-discovery/RNA-ACR-173.gtf
       412414 ../output/05.33-lncRNA-discovery/RNA-ACR-178.gtf
      1761374 total
    ../output/05.33-lncRNA-discovery/RNA-ACR-140.gtf
    ../output/05.33-lncRNA-discovery/RNA-ACR-145.gtf
    ../output/05.33-lncRNA-discovery/RNA-ACR-150.gtf
    ../output/05.33-lncRNA-discovery/RNA-ACR-173.gtf
    ../output/05.33-lncRNA-discovery/RNA-ACR-178.gtf

![gtf](http://gannet.fish.washington.edu/seashell/snaps/Monosnap_IGV_-_Session_Userssr320DesktopApul_lncRNA_igv_session.xml_2023-11-05_10-11-31.png)

![igv2](http://gannet.fish.washington.edu/seashell/snaps/Monosnap_IGV_-_Session_Userssr320DesktopApul_lncRNA_igv_session.xml_2023-11-05_12-33-58.png)

Merges all individual GTF assemblies into a single merged GTF file.

This is used to create a non-redundant set of transcripts after running
StringTie separately on multiple RNA-Seq datasets.

``` bash
/home/shared/stringtie-2.2.1.Linux_x86_64/stringtie \
--merge \
-G ../data/Amil/ncbi_dataset/data/GCF_013753865.1/genomic.gff \
-o ../output/05.33-lncRNA-discovery/stringtie_merged.gtf \
../output/05.33-lncRNA-discovery/*.gtf
```

``` bash
tail ../output/05.33-lncRNA-discovery/stringtie_merged.gtf
```

``` bash

wc -l ../output/05.33-lncRNA-discovery/stringtie_merged.gtf
head ../output/05.33-lncRNA-discovery/stringtie_merged.gtf


echo "what is possible"

grep -v '^#' ../output/05.33-lncRNA-discovery/stringtie_merged.gtf | cut -f3 | sort | uniq -c
```

    742826 ../output/05.33-lncRNA-discovery/stringtie_merged.gtf
    # /home/shared/stringtie-2.2.1.Linux_x86_64/stringtie --merge -G ../data/Amil/ncbi_dataset/data/GCF_013753865.1/genomic.gff -o ../output/05.33-lncRNA-discovery/stringtie_merged.gtf ../output/05.33-lncRNA-discovery/RNA-ACR-140.gtf ../output/05.33-lncRNA-discovery/RNA-ACR-145.gtf ../output/05.33-lncRNA-discovery/RNA-ACR-150.gtf ../output/05.33-lncRNA-discovery/RNA-ACR-173.gtf ../output/05.33-lncRNA-discovery/RNA-ACR-178.gtf
    # StringTie version 2.2.1
    NC_058066.1 StringTie   transcript  345 1190    1000    .   .   gene_id "MSTRG.1"; transcript_id "MSTRG.1.1"; 
    NC_058066.1 StringTie   exon    345 1190    1000    .   .   gene_id "MSTRG.1"; transcript_id "MSTRG.1.1"; exon_number "1"; 
    NC_058066.1 StringTie   transcript  1962    23221   1000    -   .   gene_id "MSTRG.2"; transcript_id "rna-XR_003825913.2"; gene_name "LOC114963522"; ref_gene_id "gene-LOC114963522"; 
    NC_058066.1 StringTie   exon    1962    2119    1000    -   .   gene_id "MSTRG.2"; transcript_id "rna-XR_003825913.2"; exon_number "1"; gene_name "LOC114963522"; ref_gene_id "gene-LOC114963522"; 
    NC_058066.1 StringTie   exon    18711   18775   1000    -   .   gene_id "MSTRG.2"; transcript_id "rna-XR_003825913.2"; exon_number "2"; gene_name "LOC114963522"; ref_gene_id "gene-LOC114963522"; 
    NC_058066.1 StringTie   exon    21001   21093   1000    -   .   gene_id "MSTRG.2"; transcript_id "rna-XR_003825913.2"; exon_number "3"; gene_name "LOC114963522"; ref_gene_id "gene-LOC114963522"; 
    NC_058066.1 StringTie   exon    23085   23221   1000    -   .   gene_id "MSTRG.2"; transcript_id "rna-XR_003825913.2"; exon_number "4"; gene_name "LOC114963522"; ref_gene_id "gene-LOC114963522"; 
    NC_058066.1 StringTie   transcript  3278    4416    1000    .   .   gene_id "MSTRG.3"; transcript_id "MSTRG.3.1"; 
    what is possible
     627332 exon
     115492 transcript

# 3 GFFcompare

<https://ccb.jhu.edu/software/stringtie/gffcompare.shtml>

``` bash
/home/shared/gffcompare-0.12.6.Linux_x86_64/gffcompare \
-r ../data/Amil/ncbi_dataset/data/GCF_013753865.1/genomic.gff \
-o ../output/05.33-lncRNA-discovery/gffcompare_merged \
../output/05.33-lncRNA-discovery/stringtie_merged.gtf
```

``` bash
ls ../output/05.33-lncRNA-discovery/gffcompare_merged*
```

    ../output/05.33-lncRNA-discovery/gffcompare_merged.annotated.gtf
    ../output/05.33-lncRNA-discovery/gffcompare_merged.loci
    ../output/05.33-lncRNA-discovery/gffcompare_merged.stats
    ../output/05.33-lncRNA-discovery/gffcompare_merged.stringtie_merged.gtf.refmap
    ../output/05.33-lncRNA-discovery/gffcompare_merged.stringtie_merged.gtf.tmap
    ../output/05.33-lncRNA-discovery/gffcompare_merged.tracking

``` bash
head -10 ../output/05.33-lncRNA-discovery/gffcompare_merged.annotated.gtf
```

    NC_058066.1 StringTie   transcript  9330    48114   .   +   .   transcript_id "MSTRG.5.2"; gene_id "MSTRG.5"; gene_name "LOC114963508"; xloc "XLOC_000001"; ref_gene_id "gene-LOC114963508"; cmp_ref "rna-XM_029342736.2"; class_code "j"; tss_id "TSS1";
    NC_058066.1 StringTie   exon    9330    9651    .   +   .   transcript_id "MSTRG.5.2"; gene_id "MSTRG.5"; exon_number "1";
    NC_058066.1 StringTie   exon    12840   12886   .   +   .   transcript_id "MSTRG.5.2"; gene_id "MSTRG.5"; exon_number "2";
    NC_058066.1 StringTie   exon    15503   15664   .   +   .   transcript_id "MSTRG.5.2"; gene_id "MSTRG.5"; exon_number "3";
    NC_058066.1 StringTie   exon    22159   22443   .   +   .   transcript_id "MSTRG.5.2"; gene_id "MSTRG.5"; exon_number "4";
    NC_058066.1 StringTie   exon    24771   25067   .   +   .   transcript_id "MSTRG.5.2"; gene_id "MSTRG.5"; exon_number "5";
    NC_058066.1 StringTie   exon    26653   26913   .   +   .   transcript_id "MSTRG.5.2"; gene_id "MSTRG.5"; exon_number "6";
    NC_058066.1 StringTie   exon    27072   27359   .   +   .   transcript_id "MSTRG.5.2"; gene_id "MSTRG.5"; exon_number "7";
    NC_058066.1 StringTie   exon    27659   27952   .   +   .   transcript_id "MSTRG.5.2"; gene_id "MSTRG.5"; exon_number "8";
    NC_058066.1 StringTie   exon    28248   28535   .   +   .   transcript_id "MSTRG.5.2"; gene_id "MSTRG.5"; exon_number "9";

![gff](http://gannet.fish.washington.edu/seashell/snaps/2023-11-03_09-25-24.png)

# 4 Filter

Filters the combined GTF output from GFFcompare to select only the lines
representing “transcripts” and excluding lines starting with “\#” (these
are lines in the output format from GFFcompare that don’t contain
transcript information). This step further filters for a class code of
“u”, and keep only those with lengths greater than 199 bases. The “u’
class code from the GFFcompare step is for”unknown” transcripts, that is
those that are not previously annotated in our reference GFF as protein
coding. The size filter of +200nt is a common filtering step for
isolating lncRNAs.

``` bash
awk '$3 == "transcript" && $1 !~ /^#/' \
../output/05.33-lncRNA-discovery/gffcompare_merged.annotated.gtf | grep 'class_code "u"\|class_code "x"|\class_code "i"\|class_code "y"' | awk '($5 - $4 > 199) || ($4 - $5 > 199)' > ../output/05.33-lncRNA-discovery/Apul_lncRNA_candidates.gtf
```

``` bash
wc -l ../output/05.33-lncRNA-discovery/Apul_lncRNA_candidates.gtf
head ../output/05.33-lncRNA-discovery/Apul_lncRNA_candidates.gtf
```

    16982 ../output/05.33-lncRNA-discovery/Apul_lncRNA_candidates.gtf
    NC_058066.1 StringTie   transcript  468619  469943  .   +   .   transcript_id "MSTRG.37.1"; gene_id "MSTRG.37"; xloc "XLOC_000019"; class_code "u"; tss_id "TSS42";
    NC_058066.1 StringTie   transcript  1135316 1144814 .   +   .   transcript_id "MSTRG.112.2"; gene_id "MSTRG.112"; xloc "XLOC_000056"; cmp_ref "rna-XR_003823893.2"; class_code "y"; cmp_ref_gene "LOC114952870"; tss_id "TSS98";
    NC_058066.1 StringTie   transcript  1144884 1148491 .   +   .   transcript_id "MSTRG.113.1"; gene_id "MSTRG.113"; xloc "XLOC_000057"; class_code "u"; tss_id "TSS99";
    NC_058066.1 StringTie   transcript  1153399 1165634 .   +   .   transcript_id "MSTRG.115.2"; gene_id "MSTRG.115"; xloc "XLOC_000058"; class_code "u"; tss_id "TSS100";
    NC_058066.1 StringTie   transcript  1153399 1165634 .   +   .   transcript_id "MSTRG.115.1"; gene_id "MSTRG.115"; xloc "XLOC_000058"; class_code "u"; tss_id "TSS100";
    NC_058066.1 StringTie   transcript  1153404 1165634 .   +   .   transcript_id "MSTRG.115.3"; gene_id "MSTRG.115"; xloc "XLOC_000058"; class_code "u"; tss_id "TSS100";
    NC_058066.1 StringTie   transcript  1153410 1165634 .   +   .   transcript_id "MSTRG.115.4"; gene_id "MSTRG.115"; xloc "XLOC_000058"; class_code "u"; tss_id "TSS100";
    NC_058066.1 StringTie   transcript  1154207 1155609 .   +   .   transcript_id "MSTRG.115.5"; gene_id "MSTRG.115"; xloc "XLOC_000058"; class_code "u"; tss_id "TSS101";
    NC_058066.1 StringTie   transcript  1155787 1165634 .   +   .   transcript_id "MSTRG.115.6"; gene_id "MSTRG.115"; xloc "XLOC_000058"; class_code "u"; tss_id "TSS102";
    NC_058066.1 StringTie   transcript  1222540 1225166 .   +   .   transcript_id "MSTRG.120.2"; gene_id "MSTRG.120"; xloc "XLOC_000059"; class_code "u"; tss_id "TSS103";

# 5 Bedtools

Extracts the sequence data from the `$FASTA` reference file based on the
coordinates from the filtered GTF. The resulting sequences represent
potential lncRNA candidates.

``` bash
/home/shared/bedtools2/bin/fastaFromBed \
-fi ../data/GCF_013753865.1_Amil_v2.1_genomic.fna \
-bed ../output/05.33-lncRNA-discovery/Apul_lncRNA_candidates.gtf \
-fo ../output/05.33-lncRNA-discovery/Apul_lncRNA_candidates.fasta \
-name -split
```

``` bash
fgrep -c ">" ../output/05.33-lncRNA-discovery/Apul_lncRNA_candidates.fasta
head ../output/05.33-lncRNA-discovery/Apul_lncRNA_candidates.fasta
```

# 6 CPC2

Initializes a conda environment (Anaconda) and runs CPC2, a software to
predict whether a transcript is coding or non-coding. The results are
saved to the \$OUTPUT_DIR. CPC2 uses ORF (Open Reading Frame) Analysis,
Isometric Feature Mapping (Isomap), Sequence Homology, RNA Sequence
Features, and Quality of Translation to assess coding potential and flag
any transcripts we would want to exclude using the FASTA generated in
the previous step.

``` bash
eval "$(/opt/anaconda/anaconda3/bin/conda shell.bash hook)"
python /home/shared/CPC2_standalone-1.0.1/bin/CPC2.py \
-i ../output/05.33-lncRNA-discovery/Apul_lncRNA_candidates.fasta \
-o ../output/05.33-lncRNA-discovery/Apul_CPC2
```

``` bash
wc -l head ../output/05.33-lncRNA-discovery/Apul_CPC2.txt
head ../output/05.33-lncRNA-discovery/Apul_CPC2.txt
```

    wc: head: No such file or directory
      16983 ../output/05.33-lncRNA-discovery/Apul_CPC2.txt
      16983 total
    #ID transcript_length   peptide_length  Fickett_score   pI  ORF_integrity   coding_probability  label
    transcript::NC_058066.1:468618-469943   1325    73  0.3891  11.48389835357666   1   0.0882804   noncoding
    transcript::NC_058066.1:1135315-1144814 9499    69  0.2705  10.028266716003419  1   0.0399931   noncoding
    transcript::NC_058066.1:1144883-1148491 3608    75  0.25564000000000003 9.143180274963381   1   0.0451723   noncoding
    transcript::NC_058066.1:1153398-1165634 12236   134 0.3118  9.09901943206787    1   0.295457    noncoding
    transcript::NC_058066.1:1153398-1165634 12236   134 0.3118  9.09901943206787    1   0.295457    noncoding
    transcript::NC_058066.1:1153403-1165634 12231   134 0.3118  9.09901943206787    1   0.295457    noncoding
    transcript::NC_058066.1:1153409-1165634 12225   134 0.3118  9.09901943206787    1   0.295457    noncoding
    transcript::NC_058066.1:1154206-1155609 1403    66  0.37373 8.802658271789554   1   0.0495107   noncoding
    transcript::NC_058066.1:1155786-1165634 9848    134 0.29479 9.09901943206787    1   0.241898    noncoding

\#Filter Filters the CPC2 results to get only noncoding transcripts
(using the class “noncoding” from the CPC2 results) and extracts their
IDs and matches these IDs with the sequences from the previous step to
generate a GTF of long noncoding transcripts.

Matches these IDs with the sequences from the previous step to generate
a GTF of noncoding transcripts.

``` bash
awk '$8 == "noncoding" {print $1}' ../output/05.33-lncRNA-discovery/Apul_CPC2.txt > ../output/05.33-lncRNA-discovery/Apul_noncoding_transcripts_ids.txt
```

``` bash
wc -l ../output/05.33-lncRNA-discovery/Apul_noncoding_transcripts_ids.txt
head ../output/05.33-lncRNA-discovery/Apul_noncoding_transcripts_ids.txt
```

``` bash
head ../output/05.33-lncRNA-discovery/Apul_lncRNA_candidates.fasta
fgrep -c ">" ../output/05.33-lncRNA-discovery/Apul_lncRNA_candidates.fasta
```

    >transcript::NC_058066.1:468618-469943
    taactgatcaaaacgtatcttcctacaacattaatttgacagtggcgtttctcaactgaccaatcaaaacttacatttgaaaatttggtgATGGTgcgtttacaactcgtgtatctttacgtcacacaaccatgtttgcATACTCTCTTGCaaccacgcctctcggccaatcagagcgcgcgcgtactatcttagttattttataaagatAAATACGCCCTAGGATTAGCACGCACGCTATGGTATAATTATTGATGATAACTTTGCTGGATTTACGTTTGGTTGAAGTTATCATGATATtccatcgtcgtcatcatcaacATTCTTATCGTTTATCTTCATCACAATCACCTGACACAACATGACTAAAAGCAAAGATGAAAACACTCTTACATCACCAGCCCGTGTGTGGCCATCAACGCATGCATGCGCATCACCATATCTCCTGGGTAGTGTCAGCCATGAACAGCAGTTTCGGTGTTGTTAGGTCTCgtctagtctccttcgcagccgtctttcgggacggggagcgttgcgtgacatcccgaaagacggctgcgaaggagactaggtctCGTCAAGAGTGGATCAGGTAGGAGTGTTCCTCAATCACCTTACGGTAATATCCCAGCACTGTCGGAAACCTCACCTTTAACCCAACAGCTTAAGTAATATTATTCAGCCATGTCTCGCTTACCCGGACATACTTCGTCGCTTCAATGTCATTAACAGTACTATTATTTCGGAAATGGACTTTTTGGGGAATCGTTACAGTTACAGCAACTTATTTTCAAGGAATTGTGTATTCTAATTTCCAAAGAAATTGTGGTGTTGCGTCGGTGAGACCGTAacagtgaaacatgaaaattgggtttttatcagacgagttggtaaaggtcgaattaccaccgtgaaagatttggaaagctgacttttcgagcgttagcccttcgtcagagcaagtGAATTATTCTAAAAATATACTGTGCTGGCTAGCCTGCATACATCCAATGGGAACACCCTTACCTGATCCCCTCTCGTTGTCACAGCATAGCCGTTGAACCATCTTCGGAAATGTCATAACCTCAACCTCTTTTTTGAACAAAATGTCTGTCACCACAGAAAcgacaaaataaataaatactcAAGAATTCATCGTTTAAAATAGCGATTCAAGAAGAACAGCTGTTACTGTAGTTGCGCCCACTAGCAAAGCCTCTTTTGTTGAccgctgcatgcagacgaggcttacgggtcaatacaatggaaaacgACCTGTCAGCTTCGTTTATGTGTGAAAACCGCTGTTAACTAGTACAGAATGT
    >transcript::NC_058066.1:1135315-1144814
    gtttcgagatagttgacgaggaatatatcgaaaaattaaaggacaagagcgaaaatgaaaacaagaaaaatagcccggagtggtgaaagaacgttttcaaaaaatgggcGAATGAaaaaacttgcaagcaaatttagaagagtacgagaaggatgtcctcgaccaacgattgtcgcagttttaagcattcagaaattcagtaattattTTGCCCtctgttattaacaagtaatcgcagtGGATCCTCGGAAAATTAAGGTTAATATCagttgtgttttcagaagttgctgaaattcccTTCGTCGCTGCGTGACTCTGGCAATTCctgcaacttctgaaaacacaagtgatattaatccttaattttactcgcccccatgcgattacctatactaacaaGCTACTGTGCTTAGGTGCTAGTTACCATGGGAACTATGCAGTACTCGCCCAGGATTGCTTTCATTTAATGTCCTAGCCACTATATTAGAGGACCAAAACAACAAGCGGCTGCCAAGTTGGTGGATCCTCAAAACACAGACTGACATACCAGGTATTATACAGTTGCCATGGTAACTAGCACAAAAACAGAGGCTTACGCAGACTGGGACTGCATATACAGACAGTTTCCAAAAAAACTAGCTGATTACAATTCATTAAGTGAGCCAAGTATtgagtaattattattacaactagGTTTGATGACACTATAGCGAACAGATAAACAAATGCAAGATCTTACTAACATCCCaatgataacaaaaacaataccTCAGACCGCACCATCAAACAATGAGAGCCCATTCAACATTAAAGGTACATAAATGAATTGCCGAGGCTATAATTATTATGCCAAGTAACTGGGAAGGGTGTCTTTTAATAGAATAGCAGGAGGTCAGAAATTTTTGTAGGTGCTTGCAAGTATTGTATCTTCTTTTAGAATTGTCGTACGATTCTGACAAAATGTTTTTGTGGCTGACCATTTTTGCAAGTTTGCCATTGCATCAAACCAGCATGTACAAAACACAgatcacaggtcattgttttgcCAATCCTGAAACAAACCAAATCTTTTACTAATGCTATCCCTAGGCCTTAACACTACTACCATTAGTAAAGGGTTTTGGTTTCAGTATTGGTAACACAATGACTATGAgcaatttagagtgtttcctTCCCAACTTCTGTGCATATAGTCACTGCGTAAATAAAAACATCCAGCAACATTGAGATTAAAAGTTCCCATTGATACATTATGAATAATTGTGCAACCATAGACCATCtactattagtataaatttactcgtagtctatcgtgaatctgacTGGCTATATTACtcatagactatctgctgatagtcaacagttgtgaatagccaatgaaaatcgttctatttattttgaccaatcacgaagcttcagtgcacaccGCACGctgtgtttgaaaccttgaatttgaattgccaacGTAAACACAATTatacacttttaaccatttaaactttactttacatttttatgcaatgagaccgttgtaaatttatactaaaacaattagactactcgccctcgttttctatgagcgatagtcaactcggctgggcctcgttgactatctgctcgtagaaaacgagggcgagtagtctactTGTTAACTAGAAAAATCtaagaaaaaaacataaatttcAGGAAATATTTATAGATCATGAGGTGCAGTCAGTATCAATACAGACTTAAGACATGATGACTGCGGACTACAGACTATGATTGGATGGGTGCTTACTCAGTGAAAATATGAAGTTAAGTTGGCAGACGAAAGACAACATCCCACGGTTAAAAAGCAGTGTTTGACTTTGCATgtttttacaaaacaaaatgatgTCTGTCTCAAAAATTACATGTAGATAAAcacaaatataattattataatcaagtATGAGACTAAAGAACGACAACAAACAGGTAAGAaaagaacagacatttttatGAGGGGTGATCTCCACCACAACTGTGTTTAGAATTGTGTCATGCAGAATctccaatttttctttcaagcaATCAAATTCTTTCCTCAAGAGGTCAGAATGGAAATTGAGgtaccacaaaataaattaaagaaacaGAACAGAATTTAACACAATAGAAGTGAAATAAGGACACTGAAATAGAAATTGAAGCAACAGAACAGAAATTACAAAACAGTTAACAACCATTCCAAATTCAATTCTGATGGCATGTGGCTCAAATGGGGAACCAAGTAGAACCTGGCTAAGAACATTATACTTCATGAACATTACATTTGATGACAAAGTTGCCGTGAAGGGAAAATTTATTCCTCGTTTGAAATGGCCGCTCAAAAGGAATGCACCTGTGAAGAAGGACTTTCAGCTGTTCCTCTGTAATACTGTAGCATGGTCTGCCTTGATGACAGCTTAAGACCACCTTGGCCGTTTGAGCTGAAAGAGCTTCTATGAGCACAAGGACAGTTGCGAGCTTCCCGAAGAAGATAGTGAATTTGGAGGCACAACCCCCTGTCCCAATTACTTCCCACCCATAGTTTGATACACATGTGAGAGACTAGAACGTACTGTAGTGGGATCGAATTTACCATATTTTTGGGTGTATAAGTAATCCTGATACAGTTCCTAATACCCATAAGTTTTCTGTACTGGCAACTAaaggcgcttaccatttgttaAATTGGCTGGCCCGAAAGAACTGGTACTCATAGGTGCAAATGGATGGGTCAAAGGGAGTCCCAACAACTACTCACATGATCTTGCCAGAGAAGGCggaatttctattttttttcagttatttgtATTTGTTCTCAGCCCAAATTAGCAACAATTGCAGACTACTTCAAACTCCCAGTTCCTTGCTGTTGGTTGTGATGGCTTGAGGAGCAATTCCAGATGACATGCTCAACATTCGAGTGACTAACTAACTTGTTAGCACCAAAGGAACTTATACATTTGGTAGGTTGCCCACTAAATAAAGATAACAGAAACTGTAGCTGTGAATGTCTGGGCCCTGGCAAACAAAGAGGCTTGCCATCAGATATCTGATCGCACACTTCTTGTGTTCCGCCTAGGTTTGAGAGAACTGAACAGCCCTATTTCATGAGCAGGTCCAGTTTGGCCAGAACCAGTTCTCTTAGTTCTCAGAGCTTGCATGCTCCATTTATAGACACAACAGTCTGGCCCAGCCAgatctgacaaatggtaagcgccctaagTCATCTTCcaaattagtataggtaatcgcatggcgccgagtaaaattaaggattaatatcacacgtgttttcagaagttgctgaaattgccctcgtcattgcacgactcgggcaatttcagcaacttatgaaaacaaaagtgacattaatccttaattttatgAGGATTCAttgtgattacttgttaataacatagagggcaaaattactgaatttctgaatgcttaactgaaactgtgacaatcgttggtcgaggacattgTTCttgtactcttctaaatttgcttgcaagtttctttcagtAGCACACTTTTAgaaaacgttccttcaccactccgtgttgttcttcatgttttcatatTCGCTATTgtctttaattcttcgatatattcattgttaactatctcgaaacgagacgccattgttgtacAAAAACatcttctcgattgaatgagctgttgctaggcgacctgaggaccaatcgcgagcgagtaattttgccctcttcacaaagcaaaattaagaaaaaatacacgcttcattgaccaatcagcattcagtaattttgccctttgtGTTATTAGTGACCAAAACCATTCCCTACATCCTAACTTTATGTATTGCAATCAAATAAGGTTTCCAAGAAGATTAAGGTTGAGTTGATGTTAAATTCAAAGAGAACATGTCGATATCAACCTCCTTGGATTTTAGGTACCACTGCTCATTAGGTTTTCATAACAAAATTGCTGAGTCACACAATATGGCAGAGTAATTAAAGTGACTTcaatgtttcatgttttttagtgttttgacattttgctCTGCATAAAAGACTTCCAAACATCGCCCAAAAGTTACCAAGTATGTTGAAATACCTACCCTTATGTCCAAGACAATTAGCATCGGATTtttgggattttaagataaatTGCATTGAATTCATATgatttgcatgtttttcttcaCAGGACTGATTCAACAGAATTTTTCAATTAGGTTGTACCCcaataataatatgaaaaatCTCAACTTAAAACTTCAAAAAATGGAGTCGACTTGTACACCAGTTTGACTTATACACCAGAAAATATGGTACGCTTCGAGCAATGGGTGTCATTAAACAGAGTAACAAAGTAACAAGTAAGGAGTAATGAAGTAACGGAAAACTGAGTAATGGGAAAAATTTTTGAAGGGACAAAAATTCTCAAGTCTCTTTAAAAATAATATGTTGAATCGATTATAATCGTTTTTCAGTATTGATATATAAATAGGAACTGCAAAACCAATTCACAAAGCAAGCTTGACATAGACTACAACCCATAAGCTCAGAACTAAGAGTTAAATTAAACATTGTGCCTAGCCCAGCAGACGGGGAGTGTCATTAAAATACTACATAACCGGCAGTCTGAAATTATTCTAGAATTTCATGGAAAAACTTATCCAATCTTTCTTAAAACTGTGGACATTTATTGTAAATGTCAAGGGAATTGAATTGTAATTTTGGCTTTGTGGTTGAGGAAATTTTCGAGCCTTGGTATAATAGAATCTCAATGCTTCGACACTTTGAAAACTGAAGTTGGAAACCTGTGAGAGAATATTTTGTTCAGAAATCTGATTCTAAATGacgttttgcttttgttgtgaGCGAATTCTCTTAGCTTCCGCATactatattgaaaaaaaaaaaaaaaaaaaaaacttcgaATCAGGGGCACGGTCATTCAAAAAAACTTACCAAGGAGTcaagaattcatgaataatCAGTACATTTGAAATCCACTAAAACGATCATTGAAACCTCCTGAGAGTATTTCAGCATTTTGACCTATGATCTGAAACATCAAGTGACGTTTTACATTGTTAGCATTTTGGTTTTGCAAGTTTTTCAAACACTACATAAGCAGCATCCCAAAATTCTTCCTGTTTTTCCCCTAAAATTTACATTTGAACCTGCTGGGATGTTGGTGTAGATATACCAAGCCCTTTTTGGCCTCAACAATTGCCTTGGCAAAGTAGTTCAAATAGCTATATAAACAAAACGCTTGAGGGTTTGGTGACAAAACTACACAACAGAGAAATAATTGGATGTGGTATAACACTTCAGTTGTGTGACAAAGGTTGCAGCGTTGGAAGCGACGAACACGTTGTGGGCTTAGCAATGGAAGCGACAAAGTTCTCAGTTCTCCAACTTATTTATACTGTAAGCTAAACAGCAGAAGTGTTTCCGACAAAAGATCAAAATGCTGAAATGCTCACACAAGGTTTCAGAGATCATTTAAGTGGATATTGTATTGATGTGGCTTAATGCCTGGCAGTCTAAAATCTTCACAGAATTTCCACTGAAAGattctttcaaacttttttttagtgttctttgatgtttttttttttttagaaatgaaGAGAGATTCCAAAAATTATGTTCCTTTGTAAGATAACGAAATATCCAAGCTTGGGCATACTATAAAATAGGCACTTGTCAGGTGATAGATTTCCTACACAACCTGATGAAGATAACATTAATTTCAGGATTAATTTGAGCCGCAATAGTCTTGGGCATTTAACTTGGGCCAATTTGGGtacaattttgaaaaccttggagaacgaaagaaaaaacataCACTATCTGGGTTGTTTCTACAGCGTGGTATTATGTAAACACAGAATGACATAGACCGATTTTAAACCATTTCacacattttaaaataattttgaaattttaacccTTCAACTTTTTTCCCCCGTTACTTTGGTTCTTGTCACTTCTATGACATTAATAGTCATTGCATTATgtattaaagtgcatatgaaatgaaattttttacTAACTTATTCGAAAAAGTGTTCAAAATGAGGAAGAATGGTGTTTACTTTATTGTGATACCACTATTGCTTGCCAAGTTATAGAAgaattttatgcaaattagaggacAATCATGATGccacaatgtggacacaaagtgatgtcaagtcacaaaaaattgaatatctgtgcaaatactacatctacagggttgaaattttgcaggtttgatgtactgcaagaactacacattgtgatgGTGGTTATGATGTTACCATAGCAACCTACTCGTTACCAGATCTCTACctccctaaaatgaaaaatgcctcatTTGTTGCTCCTGAGTTTAatggactttcttgtgcttgggCTGTGAAATGTCCATAtttgctcacacccactgaatgcaCAATAACTGCAAATAACCCtacttgaaggaggaaaactctgattttatcctttgaatggagaggacCTGGAGcacattgtgttgccatggaaatgtcacagtgaaCATgacatggaactttgtgatgagtgtaacaaccgTACAAAGTTTCGgttctttgcagaaaaagtcaatttttttagtttacaTCACTGTGTCcacatttttttgtcacaagtcctctaatttgcataaatcaaaatcttgaataactcaggaaccaagagtgctattaacacaacaaaataaattccATTCTTCATTCATTTGACAGCTATTTTGAATAAGGTAATAAAAagtttcatttcataggcactttaagcaaaaacgttattgcctgatgagtgtggtcacccttcgatcatacgaaacagcgttgtagcaataaaatgatgaactgaaattttgctaccattaatcattattatataaatatatatatatatataatttttctaccactgatcattattattatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatctatatctatatctatatctatatctatatctatatctatagatatctatatctatatctatatctatatatatctatatataaaagcctaagccgtgtattacactgcgataaaacactccggacatttgagaacactcgagaaatgtacaAAACACTCGCCTGGGGCTCATGTcttctacatttctctcgtgttctcaaatgtccattgtgttttatcacagtgtaatacacggcttaggcttctttattttatatatagatatatatagatatagatatagatatctatagatatagatatagatatagatatagatatatatatatatatatatatttataaaatctggatgatgtgctcaacagaaggtatgacttggcattttattactacctagtttcatgccacacaagaagtgcggcactcatcagataaaatagccaaatgcatatatatatataattttgcTACCactgatcattattatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatctatatatatatatatatatatatatatatatatatatatatatatatatagtatagaacagcaccgaactcagtaaatacatctgggcactgaagagcaaaaatgaacggTACACAATCGGCTGGAAGATCCTCGAAAAGACGAAACCGTACACAAACCTAACCGGCAAATGCCAGCTGTGTACAGCAGAAAAACACTTTATCATCACCAAACCAGAGCTGGCAACGTTGAACAAGCGCAATGAGCTAGTTTCTGCATGTAGACACCGAAGAAAGTTTATCTTGAGATACTCTATCCAGTAATTGGTTCTTTCGTGCTTTTGCGCGCGCGCTCATGTTAATGTGAATCGCCCAGGTAACGCGGGCACAAGTTTTGTAGTCTTCAACGGTTTCCTAACGGTTCTTGTAacgtatcatgttgtttttgtgtagctttcgaattattgacgacagtttgtatggaggataactgatgagtgggccaaccatgttggcctacgaaacaaatttgtattaaaatccatatggactcaaaaacagagtagttttaccacatattcaattatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatctatatctatatctatatctatatctatagatatctatatctatatctatatatatctatatataaaataaagaagcctaagccgtgtattacactttatataataacccaagtaattctcgcattctaattggttctcgcctatgatctattagaggacagacgcacagatgacgacagcgctcgattcaagttttttaaattttttgaattttgaatttgaaccaatgacaattctttgctaagcatagcaaccaatcagttcgcttcatttttttatagacataagatcacgtcagtgctgttttcgtgtctgtcaaagtggcgaaatttgaaataaaagggcattttttccgtgtattttaatttttttattatataaaacaaatagattccatgttgccgtgcgtctgttcagtaatagatcacagaggacatcaaaatgtggtaagaacatcagtgacacactcggctacgtctcgtgtgccacttttttgttcttaccacattttgacgtcatctgtgatctattactgaacagacgcacggcaacatggaatctatttgttaagcgataaaacactccggacatttgagaccACTCGAGAAATGTACAAAACACTCGCCTGGGGCTCATGTcttctacatttctctcgtgttctcaaatgtccattgtgttttatcacagtgtaatacacggcttaggcttctttatttgttagatatatatatataactataaGAGAAGCTTTTGATGATCGTCCTAAAACAAGAATCAAATAGTTGGATGAGAAGTGTTCATTGATTCCAGTAGGACTGAGGATGCCGAGTTGTAAAAcagctttttgttcttttgcagCTTTCTGTATTCCTTTTGTGTTATAGGGAAAGACCacaaacaattattgttatgttAACATTTTAGAAGGTTTCGATACatctttgtcatttttcttgttcTCAACAGCACAAAGTTGTTTGCGAAAGCAATCTACTGGTTTCCTGCCTGTTTCAATTATGTACAGTACCGCTCACAAGTAAGTCACCCCTGTTGCATGCTACGCGAACTACATAGCACACTACAGGAATCACTTGTGCTTCAGTAATCGCATAGCACGCAACGAAAAAGAAACTGCAAGTTTACTAAGAAATTGTGCCTTAGATTTGACAATTAATTGTGAATTCCGAAGAGTGGTTTCATCTCCCAAAACATTTTGCTGGCCAGTATCTGTGCTG
    >transcript::NC_058066.1:1144883-1148491
    AATAATCATACTTTCCAGCTCATAGTGAAAAATGAAACATGCGCCAACACTGATCCTAGTGAAAGGAAAATGACTTTCAGCTGGTGAGAGTCCAAAAGCATAATTTACCTAGAAAATCTCATTGAAATCGAAATAAATTATGTCACTTCTCGGCTGAAATTTGCATAATGTACTGGTACTCACAAGGCTTTGAAAATTATTATCAAATTTCAGTCAatataaaatgcagactgcagattTCAGACTcttctcattaattttgtttttggtgGCAGCATTTTGGGTTTTACTGTAAAATATGATGTTCTGTTGGCTTGCAAGGATTGTGAGTTATTTTCATTACAATCTAATCTTTTTAGGCTGTGTaatttttgtgcactgaaaCAAAACCATTTGTAAATGGGCATCTTGCACTAAaagttgaatgtttttttttttttaaatctcaaCTGCACATACTTTTTTAAACTGGCACAATTCCCTGGCTTTGCCccaaaagcaacaaaacaacaaaaggcaaggaacaaaaaaaaaaccttgaaagTTATGTTTACAGAATCAATGAATACTTTTCATTCAACTAACTCTTTCATGTTTTAATAACACCATGTTCCCAGAAGTAGCAGTTCCACCACTTCCTAAATAGGTCAAGAAGAACCCTCACTTACACTATTTGCTCTGATGACGGGCTAACAGCTGAAAGATCATCTTTGcaacttttaataataattattagtaatcTGACCTctatcaacttgtttaataccaaaCCATTGTACTTCACTTTCCCCAATGATGTGACACCACAGTTTCTCAAGAAACTAGACACTTCATTTGTTTGACAATGCTGTAGTCCAATTAAGCAAGATCACTTATCACCAACACATTACTCTGTCACCAATATTTCTGAGCAAGCCCTTGAAATTGAAACTTCAGTCCCACCCACAAATAAGTAAATTGAATCTATAAGTAATAGTGTACATTATGCAGTGCAAAAGGAGAACTGTGCTTAGAGTACAAGTATGATATACATGAGCCCCTTGATGTTCAAAGGTCTGACAGTGCTATTCAGTGGAGAAATGAATATCTATTGGATAAGTAccaccaaaacctattgagttatccagtggatggtgatttatccaatagatagtgctatccaccaTTTGAACATCAGGGACCAAGACAAAAGTACTGACTCTGTCTTTTGAAACTGGTCCAGAAATCCAGAAACCACAAATTCCATTCATGATTGCCATGCATATTTCTTGTCAGGCTTTCTTTGGTCAAAAGTAAGCTAGATAATATGCAATGTCATAAGCCTTTTGCTACAATAAAATACAAGAGAAAACTTCCCTCCACATTCTCACTAATTTGTTGGGGTTTAGAGCTTGTAATTATCATCTTTACTATCCATTATTTCTTGCTTGCTACTATGATTAGTACTGCTGTTCCTTCTTTGTCAGTTTGTTGTGCTGGTGAAATGCATTTGGGTATGAATCTACACAGTTACTGCAAAGATAAAAGATGAGCTGCACCTTcttcaaaaaaattattgtttgaaattgtGCACATGATCCTTTCTTGATGGACATTACATAACAAAGCAGACTTTTAAAGTAACATAGCATTTGTAAGGAAAGGAATTAGATGCAATATTAAAGTAGTTTGCAGTTATGCATTTGTGCGTAtttttgaaggaaaataatattgcagaataattattttctatgCCACTCTGATTTGCAAATTTGTTATGATAATTTTGCAACCATGAAAACATTTTGTATTATTACCATCTTACCATCTCACCATTACCAGGTCACATGACCATGGTTTATTAAAGTAAATGCAAAAATTCTCAATTCCTTTCCATAAAATTAACTGGCATACTAgcacaaagaacaaaaaacatcATCTTGGGACAGTACACAGTTCGCAAAGATCTTTCAATTTTTCGTAATGAAGACTGGTGATCAGCTTTAAACCTAATTCTTTGCTTATTTCTATTGACACAAGTCAAActgcttttttacactgatgaatgGCTCTGCCCGAAatatttgtacattttttaatgtaaatttaactttatttttgtttttatagtGTCACGCactgttttccaattttttgaTGCTTCAAATTTTGTATGTGCTGTCTGCTTTTACAGTGCAATGTTTGACTGGCTGCCAAGAGTAAggttttctttgctttgctttCGTATCATTCCATCAATTTCGTTGAGCTTTCTAATAGTTTGTTAGGTTTGGTTTGCTTTAAGTTTGCATCATTTGATTGATTGTTTGAGCTGTTTATTTGAATGACTACTATTCATTCGATCAGCTGGAGCATTTGAGCctttttttaagtttgtttGAGCGTTGTTTGAATCTGTTCAAAATTTGAATCACTTACTATTTCTTTGATCTTTTTACAGAGAGAACTGATATTTAAGACAAACTGTTCTGTTATGACCTAATACTCTTCACgtggaaaaacaaaagagcaaaattaatgcaCAAATGAAAAGCCAAATTTCccaaattttttctttttttttccccccACATttaataagaaagaaaaatgtttgggaTTCATGGTATTTATTTTTCCCTTGAAAGCTTATTTTTCATACTTCAAAATGATatctatttataaataaatgttgATGATTCTTTTTTCTGCCAACTAGCATTTCTCACAAAAATTAGCACTCATTTAAACTCTCCTTTTTCTGATTCTTCCAGCACAGTTTGTTGCAATCATCACCTTTCAATACATTTCAATCAGGTACTACTAATGTACAGTAGATAAGAATCAGTAGATTGAATAAAATTCAGTTTTCCCTGGACTGAACTGTAAGTTCCCTTGGTTTGAAGGCTTATGCCATTATCCATTTTCCTGACACATGACAATCATAAAATCAAAAACCTCCCAAATAAATCAAGAGCTGCTCAATATCTATCATACTTGTTTCTTGCACAACAGTTTTCCTTTTGTGATATACTTGGCTCACTTatgaaagaacaaaagaaaagatcaTAATTATAGCTACATCATCCTGCTTGCCTTGACTTTTGGATAATGGGTCTCCATCTGTCTGTTGTACCATTCTGTTAATTGAAGCAAATAAAGATTGAAAAGTTGGCATTAAAAGACGAACTTAAGCAAAAGCATTTAAGGCAGTTTTCAGAATGTTCTGAatccatttttttcatttatggcTCTGTTAACTTAATGTATAATTTTCGCAAAGGCTTTAACTAAATTTTGACACACAATTGACTCATATTGATCAGAGGGTACAAAGTTGGCCGCATAGCACTAGTGATTAAtagaatacaccttattccaaagtggcggccaataaattattcttttgtttgcatgttaattagccctcttcgccatgtataaaaaacaaaagaattttgaagcgaaaatgaggcaaagagagctaataaacatgcaaacaaaagattatttattggccgccattttggaataaggtgtatacagtagtcaataattattattagaaggTTTAGCACTTTTGTTCACCACGTTCTCATACACAGATGGATCAcattaaattcatcatc
    >transcript::NC_058066.1:1153398-1165634
    GCGAGCGTGATCTTCAGTTCGTTTCGAGGTTCAAAAAAACACATTGCCAACAACAGTGGAGTAAAGGGGTATTATTCTTGCACAGTACCTCAACAGTTCTTGCTTTTGcttcacaaacaataatatgGATTCGTACCAGttggattattattttttgggaAAGAAATGTCGGATGAAACATGAATTACGACCTACTTTGTTTACGGATCTCGTCCGTTGCAACTCAAGCGTTTTCACTGGTTTTCGTCATCAAATAACCGTTCTTGGCTTCGTCGAAATAAACCGTTGGCGAAGAAAATTGTCATGATCCATTGACAGCGTTATGGAAGCGTTATGTAAGCTAAAACAAGTGACACAGGAGGTAAGCACGTAAAGCTTAAACACATTCGATTCATCGCCATGGTTATGAACCTATTATGCgcttaagtttaagttttattgCAAAGGATGTAATCTTTTCCAAAATCTGTGGGAACTCCTCACGACAATTTGCAAGCTCAGGCAAAGCATTTTTGCAACGTTTGTCGTTTGCCAGcattcataatattatttactaTTAAAGATAATGGCAGTGGCATAAAACTGAAGTGTTTCcgaaattttcaaatgtgttgTGTTTCTCTTAAAGCTAACGTCTCACTAATGGAATGTAACACATGTTACTCTGCCCAAGCAATATGAGGAATGTTATTAATGGCTTTGAACCATGAGAGCAATTCCCATGCTGCATCAGCTGGGGCCATAGTGAGAAATGCACCACTTGTTCCATGATAATAGTTGTCAGTTATTCCATTTCGTGTTGGGatgaaagccacagttttgtcaTTGATGATGTAGAGACAAGCTGTCTTGGTCATGTTACGGGCTTGTAGAAAGTAATATGGAAGACATGCTGATTGTATCTGCTCTCTGGTGATACTTGCTGGTAACTCTCCAGACATGTGTATACTGCCAAGAACTCTGGCCTGCTGCATTTTGGTAGCTGAATCTTGGACACCAAACATCTGCCCAGAGTTGCCAGCAAATTTGTCATAGAACTTGTTGCCAATCATAATACCGGCAAGGGTGCGATGCATCCAGGTGTTGCTCAATGGCTGGCTACTGTTGTGAGTTTCAGGCGAGGAGAAATATAGTTTACTGTATGAGAGTGCAATGAAAGTGCATGCATTGCTACCTAATCGGCCAGTCAAAGTAGATTGGGAGTATTGGGGTGGAAAGTGCCAGCTTATGACTCTGTCAGATGTTGAAGGATTTATTAGTAGATGCTGAGTGACGGAATGGATTTCTACTGTTGAGGACTGAGGTTGGTTGTGAGCATGAGAGGCTTGTGCATTTGCATGTTTCGGTTGATCAATGCCCCCAATTCCTGTTAGGATATTGTCTACTTTCTGCTCCCAAACTCCAGGCCCgggttgttcaaacgatggatagcactatccaccggataaatcgctatccacaggataagtaatagcgaaaccaattattgcgatatccaatggatagtgatttatcaggtggatagcgttatccaccttttgaacaactggggccagatgtaTAATACTGCTCAGCATTTGTATCAGTGACTTCATTGTTTCCGTCAATTGTATCATCATCAGTGTCAGACTCACTTAGGTTGTCATCAGTGTCTGTATTTTGTacgtttttaaagcaaatgtcaTTTCTAATTTCGTAGCATGTATCTCCAGTTTACTCAAACAAGGTATACAGTACCTGatgctgtttccttctcttgtagggttgtgtcaagttgtttttcaagagttgGTGCAGAACATTCCCGAGAAgttccaattattttttatagtaattttgtcattgtttcacaTTATATacgattatcattattattattatcatcatcaccattgTGATACTTGTTTATCATTAGGAGATTGACCATTGCTCGTAATGGAGGTTTTAGCTTGTCTCTTAAATCAAATGAGTAGAAAGCTATTTATCGTTGTTTAATCCTGAAATTTTatacatttgtttgtttttgatagtttaATAGTGGCATTGGATAAGCTGTCCTTTTGTATTTATATGCTAAACAAAGATTAGCTAAAAAGACAATAGAAAAAAGTGTGGGCGCAACAAAAGTCCactgtttaattgttttcttttatcttaaaaacCGTGGAAATATTGTCTTGCTTTGGCTTTTTCCCTCAAGACAACAGCCAGGTATGTTATTGACTATTGTTGGTAAAACATAACGTACTTTTGGGGGAGGATTCTAGTCTGGTTCTTGCGAGAAAACAGTATCCGCGATTGTGGCGCACTTGTTGCATCAGCGTGCTCCGGGACGGTAATTGGTGGATATTCGTCCCTATCCTATCAAACAACCTGTTAAGTCTCTTGGGGATTTTCCTTCTTGAAAGTTCACTCGAGTGTTGTTGTAGCGGACAGGAGCAGAAAGCATTAGCCTCTCGGGTGGCTCTCATCAAGTGGCGAGAACAAGTCCGGCAATTTAAGTTCCCACCCACAAAACGGGCAGTTGCATCGTCAAGTCTCCTCATAAGACCTCGGCATTCTCCTCCAGCTCTggtttttttctgaaaccagTAGATCCTCCCTTCAGGTTGAGCATTTGTCCTCAAAAGTATTATCGCGTATGCGGCACAACAGGACATGATTTTATCGCGAACAAAACTTCGGTGTTCGTTCACTTTTTTCACTGCTTGCACAAATTTTAAACCACTTTGTACAAACAACGACAGGATTGGGCTGTTAGGAGGGGCATATGTTCAAAAGTTCAAGCTAGTGTTTACAGGTTAACTCTAGTTTTTCACGAGAAACCGGGGGCTattcacaaaattttgaaatagccGCCATTTTCAATCGAATTGTTGTCATGTCCAATCTTCGCGCGCCATAACTGtgcacatgcgcagacgttattcagccctgtcGATGGGTAGGGCATACTCCCACACTACATGTCAGCTCGCCCCAGATCCTGTGTGCAAGACTTACACGCCTAGCGATCATTTACCATGCACCAACCAGAAGGTTCCATTGTCCACAATAATATATTGACTTCTCATGTCGTATACTTGAACAAGTAGAGCATGAGTTTCCAGCTGTAATTGGCTGATTTTGTATATGTAATAGGACTACATGctgtccaatttggaaataattgaaTGAGAAAAATTCTGAAGACAGCCAAAATTggacgaggccgtaggccgagtccAATTTGGCAATTATTACCAATTATTTCCTAATTTCCCAATTATTTCTTAATTTCCTAATTGGATATTTCAAATTTAGACAGTTGGCAAAAGTTAATAAAAAATAGATCTGTTGGCAATATTGGATTTTGATGCAGCTATATAGATAATTAAGCAAGAGTACAAGAATGTCTCTTGACAATAACATGCACAatgcaaatacccaaatatacTCAAATTAGAACAAAATTCAGACTCAGAAAAGACTAAGAAAACAAATCTGTTACAATGATCCAGGAGCcaacaataatgttattattcTTTTCCTTCCTGTCTGAAAGGGATTGGTGTTTGAAAGCTAGCTAAAAGTACCAGAAACTGATGTATCAGTATCTCAAGATGAATCAGCATAATTTGTACAATCAATGCTCCAGCAAAGCTTATACTGTACTATATACTGTTTCAGCTGCTGGTATCATTATTTTCACCTCAAGGCTAAGGTAAAAATAATTGAGTGATTAACTTTCCTTATTGAAGCAAGCCATAGCTCAATCACTAAGATGTTATTTCCTTTTAGTGTCCTGATAACTTTGTTCTGGTGTACTGTCTTTTACACAGTATCAAATGACTTTCTGTAAATGTTGAAGTTGGTTGCCTTTAGGTTTTACCTACTACTACCTTTATCCTACTTCTTCTCTTAATTATGGCGCCCacatcataataataaattatgtcACCAATGATAATTTATTACTCTTCTTGTAAGTAAAGATACAGTAGAATCTAAATATAAAAGGCCATCATAATCACTGAGGATTATCCTGTACAGTTGTGACATTACAGAACCTAATAAAAGAGTTTTGGTTTTCACATAAATAACACAAGTCAATACattaattgaaggggtgtgtggaataaggccttaagtgacttttgatgcaatgtcaaattctcctagtcattcacaactgaatacaaggaaatttgcaaggagaatctggtaatttatcagaagtcacttaaggcttctctccaggcacccctgcaattgaaggggtgtgtggaataaggccttaagtgacttttgatgcaatgtcaaattctcctagtcgttcacaactgaatacaaggaaatttggaaggagagtctggtaatttatcagaagtcacttaaggcttttctccaggcatcccagcaatttttcttttaaatccaaCAAAAAATTTTACCAACTGAAATTCATCAACGAAAAGTAACACAAATTTAAAGCAGAAACATGCAGACTTCAAAACAGCTTCAGactaaatttaaattaaaatgcTTAAATATTTACAACAAATCATAATTTTCTTGCTCATCACAAAAAAGTGGACATCTTCATCACAAACCATCCTCACAATAATACAGTCATTACCTAGagaagacaacaacaaaaatctCAAATCAtttcacaaaacacaaaacGTTCCAATACTACACCATTCATTTGTAAGAAGGTTAGTGAGGGCATAGAAGCCACACATCACACAAAGAATTTCGTTCCCGTTACAAATCTGGAAACAGTTTATAAGGACTTAGTTCAGCCATATTTTGAATACTGTTTCCCCCCATGGGACAACTGCGACAAATAACTTAAAGATAAGATCCAAAGATTCCAGTCATGTGCTGCTAGAGTTCTTACAGTTGCTACTGTATTATGATATTCACTCCATAGACTTAATTGATTCTCTTTCTTGGGAAACACTAGATGACAGACAGCGCTATGCAAAGTCGATTTTTATGTTTGACAACATTAAATGATGGCACATCCCCAGCCTAAGAAACTCTTTTGTTAGAAGGAAGGTTGTTCAGGTTAATTaccatctaaaaaaaaagtaacaaagatATAAAGACCTGACACTACCTAAATCGTAAAGGGGATTTTtgaaaagaagttttaaatttAGTGGTGCTTATGCAGTGGAACCAGCTCTTGAATTAAACAAAACTTGTGAGTCAATCTCTTCATTTAAGAAGCTGACTACAAAATAGTTGGGTCATGACaagatatatatttttagacTAGTTAACTTTTATctcttttattgttattattattatggttattaTTACTAGCACGAAAACGAGAACAAACAGATGTATACGAGCAGGGTGTTAGAAGTGGAGCAAGGGACTTTCAGCCTATTAATGTTCACCATCACTGGAGGCATGCTGGACAAGTGCAAACACTATCACAGTAGAATCACCAAACTCATGTCTATCAAGAAAGGGAAGGATTACAGCACCACCATGGCATGGATAAGATCTAAAGTATCTTTCAGCTTGCTTACATCTGCTCTCCTCTGCCTACCAGGTTCACACACTACAAGGCGTGTCCCTCTGAACATTCAAGAGCACGACTTTGTTGTGGATAAAGAACTGGTGGGACTGGGggattaataaattattatgaactTTATTATGgcttctgtttttcttttcagattaagtgaaaaattttcataaatacaatttaattttttctatATTCTTAATTACAAAATGACAAGTCAAGTTTTCATTATAAAATCAAAGGTGTCAAACAATGtaacaatattttaaaataataggTATAATTATAGgagtttttatttaaattttttattattagcaataaagtaattggaccgagtggagtacaattcagggagtaatcactccagtaatttcaaaattggaCAAGTGCCAAGCTCGAGGCCaactttgaaattcaaatttgattttgaaaactcAAGTATTACCTCTGTCCCTACCTAACTCCAGTCCTTACCAAAAATCAATCAACAGCTTATTGATTTCATCTATATCTCCTAAAGCACGCTGCTTCAAACCATCATAATCACTTTCCAAACTCTTCTGTGCAATaactaaaaatacaaaattgttAACAACCATTACATCAACTGCATATATGAAGCACAAAATTATACATAAAACTACGTCACTTCCTGTTGGATTACAGAGTAGCTTAAAAGAACTACCTAAATATCTTTGACAACAAAATAAATCAGTTTTTAAAAGGTTAAATTAGTACAAAAATGTTTGTATAATATTTTTTAGTAAattccaactagtggtctattatcaatgctgccttctgattggttgagctactactaggctatattaTGTTATAGCCCCACTAGTTGGGAAAAGCGCCAGCCATAATTGAATGttttgacagaaaaaaaaaggattaaagtccaGCTTTAACTGCAAAAAGATGTTTTgcctcaatatttttttgagcaACTACTTGTATTTTACTacaacaattattcctctcgccctcatggcttCTGAGTAAATAGCCCATCCaaccttcggcctcatgggctattgactcagagcccaggGACAACTTAACTGACACAAACAAAGGATCCCTTTGAAATTCCAACATCAAAAAAGTGTGATTATATATGAAAAGTGTTATcatcaaaaaattgcaaaaacaacACATACATTCCTTCATGACAAAATTATTCTGCTCTAGGTGACACCATTTCCTCTCCAAATTCccaagctgaaaaaaaaagctcattattgaaaaaatacacatgaaaatgaaaacaacaaaacgtcATGCAACAACtaaatttataaaataataattattagaatagTATGCactctctcattggtcaatggGTGTGCTCAGATGAGAGTATATAGACACAGTTGTGACTTGATTGGTTGTGACTTGTTTCATGCACATTTGGTTGGCTGGTAGGAAATATGAACGCATATCCAAAAAATCTATTTCAATCAAGAAGTAAAATAAACAGCATTATCCTTCATTtgccgaatttttttttttatgagagaAGTAtcttacaaaaattaatgctaCACAGAACGTTTTTCTGTGTTGACATAGACTCTAAACACACAGGAAGTTGGAAGAACTATCAACAGTTACCAACACTGTGAACTGCGTCTgaggtttgcataactgtctcaaaCTCGGTGTTTggatgaggctatgtaaacacagaAAAAGACCTCTATTGCTTAAATTCAAACTTCCAACAACAAAACCTACCTGAGAATGAGTCTCATTTTCAATCAATTTACTCCTTGATGCATCATAAACAGCAGACAGCTTTTTAACATctgcaaggaaacaaaaaataataataaaaataataattgttaattaacctatagttcaattcaatttttcacAAGAACGcaattttacaaaaaaatttacatttcaTGTCTAGGTTTGTCCAGTAGTCcacacttctttttgtttttgttctcacttgtttcttagttcctcaataaactctaCGTCGGGttcaacaaaacgggaagccgTATTTGCAGAAGATTGTAATGAacaacaaatcttagcaataaccttgttgctaagcaactttaaaccaatcaggatcaagtaatCATCCCCTCTTGATTACTAAAAGTGCCTCATGTGATTAGgaaaaaaatgccctctgtctcagccagtCAGCCACTCAGTCATTTTTAAatgagtaaaattaaggattaatatcacgcgtgttttcagaagttgctgaaattacccgagtcgcgcatccttaattttacgaggatccattgcgattactgtaattttgccctcttcacgaagcaaaattaagaaaaaatactctcttcattgaccaatcagcattcagtaattttgtcctctatgttattaaaaaTCTAACAGGTTCAGTTGTTTCTTCTATATGCATTAAAACGTTGTTtatcattttacattttcagCAGAAccctcgaccaatcagattgctggaataAGGACATGTGACGGTCATACAGAGCGGGACAAATATTTTACTCAACTTGAAAACAGTGGATCCACTTTTCTTGCTGGCGCCAAAGCCAATCATATTACAGGATTTAGCGCACGTGACTTTTgattttgaaaggaaaacaaggaaaaacaataCATGGACTAAATGAGAAACAATGGTGTCTTCCCGAGGTAGGTGTTACACTATTGTTATATATTATgtattatgaaagaaatgttatatgcagtgcggtgtttgaaatcaaatgaagatatgatcctcgcacttgctggacaatttaagcaaatgtctcatgaacctgaaaaattcaggtgactcaacgggatttgaacccatgacctctgcgatgccggtgcagtgctctaaccaactgagctatgaagtcacacggtcatgttttcccgtgaaaggaatgtcatatgaaagaaatgttatatgcagcaagtgcgaggatcatatcttcatatttgatttcaaacaccgcactgcatataacatttctttcatagaaaaaaaattattgataaaaatTCCTAAAGAGCAAATTGTTTCACATACATGTAACTTGCAACAAGAGAAATAATACCTCTGACTTGAAGGAGATAATTCCTCTCTGCctttttcctgaaaaaaattaaaacaatacttTTTCAATACTTTCCCACAACCTAACTTGTGCAAGGCCGCTGTATTTGTTAAAATTATTTACCAACATATGCCCCTTACATAAAGTTTTACTTCCTTGACCTTTTCTTGGACTACTTTCACCCCTCACAATACAAGTGGACCTTTAAGAGGCCGATATATTAAGGTATTTTTCGTTatcatttgaatttttttcgtaaaaaCCAGTCAGATTGCGGTATATAGATCACCTGATTTTGACTGACCAATATTAAAGCGAGAAAATTACAATTGATGTTCTATACGGTTTTAAGTCTGGTTTcctcattgtttacattttctaTCTAATTTATGCATAgtccaaccaatcagattaaAGCATTTACCAATCAGGAAGCAGGAATTTTAATTGATGTACCATACGGTTTTTCACTCGTTTTCCCTGTTTCCTAATTGTTTACGTTTTCTCGACAATTTATGCATAATCGATCCAATCAGATTTGAGCATTTAGATGTGATcaaaactgaccaatcagaaagcgtgaattttgcttccttcatcggtagcaaaaaaaaaatgcaaattccagatttctcgctttctgatTGGCCAGTTTCTGggcacatgatttttttttcttacacttTTCGTATGTTTTTCCACATTTTTCTCCTCCCCTCCTTCACCCCTCCACCCCTTCTTCACCCCTTCACCCCTCCACCCTTGCGTCTTGGTCTTGGCACTAACCgtaaaccaataataataataataataataataataataactttattagcGAGTCAAGTAAAATAGAAGTTTCCCACTAAGTAAGGACAcctatctaaaaaaaaaaactagaagtACCCGTATAATCCCTATATGATCCCCTCAATAATCCCAcccacaatttaaaattaattacaatgttaAGAAAGACAAAGAGTACAGTTAATACAATTATTAGCTAAAATATGTTTAGCAAGATCTACCATCCTAATATAACGTTTTTAgttctctgaatttcctatcaATCTTAGACCAGAGCACCGGTCCTAAGTATCTGACTGAATGCTTACCATAACCCGTGGTGTTAACTCTAGGAACTACAAAATCGTTATTTCTTAAGTTATACTGATTACTTCTAAAAATAAACAACCTATAAAGATAATTTGGACATAAGCCGTTCTTAATCTTATACATTAAAATTGCAATGTCTTGTAACCTTCTATTGTATAAGGTTGGTAATTTCGCCCTTTTACTTCATGTCCTCTGTCTTCAACTGTCGTTCACTCAATGTGCGCTGTACTTACTTGCTATTTTCTTCGTATCCATTGTTCTTTCATGTTGTTTACTTCATGTCCTCTGTCTTCCCTTGCCCTTTTCTTCCTGTACACTGTTCACTTGCTTGAACACTTGACATTGGAATGAAATTTCGAAATGTCAGCAAAGCAacacattcatttttttttttgcctggaGACTGTGGAAAATTTTAACAGATGTACAGTATTGGCATGTACCAAGCCTAAACCATCTATATTCCTTGAGGTATATCATGCCTTAAATGTTGAAACTGTTTGTTTATTATTGAATTGCAGATGTGGGAATGTTGCAGCTATTTTAGAGTTGGACCAGTGTTCTCCcaaagttttagctcagcaggtaagggacaattcctgaccggtatatTTTTTATACAACTGATATAGTTtgagtaaaccttcaagaggttgcaggcggtaagaacagACTgttactgttgcttgaggcggtaaattttactggttaccgcttgataaggagaacactggttGGATGAGCATTTAAAAAGAGAATTTACAATCTTCGAAGCAGCGCCACAGGTTGGGGTTGGTTTTTACAAAAGGGGTTCATATTACTGAGAGTTCTTAGTCGGTCAAAAGAAAATCAGGGAACAGCATATTTGACTTtaagtgaaaaatgaatctttccAAAGGCCTGCGATTTGCAGAAGCAAGCACCTTTTAGTAGTAGTATTGGGGGAGGGGGCCTTGCATTTATATTTGTGCAGCTACTGCTTTTAGCACGTGATGGTATACTTTTTATGTGaatattgtttgttttggtttgtgaCTTTCAATTCTGCGTGAAGACttttaaaatagttttctttaTATAAATAGTGCCCTTCTACTTTGCTCGCCCATAGTGTCAGGAACATGATACCATGCTTTTGAACGAAGGGCTTTTCTCATCTATGGTACACTTTTATTGAATACTTCATACTCTTTATAATcatgatattaataattatttgattttattcCATAGGAGGTGAGAGGTATGCCAACTATTTCTAGGAATCCGCAGCCACATTACTTCCTTTGAATCTGCTACAGAAGTGTCCTTGGTCAACTTTTTTTGGAGATTTCCTTTCCTCTAAGCCACCGTTCAGCTATTAGGTGTGGATACTTCAGTGAACCGTGGGATCCGAGATATTGAGCATTGCACAGACGTCGAATATAGCTTGCAAGGCAATACAAATGGCTTTCGAGAGCACAAGCATAACTACATGGTGGTCttaaaagaacaataattgaTAGGCCttttgtggttttgtttttcttttttatatttagttttggaaaaagaaattcataGTTACAATTAGGAGATAACTGTATAATATACAACtacccgaaggggaggtgaatagtggtggatatatatatagtgaatagtggtggatatacatatccaccactcttcaccgaccctgagggaatagttgttttagtatttaccaaatcagatggataaaaaaacgcttcttcaatttcttcttctgaaactttcgcgaaacgacatttttctctccgttcgcaaaacagtgaatatccaaggatattccgagttacgggagccaatcagaacgcgcgaaaattgctatccactgatttggtagaTACTAAACTTGATTATTTGGGCTAACATTGTATATACCATACACTTTTATAATTGAAGtgaaaacattaatttattttacaaataaCTCAGTCATATGTTCTGCTCGTGGGAAACGTGCGACCAGGGTTACCTCCTGTTGAAAGACTAGTATCTAGTTTTGATTTGCTGGAAAGCCTGATACGCTGTATTCTTTTAATGCAAAGTGCTTAGTCTTCAATCCTTCTGTTGTTATCTGAGTCAATCATAAAACGTGTTCGAGTTTAGTATAGGAAGTGAAGCGATGCTTGAAAGACTTCTTAGGTTCTTTGGGGGAAATACTCATTTTTGAAAAATTCCCATCTCGATTCATTCTTGTTGTGAAGACTTTGAGAATAGCTAAGTGATGTCACTTGAATGGTACACAAAAAAAGCCTGAAGGGCAAGTTATTTCGCGATACGCACGCAGACGAGCAGGGCAACACGTCTCTCGAGGTGATACAGTGGTCTCGCGAGAAGGAAGTAACTTACTTTGGAGTGTACGGATTACTGGTTAAAAACTTATTTCTTGTAATAAAGGCTGTGATCGTCAGCTTGTGATTATTGCTCCAATATATAATGaacagtacatttttttt
    >transcript::NC_058066.1:1153398-1165634
    GCGAGCGTGATCTTCAGTTCGTTTCGAGGTTCAAAAAAACACATTGCCAACAACAGTGGAGTAAAGGGGTATTATTCTTGCACAGTACCTCAACAGTTCTTGCTTTTGcttcacaaacaataatatgGATTCGTACCAGttggattattattttttgggaAAGAAATGTCGGATGAAACATGAATTACGACCTACTTTGTTTACGGATCTCGTCCGTTGCAACTCAAGCGTTTTCACTGGTTTTCGTCATCAAATAACCGTTCTTGGCTTCGTCGAAATAAACCGTTGGCGAAGAAAATTGTCATGATCCATTGACAGCGTTATGGAAGCGTTATGTAAGCTAAAACAAGTGACACAGGAGGTAAGCACGTAAAGCTTAAACACATTCGATTCATCGCCATGGTTATGAACCTATTATGCgcttaagtttaagttttattgCAAAGGATGTAATCTTTTCCAAAATCTGTGGGAACTCCTCACGACAATTTGCAAGCTCAGGCAAAGCATTTTTGCAACGTTTGTCGTTTGCCAGcattcataatattatttactaTTAAAGATAATGGCAGTGGCATAAAACTGAAGTGTTTCcgaaattttcaaatgtgttgTGTTTCTCTTAAAGCTAACGTCTCACTAATGGAATGTAACACATGTTACTCTGCCCAAGCAATATGAGGAATGTTATTAATGGCTTTGAACCATGAGAGCAATTCCCATGCTGCATCAGCTGGGGCCATAGTGAGAAATGCACCACTTGTTCCATGATAATAGTTGTCAGTTATTCCATTTCGTGTTGGGatgaaagccacagttttgtcaTTGATGATGTAGAGACAAGCTGTCTTGGTCATGTTACGGGCTTGTAGAAAGTAATATGGAAGACATGCTGATTGTATCTGCTCTCTGGTGATACTTGCTGGTAACTCTCCAGACATGTGTATACTGCCAAGAACTCTGGCCTGCTGCATTTTGGTAGCTGAATCTTGGACACCAAACATCTGCCCAGAGTTGCCAGCAAATTTGTCATAGAACTTGTTGCCAATCATAATACCGGCAAGGGTGCGATGCATCCAGGTGTTGCTCAATGGCTGGCTACTGTTGTGAGTTTCAGGCGAGGAGAAATATAGTTTACTGTATGAGAGTGCAATGAAAGTGCATGCATTGCTACCTAATCGGCCAGTCAAAGTAGATTGGGAGTATTGGGGTGGAAAGTGCCAGCTTATGACTCTGTCAGATGTTGAAGGATTTATTAGTAGATGCTGAGTGACGGAATGGATTTCTACTGTTGAGGACTGAGGTTGGTTGTGAGCATGAGAGGCTTGTGCATTTGCATGTTTCGGTTGATCAATGCCCCCAATTCCTGTTAGGATATTGTCTACTTTCTGCTCCCAAACTCCAGGCCCgggttgttcaaacgatggatagcactatccaccggataaatcgctatccacaggataagtaatagcgaaaccaattattgcgatatccaatggatagtgatttatcaggtggatagcgttatccaccttttgaacaactggggccagatgtaTAATACTGCTCAGCATTTGTATCAGTGACTTCATTGTTTCCGTCAATTGTATCATCATCAGTGTCAGACTCACTTAGGTTGTCATCAGTGTCTGTATTTTGTacgtttttaaagcaaatgtcaTTTCTAATTTCGTAGCATGTATCTCCAGTTTACTCAAACAAGGTATACAGTACCTGatgctgtttccttctcttgtagggttgtgtcaagttgtttttcaagagttgGTGCAGAACATTCCCGAGAAgttccaattattttttatagtaattttgtcattgtttcacaTTATATacgattatcattattattattatcatcatcaccattgTGATACTTGTTTATCATTAGGAGATTGACCATTGCTCGTAATGGAGGTTTTAGCTTGTCTCTTAAATCAAATGAGTAGAAAGCTATTTATCGTTGTTTAATCCTGAAATTTTatacatttgtttgtttttgatagtttaATAGTGGCATTGGATAAGCTGTCCTTTTGTATTTATATGCTAAACAAAGATTAGCTAAAAAGACAATAGAAAAAAGTGTGGGCGCAACAAAAGTCCactgtttaattgttttcttttatcttaaaaacCGTGGAAATATTGTCTTGCTTTGGCTTTTTCCCTCAAGACAACAGCCAGGTATGTTATTGACTATTGTTGGTAAAACATAACGTACTTTTGGGGGAGGATTCTAGTCTGGTTCTTGCGAGAAAACAGTATCCGCGATTGTGGCGCACTTGTTGCATCAGCGTGCTCCGGGACGGTAATTGGTGGATATTCGTCCCTATCCTATCAAACAACCTGTTAAGTCTCTTGGGGATTTTCCTTCTTGAAAGTTCACTCGAGTGTTGTTGTAGCGGACAGGAGCAGAAAGCATTAGCCTCTCGGGTGGCTCTCATCAAGTGGCGAGAACAAGTCCGGCAATTTAAGTTCCCACCCACAAAACGGGCAGTTGCATCGTCAAGTCTCCTCATAAGACCTCGGCATTCTCCTCCAGCTCTggtttttttctgaaaccagTAGATCCTCCCTTCAGGTTGAGCATTTGTCCTCAAAAGTATTATCGCGTATGCGGCACAACAGGACATGATTTTATCGCGAACAAAACTTCGGTGTTCGTTCACTTTTTTCACTGCTTGCACAAATTTTAAACCACTTTGTACAAACAACGACAGGATTGGGCTGTTAGGAGGGGCATATGTTCAAAAGTTCAAGCTAGTGTTTACAGGTTAACTCTAGTTTTTCACGAGAAACCGGGGGCTattcacaaaattttgaaatagccGCCATTTTCAATCGAATTGTTGTCATGTCCAATCTTCGCGCGCCATAACTGtgcacatgcgcagacgttattcagccctgtcGATGGGTAGGGCATACTCCCACACTACATGTCAGCTCGCCCCAGATCCTGTGTGCAAGACTTACACGCCTAGCGATCATTTACCATGCACCAACCAGAAGGTTCCATTGTCCACAATAATATATTGACTTCTCATGTCGTATACTTGAACAAGTAGAGCATGAGTTTCCAGCTGTAATTGGCTGATTTTGTATATGTAATAGGACTACATGctgtccaatttggaaataattgaaTGAGAAAAATTCTGAAGACAGCCAAAATTggacgaggccgtaggccgagtccAATTTGGCAATTATTACCAATTATTTCCTAATTTCCCAATTATTTCTTAATTTCCTAATTGGATATTTCAAATTTAGACAGTTGGCAAAAGTTAATAAAAAATAGATCTGTTGGCAATATTGGATTTTGATGCAGCTATATAGATAATTAAGCAAGAGTACAAGAATGTCTCTTGACAATAACATGCACAatgcaaatacccaaatatacTCAAATTAGAACAAAATTCAGACTCAGAAAAGACTAAGAAAACAAATCTGTTACAATGATCCAGGAGCcaacaataatgttattattcTTTTCCTTCCTGTCTGAAAGGGATTGGTGTTTGAAAGCTAGCTAAAAGTACCAGAAACTGATGTATCAGTATCTCAAGATGAATCAGCATAATTTGTACAATCAATGCTCCAGCAAAGCTTATACTGTACTATATACTGTTTCAGCTGCTGGTATCATTATTTTCACCTCAAGGCTAAGGTAAAAATAATTGAGTGATTAACTTTCCTTATTGAAGCAAGCCATAGCTCAATCACTAAGATGTTATTTCCTTTTAGTGTCCTGATAACTTTGTTCTGGTGTACTGTCTTTTACACAGTATCAAATGACTTTCTGTAAATGTTGAAGTTGGTTGCCTTTAGGTTTTACCTACTACTACCTTTATCCTACTTCTTCTCTTAATTATGGCGCCCacatcataataataaattatgtcACCAATGATAATTTATTACTCTTCTTGTAAGTAAAGATACAGTAGAATCTAAATATAAAAGGCCATCATAATCACTGAGGATTATCCTGTACAGTTGTGACATTACAGAACCTAATAAAAGAGTTTTGGTTTTCACATAAATAACACAAGTCAATACattaattgaaggggtgtgtggaataaggccttaagtgacttttgatgcaatgtcaaattctcctagtcattcacaactgaatacaaggaaatttgcaaggagaatctggtaatttatcagaagtcacttaaggcttctctccaggcacccctgcaattgaaggggtgtgtggaataaggccttaagtgacttttgatgcaatgtcaaattctcctagtcgttcacaactgaatacaaggaaatttggaaggagagtctggtaatttatcagaagtcacttaaggcttttctccaggcatcccagcaatttttcttttaaatccaaCAAAAAATTTTACCAACTGAAATTCATCAACGAAAAGTAACACAAATTTAAAGCAGAAACATGCAGACTTCAAAACAGCTTCAGactaaatttaaattaaaatgcTTAAATATTTACAACAAATCATAATTTTCTTGCTCATCACAAAAAAGTGGACATCTTCATCACAAACCATCCTCACAATAATACAGTCATTACCTAGagaagacaacaacaaaaatctCAAATCAtttcacaaaacacaaaacGTTCCAATACTACACCATTCATTTGTAAGAAGGTTAGTGAGGGCATAGAAGCCACACATCACACAAAGAATTTCGTTCCCGTTACAAATCTGGAAACAGTTTATAAGGACTTAGTTCAGCCATATTTTGAATACTGTTTCCCCCCATGGGACAACTGCGACAAATAACTTAAAGATAAGATCCAAAGATTCCAGTCATGTGCTGCTAGAGTTCTTACAGTTGCTACTGTATTATGATATTCACTCCATAGACTTAATTGATTCTCTTTCTTGGGAAACACTAGATGACAGACAGCGCTATGCAAAGTCGATTTTTATGTTTGACAACATTAAATGATGGCACATCCCCAGCCTAAGAAACTCTTTTGTTAGAAGGAAGGTTGTTCAGGTTAATTaccatctaaaaaaaaagtaacaaagatATAAAGACCTGACACTACCTAAATCGTAAAGGGGATTTTtgaaaagaagttttaaatttAGTGGTGCTTATGCAGTGGAACCAGCTCTTGAATTAAACAAAACTTGTGAGTCAATCTCTTCATTTAAGAAGCTGACTACAAAATAGTTGGGTCATGACaagatatatatttttagacTAGTTAACTTTTATctcttttattgttattattattatggttattaTTACTAGCACGAAAACGAGAACAAACAGATGTATACGAGCAGGGTGTTAGAAGTGGAGCAAGGGACTTTCAGCCTATTAATGTTCACCATCACTGGAGGCATGCTGGACAAGTGCAAACACTATCACAGTAGAATCACCAAACTCATGTCTATCAAGAAAGGGAAGGATTACAGCACCACCATGGCATGGATAAGATCTAAAGTATCTTTCAGCTTGCTTACATCTGCTCTCCTCTGCCTACCAGGTTCACACACTACAAGGCGTGTCCCTCTGAACATTCAAGAGCACGACTTTGTTGTGGATAAAGAACTGGTGGGACTGGGggattaataaattattatgaactTTATTATGgcttctgtttttcttttcagattaagtgaaaaattttcataaatacaatttaattttttctatATTCTTAATTACAAAATGACAAGTCAAGTTTTCATTATAAAATCAAAGGTGTCAAACAATGtaacaatattttaaaataataggTATAATTATAGgagtttttatttaaattttttattattagcaataaagtaattggaccgagtggagtacaattcagggagtaatcactccagtaatttcaaaattggaCAAGTGCCAAGCTCGAGGCCaactttgaaattcaaatttgattttgaaaactcAAGTATTACCTCTGTCCCTACCTAACTCCAGTCCTTACCAAAAATCAATCAACAGCTTATTGATTTCATCTATATCTCCTAAAGCACGCTGCTTCAAACCATCATAATCACTTTCCAAACTCTTCTGTGCAATaactaaaaatacaaaattgttAACAACCATTACATCAACTGCATATATGAAGCACAAAATTATACATAAAACTACGTCACTTCCTGTTGGATTACAGAGTAGCTTAAAAGAACTACCTAAATATCTTTGACAACAAAATAAATCAGTTTTTAAAAGGTTAAATTAGTACAAAAATGTTTGTATAATATTTTTTAGTAAattccaactagtggtctattatcaatgctgccttctgattggttgagctactactaggctatattaTGTTATAGCCCCACTAGTTGGGAAAAGCGCCAGCCATAATTGAATGttttgacagaaaaaaaaaggattaaagtccaGCTTTAACTGCAAAAAGATGTTTTgcctcaatatttttttgagcaACTACTTGTATTTTACTacaacaattattcctctcgccctcatggcttCTGAGTAAATAGCCCATCCaaccttcggcctcatgggctattgactcagagcccaggGACAACTTAACTGACACAAACAAAGGATCCCTTTGAAATTCCAACATCAAAAAAGTGTGATTATATATGAAAAGTGTTATcatcaaaaaattgcaaaaacaacACATACATTCCTTCATGACAAAATTATTCTGCTCTAGGTGACACCATTTCCTCTCCAAATTCccaagctgaaaaaaaaagctcattattgaaaaaatacacatgaaaatgaaaacaacaaaacgtcATGCAACAACtaaatttataaaataataattattagaatagTATGCactctctcattggtcaatggGTGTGCTCAGATGAGAGTATATAGACACAGTTGTGACTTGATTGGTTGTGACTTGTTTCATGCACATTTGGTTGGCTGGTAGGAAATATGAACGCATATCCAAAAAATCTATTTCAATCAAGAAGTAAAATAAACAGCATTATCCTTCATTtgccgaatttttttttttatgagagaAGTAtcttacaaaaattaatgctaCACAGAACGTTTTTCTGTGTTGACATAGACTCTAAACACACAGGAAGTTGGAAGAACTATCAACAGTTACCAACACTGTGAACTGCGTCTgaggtttgcataactgtctcaaaCTCGGTGTTTggatgaggctatgtaaacacagaAAAAGACCTCTATTGCTTAAATTCAAACTTCCAACAACAAAACCTACCTGAGAATGAGTCTCATTTTCAATCAATTTACTCCTTGATGCATCATAAACAGCAGACAGCTTTTTAACATctgcaaggaaacaaaaaataataataaaaataataattgttaattaacctatagttcaattcaatttttcacAAGAACGcaattttacaaaaaaatttacatttcaTGTCTAGGTTTGTCCAGTAGTCcacacttctttttgtttttgttctcacttgtttcttagttcctcaataaactctaCGTCGGGttcaacaaaacgggaagccgTATTTGCAGAAGATTGTAATGAacaacaaatcttagcaataaccttgttgctaagcaactttaaaccaatcaggatcaagtaatCATCCCCTCTTGATTACTAAAAGTGCCTCATGTGATTAGgaaaaaaatgccctctgtctcagccagtCAGCCACTCAGTCATTTTTAAatgagtaaaattaaggattaatatcacgcgtgttttcagaagttgctgaaattacccgagtcgcgcatccttaattttacgaggatccattgcgattactgtaattttgccctcttcacgaagcaaaattaagaaaaaatactctcttcattgaccaatcagcattcagtaattttgtcctctatgttattaaaaaTCTAACAGGTTCAGTTGTTTCTTCTATATGCATTAAAACGTTGTTtatcattttacattttcagCAGAAccctcgaccaatcagattgctggaataAGGACATGTGACGGTCATACAGAGCGGGACAAATATTTTACTCAACTTGAAAACAGTGGATCCACTTTTCTTGCTGGCGCCAAAGCCAATCATATTACAGGATTTAGCGCACGTGACTTTTgattttgaaaggaaaacaaggaaaaacaataCATGGACTAAATGAGAAACAATGGTGTCTTCCCGAGGTAGGTGTTACACTATTGTTATATATTATgtattatgaaagaaatgttatatgcagtgcggtgtttgaaatcaaatgaagatatgatcctcgcacttgctggacaatttaagcaaatgtctcatgaacctgaaaaattcaggtgactcaacgggatttgaacccatgacctctgcgatgccggtgcagtgctctaaccaactgagctatgaagtcacacggtcatgttttcccgtgaaaggaatgtcatatgaaagaaatgttatatgcagcaagtgcgaggatcatatcttcatatttgatttcaaacaccgcactgcatataacatttctttcatagaaaaaaaattattgataaaaatTCCTAAAGAGCAAATTGTTTCACATACATGTAACTTGCAACAAGAGAAATAATACCTCTGACTTGAAGGAGATAATTCCTCTCTGCctttttcctgaaaaaaattaaaacaatacttTTTCAATACTTTCCCACAACCTAACTTGTGCAAGGCCGCTGTATTTGTTAAAATTATTTACCAACATATGCCCCTTACATAAAGTTTTACTTCCTTGACCTTTTCTTGGACTACTTTCACCCCTCACAATACAAGTGGACCTTTAAGAGGCCGATATATTAAGGTATTTTTCGTTatcatttgaatttttttcgtaaaaaCCAGTCAGATTGCGGTATATAGATCACCTGATTTTGACTGACCAATATTAAAGCGAGAAAATTACAATTGATGTTCTATACGGTTTTAAGTCTGGTTTcctcattgtttacattttctaTCTAATTTATGCATAgtccaaccaatcagattaaAGCATTTACCAATCAGGAAGCAGGAATTTTAATTGATGTACCATACGGTTTTTCACTCGTTTTCCCTGTTTCCTAATTGTTTACGTTTTCTCGACAATTTATGCATAATCGATCCAATCAGATTTGAGCATTTAGATGTGATcaaaactgaccaatcagaaagcgtgaattttgcttccttcatcggtagcaaaaaaaaaatgcaaattccagatttctcgctttctgatTGGCCAGTTTCTGggcacatgatttttttttcttacacttTTCGTATGTTTTTCCACATTTTTCTCCTCCCCTCCTTCACCCCTCCACCCCTTCTTCACCCCTTCACCCCTCCACCCTTGCGTCTTGGTCTTGGCACTAACCgtaaaccaataataataataataataataataataataactttattagcGAGTCAAGTAAAATAGAAGTTTCCCACTAAGTAAGGACAcctatctaaaaaaaaaaactagaagtACCCGTATAATCCCTATATGATCCCCTCAATAATCCCAcccacaatttaaaattaattacaatgttaAGAAAGACAAAGAGTACAGTTAATACAATTATTAGCTAAAATATGTTTAGCAAGATCTACCATCCTAATATAACGTTTTTAgttctctgaatttcctatcaATCTTAGACCAGAGCACCGGTCCTAAGTATCTGACTGAATGCTTACCATAACCCGTGGTGTTAACTCTAGGAACTACAAAATCGTTATTTCTTAAGTTATACTGATTACTTCTAAAAATAAACAACCTATAAAGATAATTTGGACATAAGCCGTTCTTAATCTTATACATTAAAATTGCAATGTCTTGTAACCTTCTATTGTATAAGGTTGGTAATTTCGCCCTTTTACTTCATGTCCTCTGTCTTCAACTGTCGTTCACTCAATGTGCGCTGTACTTACTTGCTATTTTCTTCGTATCCATTGTTCTTTCATGTTGTTTACTTCATGTCCTCTGTCTTCCCTTGCCCTTTTCTTCCTGTACACTGTTCACTTGCTTGAACACTTGACATTGGAATGAAATTTCGAAATGTCAGCAAAGCAacacattcatttttttttttgcctggaGACTGTGGAAAATTTTAACAGATGTACAGTATTGGCATGTACCAAGCCTAAACCATCTATATTCCTTGAGGTATATCATGCCTTAAATGTTGAAACTGTTTGTTTATTATTGAATTGCAGATGTGGGAATGTTGCAGCTATTTTAGAGTTGGACCAGTGTTCTCCcaaagttttagctcagcaggtaagggacaattcctgaccggtatatTTTTTATACAACTGATATAGTTtgagtaaaccttcaagaggttgcaggcggtaagaacagACTgttactgttgcttgaggcggtaaattttactggttaccgcttgataaggagaacactggttGGATGAGCATTTAAAAAGAGAATTTACAATCTTCGAAGCAGCGCCACAGGTTGGGGTTGGTTTTTACAAAAGGGGTTCATATTACTGAGAGTTCTTAGTCGGTCAAAAGAAAATCAGGGAACAGCATATTTGACTTtaagtgaaaaatgaatctttccAAAGGCCTGCGATTTGCAGAAGCAAGCACCTTTTAGTAGTAGTATTGGGGGAGGGGGCCTTGCATTTATATTTGTGCAGCTACTGCTTTTAGCACGTGATGGTATACTTTTTATGTGaatattgtttgttttggtttgtgaCTTTCAATTCTGCGTGAAGACttttaaaatagttttctttaTATAAATAGTGCCCTTCTACTTTGCTCGCCCATAGTGTCAGGAACATGATACCATGCTTTTGAACGAAGGGCTTTTCTCATCTATGGTACACTTTTATTGAATACTTCATACTCTTTATAATcatgatattaataattatttgattttattcCATAGGAGGTGAGAGGTATGCCAACTATTTCTAGGAATCCGCAGCCACATTACTTCCTTTGAATCTGCTACAGAAGTGTCCTTGGTCAACTTTTTTTGGAGATTTCCTTTCCTCTAAGCCACCGTTCAGCTATTAGGTGTGGATACTTCAGTGAACCGTGGGATCCGAGATATTGAGCATTGCACAGACGTCGAATATAGCTTGCAAGGCAATACAAATGGCTTTCGAGAGCACAAGCATAACTACATGGTGGTCttaaaagaacaataattgaTAGGCCttttgtggttttgtttttcttttttatatttagttttggaaaaagaaattcataGTTACAATTAGGAGATAACTGTATAATATACAACtacccgaaggggaggtgaatagtggtggatatatatatagtgaatagtggtggatatacatatccaccactcttcaccgaccctgagggaatagttgttttagtatttaccaaatcagatggataaaaaaacgcttcttcaatttcttcttctgaaactttcgcgaaacgacatttttctctccgttcgcaaaacagtgaatatccaaggatattccgagttacgggagccaatcagaacgcgcgaaaattgctatccactgatttggtagaTACTAAACTTGATTATTTGGGCTAACATTGTATATACCATACACTTTTATAATTGAAGtgaaaacattaatttattttacaaataaCTCAGTCATATGTTCTGCTCGTGGGAAACGTGCGACCAGGGTTACCTCCTGTTGAAAGACTAGTATCTAGTTTTGATTTGCTGGAAAGCCTGATACGCTGTATTCTTTTAATGCAAAGTGCTTAGTCTTCAATCCTTCTGTTGTTATCTGAGTCAATCATAAAACGTGTTCGAGTTTAGTATAGGAAGTGAAGCGATGCTTGAAAGACTTCTTAGGTTCTTTGGGGGAAATACTCATTTTTGAAAAATTCCCATCTCGATTCATTCTTGTTGTGAAGACTTTGAGAATAGCTAAGTGATGTCACTTGAATGGTACACAAAAAAAGCCTGAAGGGCAAGTTATTTCGCGATACGCACGCAGACGAGCAGGGCAACACGTCTCTCGAGGTGATACAGTGGTCTCGCGAGAAGGAAGTAACTTACTTTGGAGTGTACGGATTACTGGTTAAAAACTTATTTCTTGTAATAAAGGCTGTGATCGTCAGCTTGTGATTATTGCTCCAATATATAATGaacagtacatttttttt
    16982

17098 original number

# 7 subsetting fasta

``` bash
/home/shared/samtools-1.12/samtools faidx ../output/05.33-lncRNA-discovery/Apul_lncRNA_candidates.fasta \
-r ../output/05.33-lncRNA-discovery/Apul_noncoding_transcripts_ids.txt > ../output/05.33-lncRNA-discovery/Apul_lncRNA.fasta
```

``` bash
fgrep -c ">" ../output/05.33-lncRNA-discovery/Apul_lncRNA.fasta
fgrep ">" ../output/05.33-lncRNA-discovery/Apul_lncRNA.fasta | head -5

head ../output/05.33-lncRNA-discovery/Apul_lncRNA.fasta
```

    16206
    >transcript::NC_058066.1:468618-469943
    >transcript::NC_058066.1:1135315-1144814
    >transcript::NC_058066.1:1144883-1148491
    >transcript::NC_058066.1:1153398-1165634
    >transcript::NC_058066.1:1153398-1165634
    >transcript::NC_058066.1:468618-469943
    taactgatcaaaacgtatcttcctacaacattaatttgacagtggcgtttctcaactgac
    caatcaaaacttacatttgaaaatttggtgATGGTgcgtttacaactcgtgtatctttac
    gtcacacaaccatgtttgcATACTCTCTTGCaaccacgcctctcggccaatcagagcgcg
    cgcgtactatcttagttattttataaagatAAATACGCCCTAGGATTAGCACGCACGCTA
    TGGTATAATTATTGATGATAACTTTGCTGGATTTACGTTTGGTTGAAGTTATCATGATAT
    tccatcgtcgtcatcatcaacATTCTTATCGTTTATCTTCATCACAATCACCTGACACAA
    CATGACTAAAAGCAAAGATGAAAACACTCTTACATCACCAGCCCGTGTGTGGCCATCAAC
    GCATGCATGCGCATCACCATATCTCCTGGGTAGTGTCAGCCATGAACAGCAGTTTCGGTG
    TTGTTAGGTCTCgtctagtctccttcgcagccgtctttcgggacggggagcgttgcgtga

# 8 Getting genome feature track

``` python
# Open the input file and the output file
with open('../output/05.33-lncRNA-discovery/Apul_noncoding_transcripts_ids.txt', 'r') as infile, open('../output/05.33-lncRNA-discovery/Apul_lncRNA.bed', 'w') as outfile:
    # Process each line in the input file
    for line in infile:
        # Remove 'transcript::' and then split the line by ':' to extract the relevant parts
        parts = line.strip().replace('transcript::', '').split(':')
        chromosome = parts[0]
        # Split the position part by '-' to get start and end positions
        start, end = parts[1].split('-')
        
        # BED format requires the start position to be 0-based
        # Convert the start position to 0-based by subtracting 1
        start = str(int(start) - 1)
        
        # Write the chromosome, start, and end positions to the output file
        # Separate each field with a tab character
        outfile.write(f'{chromosome}\t{start}\t{end}\n')

# After running this script, 'output.bed' will contain the converted data in BED format.
```

``` bash
head ../output/05.33-lncRNA-discovery/Apul_lncRNA.bed 
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
