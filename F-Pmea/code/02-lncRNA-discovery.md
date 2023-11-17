02-Pmea-lncRNA-discovery
================
Steven Roberts
17 November, 2023

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

wget -r \
--no-directories --no-parent \
-P ../data/fastq/ \
-A "*fastq.gz" https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_meandrina/trimmed/
```

``` bash
ls ../data/fastq/
```

    RNA-POC-47-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    RNA-POC-47-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    RNA-POC-48-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    RNA-POC-48-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    RNA-POC-50-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    RNA-POC-50-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    RNA-POC-53-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    RNA-POC-53-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    RNA-POC-57-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    RNA-POC-57-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz

## 1.2 Genome

``` bash
cd ../data

curl -O https://owl.fish.washington.edu/halfshell/genomic-databank/Pocillopora_meandrina_HIv1.assembly.fasta
```

``` bash
cd ../data

curl -O https://owl.fish.washington.edu/halfshell/genomic-databank/Pocillopora_meandrina_HIv1.genes.gff3
```

``` bash
head -3 ../data/Pocillopora_meandrina_HIv1.genes.gff3

grep -v '^#' ../data/Pocillopora_meandrina_HIv1.genes.gff3 | cut -f3 | sort | uniq -c
```

## 1.3 HiSat

``` bash
/home/shared/hisat2-2.2.1/hisat2-build \
../data/Pocillopora_meandrina_HIv1.assembly.fasta \
../output/02-lncRNA-discovery/Pocillopora_meandrina_HIv1.assembly.index \
-p 40 \
../data/Pocillopora_meandrina_HIv1.genes.gff3 \
2> ../output/02-lncRNA-discovery/hisat2-build_stats.txt
```

``` bash
find ../data/fastq/*gz
```

    ../data/fastq/RNA-POC-47-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POC-47-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POC-48-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POC-48-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POC-50-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POC-50-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POC-53-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POC-53-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POC-57-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POC-57-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz

``` bash
find ../data/fastq/*R2_001.fastp-trim.20230519.fastq.gz \
| xargs basename -s -S1-TP2_R2_001.fastp-trim.20230519.fastq.gz | xargs -I{} \
echo {}
```

    RNA-POC-47
    RNA-POC-48
    RNA-POC-50
    RNA-POC-53
    RNA-POC-57

``` bash
find ../data/fastq/*R2_001.fastp-trim.20230519.fastq.gz \
| xargs basename -s -S1-TP2_R2_001.fastp-trim.20230519.fastq.gz | xargs -I{} \
/home/shared/hisat2-2.2.1/hisat2 \
-x ../output/02-lncRNA-discovery/Pocillopora_meandrina_HIv1.assembly.index \
-p 20 \
-1 ../data/fastq/{}-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz \
-2 ../data/fastq/{}-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz \
-S ../output/02-lncRNA-discovery/{}.sam \
2> ../output/02-lncRNA-discovery/hisat.out
```

``` bash
cat ../output/02-lncRNA-discovery/hisat.out
```

## 1.4 convert to bams

``` bash
for samfile in ../output/02-lncRNA-discovery/*.sam; do
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

![igv](http://gannet.fish.washington.edu/seashell/snaps/Monosnap_IGV_2023-11-05_09-57-16.png)

``` bash
ls ../output/02-lncRNA-discovery/*xml
```

# 2 StringTie

StringTie uses the sorted BAM files to assemble transcripts for each
sample, outputting them as GTF (Gene Transfer Format) files. And then
merges all individual GTF assemblies into a single merged GTF file. This
step extracts transcript information and merges GTFs from all samples–an
important step in creating a canonical list of lncRNAs across all
samples included in the pipeline.

``` bash
find ../output/02-lncRNA-discovery/*sorted.bam \
| xargs basename -s .sorted.bam | xargs -I{} \
/home/shared/stringtie-2.2.1.Linux_x86_64/stringtie \
-p 12 \
-G ../data/Pocillopora_meandrina_HIv1.genes.gff3 \
-o ../output/02-lncRNA-discovery/{}.gtf \
../output/02-lncRNA-discovery/{}.sorted.bam
```

``` bash
wc -l ../output/02-lncRNA-discovery/RNA*.gtf
head ../output/02-lncRNA-discovery/RNA*.gtf
```

       316036 ../output/02-lncRNA-discovery/RNA-POC-47.gtf
       328549 ../output/02-lncRNA-discovery/RNA-POC-48.gtf
       351599 ../output/02-lncRNA-discovery/RNA-POC-50.gtf
       358350 ../output/02-lncRNA-discovery/RNA-POC-53.gtf
       264594 ../output/02-lncRNA-discovery/RNA-POC-57.gtf
      1619128 total
    ==> ../output/02-lncRNA-discovery/RNA-POC-47.gtf <==
    # /home/shared/stringtie-2.2.1.Linux_x86_64/stringtie -p 12 -G ../data/Pocillopora_meandrina_HIv1.genes.gff3 -o ../output/02-lncRNA-discovery/RNA-POC-47.gtf ../output/02-lncRNA-discovery/RNA-POC-47.sorted.bam
    # StringTie version 2.2.1
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  10771   23652   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "46.000969"; FPKM "5.151896"; TPM "10.378988";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    10771   11117   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "2.982095";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    12784   12875   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "2"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "19.432919";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    13540   13643   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "3"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "23.592274";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    14320   14392   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "4"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "25.351660";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    14677   14752   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "5"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "30.846117";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    15026   15131   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "6"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "25.737658";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    15460   15575   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "7"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "36.226398";

    ==> ../output/02-lncRNA-discovery/RNA-POC-48.gtf <==
    # /home/shared/stringtie-2.2.1.Linux_x86_64/stringtie -p 12 -G ../data/Pocillopora_meandrina_HIv1.genes.gff3 -o ../output/02-lncRNA-discovery/RNA-POC-48.gtf ../output/02-lncRNA-discovery/RNA-POC-48.sorted.bam
    # StringTie version 2.2.1
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  168567  168771  1000    .   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "6.292683"; FPKM "0.779261"; TPM "1.496312";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    168567  168771  1000    .   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; cov "6.292683";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  171058  171339  1000    .   .   gene_id "STRG.2"; transcript_id "STRG.2.1"; cov "9.698582"; FPKM "1.201034"; TPM "2.306187";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    171058  171339  1000    .   .   gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "1"; cov "9.698582";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  173025  175585  1000    -   .   gene_id "STRG.3"; transcript_id "STRG.3.1"; cov "10.911830"; FPKM "1.351278"; TPM "2.594680";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    173025  175585  1000    -   .   gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "1"; cov "10.911830";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  174509  176444  1000    -   .   gene_id "STRG.3"; transcript_id "STRG.3.2"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20918.t1"; cov "3.123374"; FPKM "0.386786"; TPM "0.742694";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    174509  175333  1000    -   .   gene_id "STRG.3"; transcript_id "STRG.3.2"; exon_number "1"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20918.t1"; cov "3.302791";

    ==> ../output/02-lncRNA-discovery/RNA-POC-50.gtf <==
    # /home/shared/stringtie-2.2.1.Linux_x86_64/stringtie -p 12 -G ../data/Pocillopora_meandrina_HIv1.genes.gff3 -o ../output/02-lncRNA-discovery/RNA-POC-50.gtf ../output/02-lncRNA-discovery/RNA-POC-50.sorted.bam
    # StringTie version 2.2.1
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  10771   23652   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "37.975914"; FPKM "4.307361"; TPM "8.386994";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    10771   11117   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "17.452450";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    12784   12875   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "2"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "19.152174";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    13540   13643   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "3"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "25.913462";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    14320   14392   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "4"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "18.417809";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    14677   14752   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "5"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "19.087801";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    15026   15131   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "6"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "22.332779";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    15460   15575   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "7"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "26.038803";

    ==> ../output/02-lncRNA-discovery/RNA-POC-53.gtf <==
    # /home/shared/stringtie-2.2.1.Linux_x86_64/stringtie -p 12 -G ../data/Pocillopora_meandrina_HIv1.genes.gff3 -o ../output/02-lncRNA-discovery/RNA-POC-53.gtf ../output/02-lncRNA-discovery/RNA-POC-53.sorted.bam
    # StringTie version 2.2.1
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  223131  223665  1000    .   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "7.887850"; FPKM "0.951051"; TPM "1.908362";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    223131  223665  1000    .   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; cov "7.887850";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  10771   23652   1000    +   .   gene_id "STRG.2"; transcript_id "STRG.2.1"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "73.107300"; FPKM "8.814664"; TPM "17.687349";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    10771   11117   1000    +   .   gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "1"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "9.841498";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    12784   12875   1000    +   .   gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "2"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "26.038044";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    13540   13643   1000    +   .   gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "3"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "20.379808";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    14320   14392   1000    +   .   gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "4"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "22.808220";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    14677   14752   1000    +   .   gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "5"; reference_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; cov "32.276318";

    ==> ../output/02-lncRNA-discovery/RNA-POC-57.gtf <==
    # /home/shared/stringtie-2.2.1.Linux_x86_64/stringtie -p 12 -G ../data/Pocillopora_meandrina_HIv1.genes.gff3 -o ../output/02-lncRNA-discovery/RNA-POC-57.gtf ../output/02-lncRNA-discovery/RNA-POC-57.sorted.bam
    # StringTie version 2.2.1
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  10866   15555   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "9.660200"; FPKM "1.752686"; TPM "3.484353";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    10866   11117   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; cov "7.884921";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    12784   12875   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "2"; cov "11.532609";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    13540   13643   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "3"; cov "8.730769";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    14320   14392   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "4"; cov "19.746574";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    14677   14752   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "5"; cov "18.355263";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    15026   15131   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "6"; cov "6.924528";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    15460   15555   1000    +   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "7"; cov "2.000000";

``` bash
head ../data/Pocillopora_meandrina_HIv1.genes.gff3
```

Merges all individual GTF assemblies into a single merged GTF file.

This is used to create a non-redundant set of transcripts after running
StringTie separately on multiple RNA-Seq datasets.

``` bash
/home/shared/stringtie-2.2.1.Linux_x86_64/stringtie \
--merge \
-G ../data/Pocillopora_meandrina_HIv1.genes.gff3 \
-o ../output/02-lncRNA-discovery/stringtie_merged.gtf \
../output/02-lncRNA-discovery/*.gtf
```

``` bash

wc -l ../output/02-lncRNA-discovery/stringtie_merged.gtf
head ../output/02-lncRNA-discovery/stringtie_merged.gtf


echo "what is possible"

grep -v '^#' ../output/02-lncRNA-discovery/stringtie_merged.gtf | cut -f3 | sort | uniq -c
```

    559644 ../output/02-lncRNA-discovery/stringtie_merged.gtf
    # /home/shared/stringtie-2.2.1.Linux_x86_64/stringtie --merge -G ../data/Pocillopora_meandrina_HIv1.genes.gff3 -o ../output/02-lncRNA-discovery/stringtie_merged.gtf ../output/02-lncRNA-discovery/RNA-POC-47.gtf ../output/02-lncRNA-discovery/RNA-POC-48.gtf ../output/02-lncRNA-discovery/RNA-POC-50.gtf ../output/02-lncRNA-discovery/RNA-POC-53.gtf ../output/02-lncRNA-discovery/RNA-POC-57.gtf
    # StringTie version 2.2.1
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  10771   23776   1000    +   .   gene_id "MSTRG.1"; transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; 
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    10771   11117   1000    +   .   gene_id "MSTRG.1"; transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; exon_number "1"; 
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    12784   12875   1000    +   .   gene_id "MSTRG.1"; transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; exon_number "2"; 
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    13540   13643   1000    +   .   gene_id "MSTRG.1"; transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; exon_number "3"; 
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    14320   14392   1000    +   .   gene_id "MSTRG.1"; transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; exon_number "4"; 
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    14677   14752   1000    +   .   gene_id "MSTRG.1"; transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; exon_number "5"; 
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    15026   15131   1000    +   .   gene_id "MSTRG.1"; transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; exon_number "6"; 
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    15460   15575   1000    +   .   gene_id "MSTRG.1"; transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; exon_number "7"; 
    what is possible
     482050 exon
      77592 transcript

# 3 GFFcompare

<https://ccb.jhu.edu/software/stringtie/gffcompare.shtml>

``` bash
/home/shared/gffcompare-0.12.6.Linux_x86_64/gffcompare \
-r ../data/Pocillopora_meandrina_HIv1.genes.gff3 \
-o ../output/02-lncRNA-discovery/gffcompare_merged \
../output/02-lncRNA-discovery/stringtie_merged.gtf
```

``` bash
ls ../output/02-lncRNA-discovery/gffcompare_merged*
```

    ../output/02-lncRNA-discovery/gffcompare_merged.annotated.gtf
    ../output/02-lncRNA-discovery/gffcompare_merged.loci
    ../output/02-lncRNA-discovery/gffcompare_merged.stats
    ../output/02-lncRNA-discovery/gffcompare_merged.stringtie_merged.gtf.refmap
    ../output/02-lncRNA-discovery/gffcompare_merged.stringtie_merged.gtf.tmap
    ../output/02-lncRNA-discovery/gffcompare_merged.tracking

``` bash
cat ../output/02-lncRNA-discovery/gffcompare_merged.stats
```

    # gffcompare v0.12.6 | Command line was:
    #/home/shared/gffcompare-0.12.6.Linux_x86_64/gffcompare -r ../data/Pocillopora_meandrina_HIv1.genes.gff3 -o ../output/02-lncRNA-discovery/gffcompare_merged ../output/02-lncRNA-discovery/stringtie_merged.gtf
    #

    #= Summary for dataset: ../output/02-lncRNA-discovery/stringtie_merged.gtf 
    #     Query mRNAs :   77592 in   47812 loci  (57737 multi-exon transcripts)
    #            (11841 multi-transcript loci, ~1.6 transcripts per locus)
    # Reference mRNAs :   31840 in   31765 loci  (23907 multi-exon)
    # Super-loci w/ reference transcripts:    29451
    #-----------------| Sensitivity | Precision  |
            Base level:   100.0     |    60.5    |
            Exon level:    95.8     |    69.0    |
          Intron level:   100.0     |    80.1    |
    Intron chain level:   100.0     |    41.4    |
      Transcript level:    99.2     |    40.7    |
           Locus level:    99.2     |    61.3    |

         Matching intron chains:   23907
           Matching transcripts:   31588
                  Matching loci:   31525

              Missed exons:       0/208535  (  0.0%)
               Novel exons:   51995/292394  ( 17.8%)
            Missed introns:       5/176695  (  0.0%)
             Novel introns:   24810/220568  ( 11.2%)
               Missed loci:       0/31765   (  0.0%)
                Novel loci:   18361/47812   ( 38.4%)

     Total union super-loci across all input datasets: 47812 
    77592 out of 77592 consensus transcripts written in ../output/02-lncRNA-discovery/gffcompare_merged.annotated.gtf (0 discarded as redundant)

``` bash
head -10 ../output/02-lncRNA-discovery/gffcompare_merged.annotated.gtf
```

    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  10771   23776   .   +   .   transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; gene_id "MSTRG.1"; xloc "XLOC_000001"; cmp_ref "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; class_code "="; tss_id "TSS1";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    10771   11117   .   +   .   transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; gene_id "MSTRG.1"; exon_number "1";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    12784   12875   .   +   .   transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; gene_id "MSTRG.1"; exon_number "2";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    13540   13643   .   +   .   transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; gene_id "MSTRG.1"; exon_number "3";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    14320   14392   .   +   .   transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; gene_id "MSTRG.1"; exon_number "4";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    14677   14752   .   +   .   transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; gene_id "MSTRG.1"; exon_number "5";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    15026   15131   .   +   .   transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; gene_id "MSTRG.1"; exon_number "6";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    15460   15575   .   +   .   transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; gene_id "MSTRG.1"; exon_number "7";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    17213   17287   .   +   .   transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; gene_id "MSTRG.1"; exon_number "8";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   exon    17878   17995   .   +   .   transcript_id "Pocillopora_meandrina_HIv1___RNAseq.g20902.t1"; gene_id "MSTRG.1"; exon_number "9";

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
../output/02-lncRNA-discovery/gffcompare_merged.annotated.gtf | grep 'class_code "u"\|class_code "x"|\class_code "i"\|class_code "y"' | awk '($5 - $4 > 199) || ($4 - $5 > 199)' > ../output/02-lncRNA-discovery/Pmea_lncRNA_candidates.gtf
```

``` bash
wc -l ../output/02-lncRNA-discovery/Pmea_lncRNA_candidates.gtf
head ../output/02-lncRNA-discovery/Pmea_lncRNA_candidates.gtf
```

    14307 ../output/02-lncRNA-discovery/Pmea_lncRNA_candidates.gtf
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  164392  165433  .   +   .   transcript_id "MSTRG.18.1"; gene_id "MSTRG.18"; xloc "XLOC_000010"; class_code "u"; tss_id "TSS17";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  164763  165433  .   +   .   transcript_id "MSTRG.18.2"; gene_id "MSTRG.18"; xloc "XLOC_000010"; class_code "u"; tss_id "TSS18";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  245810  248612  .   +   .   transcript_id "MSTRG.29.1"; gene_id "MSTRG.29"; xloc "XLOC_000017"; class_code "u"; tss_id "TSS26";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  540345  541981  .   +   .   transcript_id "MSTRG.71.3"; gene_id "MSTRG.71"; xloc "XLOC_000039"; class_code "u"; tss_id "TSS54";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  551804  553744  .   +   .   transcript_id "MSTRG.76.1"; gene_id "MSTRG.76"; xloc "XLOC_000040"; class_code "u"; tss_id "TSS55";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  608054  611735  .   +   .   transcript_id "MSTRG.89.1"; gene_id "MSTRG.89"; xloc "XLOC_000045"; class_code "u"; tss_id "TSS69";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  608120  611735  .   +   .   transcript_id "MSTRG.89.2"; gene_id "MSTRG.89"; xloc "XLOC_000045"; class_code "u"; tss_id "TSS69";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  849629  852046  .   +   .   transcript_id "MSTRG.129.1"; gene_id "MSTRG.129"; xloc "XLOC_000063"; class_code "u"; tss_id "TSS91";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  850075  852046  .   +   .   transcript_id "MSTRG.129.2"; gene_id "MSTRG.129"; xloc "XLOC_000063"; class_code "u"; tss_id "TSS92";
    Pocillopora_meandrina_HIv1___Sc0000000  StringTie   transcript  958614  960727  .   +   .   transcript_id "MSTRG.142.1"; gene_id "MSTRG.142"; xloc "XLOC_000069"; class_code "u"; tss_id "TSS105";

# 5 Bedtools

Extracts the sequence data from the `$FASTA` reference file based on the
coordinates from the filtered GTF. The resulting sequences represent
potential lncRNA candidates.

``` bash
/home/shared/bedtools2/bin/fastaFromBed \
-fi ../data/Pocillopora_meandrina_HIv1.assembly.fasta \
-bed ../output/02-lncRNA-discovery/Pmea_lncRNA_candidates.gtf \
-fo ../output/02-lncRNA-discovery/Pmea_lncRNA_candidates.fasta \
-name -split
```

``` bash
fgrep -c ">" ../output/02-lncRNA-discovery/Pmea_lncRNA_candidates.fasta
head ../output/02-lncRNA-discovery/Pmea_lncRNA_candidates.fasta
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
-i ../output/02-lncRNA-discovery/Pmea_lncRNA_candidates.fasta \
-o ../output/02-lncRNA-discovery/Pmea_CPC2
```

``` bash
wc -l ../output/02-lncRNA-discovery/Pmea_CPC2.txt
head ../output/02-lncRNA-discovery/Pmea_CPC2.txt
```

    14308 ../output/02-lncRNA-discovery/Pmea_CPC2.txt
    #ID transcript_length   peptide_length  Fickett_score   pI  ORF_integrity   coding_probability  label
    transcript::Pocillopora_meandrina_HIv1___Sc0000000:164391-165433    1042    38  0.36996999999999997 7.97959041595459    1   0.010528    noncoding
    transcript::Pocillopora_meandrina_HIv1___Sc0000000:164762-165433    671 35  0.40464999999999995 10.916125297546387  1   0.0237837   noncoding
    transcript::Pocillopora_meandrina_HIv1___Sc0000000:245809-248612    2803    69  0.32791000000000003 10.203878593444824  1   0.0463183   noncoding
    transcript::Pocillopora_meandrina_HIv1___Sc0000000:540344-541981    1637    47  0.31204 5.33748836517334    1   0.0201418   noncoding
    transcript::Pocillopora_meandrina_HIv1___Sc0000000:551803-553744    1941    61  0.33813000000000004 5.948223304748535   1   0.034052    noncoding
    transcript::Pocillopora_meandrina_HIv1___Sc0000000:608053-611735    3682    63  0.25289 9.650481986999512   1   0.0372955   noncoding
    transcript::Pocillopora_meandrina_HIv1___Sc0000000:608119-611735    3616    63  0.25289 9.650481986999512   1   0.0372955   noncoding
    transcript::Pocillopora_meandrina_HIv1___Sc0000000:849628-852046    2418    58  0.25349000000000005 9.87547664642334    1   0.0347956   noncoding
    transcript::Pocillopora_meandrina_HIv1___Sc0000000:850074-852046    1972    58  0.27174000000000004 9.87547664642334    1   0.0296687   noncoding

\#Filter Filters the CPC2 results to get only noncoding transcripts
(using the class “noncoding” from the CPC2 results) and extracts their
IDs and matches these IDs with the sequences from the previous step to
generate a GTF of long noncoding transcripts.

Matches these IDs with the sequences from the previous step to generate
a GTF of noncoding transcripts.

``` bash
awk '$8 == "noncoding" {print $1}' ../output/02-lncRNA-discovery/Pmea_CPC2.txt > ../output/02-lncRNA-discovery/Pmea_noncoding_transcripts_ids.txt
```

``` bash
wc -l ../output/02-lncRNA-discovery/Pmea_noncoding_transcripts_ids.txt
head ../output/02-lncRNA-discovery/Pmea_noncoding_transcripts_ids.txt
```

``` bash
head ../output/02-lncRNA-discovery/Pmea_lncRNA_candidates.fasta
fgrep -c ">" ../output/02-lncRNA-discovery/Pmea_lncRNA_candidates.fasta
```

# 7 subsetting fasta

``` bash
/home/shared/samtools-1.12/samtools faidx ../output/02-lncRNA-discovery/Pmea_lncRNA_candidates.fasta \
-r ../output/02-lncRNA-discovery/Pmea_noncoding_transcripts_ids.txt > ../output/02-lncRNA-discovery/Pmea_lncRNA.fasta
```

ddd

``` bash
fgrep -c ">" ../output/02-lncRNA-discovery/Pmea_lncRNA.fasta
fgrep ">" ../output/02-lncRNA-discovery/Pmea_lncRNA.fasta | head -5

head ../output/02-lncRNA-discovery/Pmea_lncRNA.fasta
```

# 8 Getting genome feature track

``` python
# Open the input file and the output file
with open('../output/02-lncRNA-discovery/Pmea_noncoding_transcripts_ids.txt', 'r') as infile, open('../output/02-lncRNA-discovery/Pmea_lncRNA.bed', 'w') as outfile:
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
head ../output/02-lncRNA-discovery/Pmea_lncRNA.bed 
```
