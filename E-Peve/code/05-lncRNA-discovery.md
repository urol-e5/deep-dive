05-Peve-lncRNA-discovery
================
Steven Roberts
06 November, 2023

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
-A "*fastq.gz" https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_evermanni/trimmed/
```

``` bash
ls ../data/fastq/
```

    RNA-POR-71-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    RNA-POR-71-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    RNA-POR-73-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    RNA-POR-73-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    RNA-POR-76-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    RNA-POR-76-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    RNA-POR-79-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    RNA-POR-79-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    RNA-POR-82-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    RNA-POR-82-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz

## 1.2 Genome

``` bash
cd ../data

curl -O https://gannet.fish.washington.edu/seashell/snaps/Porites_evermanni_v1.fa
```

``` bash
cd ../data

curl -O https://gannet.fish.washington.edu/seashell/snaps/Porites_evermanni_v1.annot.gff
```

## 1.3 HiSat

``` bash
/home/shared/hisat2-2.2.1/hisat2_extract_exons.py \
../data/Porites_evermanni_v1.annot.gff \
> ../output/05-lncRNA-discovery/m_exon.tab
```

``` bash
head ../output/05-lncRNA-discovery/m_exon.tab
wc -l ../output/05-lncRNA-discovery/m_exon.tab
```

    0 ../output/05-lncRNA-discovery/m_exon.tab

``` bash

/home/shared/hisat2-2.2.1/hisat2_extract_splice_sites.py \
../data/Porites_evermanni_v1.annot.gff \
> ../output/05-lncRNA-discovery/m_splice_sites.tab
```

``` bash
wc -l ../output/05-lncRNA-discovery/m_splice_sites.tab
```

    0 ../output/05-lncRNA-discovery/m_splice_sites.tab

``` bash
/home/shared/hisat2-2.2.1/hisat2-build \
../data/Porites_evermanni_v1.fa \
../data/Porites_evermanni_v1.index \
--exon ../output/05-lncRNA-discovery/m_exon.tab \
--ss ../output/05-lncRNA-discovery/m_splice_sites.tab \
-p 40 \
../data/Porites_evermanni_v1.annot.gff \
2> ../output/05-lncRNA-discovery/hisat2-build_stats.txt
```

``` bash
find ../data/fastq/*gz
```

    ../data/fastq/RNA-POR-71-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POR-71-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POR-73-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POR-73-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POR-76-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POR-76-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POR-79-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POR-79-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POR-82-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    ../data/fastq/RNA-POR-82-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz

``` bash
find ../data/fastq/*R2_001.fastp-trim.20230519.fastq.gz \
| xargs basename -s -S1-TP2_R2_001.fastp-trim.20230519.fastq.gz | xargs -I{} \
echo {}
```

    RNA-POR-71
    RNA-POR-73
    RNA-POR-76
    RNA-POR-79
    RNA-POR-82

``` bash
find ../data/fastq/*R2_001.fastp-trim.20230519.fastq.gz \
| xargs basename -s -S1-TP2_R2_001.fastp-trim.20230519.fastq.gz | xargs -I{} \
/home/shared/hisat2-2.2.1/hisat2 \
-x ../data/Porites_evermanni_v1.index \
-p 20 \
-1 ../data/fastq/{}-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz \
-2 ../data/fastq/{}-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz \
-S ../output/05-lncRNA-discovery/{}.sam \
2> ../output/05-lncRNA-discovery/hisat.out
```

``` bash
cat ../output/05-lncRNA-discovery/hisat.out
```

    50831351 reads; of these:
      50831351 (100.00%) were paired; of these:
        11712763 (23.04%) aligned concordantly 0 times
        30209321 (59.43%) aligned concordantly exactly 1 time
        8909267 (17.53%) aligned concordantly >1 times
        ----
        11712763 pairs aligned concordantly 0 times; of these:
          2145413 (18.32%) aligned discordantly 1 time
        ----
        9567350 pairs aligned 0 times concordantly or discordantly; of these:
          19134700 mates make up the pairs; of these:
            13827000 (72.26%) aligned 0 times
            3449026 (18.02%) aligned exactly 1 time
            1858674 (9.71%) aligned >1 times
    86.40% overall alignment rate
    51385213 reads; of these:
      51385213 (100.00%) were paired; of these:
        11914070 (23.19%) aligned concordantly 0 times
        28680865 (55.82%) aligned concordantly exactly 1 time
        10790278 (21.00%) aligned concordantly >1 times
        ----
        11914070 pairs aligned concordantly 0 times; of these:
          2306667 (19.36%) aligned discordantly 1 time
        ----
        9607403 pairs aligned 0 times concordantly or discordantly; of these:
          19214806 mates make up the pairs; of these:
            13255152 (68.98%) aligned 0 times
            3430116 (17.85%) aligned exactly 1 time
            2529538 (13.16%) aligned >1 times
    87.10% overall alignment rate
    49828147 reads; of these:
      49828147 (100.00%) were paired; of these:
        11096891 (22.27%) aligned concordantly 0 times
        29566240 (59.34%) aligned concordantly exactly 1 time
        9165016 (18.39%) aligned concordantly >1 times
        ----
        11096891 pairs aligned concordantly 0 times; of these:
          1929423 (17.39%) aligned discordantly 1 time
        ----
        9167468 pairs aligned 0 times concordantly or discordantly; of these:
          18334936 mates make up the pairs; of these:
            13987351 (76.29%) aligned 0 times
            2687611 (14.66%) aligned exactly 1 time
            1659974 (9.05%) aligned >1 times
    85.96% overall alignment rate
    49976568 reads; of these:
      49976568 (100.00%) were paired; of these:
        10514435 (21.04%) aligned concordantly 0 times
        27932980 (55.89%) aligned concordantly exactly 1 time
        11529153 (23.07%) aligned concordantly >1 times
        ----
        10514435 pairs aligned concordantly 0 times; of these:
          2532283 (24.08%) aligned discordantly 1 time
        ----
        7982152 pairs aligned 0 times concordantly or discordantly; of these:
          15964304 mates make up the pairs; of these:
            10928968 (68.46%) aligned 0 times
            2662474 (16.68%) aligned exactly 1 time
            2372862 (14.86%) aligned >1 times
    89.07% overall alignment rate
    48908730 reads; of these:
      48908730 (100.00%) were paired; of these:
        12403366 (25.36%) aligned concordantly 0 times
        28498492 (58.27%) aligned concordantly exactly 1 time
        8006872 (16.37%) aligned concordantly >1 times
        ----
        12403366 pairs aligned concordantly 0 times; of these:
          2215956 (17.87%) aligned discordantly 1 time
        ----
        10187410 pairs aligned 0 times concordantly or discordantly; of these:
          20374820 mates make up the pairs; of these:
            15478680 (75.97%) aligned 0 times
            3328833 (16.34%) aligned exactly 1 time
            1567307 (7.69%) aligned >1 times
    84.18% overall alignment rate

## 1.4 convert to bams

``` bash
for samfile in ../output/05-lncRNA-discovery/*.sam; do
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
ls ../output/05.32-lncRNA-discovery/*xml
```

# 2 StringTie

StringTie uses the sorted BAM files to assemble transcripts for each
sample, outputting them as GTF (Gene Transfer Format) files. And then
merges all individual GTF assemblies into a single merged GTF file. This
step extracts transcript information and merges GTFs from all samples–an
important step in creating a canonical list of lncRNAs across all
samples included in the pipeline.

``` bash
find ../output/05-lncRNA-discovery/*sorted.bam \
| xargs basename -s .sorted.bam | xargs -I{} \
/home/shared/stringtie-2.2.1.Linux_x86_64/stringtie \
-p 12 \
-G ../data/Porites_evermanni_v1.annot.gff \
-o ../output/05-lncRNA-discovery/{}.gtf \
../output/05-lncRNA-discovery/{}.sorted.bam
```

``` bash
wc -l ../output/05-lncRNA-discovery/RNA*.gtf
ls ../output/05-lncRNA-discovery/RNA*.gtf
```

       371509 ../output/05-lncRNA-discovery/RNA-POR-71.gtf
       269300 ../output/05-lncRNA-discovery/RNA-POR-73.gtf
       308422 ../output/05-lncRNA-discovery/RNA-POR-76.gtf
       297750 ../output/05-lncRNA-discovery/RNA-POR-79.gtf
       387644 ../output/05-lncRNA-discovery/RNA-POR-82.gtf
      1634625 total
    ../output/05-lncRNA-discovery/RNA-POR-71.gtf
    ../output/05-lncRNA-discovery/RNA-POR-73.gtf
    ../output/05-lncRNA-discovery/RNA-POR-76.gtf
    ../output/05-lncRNA-discovery/RNA-POR-79.gtf
    ../output/05-lncRNA-discovery/RNA-POR-82.gtf

![gtf](http://gannet.fish.washington.edu/seashell/snaps/Monosnap_IGV_-_Session_Userssr320DesktopApul_lncRNA_igv_session.xml_2023-11-05_10-11-31.png)

Merges all individual GTF assemblies into a single merged GTF file.

This is used to create a non-redundant set of transcripts after running
StringTie separately on multiple RNA-Seq datasets.

``` bash
/home/shared/stringtie-2.2.1.Linux_x86_64/stringtie \
--merge \
-G ../data/Porites_evermanni_v1.annot.gff \
-o ../output/05-lncRNA-discovery/stringtie_merged.gtf \
../output/05-lncRNA-discovery/*.gtf
```

``` bash

wc -l ../output/05-lncRNA-discovery/stringtie_merged.gtf
head ../output/05-lncRNA-discovery/stringtie_merged.gtf


echo "what is possible"

grep -v '^#' ../output/05-lncRNA-discovery/stringtie_merged.gtf | cut -f3 | sort | uniq
```

    542182 ../output/05-lncRNA-discovery/stringtie_merged.gtf
    # /home/shared/stringtie-2.2.1.Linux_x86_64/stringtie --merge -G ../data/Porites_evermanni_v1.annot.gff -o ../output/05-lncRNA-discovery/stringtie_merged.gtf ../output/05-lncRNA-discovery/RNA-POR-71.gtf ../output/05-lncRNA-discovery/RNA-POR-73.gtf ../output/05-lncRNA-discovery/RNA-POR-76.gtf ../output/05-lncRNA-discovery/RNA-POR-79.gtf ../output/05-lncRNA-discovery/RNA-POR-82.gtf
    # StringTie version 2.2.1
    Porites_evermani_scaffold_1 StringTie   transcript  2566    6331    1000    -   .   gene_id "MSTRG.1"; transcript_id "MSTRG.1.1"; 
    Porites_evermani_scaffold_1 StringTie   exon    2566    3444    1000    -   .   gene_id "MSTRG.1"; transcript_id "MSTRG.1.1"; exon_number "1"; 
    Porites_evermani_scaffold_1 StringTie   exon    4284    4478    1000    -   .   gene_id "MSTRG.1"; transcript_id "MSTRG.1.1"; exon_number "2"; 
    Porites_evermani_scaffold_1 StringTie   exon    4726    4795    1000    -   .   gene_id "MSTRG.1"; transcript_id "MSTRG.1.1"; exon_number "3"; 
    Porites_evermani_scaffold_1 StringTie   exon    5415    5551    1000    -   .   gene_id "MSTRG.1"; transcript_id "MSTRG.1.1"; exon_number "4"; 
    Porites_evermani_scaffold_1 StringTie   exon    6274    6331    1000    -   .   gene_id "MSTRG.1"; transcript_id "MSTRG.1.1"; exon_number "5"; 
    Porites_evermani_scaffold_1 StringTie   transcript  2566    6331    1000    -   .   gene_id "MSTRG.1"; transcript_id "MSTRG.1.2"; 
    Porites_evermani_scaffold_1 StringTie   exon    2566    3444    1000    -   .   gene_id "MSTRG.1"; transcript_id "MSTRG.1.2"; exon_number "1"; 
    what is possible
    exon
    transcript

# 3 GFFcompare

<https://ccb.jhu.edu/software/stringtie/gffcompare.shtml>

``` bash
/home/shared/gffcompare-0.12.6.Linux_x86_64/gffcompare \
-r ../data/Porites_evermanni_v1.annot.gff \
-o ../output/05-lncRNA-discovery/gffcompare_merged \
../output/05-lncRNA-discovery/stringtie_merged.gtf
```

``` bash
ls ../output/05-lncRNA-discovery/gffcompare_merged*
```

    ../output/05-lncRNA-discovery/gffcompare_merged.annotated.gtf
    ../output/05-lncRNA-discovery/gffcompare_merged.loci
    ../output/05-lncRNA-discovery/gffcompare_merged.stats
    ../output/05-lncRNA-discovery/gffcompare_merged.stringtie_merged.gtf.refmap
    ../output/05-lncRNA-discovery/gffcompare_merged.stringtie_merged.gtf.tmap
    ../output/05-lncRNA-discovery/gffcompare_merged.tracking

``` bash
cat ../output/05-lncRNA-discovery/gffcompare_merged.stats
```

    # gffcompare v0.12.6 | Command line was:
    #/home/shared/gffcompare-0.12.6.Linux_x86_64/gffcompare -r ../data/Porites_evermanni_v1.annot.gff -o ../output/05-lncRNA-discovery/gffcompare_merged ../output/05-lncRNA-discovery/stringtie_merged.gtf
    #

    #= Summary for dataset: ../output/05-lncRNA-discovery/stringtie_merged.gtf 
    #     Query mRNAs :   74686 in   50813 loci  (59310 multi-exon transcripts)
    #            (11697 multi-transcript loci, ~1.5 transcripts per locus)
    # Reference mRNAs :   40389 in   40240 loci  (32889 multi-exon)
    # Super-loci w/ reference transcripts:    38915
    #-----------------| Sensitivity | Precision  |
            Base level:   100.0     |    74.6    |
            Exon level:    98.0     |    76.4    |
          Intron level:   100.0     |    82.0    |
    Intron chain level:   100.0     |    55.5    |
      Transcript level:    99.7     |    53.9    |
           Locus level:    99.7     |    76.4    |

         Matching intron chains:   32889
           Matching transcripts:   40267
                  Matching loci:   40119

              Missed exons:       0/235837  (  0.0%)
               Novel exons:   45306/305701  ( 14.8%)
            Missed introns:      81/195463  (  0.0%)
             Novel introns:   22911/238378  (  9.6%)
               Missed loci:       0/40240   (  0.0%)
                Novel loci:   11898/50813   ( 23.4%)

     Total union super-loci across all input datasets: 50813 
    74686 out of 74686 consensus transcripts written in ../output/05-lncRNA-discovery/gffcompare_merged.annotated.gtf (0 discarded as redundant)

``` bash
head -10 ../output/05-lncRNA-discovery/gffcompare_merged.annotated.gtf
```

    Porites_evermani_scaffold_1 StringTie   transcript  7720    11436   .   +   .   transcript_id "MSTRG.2.1"; gene_id "MSTRG.2"; xloc "XLOC_000001"; cmp_ref "Peve_00000039"; class_code "j"; tss_id "TSS1";
    Porites_evermani_scaffold_1 StringTie   exon    7720    7745    .   +   .   transcript_id "MSTRG.2.1"; gene_id "MSTRG.2"; exon_number "1";
    Porites_evermani_scaffold_1 StringTie   exon    8484    8620    .   +   .   transcript_id "MSTRG.2.1"; gene_id "MSTRG.2"; exon_number "2";
    Porites_evermani_scaffold_1 StringTie   exon    9236    9305    .   +   .   transcript_id "MSTRG.2.1"; gene_id "MSTRG.2"; exon_number "3";
    Porites_evermani_scaffold_1 StringTie   exon    9554    9748    .   +   .   transcript_id "MSTRG.2.1"; gene_id "MSTRG.2"; exon_number "4";
    Porites_evermani_scaffold_1 StringTie   exon    10689   11436   .   +   .   transcript_id "MSTRG.2.1"; gene_id "MSTRG.2"; exon_number "5";
    Porites_evermani_scaffold_1 StringTie   transcript  8492    11436   .   +   .   transcript_id "Peve_00000039"; gene_id "MSTRG.2"; xloc "XLOC_000001"; cmp_ref "Peve_00000039"; class_code "="; tss_id "TSS2";
    Porites_evermani_scaffold_1 StringTie   exon    8492    8620    .   +   .   transcript_id "Peve_00000039"; gene_id "MSTRG.2"; exon_number "1";
    Porites_evermani_scaffold_1 StringTie   exon    9236    9320    .   +   .   transcript_id "Peve_00000039"; gene_id "MSTRG.2"; exon_number "2";
    Porites_evermani_scaffold_1 StringTie   exon    9554    9748    .   +   .   transcript_id "Peve_00000039"; gene_id "MSTRG.2"; exon_number "3";

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
../output/05-lncRNA-discovery/gffcompare_merged.annotated.gtf | grep 'class_code "u"\|class_code "x"|\class_code "i"\|class_code "y"' | awk '($5 - $4 > 199) || ($4 - $5 > 199)' > ../output/05-lncRNA-discovery/Peve_lncRNA_candidates.gtf
```

``` bash
wc -l ../output/05-lncRNA-discovery/Peve_lncRNA_candidates.gtf
head ../output/05-lncRNA-discovery/Peve_lncRNA_candidates.gtf
```

    8105 ../output/05-lncRNA-discovery/Peve_lncRNA_candidates.gtf
    Porites_evermani_scaffold_1 StringTie   transcript  422644  423512  .   -   .   transcript_id "MSTRG.27.1"; gene_id "MSTRG.27"; xloc "XLOC_000078"; class_code "u"; tss_id "TSS109";
    Porites_evermani_scaffold_1 StringTie   transcript  1084868 1089422 .   -   .   transcript_id "MSTRG.67.1"; gene_id "MSTRG.67"; xloc "XLOC_000104"; class_code "u"; tss_id "TSS140";
    Porites_evermani_scaffold_1 StringTie   transcript  372246  372449  .   .   .   transcript_id "MSTRG.24.1"; gene_id "MSTRG.24"; xloc "XLOC_000124"; class_code "u"; tss_id "TSS166";
    Porites_evermani_scaffold_1 StringTie   transcript  683879  684280  .   .   .   transcript_id "MSTRG.46.1"; gene_id "MSTRG.46"; xloc "XLOC_000126"; class_code "u"; tss_id "TSS168";
    Porites_evermani_scaffold_1 StringTie   transcript  1202045 1202328 .   .   .   transcript_id "MSTRG.70.1"; gene_id "MSTRG.70"; xloc "XLOC_000132"; class_code "u"; tss_id "TSS174";
    Porites_evermani_scaffold_1 StringTie   transcript  1477782 1478077 .   .   .   transcript_id "MSTRG.81.1"; gene_id "MSTRG.81"; xloc "XLOC_000133"; class_code "u"; tss_id "TSS175";
    Porites_evermani_scaffold_10    StringTie   transcript  260428  263111  .   +   .   transcript_id "MSTRG.120.1"; gene_id "MSTRG.120"; xloc "XLOC_000142"; class_code "u"; tss_id "TSS186";
    Porites_evermani_scaffold_10    StringTie   transcript  260441  263111  .   +   .   transcript_id "MSTRG.120.2"; gene_id "MSTRG.120"; xloc "XLOC_000142"; class_code "u"; tss_id "TSS186";
    Porites_evermani_scaffold_10    StringTie   transcript  890034  892111  .   +   .   transcript_id "MSTRG.154.1"; gene_id "MSTRG.154"; xloc "XLOC_000164"; class_code "u"; tss_id "TSS212";
    Porites_evermani_scaffold_10    StringTie   transcript  1071613 1085055 .   +   .   transcript_id "MSTRG.166.1"; gene_id "MSTRG.166"; xloc "XLOC_000170"; class_code "u"; tss_id "TSS218";

# 5 Bedtools

Extracts the sequence data from the `$FASTA` reference file based on the
coordinates from the filtered GTF. The resulting sequences represent
potential lncRNA candidates.

``` bash
/home/shared/bedtools2/bin/fastaFromBed \
-fi ../data/Porites_evermanni_v1.fa \
-bed ../output/05-lncRNA-discovery/Peve_lncRNA_candidates.gtf \
-fo ../output/05-lncRNA-discovery/Peve_lncRNA_candidates.fasta \
-name -split
```

``` bash
fgrep -c ">" ../output/05-lncRNA-discovery/Peve_lncRNA_candidates.fasta
head ../output/05-lncRNA-discovery/Peve_lncRNA_candidates.fasta
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
-i ../output/05-lncRNA-discovery/Peve_lncRNA_candidates.fasta \
-o ../output/05-lncRNA-discovery/Peve_CPC2
```

``` bash
wc -l ../output/05-lncRNA-discovery/Peve_CPC2.txt
head ../output/05-lncRNA-discovery/Peve_CPC2.txt
```

    8106 ../output/05-lncRNA-discovery/Peve_CPC2.txt
    #ID transcript_length   peptide_length  Fickett_score   pI  ORF_integrity   coding_probability  label
    transcript::Porites_evermani_scaffold_1:422643-423512   869 30  0.36561 4.0500284194946286  1   0.0162382   noncoding
    transcript::Porites_evermani_scaffold_1:1084867-1089422 4555    111 0.28259 6.533494377136231   1   0.214842    noncoding
    transcript::Porites_evermani_scaffold_1:372245-372449   204 15  0.26338999999999996 4.0500284194946286  -1  0.0543788   noncoding
    transcript::Porites_evermani_scaffold_1:683878-684280   402 47  0.40454999999999997 9.512841987609864   1   0.0312414   noncoding
    transcript::Porites_evermani_scaffold_1:1202044-1202328 284 68  0.37632 11.633916282653807  1   0.0619633   noncoding
    transcript::Porites_evermani_scaffold_1:1477781-1478077 296 28  0.42535 11.999967765808105  -1  0.158961    noncoding
    transcript::Porites_evermani_scaffold_10:260427-263111  2684    82  0.33597000000000005 5.269224739074706   1   0.131892    noncoding
    transcript::Porites_evermani_scaffold_10:260440-263111  2671    82  0.33597000000000005 5.269224739074706   1   0.131892    noncoding
    transcript::Porites_evermani_scaffold_10:890033-892111  2078    62  0.31204 6.884190940856934   1   0.0264183   noncoding

\#Filter Filters the CPC2 results to get only noncoding transcripts
(using the class “noncoding” from the CPC2 results) and extracts their
IDs and matches these IDs with the sequences from the previous step to
generate a GTF of long noncoding transcripts.

Matches these IDs with the sequences from the previous step to generate
a GTF of noncoding transcripts.

``` bash
awk '$8 == "noncoding" {print $1}' ../output/05-lncRNA-discovery/Peve_CPC2.txt > ../output/05-lncRNA-discovery/Peve_noncoding_transcripts_ids.txt
```

``` bash
wc -l ../output/05-lncRNA-discovery/Peve_noncoding_transcripts_ids.txt
head ../output/05-lncRNA-discovery/Peve_noncoding_transcripts_ids.txt
```

``` bash
head ../output/05-lncRNA-discovery/Apul_lncRNA_candidates.fasta
fgrep -c ">" ../output/05-lncRNA-discovery/Apul_lncRNA_candidates.fasta
```

# 7 subsetting fasta

``` bash
/home/shared/samtools-1.12/samtools faidx ../output/05.32-lncRNA-discovery/Apul_lncRNA_candidates.fasta \
-r ../output/05.32-lncRNA-discovery/Apul_noncoding_transcripts_ids.txt > ../output/05.32-lncRNA-discovery/Apul_lncRNA.fasta
```

``` bash
fgrep -c ">" ../output/05.32-lncRNA-discovery/Apul_lncRNA.fasta
fgrep ">" ../output/05.32-lncRNA-discovery/Apul_lncRNA.fasta | head -5

head ../output/05.32-lncRNA-discovery/Apul_lncRNA.fasta
```

# 8 Getting genome feature track

``` python
# Open the input file and the output file
with open('../output/05.32-lncRNA-discovery/Apul_noncoding_transcripts_ids.txt', 'r') as infile, open('../output/05.32-lncRNA-discovery/Apul_lncRNA.bed', 'w') as outfile:
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
head ../output/05.32-lncRNA-discovery/Apul_lncRNA.bed 
```
