02 Peve lncRNA align
================
Steven Roberts
07 September, 2023

- <a href="#001-analysis" id="toc-001-analysis">0.0.1 Analysis:</a>
- <a href="#002-suggestions" id="toc-002-suggestions">0.0.2
  Suggestions:</a>
- <a href="#1-blast-to-srna-genome-proteome"
  id="toc-1-blast-to-srna-genome-proteome">1 Blast to sRNA, genome,
  proteome….</a>
  - <a href="#11-lncrna-fasta" id="toc-11-lncrna-fasta">1.1 lncRNA fasta</a>
- <a href="#2-make-some-blastdb" id="toc-2-make-some-blastdb">2 make some
  blastdb</a>
  - <a href="#21-srna" id="toc-21-srna">2.1 sRNA</a>
  - <a href="#22-genome" id="toc-22-genome">2.2 Genome</a>
  - <a href="#23-proteins" id="toc-23-proteins">2.3 Proteins</a>
  - <a href="#24-genes" id="toc-24-genes">2.4 Genes</a>
- <a href="#3-blast-comparison" id="toc-3-blast-comparison">3 Blast
  comparison</a>
- <a href="#4-result---database-srna" id="toc-4-result---database-srna">4
  Result - database: sRNA</a>
- <a href="#5-result--database-genome"
  id="toc-5-result--database-genome">5 Result -database: genome</a>
- <a href="#6-result---database-proteome"
  id="toc-6-result---database-proteome">6 Result - database: proteome</a>

\#TLDR

The code and its output suggest that you are running BLAST comparisons
on a set of long non-coding RNAs (lncRNAs) against different
databases—small RNAs (sRNAs), genome, and proteome—and observing
different patterns of hits.

1.  **sRNA Database**:
    - In the first run, you used a stringent e-value (1E-40) and got 0
      hits.
    - In the second run without e-value restriction, you got 1,798,786
      hits.
2.  **Genome Database**:
    - Without an e-value restriction, you got 14,547,262 hits.
    - With an e-value restriction of 1E-40, the hits reduced to 26,669.
3.  **Proteome Database**:
    - You got 351,260 hits when querying against the proteome database.

### 0.0.1 Analysis:

- **sRNA Database**: The 0 hits in the first run suggest that no lncRNAs
  had significant similarity to sRNAs under the stringent conditions.
  However, when the stringency was relaxed, a high number of hits were
  observed. This suggests that there is some level of sequence
  similarity between lncRNAs and sRNAs, but it may not be functionally
  relevant (at least under stringent e-value settings).

- **Genome Database**: As expected, you see a vast number of hits when
  the e-value is not restricted. It narrows down to a more manageable
  number when stringency is increased, but still, a large number of hits
  are observed. This could be due to the fact that the lncRNAs are part
  of the genome, and thus, self-matches and paralogous sequences might
  be increasing the hits.

- **Proteome Database**: The number of hits suggests that some of these
  lncRNAs might have regions that could be translated into protein
  sequences or resemble known proteins, although lncRNAs are generally
  not translated.

### 0.0.2 Suggestions:

1.  **Analyze Overlap**: You could analyze the overlapping hits between
    these databases to see if any lncRNAs show hits in all three
    databases.

2.  **Functional Annotation**: Use other bioinformatic tools to predict
    the functional roles of the lncRNAs that show significant hits.

3.  **Alignment Visualization**: You might want to visualize some of
    these alignments to better understand the areas of similarity.

4.  **Statistical Significance**: You may also apply statistical tests
    to see if the number of hits in any of these databases is
    significantly higher than expected by chance.

5.  **Investigate Parameters**: Review BLAST parameters like
    `-max_target_seqs` and `-max_hsps` based on your research questions,
    as these can impact your results.

Remember, the number of hits alone doesn’t convey functional or
biological significance; further investigation is needed.

# 1 Blast to sRNA, genome, proteome….

## 1.1 lncRNA fasta

``` bash
head ../../DEF-cross-species/data/peve_bedtools_lncRNAs.fasta
```

    ## >::Porites_evermani_scaffold_1:422643-423512
    ## GGCAAAGCCACAATCCATGATAAATGAGGGCATAAGCCCGAGGAAGAAGAACTCATAGATCTTGTCCTAATCCCTTTTGGGAGAGCACATTTTTTTCTTTATGCACTCACTGTGGATAAACACTCAATGGATTTTAGAGAAAAGGTGGACTGTAAGCAGTATAATACCTCCTGGAATTTGCCAGTAACTGAAAGAATCTAATCTGAAAAGTCTCTAGGTGTATACTGGGCAACATAGTAATAATTAGTGCATTTTAGAAGATAGGAATGAACGTGGCACATGTATCAAAATTGCATCTTACAAAATAGTTATATTATCAAGAAATCCCTTGCTAAATAATAATTGTAAATCAAGACATGTTTCTAGAACTGGAAACTCCTGGAAATCTGTATGCATCATGTCAATCACTGTACCAATGTTCCCTCATAGAGGGTTTTATATACTAGAAACACTTTGTGAGACTTAAATGTGTTATGCAATTCCAGGTTTGACTACTGTTCCAGGGTCTCAGTCAGTTAAGATGCAGGCATAAGTGAGAGACAGCCCTCCATTCCTTCCTCTCCCTTTTTTTTGGCACTTGGTTTTCTATTTGCTGAATGCCAGTAACTGAGAGCCTAGAACAACCTAGGTTTATGAAAAGTTTCAAGTAATGCTGATACAACTGTGGGAGTTTTGAACCCAGGAGTCACTTCAAAAGTAGGTTGAGTTTGATCGTCCAGGTGAACGTAGTCCTGAATAGGACTGTTGTTGTTGACAGTGACTGACGTTTCGACAACCTGTGCGGTAGTCATTTTCCTATTCAGGATTCATTCACCCGAACATACTTTTGGAATGCTGATACAGTTAACTATCCTTACACACAAACAATG
    ## >::Porites_evermani_scaffold_1:1084867-1089422
    ## AGGGATTAGCTGTTTTTAATTTTTGTGACTTCCAGAGAGTACCCAGTACCCAGCATTGATAATATTTTCGTTTTTATTGAGTATGTGCAATAGTCTAAGATAGAAATACATACTTTCAAACAATACGGTGTGAATACAGTGTGCATGCCCTATGCAAACCAGTATTTCAGTGTATACCATTCTGTTTCTGAACTAGGGGAATGAAACAAGTTGTAACAATTTTGAACAGACACACTTGTGTAGCGAATTTAAGTTAGAGAATATTTACTTCGGAATAAACAATTAAAAAAATAGCAAACATTAGAACCCACAAAAATTTCAGCAAAATAAAATAATAACACAGCAATTAACAATTTTTCATCCTCCAGTCAACTTGTAATTGGCTAATTAAAATGAACAAAACTACCTTTATGAGCATCATGGACACAAATATAGTACTACATGTAACAACAGTATTCTTTTCAATCTCACTGAATAGTGGCCACAAGAATAATAATTATTTTGCTCTATATTCAACGCCTAGCACACTTATTTCAATGCTGTAGTTAGACTGCATATTGCACATCTGATTTGAGAAAATCAACTGGCCCTCATACGGCATCTTTCCTGCTTTTTTCTTTTTCTAGGATATTCACAAGAATGATCCCCATTCCCTTTCAAAAAAGGGAAGGGAAGCAGAGTATAAGCGAAAGACTAGAGGGGAAACAGAGGTTTGAGGCCTTCCTCAGAGCAGCACCGAAAAGTGTGGGTGTTTTTATAAAACTTCATGAAGCATATGAGGAAATGACTAAAAGCCCACTCTTTAAGAACATTGAACAGGAGATCCAACCACTTGATCTTGACGTTCGGTCACAATCAAATTTAGCAGTTCTCCTTTATTCCACTGGCAAGGTACTAATAAAGAACAAAGTAGATAAAAATAATAATTAATACTGTGAAATTCTAAAATGCCCTGGGCCGAGTTGTTCAAAGCTGGGTTAAGATAACCCAGGGTTATTTGAAGAGATTTGAATTCAGATTTGAAAGCTTAAAAAGCATTTCGGTTTAAGTTCTTTGTGTTGACAAGTTGATGATTGGAAGCTCTAAAAATAACAGAGAAAATTTCCGAGAAAATACTTTAGAACACAAGAACAAGAAACCCGGGTTAAATTTAACCCCAGGTTAAGCGCTAACCAGCCTTTGAAGAACTGGGCCCTGGTTACTTTCTACTTAAATAAGCAAAAGGGGCCACTACTTTAAGGTTGTTGTTGTCTTGGGGTAATCATTCCTTTTGGATGGTTAAAAAAAACTTAGTAGTTTTTTTTGAGGGTGGCATTTATTGGTACTTATGCTTCCATGTGTGGCATTTGTTCGAGGGTGGCATTTAATCGAATAATAATAGTTCACAGTTTTAGTCGAGACAGTTTCGCACATATGGTACCATAGTTCACAGATTATGACACATTGAGGCTTTTATTAAATTGATTTGTTACTATGAATGGTGCTAGCTGTGTTTTGGAAAAAAACTTGTTATTATAGTGTAGTTAGAATAAACTTCATCTGGTGTGACTTTAATCCTTTCCTGAACTGAAATTATTTCATCGGAAAGTGTATGAATTTATCACTTTATTATTATCCATGCTGTAAGGATGTGAAACTAAGTCTAAGGAAAAATAATGTATGTCACACCTTCTTTAATTTTGTGTACAATGTTGTTACAGATCTTGGATTATATTGGAATGAAAAACCGGTACTGGGCAGACATACTTCTATTTATAAAAGAAGTTATGGTCCCGAAGTACCAAGCTCTGCTGGTTGGTAACATTGAGAAAATGCAGAACATGTCCTCCAATCAGCAATTCAGTGGCATGGAAGATTCTTGCCTTCAGGCAATTGAACAGATTGTTCGAACTGCCATGGAATACCCAAATAAGGTAAGCGTTTGTATTTTTATTAATCACTATTACATGTATTATTTTAATAAGACATTTACCGGATGGTACATGTAGGTAGTGGACCATGAAGCCCCGGGCGGGGGGAAGGAGTTGGCTAAGGTACATGTAGCCTTCCCCAATGTATGCTCAGTTGTATCTTTAAAGCCTTATTCATGAATAAAAGCACACTAATATTTTAAGATACCAACAAAAATCCAAAGCCCAAGGGTTTAAATTGGGCATCAATTTAGATAAAGATTTCAAGAGATGAAGAGCTAAACTGGTCACAATAGACTGTAATATCCTTAATCTCCTTTGTCTTCCTTCACAGGCCAAAGACTACCTCAAACCTAAAAGTGCAACAACAAGTTGCCTTTGTGGTGAAGTGGTATCCTTCAATGTTCCCCAAGCATTCCAACATGCGGGTATTACTGTTTCTTCTACAATCTTCAGACAGTATATAAAAAAAAAATGAGTTGGGCTGTACTGGAAAAAACATGCACTTACCTAGTGGGAAAAATGTAACTTGCGTCCACCTTAGCCGAGGAAAACTGCCACAAGAAATTCTGGAAATCCTAAACAAAGGTAATTTGATTTTTATTAAGTATTGACGCATTTAACTTAAAGATAAAACTTTGGCATTCTTCTGGCATACCAAAATTTACATCAAACCATTAAACCTTCAATGCATTTACAGCCACAAGCAGCTCTACTACTGTTACCACTTCCTCAGAATCTCTTGATGTTGCACCAGTAATGGAAGCACAAAAGTCAGGTATATAGTGCTATTTTAATTAATTAAACAAAGAAAACAGTGTGGAGACAGAAACGATCCTTCCTTTTCGTATATGTTGTAGTTAATTAATTTTTTTTTTATCTCGGGTAATTTTTATTTTTCCTTTCATAAGAATACATTACCATATCCAAAAACAAAAGAAAAAAAAAGTACCCTGGGCTAAAAAATTAACTACAATAGACACTCTCAACCAATGAAAGTGAAGTCAAGGTGGGTGGTGGACCAGATGAGACAAACACCCAAACATAGAACATGTGGAAATACATGTAACAACAATCAAATCTATTATGAGATCACATTGTTTGTTTAGGATGAAACTAGATGATAAGCAATCTACAAAAATTAGTGCTAGACTAAATTATACATGTATACATGTATAAGAGGTAATGTGCTCAAGGCTGACATTTAATATTTATTATTATTACTAATATTATTATTTGGTAGATAAAATGCATTTAGTTAAGTTAAGTTAAACAGGAAATTAAATGGGCAGTCCAATCTAGTATTTAAAGTAGAATTAAAGATGACAAATTTGGCATTCTTCTGGTATATCAGAATTTATCTTCTGTCTAATAATAACCCTTTTCACCCTGAAATTGTGTTTGCAGCCACCAGTAGCTCTACTACTGTAAAAACCTCCTCAGAATCTCTTAATGCTGCAGAGGTGGTGGAAGTACAAATGTCAGGTAAGAGCTATTCCCCGCTTATACTCTACATATAATAACAAAGATTGCATGTTTGTATCTTAATAATTATTTTTTTACTATTAATAGGCCAGTTCCCCCGGGCCTCTGTTTCAAAACTAGGGTAGGTGCTCAGCCTTTGGTATGGAGATCACTTTTTCATTCTCATGCAAATAAAACTCATTTTCGCAATAAAGGTTGTGCACCCAGCTTCATTTTGAAAGTGAGGGTTTTTGGAGCTCGCAGTTGGCCTATTTTAATGAGAGACAAGATTGTGTGATTATTATGTCACTGCAACCATCCACTGGCATGGCATTGCATGGTCTCTTTTCCTATAACCCTTACTAATTACAGTGGAGTCTCTATTAAAGGGGACACCCTCGGGACCAAGGCAAGAGGAGGTTGGGTTTTTAGTTAATATTGATAAAGGCATAAAATGTTTTCCTTTCATTTTGCCTTAAATCTGCTGTTGCCATAATTTTAAGCAGCTTGATAAAGCATTGCAAAATCATGAATAACAAACGTGTATTTCTGCGATGTTGCATGTTGAATTTCCAAAACTAGTACAATACAAATGATAATTGATGTACGTGAGATAATTCATGATTGACAAGCTGTTGTAATTATGACAAATGTACCCTTAAACTTATCAACAGGTTTCGTGGTCAGTCACTTTTTGTGTGCTAAGTTATCCCCTGAATGGTAAAGAGGTTAAACACAGGTTTTTCCTTTCTAAAAATAGAGGTGTCCCTTCAAGAGAGGTAACAAATAGTAATACAACAAAATATTTACATGTACCTCTCCAATCTTTACCCAGCTGCTGTTGCTTCTAATAATACTGTAGAAGAATCTAACATTGATGCTGAACCTGCAATGGTGGAAAACCAAGTGTCAGGTAAGTAAACAAATTGCTTGTTATAATAGGCCGACATTTTACAGTTATGACTTATGGATGGAAGCGAGAGTGGAGGTGACCTTATTTTGATACAAACCTCCTTCTTTGCCATGAAAATTGTCCTTCAAAAATAATAGTTAGCATAAGAACAACTTGATTCAACTTGAATCCAGTCAATTGTATGAAAGCTACGTAAACATACTCACTTAAGAACAAGGACGCTTTTCCTGGTGAGCCTTTCAGACGCCTTTTTTTTTTTAATTTCAAC
    ## >::Porites_evermani_scaffold_1:372245-372449
    ## GCCTTTTCTAGTTCCAGGTTTTAGTCTTTTTCAATAACGTTGGTTGTAAATTTTGTTTTTCCAACCTTTTTAACACTTAGAGTCTATTTGTAAGCCATTTTTATATTGTAAGGCAGGTTTTTTATCTCTGTGAAGCAGAAAACAGGGCATTATTTTTATGGATATACAGTGGAACCTCTCTAATACAGACACCGAAGGGACAGA
    ## >::Porites_evermani_scaffold_1:683878-684280
    ## TCTTGACTTTTACTTTTCGCTTTCTCTCCTCCCTTCTTTTTTCGCTTTTCTTGCCTCATTTTTTTTTCTCTTGCTGGGCATTTAGTAGGCTTCATTTTGGTGGGAAGAGTTTTTAGGAAAGCTTTTAGGATCTTAGGATTAGGTGAAAGGAAAGGTAGGTGGGTAATGGAACAAGATTTTCATGGAGATTTTCAGGTCCTTGTCACGTGGTTTTTTGCTTCTTTCTCCGGTGTCCTTGACTGAATTGTGCTCATTCTGGTATGGTTTGAAAGATCTCTTCACTCTGCACAAGTTAGCGAAGAAAGTTGTCCTTGACCGTTAAAACTGATGACGTCACAAAGGGTAGAAAGGACCTGGATCCGCACGGGCGGTTACGGGCGGTTCAGGGGCGAATGGGTTAAG
    ## >::Porites_evermani_scaffold_1:1202044-1202328
    ## GGGAGTAGCCTATGGATGAAAACATTTTGAGAACATGGTCAGGCAATGGTTTCGACTCCCTTTGGTCATAGCCTGCTCCAGGCGTTCTGATTGTGGAGCGTGGCGGCGTCGCTGTTTTTCCCGTCCCCACGATCTGAACGCCTGGAACAGGCTACTTTCGGCATTGCTAAACTTCTTACCCACAATTCGCGTTCCGTTTGTTTTTGTTGCTGTTGGCGGTTTTGTTGTTTTTTGTTGATGATGTTGTTGTTTTTTTAGTTCAGAGTGTTTCTCGACGTCTAGCG

# 2 make some blastdb

## 2.1 sRNA

``` bash
ls ../../DEF-cross-species/data/blast/peve_sRNA*
```

    ## ../../DEF-cross-species/data/blast/peve_sRNA.ndb
    ## ../../DEF-cross-species/data/blast/peve_sRNA.nhr
    ## ../../DEF-cross-species/data/blast/peve_sRNA.nin
    ## ../../DEF-cross-species/data/blast/peve_sRNA.not
    ## ../../DEF-cross-species/data/blast/peve_sRNA.nsq
    ## ../../DEF-cross-species/data/blast/peve_sRNA.ntf
    ## ../../DEF-cross-species/data/blast/peve_sRNA.nto

## 2.2 Genome

<https://gannet.fish.washington.edu/seashell/snaps/Porites_evermanni_v1.fa>

``` bash
cd ../data/
curl -O https://gannet.fish.washington.edu/seashell/snaps/Porites_evermanni_v1.fa
```

``` bash
head ../data/Porites_evermanni_v1.fa
```

    ## >Porites_evermani_scaffold_1
    ## GGCGGGGGGGGGGGGGGGGGGGGTACTCCCATACATTACCTATACGGGTATGTGCCGCCC
    ## AAAAGGGGCCGTGATTTTGAAGCTCCTGATTTAGAACGGGGTATCCATTTCAGAGGCGTT
    ## TTCTAGAACGGGGTGTAATATTTCGAACGCACGAAAGCTCCACTTTTGTGTAAGCAGCCA
    ## TTTGAAATTATTCAAGGACAGATTGCTTTTAAAAATACGGTTCAGCGCGTTAACAAGCAA
    ## ACCGTTGTACTCTTGTTGCACCCTAGAACGGTGTATAAAAAATTGGCCCATTTCTAGAAC
    ## GGGGTATCAGTTTTAGGGAGAATTCTAGAACGGGGTATAAAAAATTGGCCCTTTTCTGAA
    ## CGGGGCATCAATGTTAGGGGAAATTTTTTCCAGAACGGGGTGCCAATTTGGAGTCCCGGG
    ## CGGCACATACCCACCCAAAAAATACCCAAGTGCCCCCCCGGGGTCTAAACCCACATATTC
    ## TTCACACTGTTCACAATTTACCTCTTTTGGCTCTTCTAAGGAGAGCTCATCTAAATATTG

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/makeblastdb \
-in ../data/Porites_evermanni_v1.fa \
-dbtype nucl \
-out ../data/blast/pmea_genome
```

``` bash
ls ../data/blast/pmea_genome*
```

    ## ../data/blast/pmea_genome.ndb
    ## ../data/blast/pmea_genome.nhr
    ## ../data/blast/pmea_genome.nin
    ## ../data/blast/pmea_genome.not
    ## ../data/blast/pmea_genome.nsq
    ## ../data/blast/pmea_genome.ntf
    ## ../data/blast/pmea_genome.nto

## 2.3 Proteins

<https://gannet.fish.washington.edu/seashell/snaps/Porites_evermanni_v1.annot.pep.fa>

``` bash
cd ../data/

curl -O https://gannet.fish.washington.edu/seashell/snaps/Porites_evermanni_v1.annot.pep.fa
```

``` bash
head ../data/Porites_evermanni_v1.annot.pep.fa
```

    ## >Peve_00000001 assembled CDS
    ## KILADAGAVISGAQLGLGILTTILNTLGSISRKVAIGVDNESGYKWRAVNIYFYSGTSNR
    ## VLPHDVSSGTALLYGARKTAGPVARGAVGVLTYYIPHIDKTLAVMYSVPFDYNWYSNWFD
    ## VWLYSGKRRANYNLWYRMYYDNPFKGDNYWHERDLGSGLRARGAMSNSGQATIEIHILKQ
    ## *
    ## >Peve_00000002 assembled CDS
    ## MPSVSRSASPPSSPMSPPRKKTRFSLKLKSKKKASSESDPDKREKLYPIFSPPTAQNASD
    ## DETCEPPKPPPPYAGLRNHGNICYANAVVQVLRHCPGILESVDELDSLVKEKNYCATEGK
    ## SFDTRTDGKENGSEDLENCCDQKHVVCELRELFSQMEELEQDYNDNKENCEHLVQRRNTL
    ## VLAAKPLDFMQTFRGKNPLFEDNLQHDAQEFLCSLLVDLQDTETEIKKKREDPRKMYIGD

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/makeblastdb \
-in ../data/Porites_evermanni_v1.annot.pep.fa \
-dbtype prot \
-out ../data/blast/pmea_proteome
```

``` bash
ls ../data/blast/pmea_proteome*
```

    ## ../data/blast/pmea_proteome.pdb
    ## ../data/blast/pmea_proteome.phr
    ## ../data/blast/pmea_proteome.pin
    ## ../data/blast/pmea_proteome.pot
    ## ../data/blast/pmea_proteome.psq
    ## ../data/blast/pmea_proteome.ptf
    ## ../data/blast/pmea_proteome.pto

## 2.4 Genes

???

# 3 Blast comparison

# 4 Result - database: sRNA

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/blastn \
-task blastn \
-query ../../DEF-cross-species/data/peve_bedtools_lncRNAs.fasta \
-db ../../DEF-cross-species/data/blast/peve_sRNA \
-out ../output/02-Peve-lncRNA-align/lncRNA_sRNA_blastn.tab \
-evalue 1E-40 \
-num_threads 40 \
-max_target_seqs 5 \
-max_hsps 1 \
-outfmt 6
```

``` bash

echo "Number of hits?"
wc -l ../output/02-Peve-lncRNA-align/lncRNA_sRNA_blastn.tab


echo "File header"
head ../output/02-Peve-lncRNA-align/lncRNA_sRNA_blastn.tab
```

    ## Number of hits?
    ## 0 ../output/02-Peve-lncRNA-align/lncRNA_sRNA_blastn.tab
    ## File header

(nothing)

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/blastn \
-task blastn \
-query ../../DEF-cross-species/data/peve_bedtools_lncRNAs.fasta \
-db ../../DEF-cross-species/data/blast/peve_sRNA \
-out ../output/02-Peve-lncRNA-align/lncRNA_sRNA_blastn02.tab \
-num_threads 20 \
-outfmt 6
```

``` bash

echo "Number of hits?"
wc -l ../output/02-Peve-lncRNA-align/lncRNA_sRNA_blastn02.tab


echo "File header"
head ../output/02-Peve-lncRNA-align/lncRNA_sRNA_blastn02.tab
```

    ## Number of hits?
    ## 1798786 ../output/02-Peve-lncRNA-align/lncRNA_sRNA_blastn02.tab
    ## File header
    ## ::Porites_evermani_scaffold_1:422643-423512  1830741-1   100.000 35  0   0   726 760 35  1   3.58e-10    64.4
    ## ::Porites_evermani_scaffold_1:422643-423512  1848065-1   100.000 35  0   0   748 782 35  1   3.58e-10    64.4
    ## ::Porites_evermani_scaffold_1:422643-423512  907001-1    100.000 35  0   0   261 295 1   35  3.58e-10    64.4
    ## ::Porites_evermani_scaffold_1:422643-423512  874194-1    100.000 35  0   0   713 747 1   35  3.58e-10    64.4
    ## ::Porites_evermani_scaffold_1:422643-423512  439663-3    100.000 34  0   0   753 786 1   34  1.25e-09    62.6
    ## ::Porites_evermani_scaffold_1:422643-423512  2530059-1   97.143  35  1   0   723 757 1   35  1.52e-08    59.9
    ## ::Porites_evermani_scaffold_1:422643-423512  676574-2    97.143  35  1   0   734 768 1   35  1.52e-08    59.9
    ## ::Porites_evermani_scaffold_1:422643-423512  510681-2    97.143  35  1   0   753 787 1   35  1.52e-08    59.9
    ## ::Porites_evermani_scaffold_1:422643-423512  535611-2    97.143  35  1   0   737 771 1   35  1.52e-08    59.9
    ## ::Porites_evermani_scaffold_1:422643-423512  342166-4    97.143  35  1   0   752 786 35  1   1.52e-08    59.9

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/blastn \
-task blastn \
-query ../../DEF-cross-species/data/peve_bedtools_lncRNAs.fasta \
-db ../../DEF-cross-species/data/blast/peve_sRNA \
-out ../output/02-Peve-lncRNA-align/lncRNA_sRNA_blastn03.tab \
-evalue 1E-08 \
-num_threads 20 \
-outfmt 6
```

``` bash
echo "Number of hits?"
wc -l ../output/02-Peve-lncRNA-align/lncRNA_sRNA_blastn03.tab

echo "File header"
head ../output/02-Peve-lncRNA-align/lncRNA_sRNA_blastn03.tab
```

    ## Number of hits?
    ## 29152 ../output/02-Peve-lncRNA-align/lncRNA_sRNA_blastn03.tab
    ## File header
    ## ::Porites_evermani_scaffold_1:422643-423512  1830741-1   100.000 35  0   0   726 760 35  1   3.58e-10    64.4
    ## ::Porites_evermani_scaffold_1:422643-423512  1848065-1   100.000 35  0   0   748 782 35  1   3.58e-10    64.4
    ## ::Porites_evermani_scaffold_1:422643-423512  907001-1    100.000 35  0   0   261 295 1   35  3.58e-10    64.4
    ## ::Porites_evermani_scaffold_1:422643-423512  874194-1    100.000 35  0   0   713 747 1   35  3.58e-10    64.4
    ## ::Porites_evermani_scaffold_1:422643-423512  439663-3    100.000 34  0   0   753 786 1   34  1.25e-09    62.6
    ## ::Porites_evermani_scaffold_1:1084867-1089422    2127498-1   100.000 35  0   0   962 996 35  1   7.70e-10    64.4
    ## ::Porites_evermani_scaffold_1:683878-684280  1639214-1   100.000 35  0   0   222 256 1   35  1.60e-10    64.4
    ## ::Porites_evermani_scaffold_1:683878-684280  1195174-1   100.000 35  0   0   248 282 1   35  1.60e-10    64.4
    ## ::Porites_evermani_scaffold_1:683878-684280  1173775-1   100.000 35  0   0   244 278 1   35  1.60e-10    64.4
    ## ::Porites_evermani_scaffold_1:683878-684280  903402-1    100.000 35  0   0   262 296 1   35  1.60e-10    64.4

# 5 Result -database: genome

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/blastn \
-task blastn \
-query ../../DEF-cross-species/data/peve_bedtools_lncRNAs.fasta \
-db ../data/blast/pmea_genome \
-out ../output/02-Peve-lncRNA-align/lncRNA_genome_blastn.tab \
-num_threads 20 \
-outfmt 6
```

``` bash

echo "Number of hits?"
wc -l ../output/02-Peve-lncRNA-align/lncRNA_genome_blastn.tab


echo "File header"
head ../output/02-Peve-lncRNA-align/lncRNA_genome_blastn.tab
```

    ## Number of hits?
    ## 14547262 ../output/02-Peve-lncRNA-align/lncRNA_genome_blastn.tab
    ## File header
    ## ::Porites_evermani_scaffold_1:422643-423512  Porites_evermani_scaffold_1 100.000 869 0   0   1   869 422644  423512  0.0 1568
    ## ::Porites_evermani_scaffold_1:422643-423512  Porites_evermani_scaffold_1 94.958  119 6   0   663 781 1690031 1690149 7.22e-46    188
    ## ::Porites_evermani_scaffold_1:422643-423512  Porites_evermani_scaffold_1 90.196  102 10  0   694 795 1764562 1764461 3.28e-31    140
    ## ::Porites_evermani_scaffold_1:422643-423512  Porites_evermani_scaffold_1 81.690  142 19  3   657 794 967689  967827  4.87e-29    132
    ## ::Porites_evermani_scaffold_1:422643-423512  Porites_evermani_scaffold_1 84.337  83  9   2   714 795 642727  642806  6.34e-15    86.9
    ## ::Porites_evermani_scaffold_1:422643-423512  Porites_evermani_scaffold_1 87.879  66  5   1   695 760 766442  766504  7.72e-14    82.4
    ## ::Porites_evermani_scaffold_1:422643-423512  Porites_evermani_scaffold_1 87.692  65  5   1   731 795 1340298 1340237 2.69e-13    80.6
    ## ::Porites_evermani_scaffold_1:422643-423512  Porites_evermani_scaffold_1 83.333  72  12  0   663 734 1340494 1340423 3.28e-12    77.0
    ## ::Porites_evermani_scaffold_1:422643-423512  Porites_evermani_scaffold_1 89.091  55  6   0   663 717 1269951 1270005 4.00e-11    73.4
    ## ::Porites_evermani_scaffold_1:422643-423512  Porites_evermani_scaffold_1 80.645  62  9   1   730 791 1333238 1333180 3.07e-06    57.2

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/blastn \
-task blastn \
-query ../../DEF-cross-species/data/peve_bedtools_lncRNAs.fasta \
-db ../data/blast/pmea_genome \
-out ../output/02-Peve-lncRNA-align/lncRNA_genome_blastn02.tab \
-evalue 1E-40 \
-num_threads 20 \
-max_target_seqs 5 \
-max_hsps 1 \
-outfmt 6
```

``` bash

echo "Number of hits?"
wc -l ../output/02-Peve-lncRNA-align/lncRNA_genome_blastn02.tab


echo "File header"
head ../output/02-Peve-lncRNA-align/lncRNA_genome_blastn02.tab
```

    ## Number of hits?
    ## 26669 ../output/02-Peve-lncRNA-align/lncRNA_genome_blastn02.tab
    ## File header
    ## ::Porites_evermani_scaffold_1:422643-423512  Porites_evermani_scaffold_1 100.000 869 0   0   1   869 422644  423512  0.0 1568
    ## ::Porites_evermani_scaffold_1:422643-423512  Porites_evermani_scaffold_1058  95.035  141 7   0   655 795 26887   27027   3.50e-56    224
    ## ::Porites_evermani_scaffold_1:422643-423512  Porites_evermani_scaffold_940   94.366  142 8   0   654 795 11486   11627   1.22e-55    221
    ## ::Porites_evermani_scaffold_1:422643-423512  Porites_evermani_scaffold_546   94.444  144 7   1   653 795 161133  160990  1.22e-55    221
    ## ::Porites_evermani_scaffold_1:422643-423512  Porites_evermani_scaffold_2450  94.964  139 7   0   657 795 5502    5640    4.26e-55    220
    ## ::Porites_evermani_scaffold_1:1084867-1089422    Porites_evermani_scaffold_1 100.000 4555    0   0   1   4555    1084868 1089422 0.0 8215
    ## ::Porites_evermani_scaffold_1:1084867-1089422    Porites_evermani_scaffold_1840  85.342  921 100 10  16  932 36761   35872   0.0 1040
    ## ::Porites_evermani_scaffold_1:1084867-1089422    Porites_evermani_scaffold_1041  90.661  257 21  2   942 1197    52891   52637   1.88e-93    350
    ## ::Porites_evermani_scaffold_1:1084867-1089422    Porites_evermani_scaffold_705   91.600  250 16  4   955 1201    50079   49832   8.01e-92    343
    ## ::Porites_evermani_scaffold_1:1084867-1089422    Porites_evermani_scaffold_3651  91.129  248 19  2   953 1199    5737    5982    2.80e-91    343

# 6 Result - database: proteome

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/blastx \
-query ../../DEF-cross-species/data/peve_bedtools_lncRNAs.fasta \
-db ../data/blast/pmea_proteome \
-out ../output/02-Peve-lncRNA-align/lncRNA_proteome_blastx.tab \
-num_threads 20 \
-outfmt 6
```

``` bash

echo "Number of hits?"
wc -l ../output/02-Peve-lncRNA-align/lncRNA_proteome_blastx.tab


echo "File header"
head ../output/02-Peve-lncRNA-align/lncRNA_proteome_blastx.tab
```

    ## Number of hits?
    ## 351260 ../output/02-Peve-lncRNA-align/lncRNA_proteome_blastx.tab
    ## File header
    ## ::Porites_evermani_scaffold_1:422643-423512  Peve_00020583   75.000  40  7   2   791 672 1   37  2.36e-04    40.0
    ## ::Porites_evermani_scaffold_1:422643-423512  Peve_00044132   66.667  27  9   0   740 660 143 169 3.17e-04    41.2
    ## ::Porites_evermani_scaffold_1:422643-423512  Peve_00030450   60.606  33  12  1   740 642 926 957 4.05e-04    42.4
    ## ::Porites_evermani_scaffold_1:422643-423512  Peve_00029860   82.609  23  3   1   795 727 113 134 0.002   40.0
    ## ::Porites_evermani_scaffold_1:422643-423512  Peve_00012851   56.410  39  12  1   791 675 1   34  0.007   37.0
    ## ::Porites_evermani_scaffold_1:422643-423512  Peve_00000808   50.000  32  16  0   740 645 17  48  0.074   33.9
    ## ::Porites_evermani_scaffold_1:422643-423512  Peve_00043711   60.870  23  9   0   740 672 95  117 0.15    32.7
    ## ::Porites_evermani_scaffold_1:422643-423512  Peve_00014361   88.235  17  2   0   663 713 112 128 0.20    33.5
    ## ::Porites_evermani_scaffold_1:422643-423512  Peve_00010764   77.419  31  7   0   794 702 27  57  0.74    30.4
    ## ::Porites_evermani_scaffold_1:422643-423512  Peve_00019544   69.048  42  13  0   687 812 698 739 0.86    32.0

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/blastx \
-query ../../DEF-cross-species/data/peve_bedtools_lncRNAs.fasta \
-db ../data/blast/pmea_proteome \
-out ../output/02-Peve-lncRNA-align/lncRNA_proteome_blastx02.tab \
-evalue 1E-40 \
-num_threads 20 \
-max_target_seqs 5 \
-max_hsps 1 \
-outfmt 6
```

``` bash

echo "Number of hits?"
wc -l ../output/02-Peve-lncRNA-align/lncRNA_proteome_blastx02.tab


echo "File header"
head ../output/02-Peve-lncRNA-align/lncRNA_proteome_blastx02.tab
```

    ## Number of hits?
    ## 1684 ../output/02-Peve-lncRNA-align/lncRNA_proteome_blastx02.tab
    ## File header
    ## ::Porites_evermani_scaffold_10:260427-263111 Peve_00036837   51.678  149 68  3   1525    1965    1   147 5.36e-44    156
    ## ::Porites_evermani_scaffold_10:260440-263111 Peve_00036837   51.678  149 68  3   1512    1952    1   147 5.30e-44    156
    ## ::Porites_evermani_scaffold_10:25818-26461   Peve_00041908   52.577  97  42  1   6   296 76  168 1.79e-50    87.0
    ## ::Porites_evermani_scaffold_1011:155450-155947   Peve_00036185   97.203  143 4   0   452 24  150 292 1.15e-59    187
    ## ::Porites_evermani_scaffold_1013:148823-153034   Peve_00015807   58.736  269 47  5   3   809 59  263 2.61e-73    245
    ## ::Porites_evermani_scaffold_1013:148823-153034   Peve_00025583   57.249  269 51  3   3   809 58  262 1.90e-72    243
    ## ::Porites_evermani_scaffold_1013:148823-153034   Peve_00005800   57.249  269 51  3   3   809 59  263 6.57e-72    241
    ## ::Porites_evermani_scaffold_1013:148823-153034   Peve_00004106   58.736  269 47  3   3   809 60  264 4.23e-71    239
    ## ::Porites_evermani_scaffold_1013:148823-153034   Peve_00015805   55.390  269 56  3   3   809 59  263 7.84e-70    235
    ## ::Porites_evermani_scaffold_1016:8862-9109   Peve_00015308   87.500  80  10  0   2   241 1125    1204    2.38e-42    146
