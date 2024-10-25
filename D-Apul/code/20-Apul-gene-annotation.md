20-Apul Annotation
================
Steven Roberts
26 August, 2024

- <a href="#1-protein" id="toc-1-protein">1 Protein</a>
- <a href="#2-download-swiss-prot-information"
  id="toc-2-download-swiss-prot-information">2 <strong>Download Swiss-Prot
  Information</strong></a>
- <a href="#3-join-blast-with-go-info"
  id="toc-3-join-blast-with-go-info">3 <strong>Join blast with GO
  info</strong></a>
- <a href="#4-rnas" id="toc-4-rnas">4 RNAs</a>
- <a href="#5-download-swiss-prot-information"
  id="toc-5-download-swiss-prot-information">5 <strong>Download Swiss-Prot
  Information</strong></a>
- <a href="#6-join-blast-with-go-info"
  id="toc-6-join-blast-with-go-info">6 <strong>Join blast with GO
  info</strong></a>

Apul

# 1 Protein

``` bash
cd ../data

curl -O https://gannet.fish.washington.edu/seashell/snaps/GCF_013753865.1_Amil_v2.1.protein.faa
```

``` bash
head ../data/GCF_013753865.1_Amil_v2.1.protein.faa
```

``` bash
ls ../../data/blast_dbs
```

``` bash
cd ../../data/blast_dbs
curl -O https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
mv uniprot_sprot.fasta.gz uniprot_sprot_r2024_04.fasta.gz
gunzip -k uniprot_sprot_r2024_04.fasta.gz
```

``` bash
head ../../data/blast_dbs/uniprot_sprot_r2024_04.fasta
echo "Number of Sequences"
grep -c ">" ../../data/blast_dbs/uniprot_sprot_r2024_04.fasta
```

``` bash
/home/shared/ncbi-blast-2.15.0+/bin/makeblastdb \
-in ../../data/blast_dbs/uniprot_sprot_r2024_04.fasta \
-dbtype prot \
-out ../../data/blast_dbs/uniprot_sprot_r2024_04
```

``` bash
fasta="../data/GCF_013753865.1_Amil_v2.1.protein.faa"

/home/shared/ncbi-blast-2.15.0+/bin/blastp \
-query $fasta \
-db ../../data/blast_dbs/uniprot_sprot_r2024_04 \
-out ../output/20-Apul-gene-annotation/blastp_out.tab \
-evalue 1E-05 \
-num_threads 48 \
-max_target_seqs 1 \
-max_hsps 1 \
-outfmt 6
```

``` bash
wc -l ../output/20-Apul-gene-annotation/blastp_out.tab
```

``` bash
tr '|' '\t' < ../output/20-Apul-gene-annotation/blastp_out.tab \
> ../output/20-Apul-gene-annotation/blastp_out_sep.tab

head -1 ../output/20-Apul-gene-annotation/blastp_out_sep.tab
```

![](http://gannet.fish.washington.edu/seashell/snaps/Monosnap__in_UniProtKB_search_571864__UniProt_2024-08-23_16-18-03.png)

# 2 **Download Swiss-Prot Information**

    sr320@raven:/home/shared/8TB_HDD_01/sr320/github/deep-dive/data$ curl -H "Accept: text/plain; format=tsv" "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Cgo_p%2Cgo_c%2Cgo%2Cgo_f%2Cgo_id&format=tsv&query=%28*%29+AND+%28reviewed%3Atrue%29" -o SwissProt-Annot-GO_082324.tsv

``` bash
wc -l ../../data/SwissProt-Annot-GO_082324.tsv 
```

# 3 **Join blast with GO info**

``` r
bltabl <- read.csv("../output/20-Apul-gene-annotation/blastp_out_sep.tab", sep = '\t', header = FALSE)

spgo <- read.csv("../../data/SwissProt-Annot-GO_082324.tsv", sep = '\t', header = TRUE)
```

``` r
annot_tab <- left_join(bltabl, spgo, by = c("V3" = "Entry")) %>%
  select(
    query = V1,
    blast_hit = V3,
    evalue = V13,
    ProteinNames = Protein.names,
    BiologicalProcess = Gene.Ontology..biological.process.,
    GeneOntologyIDs = Gene.Ontology.IDs
  )
```

``` r
head(annot_tab)
```

``` r
write.table(annot_tab, 
            file = "../output/20-Apul-gene-annotation/Apul-protein-GO.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
```

``` bash
head ../output/20-Apul-gene-annotation/Apul-protein-GO.tsv
```

# 4 RNAs

<https://gannet.fish.washington.edu/acropora/E5-deep-dive/Transcripts/Apul_GCF_013753865.1_rna.fna>

``` bash
cd ../data

curl -O https://gannet.fish.washington.edu/acropora/E5-deep-dive/Transcripts/Apul_GCF_013753865.1_rna.fna
```

``` bash
head ../data/Apul_GCF_013753865.1_rna.fna
```

``` bash
fasta="../data/Apul_GCF_013753865.1_rna.fna"

/home/shared/ncbi-blast-2.15.0+/bin/blastx \
-query $fasta \
-db ../../data/blast_dbs/uniprot_sprot_r2024_04 \
-out ../output/20-Apul-gene-annotation/blastx_out.tab \
-evalue 1E-05 \
-num_threads 48 \
-max_target_seqs 1 \
-max_hsps 1 \
-outfmt 6
```

``` bash
wc -l ../output/20-Apul-gene-annotation/blastx_out.tab
```

``` bash
tr '|' '\t' < ../output/20-Apul-gene-annotation/blastx_out.tab \
> ../output/20-Apul-gene-annotation/blastx_out_sep.tab

head -1 ../output/20-Apul-gene-annotation/blastx_out_sep.tab
```

![](http://gannet.fish.washington.edu/seashell/snaps/Monosnap__in_UniProtKB_search_571864__UniProt_2024-08-23_16-18-03.png)

# 5 **Download Swiss-Prot Information**

    sr320@raven:/home/shared/8TB_HDD_01/sr320/github/deep-dive/data$ curl -H "Accept: text/plain; format=tsv" "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Cgo_p%2Cgo_c%2Cgo%2Cgo_f%2Cgo_id&format=tsv&query=%28*%29+AND+%28reviewed%3Atrue%29" -o SwissProt-Annot-GO_082324.tsv

``` bash
wc -l ../../data/SwissProt-Annot-GO_082324.tsv 
```

# 6 **Join blast with GO info**

``` r
bltabl <- read.csv("../output/20-Apul-gene-annotation/blastx_out_sep.tab", sep = '\t', header = FALSE)

spgo <- read.csv("../../data/SwissProt-Annot-GO_082324.tsv", sep = '\t', header = TRUE)
```

``` r
annot_tab <- left_join(bltabl, spgo, by = c("V3" = "Entry")) %>%
  select(
    query = V1,
    blast_hit = V3,
    evalue = V13,
    ProteinNames = Protein.names,
    BiologicalProcess = Gene.Ontology..biological.process.,
    GeneOntologyIDs = Gene.Ontology.IDs
  )
```

``` r
head(annot_tab)
```

``` r
write.table(annot_tab, 
            file = "../output/20-Apul-gene-annotation/Apul-rna-GO.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
```

``` bash
head ../output/20-Apul-gene-annotation/Apul-rna-GO.tsv
```
