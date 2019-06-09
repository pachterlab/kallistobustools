---
layout: page
title: "Tutorial: Species Mixing"
---

{% include JB/setup %}

This page provides instructions for how to pre-process 1k 1:1 mixture of fresh frozen human (HEK293T) and mouse (NIH3T3) cells (v3 chemistry) from [the 10x website](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/1k_hgmm_v3). Details for each of the steps are expanded on the [explanation page](getting_started_explained.md).

__Note:__ for the instructions, command line arguments are preceeded by`$`. For example, if you see `$ cd my_folder` then type `cd my_folder`. 

#### 0. Download and install software
Obtain ```kallisto``` from the [__kallisto__ installation page](https://pachterlab.github.io/kallisto/download), and ```bustools``` from the [bustools installation page](https://github.com/BUStools/bustools).

#### 1. Download materials
Prepare a folder:
```
$ mkdir kallisto_bustools_species_mixing/; cd kallisto_bustools_species_mixing/
```
Download the following files:

- Mouse transcriptome `Mus_musculus.GRCm38.cdna.all.fa.gz`
- Human transcriptome `Homo_sapiens.GRCh38.cdna.all.fa.gz`
- 10x Chromium v3 chemistry barcode whitelist `10xv3_whitelist.txt`
- Mouse Transcripts to Genes map
- Human Transcripts to Genes map
- FastQ files from the 10x website

```
$ wget ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
$ wget ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
$ wget https://github.com/bustools/getting_started/releases/download/getting_started/10xv3_whitelist.txt
$ wget https://github.com/bustools/getting_started/releases/download/getting_started/GRCm38_transcripts_to_genes.txt
$ wget https://github.com/bustools/getting_started/releases/download/getting_started/GRCh38_transcripts_to_genes.txt
$ wget http://cf.10xgenomics.com/samples/cell-exp/3.0.2/1k_hgmm_v3/1k_hgmm_v3_fastqs.tar
$ tar -xvf 1k_hgmm_v3_fastqs.tar
```
#### 2. Build Index
Build the species index (alternatively download a pre-built index from the [kallisto transcriptome indices](https://github.com/pachterlab/kallisto-transcriptome-indices) page):
```
$ gunzip Mus_musculus.GRCm38.cdna.all.fa.gz
$ gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
$ cat Mus_musculus.GRCm38.cdna.all.fa Homo_sapiens.GRCh38.cdna.all.fa > GRCm38_GRCh38.fa
$ kallisto index -i GRCm38_GRCh38.idx -k 31 GRCm38_GRCh38.fa
```

#### 3. Run kallisto
Pseudoalign the reads:
```
$ kallisto bus -i GRCm38_GRCh38.idx -o bus_output/ -x 10xv3 -t 4 \
  1k_hgmm_v3_fastqs/1k_hgmm_v3_S1_L001_R1_001.fastq.gz 1k_hgmm_v3_fastqs/1k_hgmm_v3_S1_L001_R2_001.fastq.gz \
  1k_hgmm_v3_fastqs/1k_hgmm_v3_S1_L002_R1_001.fastq.gz 1k_hgmm_v3_fastqs/1k_hgmm_v3_S1_L002_R2_001.fastq.gz \
  1k_hgmm_v3_fastqs/1k_hgmm_v3_S1_L003_R1_001.fastq.gz 1k_hgmm_v3_fastqs/1k_hgmm_v3_S1_L003_R2_001.fastq.gz \
  1k_hgmm_v3_fastqs/1k_hgmm_v3_S1_L004_R1_001.fastq.gz 1k_hgmm_v3_fastqs/1k_hgmm_v3_S1_L004_R2_001.fastq.gz
```
#### 4. Run bustools
Correct, sort, and count the bus file. This creates the gene count matrix:
```
$ cd bus_output/
$ mkdir genecount/
$ bustools correct -w ../10xv3_whitelist.txt -p output.bus | bustools sort -t 4 -p | bustools count -o genecount/gene -g ../transcripts_to_genes.txt -e matrix.ec -t transcripts.txt --genecounts output.correct.sort.bus
```

#### 5. Load count matrices into notebook
