---
layout: page
title: "Getting Started"
group: navigation
---

{% include JB/setup %}

This page provides instructions for how to pre-process [mouse retinal cells SRR8599150](https://www.ncbi.nlm.nih.gov/sra/?term=SRR8599150) from [this paper](https://doi.org/10.1016/j.immuni.2019.02.007) using the __kallisto &#124; bustools workflow__. Details for each of the steps are expanded on the [explanation page](getting_started_explained.md).

__Note:__ for the instructions, command line arguments are preceeded by`$`. For example, if you see `$ cd my_folder` then type `cd my_folder`. 

#### 0. Download and install software
Obtain ```kallisto``` from the [__kallisto__ installation page](https://pachterlab.github.io/kallisto/download), and ```bustools``` from the [bustools installation page](https://github.com/BUStools/bustools).

#### 1. Download materials
Prepare a folder:
```
$ mkdir kallisto_bustools_getting_started/; cd kallisto_bustools_getting_started/
```
Download the following files:

- Human cDNA Transcripts `cDNA.correct_header.fa.gz`
- Human introns Transcripts `introns.correct_header.fa.gz`
- cDNA Transcripts to Capture `cDNA_transcripts.to_capture.txt.gz`
- Introns Transcripts to Capture `introns_transcripts.to_capture.txt.gz`
- cDNA/introns Transcripts to Genes map `cDNA_introns.t2g.txt.gz`
- 10x Chromium v2 chemistry barcode whitelist `10xv2_whitelist.txt`
- [SRR6470906 BAM File 1](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR6470906)
- [SRR6470907 BAM File 2](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR6470907)

```
$ wget https://github.com/BUStools/getting_started/releases/download/velocity_tutorial/cDNA.correct_header.fa.gz
$ wget https://github.com/BUStools/getting_started/releases/download/velocity_tutorial/introns.correct_header.fa.gz
$ wget https://github.com/BUStools/getting_started/releases/download/velocity_tutorial/cDNA_transcripts.to_capture.txt.gz
$ wget https://github.com/BUStools/getting_started/releases/download/velocity_tutorial/introns_transcripts.to_capture.txt.gz
$ wget https://github.com/BUStools/getting_started/releases/download/velocity_tutorial/cDNA_introns.t2g.txt.gz
$ wget https://github.com/BUStools/getting_started/releases/download/velocity_tutorial/10xv2_whitelist.txt
$ wget https://sra-download.ncbi.nlm.nih.gov/traces/sra57/SRZ/006470/SRR6470906/10X_17_029.bam
$ wget https://sra-download.ncbi.nlm.nih.gov/traces/sra57/SRZ/006470/SRR6470907/10X_17_028.bam
```
#### 2. Build Index
Build the species index (alternatively download a pre-built index from the [kallisto transcriptome indices](https://github.com/pachterlab/kallisto-transcriptome-indices) page):
```
$ gunzip Mus_musculus.GRCm38.cdna.all.fa.gz
$ kallisto index -i Mus_musculus.GRCm38.cdna.all.idx -k 31 Mus_musculus.GRCm38.cdna.all.fa
```

#### 3. Run kallisto
Pseudoalign the reads:
```
$ kallisto bus -i Mus_musculus.GRCm38.cdna.all.idx -o bus_output/ -x 10xv2 -t 10 SRR8599150_S1_L001_R1_001.fastq.gz SRR8599150_S1_L001_R2_001.fastq.gz
```
#### 4. Run bustools
Correct, sort, and count the bus file. This creates the gene count matrix:
```
$ cd bus_output/
$ mkdir genecount/ tmp/
$ bustools correct -w ../10xv2_whitelist.txt -p output.bus | bustools sort -T tmp/ -t 10 -p - | bustools count -o genecount/genes -g ../transcripts_to_genes.txt -e matrix.ec -t transcripts.txt --genecounts -
```

#### 5. Load count matrices into notebook
