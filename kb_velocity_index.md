---
layout: page
title: "Building a cDNA and intron index with kb"
---

{% include JB/setup %}

This tutorial provides instructions for how to generate a RNA velocity index to use with __kallisto &#124; bustools__ using `kb`.

__Note:__ for the instructions, command line arguments are preceeded by`$`. For example, if you see `$ cd my_folder` then type `cd my_folder`.

#### 0. Download materials
Prepare a folder
```
$ mkdir velocity_getting_started/; cd velocity_getting_started/
```

Download the genomic (DNA) FASTA and GTF annotations for your desired organism from the database of your choice. This tutorial uses human reference files downloaded from [Ensembl](https://uswest.ensembl.org/info/data/ftp/index.html).
```
$ wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
$ wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
```
Extract the files
```
$ gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
$ gunzip Homo_sapiens.GRCh38.98.gtf.gz
```

#### 1. Build the index
`kb` automatically splits the genome into a cDNA and intron FASTA file and uses these to build a kallisto index. This method is based on [La Manno, et al. 2019](https://doi.org/10.1038/s41586-018-0414-6).
```
$ kb ref -i index.idx -g transcripts_to_genes.txt -f1 cdna.fa -f2 intron.fa -c1 cdna_transcripts_to_capture.txt -c2 intron_transcripts_to_capture.txt --workflow lamanno Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.98.gtf
```
__Note__: The `--lamanno` option has been deprecated since version 0.24.5.

#### 2. Align your reads and generate a velocity matrix
See [this](kb_velocity_tutorial.html) tutorial for how to proceed.
