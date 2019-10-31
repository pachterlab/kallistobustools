---
layout: page
title: "Constructing a transcriptome index with kb"
---

{% include JB/setup %}

This tutorial provides instructions for how to generate a transcriptome index to use with __kallisto &#124; bustools__ using `kb`.

__Note:__ for the instructions, command line arguments are preceeded by`$`. For example, if you see `$ cd my_folder` then type `cd my_folder`.

#### 0. Download materials
Prepare a folder
```
$ mkdir transcriptome_index/; cd transcriptome_index/
```

Download the genomic (DNA) FASTA and GTF annotations for your desired organism from the database of your choice. This tutorial uses mouse reference files downloaded from [Ensembl](https://uswest.ensembl.org/info/data/ftp/index.html).
```
$ wget ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
$ wget ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz
```

#### 1. Build the index
`kb` automatically splits the genome into a cDNA FASTA file and uses that to build a kallisto index.
```
$ kb ref -i transcriptome.idx -g transcripts_to_genes.txt -f1 cdna.fa Mus_musculus.GRCm38.dna.primary_assembly.fa.gz Mus_musculus.GRCm38.98.gtf.gz
```

#### 2. Align your reads and generate a count matrix
See [this](kb_getting_started.html) tutorial for how to proceed.
