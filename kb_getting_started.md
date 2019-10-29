---
layout: page
title: "Getting Started"
group: navigation
---

{% include JB/setup %}

`kb` is a python package for rapidly pre-processing single-cell RNA-seq data. It is a wrapper for the methods described in [Melsted, Booeshaghi et al., 2019](https://doi.org/10.1101/673285). The goal of the wrapper is to simplify downloading and running of the kallisto and bustools programs. It was inspired by Sten Linnarsson’s [loompy fromfq command](http://linnarssonlab.org/loompy/kallisto/index.html).


The kb program consists of two parts:

The `kb count` command runs the kallisto and bustools programs. It can be used for pre-processing of data from a variety of single-cell RNA-seq technologies, and for a number of different workflows (e.g. production of gene count matrices, RNA velocity analyses, etc.). The output can be saved in a variety of formats including mix and loom. Tutorials are provided below.

The `kb ref` command builds or downloads a species-specific index for pseudoalignment of reads. This command must be run prior to `kb count`, and it runs the `kallisto index`.

#### 0. Install kb
Install `kb` from PyPi with `pip`:
```
pip install kb-python
```

#### 1. Download the materials
Prepare a folder:
```
$ mkdir kallisto_bustools_getting_started/; cd kallisto_bustools_getting_started/
```
Download the following files:
- Read 1 fastq file `SRR8599150_S1_L001_R1_001.fastq.gz`
- Read 2 fastq file `SRR8599150_S1_L001_R2_001.fastq.gz`

```
$ wget https://github.com/bustools/getting_started/releases/download/getting_started/SRR8599150_S1_L001_R1_001.fastq.gz
$ wget https://github.com/bustools/getting_started/releases/download/getting_started/SRR8599150_S1_L001_R2_001.fastq.gz
```

#### 2. Download the index
Download the pre-built mouse index using `kb`.
```
$ kb ref -d mouse -i index.idx -g transcripts_to_genes.txt
```

#### 3. Generate count matrices
The following command will
1. Pseudoalign the reads into a BUS file.
2. Correct, sort, and count the BUS file into a gene count matrix.
```
$ kb count -i index.idx -g transcripts_to_genes.txt -x 10xv2 -o output SRR8599150_S1_L001_R1_001.fastq.gz SRR8599150_S1_L001_R2_001.fastq.gz
```

#### 4. Load the count matrices into a notebook
See [this python notebook](https://github.com/BUStools/getting_started/blob/master/getting_started.ipynb) for how to load the count matrices into [ScanPy](https://scanpy.readthedocs.io/en/latest/index.html) for analysis.

<!-- #### Tutorials
- [Count matrices with `kb`](kb_count_matrix_tutorial.html)
- [Velocity matrices with `kb`](kb_velocity_matrix_tutorial.html) -->
