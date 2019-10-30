---
layout: page
title: "Getting Started"
group: navigation
---

{% include JB/setup %}

<p align="center">
  <a href="secret.html">
    <img src="assets/secret_tsne.jpg" width="70%">
  </a>
</p>

This page provides instructions for how to pre-process the [mouse retinal cells SRR8599150](https://www.ncbi.nlm.nih.gov/sra/?term=SRR8599150) dataset from [Koren et al., 2019](https://doi.org/10.1016/j.immuni.2019.02.007) using the __kallisto &#124; bustools workflow__.

__Note:__ command line arguments are preceeded by`$`. For example, if you see `$ cd my_folder` then type `cd my_folder`.


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
