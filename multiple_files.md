---
layout: page
title: "Processing Multiple Lanes at Once"
---

{% include JB/setup %}

<!--- <p align="center">
  <a href="secret.html">
    <img src="assets/" width="70%">
  </a>
</p> -->

This page provides instructions for how to pre-process the [mouse T cells SRR8206317](https://www.ncbi.nlm.nih.gov/sra/?term=SRR8206317) dataset from [Miller & Sen et al., 2019](https://doi.org/10.1038/s41590-019-0312-6) using the __kallisto &#124; bustools workflow__. 

__Note:__ command line arguments are preceeded by`$`. For example, if you see `$ cd my_folder` then type `cd my_folder`. 

#### 0. Download and install the software
Obtain ```kallisto``` from the [__kallisto__ installation page](https://pachterlab.github.io/kallisto/download), and ```bustools``` from the [bustools installation page](https://bustools.github.io/download). A video tutorial for how to install the software can be viewed [here](https://youtu.be/thvtp7Ik6ts). Download and install ```bamtofastq``` from [here](https://support.10xgenomics.com/docs/bamtofastq) to generate the original FASTQ files from the BAM files provided by the authors. For a brief tutorial on how to install ```bamtofastq``` please see [this page](install_bamtofastq.html)

__Note:__ this dataset is v2 chemistry. If you would like to process v3 chemistry then you would use the [10xv3 whitelist](https://github.com/BUStools/getting_started/releases).

#### 1. Download the materials
Prepare a folder:
```
$ mkdir kallisto_bustools_multiple_lanes/; cd kallisto_bustools_multiple_lanes/
```
Download the following files:

- Mouse transcriptome `Mus_musculus.GRCm38.cdna.all.fa.gz`
- 10x Chromium v2 chemistry barcode whitelist `10xv2_whitelist.txt`
- Transcripts to Genes map
- Bam File `d10_Tet_possorted_genome_bam.bam`

Since the FASTQ files were not made available, we will download the BAM file from the European Nucleotide Archive and then use the 10x Genomics utility `bamtofastq`.

```
$ wget ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
$ wget https://github.com/bustools/getting_started/releases/download/getting_started/10xv2_whitelist.txt
$ wget https://github.com/bustools/getting_started/releases/download/getting_started/transcripts_to_genes.txt
$ wget ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR820/SRR8206317/d10_Tet_possorted_genome_bam.bam
```
#### 2. Build an index
Build the species index (alternatively download a pre-built index from the [kallisto transcriptome indices](https://github.com/pachterlab/kallisto-transcriptome-indices) page):
```
$ gunzip Mus_musculus.GRCm38.cdna.all.fa.gz
$ kallisto index -i Mus_musculus.GRCm38.cdna.all.idx -k 31 Mus_musculus.GRCm38.cdna.all.fa
```

#### 3. Generate the FASTQs from the BAM file
Use the `bamtofastq` utility to generate the FASTqs.
```
$ bamtofastq --reads-per-fastq=500000000 d10_Tet_possorted_genome_bam.bam ./fastqs
```

#### 3. Run kallisto
Pseudoalign the reads:
```
$ kallisto bus -i Mus_musculus.GRCm38.cdna.all.idx -o bus_output/ -x 10xv2 -t 4 \
bamtofastq_S1_L001_R1_001.fastq.gz \
bamtofastq_S1_L001_R2_001.fastq.gz \
bamtofastq_S1_L002_R1_001.fastq.gz \
bamtofastq_S1_L002_R2_001.fastq.gz \
bamtofastq_S1_L003_R1_001.fastq.gz \
bamtofastq_S1_L003_R2_001.fastq.gz \
bamtofastq_S1_L004_R1_001.fastq.gz \
bamtofastq_S1_L004_R2_001.fastq.gz 
```
**Note:** The `\` is to indicate the continuation of a line and will not affect the running of the program.

#### 4. Run bustools
Correct, sort, and count the bus file. This creates the gene count matrix:
```
$ cd bus_output/
$ mkdir genecount/ tmp/
$ bustools correct -w ../10xv2_whitelist.txt -p output.bus | bustools sort -T tmp/ -t 4 -p - | bustools count -o genecount/genes -g ../transcripts_to_genes.txt -e matrix.ec -t transcripts.txt --genecounts -
```

#### 5. Load the count matrices into a notebook [from getting started tutorial]
See [this python notebook](https://github.com/BUStools/getting_started/blob/master/getting_started.ipynb) for how to load the count matrices into [ScanPy](https://scanpy.readthedocs.io/en/latest/index.html) for analysis.
