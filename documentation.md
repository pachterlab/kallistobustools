---
layout: page
title: "Documentation"
group: navigation
---

{% include JB/setup %}

#### kallisto and bustools manuals

The kallisto &#124; bustools single-cell RNA-seq workflow requires two programs: `kallisto` and `bustools`

The `kallisto` manual is available at: [https://pachterlab.github.io/kallisto/download](https://pachterlab.github.io/kallisto/download)

The `bustools` manual is available at: [https://bustools.github.io/manual](https://bustools.github.io/manual)

#### Overview of the `kallisto | bustools` workflow

For a hands-on introduction to the workflow see the [getting started](https://www.kallistobus.tools/getting_started) page and accompanying [video](https://www.youtube.com/watch?v=hWxnL86sak8).

#### Description of associated files: 

1. A pair of read files of the form read 1 and read 2 with `.fastq` or `.fastq.gz` extensions are required. For example, `SRR8599150_S1_L001_R1_001.fastq.gz` and `SRR8599150_S1_L001_R2_001.fastq.gz` are the files used in the [getting started](https://www.kallistobus.tools/getting_started) tutorial.

2. A set of target sequences, typically a reference transcriptome is needed for pseudoalignment. Species transcriptomes can be downloaded from the [Ensembl database page](https://uswest.ensembl.org/info/data/ftp/index.html). They are usually around 100MB in size. For example the mouse transcriptome downloaded from Ensembl is named `Mus_musculus.GRCm38.cdna.all.fa.gz`

3. A kallisto index must be constructed from the reference transcriptome (usually `.idx` extension). Indices can be downloaded (as long as they match the reference transcriptome), or built with the `kallisto index` command (for details see the [kallisto manual](https://pachterlab.github.io/kallisto/manual)). Standard indices can usually be built on a laptop with 8Gb of RAM in 10--30 min. depending on reference transcriptome size and hardware specifications. The resulting index file will be under 4GB in size. Pre built indices are available for the human transcriptome as well as many model organisms from the [kallisto transcriptome indices website](https://github.com/pachterlab/kallisto-transcriptome-indices).

4. A barcode whitelist with `.txt` extension must be input for barcode error correction. For example, a subset of barcodes from the 10x  Genomics v2 chemistry whitelist may look as follows:
```
AAACCTGAGAAACCAT
AAACCTGAGAAACCGC
AAACCTGAGAAACCTA
AAACCTGAGAAACGAG
AAACCTGAGAAACGCC
```

5. A transcript-to-gene file `.tsv` file is required. This is a tsv (tab separated value) file containing a mapping between the transcript Ensembl id and the gene Ensembl id. It may also have the gene name, but that is not required. Care must be taken to match the exact names in the reference transcriptome (e.g. Ensembl ID versions may or may not have been included). The `t2g.py` script can produce the transcript-to-gene file; a bustools command will be available shortly. For example, to make the transcript-to-gene file for the [getting started](https://www.kallistobus.tools/getting_started) tutorial withthe mouse `Mus_musculus.GRCm38.96.gtf` GTF file you'd use the command:
```
./t2g.py --use_version < Mus_musculus.GRCm38.96.gtf > transcripts_to_genes.txt
```
This will create the file `transcripts_to_genes.txt`. Some examples of such files are provided below: 
Gene names, no Ensembl version:
```
ENSMUST00000162897      ENSMUSG00000051951      Xkr4
ENSMUST00000159265      ENSMUSG00000051951      Xkr4
ENSMUST00000161581      ENSMUSG00000089699      Gm1992
ENSMUST00000194643      ENSMUSG00000102343      Gm37381
```
No gene names and no Ensembl version:
```
ENSMUST00000162897      ENSMUSG00000051951     
ENSMUST00000159265      ENSMUSG00000051951
ENSMUST00000161581      ENSMUSG00000089699
ENSMUST00000194643      ENSMUSG00000102343
```
No gene names and with Ensembl version:
```
ENSMUST00000162897.1      ENSMUSG00000051951     
ENSMUST00000159265.2      ENSMUSG00000051951
ENSMUST00000161581.1      ENSMUSG00000089699
ENSMUST00000194643.1      ENSMUSG00000102343
```

6. Running `kallisto` and `bustools`: once the read 1 and read 2 fastqs, thte kallisto index, and the transcript-to-gene files are ready, the count matrix can be generated with just a few commands; See the [getting started](https://www.kallistobus.tools/getting_started) tutorial. A detailed description of each file you should have at each step is provided at the [Getting Started Explained](https://www.kallistobus.tools/getting_started_explained.html) page.

7. Enjoy your pre-processed data! The [tutorials page](https://www.kallistobus.tools/tutorials) contains examples of how to parse BUS files and/or count matrices for downstream analysis.
