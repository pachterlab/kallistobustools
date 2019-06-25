---
layout: page
title: "Documentation"
group: navigation
---

{% include JB/setup %}

# kallisto and bustools manuals

The kallisto bus single cell workflow requires two programs: `kallisto` and `bustools`

To `kallisto` manual is available at: [https://pachterlab.github.io/kallisto/download](https://pachterlab.github.io/kallisto/download)

The `bustools` manual is available at: [https://bustools.github.io/manual](https://bustools.github.io/manual)

# Overview of the `kallisto | bus` workflow

For a complete walkthrough of the workflow see the [getting started](https://www.kallistobus.tools/getting_started) page. 
A video tutorial describing how to go from downloading the fastq file to generating the count matrix is [here](https://www.youtube.com/watch?v=hWxnL86sak8).

<p align="center">
<img src="https://user-images.githubusercontent.com/10369156/58990086-06c81500-879a-11e9-886b-7e4a690c5862.png" width="50%" />
</p>





# Description of files used 


###  Read 1 and read 2 `.fastq` or `.fastq.gz` files

For example, `SRR8599150_S1_L001_R1_001.fastq.gz` and `SRR8599150_S1_L001_R2_001.fastq.gz` are the files we use in the [getting started](https://www.kallistobus.tools/getting_started) tutorial.

### Barcode whitelist (usually `.txt`)

For example, in 10x v2 chemistry there are about  737,000 barcodes. They are provided in a `.txt` file that comes with cell ranger. One of the 10x v2 chemistry whitelists is called  `737K-august-2016.txt`, and this is what it looks like.
```
AAACCTGAGAAACCAT
AAACCTGAGAAACCGC
AAACCTGAGAAACCTA
AAACCTGAGAAACGAG
AAACCTGAGAAACGCC
```


### Species transcriptome (usually a `cdna.all.fa.gz`)

Species transcriptomes can be downloaded from the Ensembl website database page at [https://uswest.ensembl.org/info/data/ftp/index.html](https://uswest.ensembl.org/info/data/ftp/index.html). They are usually around 100MB in size. For example the mouse transcriptome downloaded from ensembl is named `Mus_musculus.GRCm38.cdna.all.fa.gz`


### Species kallisto index (usually `.idx`)

To process the sample with kallisto, first you need to use kallisto to build an index. For a given transcriptome and `kmer` size the transcriptome index needs to be built only once. The default `kmer` size is 31, and you don't need to change that. In order to build the index, first you need to download from Ensembl the appropriate transcriptome, described above.

You should be able to build the index on a laptop with 8GB RAM, it usually takes 15-30 min if you have enough memory. The resulting index file will be about 2GB. Build species indeices is done using the using  command `kallisto index`, see the [kallisto manual](https://pachterlab.github.io/kallisto/download) to learn more. Building the index only needs to be done once for each species. 

Pre built indices are available from the human transcriptome and many model organism transcriptomes are available from the [kallisto transcriptome indices](https://github.com/pachterlab/kallisto-transcriptome-indices) page.

###  Transcript to gene file `.tsv` file

This is a tsv (tab separated value) file containing a mapping between the transcript ensembl id and the gene ensembl id. It may also have the gene name, but that is not required. Be careful about ensembl id versions, they may or may not be present on the kallisto index depending on the file you used to generate the transcriptome. If you are getting empty matrices after running `bustools count` this is a likely reason (so check your file and remove ensembl id versions if needed.)

We offer the `t2g.py` script that you can use to generate this file from the organism GTF file, but if you have a transcript to gene map you can use that and skip dealing with GTF files. For example, if you're following the [getting started](https://www.kallistobus.tools/getting_started) tutorial then to make the transcript to gene map using `t2g.py` script to parse the mouse `Mus_musculus.GRCm38.96.gtf` GTF file you'd use the command:
```
./t2g.py --use_version < Mus_musculus.GRCm38.96.gtf > transcripts_to_genes.txt
```
This would create a file `transcripts_to_genes.txt`. **It is always a good idea to inspect the files you generate in this way. We provide a few examples of how they could look like below**

**Example excerpt of a mouse transcript to gene file including gene names, but no ensembl version**
```
ENSMUST00000162897      ENSMUSG00000051951      Xkr4
ENSMUST00000159265      ENSMUSG00000051951      Xkr4
ENSMUST00000161581      ENSMUSG00000089699      Gm1992
ENSMUST00000194643      ENSMUSG00000102343      Gm37381
```

**Example excerpt of a mouse transcript to gene file with no gene names and no ensembl version**
```
ENSMUST00000162897      ENSMUSG00000051951     
ENSMUST00000159265      ENSMUSG00000051951
ENSMUST00000161581      ENSMUSG00000089699
ENSMUST00000194643      ENSMUSG00000102343
```

**Example excerpt of a mouse transcript to gene file with no gene names and including ensembl version**
```
ENSMUST00000162897.1      ENSMUSG00000051951     
ENSMUST00000159265.2      ENSMUSG00000051951
ENSMUST00000161581.1      ENSMUSG00000089699
ENSMUST00000194643.1      ENSMUSG00000102343
```

**Example excerpt of an *Arabidopsis thaliana* transcript to gene file with gene names and including ensembl version**
```
AT1G01010.1     AT1G01010       NAC001
AT1G01020.1     AT1G01020       ARV1
AT1G01020.2     AT1G01020       ARV1
AT1G01030.1     AT1G01030       NGA3
```

**Example excerpt of a *Caenorhabditis elegans* transcript to gene file without gene names and including ensembl version**
```
Y74C9A.3        WBGene00022277
Y74C9A.2a.1     WBGene00022276
Y74C9A.2a.2     WBGene00022276
Y74C9A.2a.3     WBGene00022276
Y74C9A.2a.4     WBGene00022276
Y74C9A.2a.5     WBGene00022276
Y74C9A.2b       WBGene00022276
Y74C9A.4b       WBGene00022278
```

# Running `kallisto` and `bustools` in 3 easy steps

If you have the ** read 1 and read 2 fastqs**, ** kallisto index** and **transcript to gene** files ready, then you can generate the count matrix with just a few commands. These examples are from the [getting started](https://www.kallistobus.tools/getting_started) tutorial. A detailed description of each file you should have at each step is given at the [Getting Started Explained](https://www.kallistobus.tools/getting_started_explained.html) page.


## 1) Run `kallisto bus` to quantify the data
```
kallisto bus -i Mus_musculus.GRCm38.cdna.all.idx -o bus_output/ -x 10xv2 -t 10 SRR8599150_S1_L001_R1_001.fastq.gz SRR8599150_S1_L001_R2_001.fastq.gz
```

## 2) Correct and sort the bus file with `bustools correct` and `bustools sort`
```
cd bus_output/
bustools correct -w ../10xv2_whitelist.txt -o output.correct.bus output.bus
bustools sort -t 4 -o output.correct.sort.bus output.correct.bus
mkdir eqclass
mkdir genecount
```

## 3) Produce the gene count matrix with bustools
```
bustools count -o genecount/gene -g ../transcripts_to_genes.txt -e matrix.ec -t transcripts.txt --genecounts output.correct.sort.bus
```

Enjoy your freshly processed data!

3 optional) Or produce the transcript compatibility count matrix with bustools
```
bustools count -o eqclass/tcc -g ../transcripts_to_genes.txt -e matrix.ec -t transcripts.txt output.correct.sort.bus
```
