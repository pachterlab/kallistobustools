---
layout: page
title: "Getting Started"
group: navigation
---

{% include JB/setup %}

In this section we will walk through how to download a dataset from the SRA, how to process it with __kallisto__ and __bustools__, and how to load the data into a jupyter notebook for downstream processing. 

Before we begin, make sure that you have 
1. downloaded and installed __kallisto__ from the [__kallisto__ installation page](https://pachterlab.github.io/kallisto/download), and have
2. downloaded and installed __bustools__ from the [__bustools__ repository](https://github.com/BUStools/bustools).

*Note: for this tutorial, command line arguments are everything after the '$' and user input is required anywhere you see <user_input>. So if you see `$ cd <my_folder>` then you would type `cd folder`.*

## 0. Make sure that __kallisto__ and __bustools__ are installed correctly
Open up your terminal and run the following commands. This is the expected output. 

```
$ kallisto
kallisto 0.45.1

Usage: kallisto <CMD> [arguments] ..

Where <CMD> can be one of:

    index         Builds a kallisto index
    quant         Runs the quantification algorithm
    bus           Generate BUS files for single-cell data
    pseudo        Runs the pseudoalignment step
    merge         Merges several batch runs
    h5dump        Converts HDF5-formatted results to plaintext
    inspect       Inspects and gives information about an index
    version       Prints version information
    cite          Prints citation information

Running kallisto <CMD> without arguments prints usage information for <CMD>
```
and 
```
$ bustools
Usage: bustools <CMD> [arguments] ..

Where <CMD> can be one of:

sort            Sort bus file by barcodes and UMI
text            Output as tab separated text file
merge           Merge bus files from same experiment
correct         Error correct bus files
count           Generate count matrices from bus file
capture         Capture reads mapping to a transcript capture list

Running bustools <CMD> without arguments prints usage information for <CMD>
```

If you don't see this, then you have either not installed the programs correctly, or you have not told your terminal to "point" to the program so that you can use it. See *insert here* for how to correct this.

## 1. Downloading a dataset & index
### Index
The kallisto | bustools workflow uses standard ensembl transcriptome fasta file to build an index. This index makes it easy (and fast!) to pseudoalign RNA sequencing reads. Navigate to the ensembl website:
`http://uswest.ensembl.org/` and select your species of interest. For getting started, select `Mouse (Mus Musculus)`.

Once on this page `http://uswest.ensembl.org/Mus_musculus/Info/Index` select `Download Fasta` under the __Gene annotation__ section. This will take you to `ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/`. Select `cdna`. Right-click on `Mus_musculus.GRCm38.cdna.all.fa.gz` and select `Copy Link Address`.

On your terminal navigate to a folder where you want to download your index and data.

```$ cd <destination_folder/>```

Then download your index where the link is one that you copied.

```
[Linux/Mac]$ wget ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
[Windows]$ curl ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
```

### Dataset
Steps to Download the data 


## 2. Pre-processing with __kallisto bus__
### Build the index
We first need to decompress (unzip) the reference fasta file we downloaded.

```$ gunzip Mus_musculus.GRCm38.cdna.all.fa.gz```

Now we can build the kallisto index. I recommend naming the index `Mus_musculus.GRCm38.cdna.all.idx` and using a kmer size of `31`. Note that a kmer size of 31 is default, and always must be odd.

```$ kallisto index Mus_musculus.GRCm38.cdna.all.fa -i <name_of_index.idx> -k <kmer_size> Mus_musculus.GRCm38.cdna.all.fa```

### Pseudoaligning reads



Other useful tutorial notebooks on the __BUStools__ repository include the [10x_hgmm_100 notebook](https://github.com/BUStools/BUS_notebooks_python/blob/master/dataset-notebooks/10x_hgmm_100_python/10x_hgmm_100.ipynb) which details the analysis of a small, and therefore easily downloadable dataset. Links to other tutorial notebooks are posted on the [__BUStools__ python notebook website](https://github.com/BUStools/BUS_notebooks_python) and the [__BUStools__ R notebook website](https://github.com/BUStools/BUS_notebooks_R).
