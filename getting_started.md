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

__Note: for this tutorial, command line arguments are everything after the `$` and user input is required anywhere you see `<user_input>`. So if you see `$ cd <my_folder>` then you would type `cd folder`.__

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

```
$ cd <destination_folder/>
```

Then download your index where the link is one that you copied.

```
[Linux/Mac]$ wget ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
[Windows]$ curl ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
```

Now do the exact same as above but instead of clicking `Download Fasta` click `Download GTF`. Right-click on `Mus_musculus.GRCm38.96.gtf` select `Copy Link Address` and download this file on your terminal.

### Dataset
Steps to Download the data 

__Summary:__ Download the index from ensembl, download the data. Type `$ ls -1` and you should see

```
$ ls -1
Mus_musculus.GRCm38.96.gtf.gz
Mus_musculus.GRCm38.cdna.all.fa.gz
SRR8599150_S1_L001_I1_001.fastq.gz
SRR8599150_S1_L001_R1_001.fastq.gz
SRR8599150_S1_L001_R2_001.fastq.gz
```

## 2. Setting up the index and gene map
### Build the index
We first need to decompress (unzip) the reference fasta file we downloaded.

```
$ gunzip Mus_musculus.GRCm38.cdna.all.fa.gz
```

Now we can build the kallisto index. I recommend naming the index `Mus_musculus.GRCm38.cdna.all.idx` and using a kmer size of `31`. Note that a kmer size of 31 is default, and always must be odd.

```
$ kallisto index -i <index_name.idx> -k <kmer_size> Mus_musculus.GRCm38.cdna.all.fa

[build] loading fasta file Mus_musculus.GRCm38.cdna.all.fa
[build] k-mer length: 31
[build] warning: clipped off poly-A tail (longer than 10)
        from 641 target sequences
[build] warning: replaced 3 non-ACGUT characters in the input sequence
        with pseudorandom nucleotides
[build] counting k-mers ... done.
[build] building target de Bruijn graph ...  done
[build] creating equivalence classes ...  done
[build] target de Bruijn graph has 734746 contigs and contains 100614952 k-mers

```
### Building the `transcripts_to_genes.txt` map
Insert steps here

## 3. Pseudoaligning reads with __kallisto bus__
First we take note of the technology that was used to generate our library. The 10x Chromium V2 chemistry was used to generate the data we downloaded above. Now we pseudo align the reads

```
$ kallisto bus -i <your_index.idx> -o <bus_output_folder/> -x 10xv2 -t 10 SRR8599150_S1_L001_R1_001.fastq.gz SRR8599150_S1_L001_R2_001.fastq.gz

[index] k-mer length: 31
[index] number of targets: 118,489
[index] number of k-mers: 100,614,952
[index] number of equivalence classes: 433,624
[quant] will process sample 1: SRR8599150_S1_L001_R1_001.fastq.gz
                               SRR8599150_S1_L001_R2_001.fastq.gz
[quant] finding pseudoalignments for the reads ... done
[quant] processed 8,860,361 reads, 3,431,849 reads pseudoaligned
```
__Note:__ The order of the `.fastq` files is important, `R1` comes first then `R2` goes second. Please see the [__Tutorials__] page.

## Processing BUS file with __bustools__

Other useful tutorial notebooks on the __BUStools__ repository include the [10x_hgmm_100 notebook](https://github.com/BUStools/BUS_notebooks_python/blob/master/dataset-notebooks/10x_hgmm_100_python/10x_hgmm_100.ipynb) which details the analysis of a small, and therefore easily downloadable dataset. Links to other tutorial notebooks are posted on the [__BUStools__ python notebook website](https://github.com/BUStools/BUS_notebooks_python) and the [__BUStools__ R notebook website](https://github.com/BUStools/BUS_notebooks_R).
