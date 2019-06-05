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
```http://uswest.ensembl.org/``` and select your species of interest. For getting started, select ```Mouse (Mus Musculus)```.

Once on this page ```http://uswest.ensembl.org/Mus_musculus/Info/Index``` select `Download Fasta` under the __Gene annotation__ section. This will take you to ```ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/```. Select `cdna`. Right-click on ```Mus_musculus.GRCm38.cdna.all.fa.gz``` and select ```Copy Link Address```.

On your terminal make a folder where you want to download your index and data.

```
$ mkdir kallisto_bustools_getting_started
$ cd kallisto_bustools_getting_started
```

Then download the fasta reference using the link is one that you copied.

```
[Linux/Mac]$ wget ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
[Windows]$ curl ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
```

Now download the GTF file. Do the exact same as above but instead of clicking ```Download Fasta``` click ```Download GTF``` under the __Gene annotation__ section. Right-click on ```Mus_musculus.GRCm38.96.gtf``` select ```Copy Link Address``` and download this file on your terminal.

### Dataset
Steps to download the data 

### Barcode whitelist
Steps to download the barcode whitelist

__tl;dr/Summary:__ Download the transcriptome reference and GTF file from ensembl, download the data, and download the barcode whitelist. Type `$ ls -1` and you should see

```
$ ls -1
whitelist.txt
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

Now we can build the kallisto index. I recommend naming the index ```Mus_musculus.GRCm38.cdna.all.idx``` and using a kmer size of `31`. Note that a kmer size of 31 is default, and always must be odd.

```
$ kallisto index -i Mus_musculus.GRCm38.cdna.all.idx -k 31 Mus_musculus.GRCm38.cdna.all.fa

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

## 3. Pseudoaligning reads with ```kallisto bus```
The 10x Chromium V2 chemistry was used to generate the data we downloaded above. The technology dictates the Barcode/UMI structure and the whitelist used for barcode error correction. We have to specify the technology in the __kallisto bus__ command and the whitelist in the __bustools__ command. Now we pseudo align the reads

```
$ kallisto bus -i Mus_musculus.GRCm38.cdna.all.idx -o bus_output/ -x 10xv2 -t 10 SRR8599150_S1_L001_R1_001.fastq.gz SRR8599150_S1_L001_R2_001.fastq.gz

[index] k-mer length: 31
[index] number of targets: 118,489
[index] number of k-mers: 100,614,952
[index] number of equivalence classes: 433,624
[quant] will process sample 1: SRR8599150_S1_L001_R1_001.fastq.gz
                               SRR8599150_S1_L001_R2_001.fastq.gz
[quant] finding pseudoalignments for the reads ... done
[quant] processed 8,860,361 reads, 3,431,849 reads pseudoaligned
```
__Note:__ For single cell sequencing you always need at least two fastq files and the order of the ```.fastq``` files is important, ```R1``` comes first then ```R2``` goes second. Please see the __Tutorials__ page if you want to know how to process more than one set of fastq files in one go.

## Processing BUS file with ```bustools```
```bustools``` allows us to go from a __BUS__ file, to a equivalence-class-UMI count matrix or a gene-UMI count matrix that can be loaded directly into python for analysis. We will use __bustools__ to do the following: 

1. Correct the barcodes: fix the barcodes that are within one hamming distance of the barcodes in the whitelist using ```whitelist.txt```
2. Sort the busfile: organize the busfile by barcode, umi, set, and multiplicity
3. Count the busfile: generate the umi count matrix using ```transcripts_to_genes.txt```

First navigate to your bus output directory. Your folder should contain the following items:
```
$ cd bus_output/
$ ls -1
matrix.ec
output.bus
run_info.json
transcripts.txt
```

Second correct the barcodes using the `whitelist.txt`. This makes a corrected bus file ```output.correct.bus```
```
$ bustools correct -w ../whitelist.txt -o output.correct.bus output.bus
Found 737280 barcodes in the whitelist
Number of hamming dist 1 barcodes = 20550336
Processed 3431849 bus records
In whitelist = 3281671
Corrected = 36927
Uncorrected = 113251
```

Third sort the busfile. This makes a sorted bus file ```output.correct.sort.bus```
```
$ bustools sort -t 4 -o output.correct.sort.bus output.correct.bus
Read in 3318598 number of busrecords
```

Fourth count the busfile using the `transcripts_to_genes.txt` to make the Equivalence Class Matrix (TCC)
```
$ bustools count -o eqcount/tcc -g ../transcripts_to_genes.txt -e matrix.ec -t transcripts.txt output.correct.sort.bus
bad counts = 0, rescued  =0, compacted = 0
```

or the Gene Count matrix
```
$ bustools count -o genecount/genes -g ../transcripts_to_genes.txt -e matrix.ec -t transcripts.txt --genecounts output.correct.sort.bus
bad counts = 0, rescued  =0, compacted = 0
```

Now you have your matrices!



Other useful tutorial notebooks on the __BUStools__ repository include the [10x_hgmm_100 notebook](https://github.com/BUStools/BUS_notebooks_python/blob/master/dataset-notebooks/10x_hgmm_100_python/10x_hgmm_100.ipynb) which details the analysis of a small, and therefore easily downloadable dataset. Links to other tutorial notebooks are posted on the [__BUStools__ python notebook website](https://github.com/BUStools/BUS_notebooks_python) and the [__BUStools__ R notebook website](https://github.com/BUStools/BUS_notebooks_R).
