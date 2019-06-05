---
layout: page
title: "Getting Started"
group: navigation
---

{% include JB/setup %}

In this section we will walk through how to download a dataset, process it with ```kallisto``` and ```bustools```, and load umi count matrices into a jupyter notebook for downstream processing. 

Before we begin, make sure that you have 
1. downloaded and installed ```kallisto``` from the [__kallisto__ installation page](https://pachterlab.github.io/kallisto/download), and have
2. downloaded and installed ```bustools``` from the [__bustools__ repository](https://github.com/BUStools/bustools).

### These are all of the commands that we will run in this tutorial
```
$ mkdir kallisto_bustools_getting_started
$ cd kallisto_bustools_getting_started
$ wget ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
$ gunzip Mus_musculus.GRCm38.cdna.all.fa.gz
$ kallisto index -i Mus_musculus.GRCm38.cdna.all.idx -k 31 Mus_musculus.GRCm38.cdna.all.fa
$ kallisto bus -i Mus_musculus.GRCm38.cdna.all.idx -o bus_output/ -x 10xv2 -t 10 SRR8599150_S1_L001_R1_001.fastq.gz SRR8599150_S1_L001_R2_001.fastq.gz
$ bustools correct -w ../whitelist.txt -o output.correct.bus output.bus
$ bustools sort -t 4 -o output.correct.sort.bus output.correct.bus
$ bustools count -o eqcount/tcc -g ../transcripts_to_genes.txt -e matrix.ec -t transcripts.txt output.correct.sort.bus
```

__Note:__ for these instructions, command line arguments are everything after the `$`. So if you see `$ cd my_folder` then you would type `cd my_folder` on your terminal.  

&nbsp;
&nbsp;
&nbsp;

## 0. Make sure that ```kallisto``` and ```bustools``` are installed correctly
Open up your terminal and run the following commands. These are the expected outputs:

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

If you don't see this then you have either (a) not installed the programs correctly, or (b) you have not told your terminal to "point" to the program so that you can use it. See *insert here* for how to correct these issues.  

&nbsp;
&nbsp;
&nbsp;

## 1. Download a reference, whitelist, and dataset
### (a) Reference
#### Transcriptome FASTA
The kallisto | bustools workflow uses a standard ensembl transcriptome fasta file reference to build an index. This index makes it easy (and fast!) to pseudoalign RNA sequencing reads. 

Navigate to the ensembl website ```http://uswest.ensembl.org/``` and select your species of interest. We will be using the ```Mouse (Mus Musculus)``` reference.

Once on ```http://uswest.ensembl.org/Mus_musculus/Info/Index``` select `Download Fasta` under the __Gene annotation__ section. This will take you to ```ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/```. Select `cdna`. Right-click on ```Mus_musculus.GRCm38.cdna.all.fa.gz``` and select ```Copy Link Address```.

On your terminal make a folder where you want to download your index and data.

```
$ mkdir kallisto_bustools_getting_started
$ cd kallisto_bustools_getting_started
```

Then download the fasta reference using the link is one that you copied.

```
$ wget ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
```

#### Genome annotation GTF
Next download the GTF file. Do the exact same as above but instead of clicking ```Download Fasta``` click ```Download GTF``` under the __Gene annotation__ section. Right-click on ```Mus_musculus.GRCm38.96.gtf``` select ```Copy Link Address``` and download this file on your terminal.

### (b) Barcode whitelist
Steps to download the barcode whitelist

### (c) Dataset
Steps to download the data 

### tl;dr/Summary
Download the transcriptome reference and GTF file from ensembl, download the barcode whitelist, and download the data. You should have the following files

```
kallisto_bustools_getting_started/
├── Mus_musculus.GRCm38.96.gtf.gz
├── Mus_musculus.GRCm38.cdna.all.fa.gz
├── SRR8599150_S1_L001_I1_001.fastq.gz
├── SRR8599150_S1_L001_R1_001.fastq.gz
├── SRR8599150_S1_L001_R2_001.fastq.gz
└── whitelist.txt

0 directories, 6 files
```  

&nbsp;
&nbsp;
&nbsp;

## 2. Build the index and gene map
### (a) Index
We first need to decompress (unzip) the reference fasta file we downloaded.

```
$ gunzip Mus_musculus.GRCm38.cdna.all.fa.gz
```

Now we can build the kallisto index. I recommend naming the index ```Mus_musculus.GRCm38.cdna.all.idx``` and using a kmer size of `31`. Note that a kmer size must always be odd and defaults to 31.

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

### (b) Gene map
As above, decompress (unzip) the genome GTF file we downloaded.
```
$ gunzip Mus_musculus.GRCm38.96.gtf.gz
```
__CHANGE!!__
Next use ```bustools genemap``` to make the gene map
```
$ bustools genemap -o transcripts_to_genes.txt Mus_musculus.GRCm38.96.gtf.gz
```

### tl;dr/Summary
Build the kallisto index from the reference fasta file, build the transcripts to genes map. You should have the following files

```
kallisto_bustools_getting_started/
├── Mus_musculus.GRCm38.96.gtf
├── Mus_musculus.GRCm38.cdna.all.fa
├── Mus_musculus.GRCm38.cdna.all.idx
├── SRR8599150_S1_L001_I1_001.fastq.gz
├── SRR8599150_S1_L001_R1_001.fastq.gz
├── SRR8599150_S1_L001_R2_001.fastq.gz
├── transcripts_to_genes.txt
└── whitelist.txt

0 directories, 8 files
```  

&nbsp;
&nbsp;
&nbsp;

## 3. Pseudoalign the reads with ```kallisto bus```
The 10x Chromium V2 technology was used to generate the data we downloaded above. The technology dictates the Barcode/UMI structure and the whitelist used for barcode error correction. We have to specify the technology in the ```kallisto bus``` command and the whitelist in the ```bustools``` command. 

### (a) Run ```kallisto bus``` to pseudoalign the reads
This will create a BUS file which will be located in ```bus_output/```.
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
__Note:__ You always need an even number of fastq files and the order of the ```.fastq``` files is important, ```R1``` comes first then ```R2``` goes second. Please see the __Tutorials__ page if you want to know how to process more than one set of fastq files in one go.

### tl;dr/Summary
Pseudoalign the single-cell RNA-seq reads using ```kallisto bus```. You should have the following files 

```
kallisto_bustools_getting_started/
├── bus_output
│   ├── matrix.ec
│   ├── output.bus
│   ├── run_info.json
│   └── transcripts.txt
├── Mus_musculus.GRCm38.96.gtf
├── Mus_musculus.GRCm38.cdna.all.fa
├── Mus_musculus.GRCm38.cdna.all.idx
├── SRR8599150_S1_L001_I1_001.fastq.gz
├── SRR8599150_S1_L001_R1_001.fastq.gz
├── SRR8599150_S1_L001_R2_001.fastq.gz
├── transcripts_to_genes.txt
└── whitelist.txt

1 directories, 12 files
```  

&nbsp;
&nbsp;
&nbsp;

## 4. Process the BUS file with ```bustools```
```bustools``` allows us to go from a __BUS__ file, to a equivalence-class-UMI count matrix or a gene-UMI count matrix that can be loaded directly into python for analysis. We will use __bustools__ to do the following: 

1. correct the barcodes: fix the barcodes that are within one hamming distance of the barcodes in the whitelist using ```whitelist.txt```,
2. sort the busfile: organize the busfile by barcode, umi, set, and multiplicity, and
3. count the busfile: generate the umi count matrix using ```transcripts_to_genes.txt```.

First navigate to your bus output directory. Your folder should contain the following items:
```
$ cd bus_output/
$ ls -1
matrix.ec
output.bus
run_info.json
transcripts.txt
```
### (a) Correct the barcodes in the busfile with ```bustools correct``` and the `whitelist.txt`
This makes a corrected bus file ```output.correct.bus```
```
$ bustools correct -w ../whitelist.txt -o output.correct.bus output.bus
Found 737280 barcodes in the whitelist
Number of hamming dist 1 barcodes = 20550336
Processed 3431849 bus records
In whitelist = 3281671
Corrected = 36927
Uncorrected = 113251
```
### (b) Sort the busfile with ```bustools sort```
This makes a sorted bus file ```output.correct.sort.bus```. This step __cannot__ be skipped. Sorting takes BUS records that are the same and "collapses them" into one BUS record with multiplicity. Note that this is different than UMI collapsing and serves the purpose of making the bus file smaller and making UMI counting more efficient.
```
$ bustools sort -t 4 -o output.correct.sort.bus output.correct.bus
Read in 3318598 number of busrecords
```

### (c) Count the UMIs in the busfile with ```bustools count``` and the ```transcripts_to_genes.txt```
For organization first make the following two folders:
```
$ mkdir eqclass
$ mkdir genecount
```

then make the Equivalence Class Matrix (TCC), 
```
$ bustools count -o eqcount/tcc -g ../transcripts_to_genes.txt -e matrix.ec -t transcripts.txt output.correct.sort.bus
bad counts = 0, rescued  =38627, compacted = 65899
```

or the Gene Count Matrix
```
$ bustools count -o genecount/genes -g ../transcripts_to_genes.txt -e matrix.ec -t transcripts.txt --genecounts output.correct.sort.bus
bad counts = 0, rescued  =0, compacted = 0
```

Now you have your matrices!

### tl;dr/Summary
Use ```bustools``` to correct the barcodes in the busfile, sort the busfile, and count the busfile. After all of these steps, your files should be structured like this:
```
kallisto_bustools_getting_started/
├── bus_output
│   ├── eqclass
│   │   ├── tcc.barcodes.txt
│   │   ├── tcc.ec.txt
│   │   └── tcc.mtx
│   ├── genecount
│   │   ├── gene.barcodes.txt
│   │   ├── gene.genes.txt
│   │   └── gene.mtx
│   ├── matrix.ec
│   ├── output.bus
│   ├── output.correct.bus
│   ├── output.correct.sort.bus
│   ├── run_info.json
│   └── transcripts.txt
├── Mus_musculus.GRCm38.96.gtf
├── Mus_musculus.GRCm38.cdna.all.fa
├── Mus_musculus.GRCm38.cdna.all.idx
├── SRR8599150_S1_L001_I1_001.fastq.gz
├── SRR8599150_S1_L001_R1_001.fastq.gz
├── SRR8599150_S1_L001_R2_001.fastq.gz
├── transcripts_to_genes.txt
└── whitelist.txt

3 directories, 20 files
```

And now we can load the data into python.  

&nbsp;
&nbsp;
&nbsp;

## 5. Load and analyze the matrices with python/R
Analysis of matrices with python: python-notebook
Analysis of matrices with R: R-notebook

Other useful tutorial notebooks on the __BUStools__ repository include the [10x_hgmm_100 notebook](https://github.com/BUStools/BUS_notebooks_python/blob/master/dataset-notebooks/10x_hgmm_100_python/10x_hgmm_100.ipynb) which details the analysis of a small, and therefore easily downloadable dataset. Links to other tutorial notebooks are posted on the [__BUStools__ python notebook website](https://github.com/BUStools/BUS_notebooks_python) and the [__BUStools__ R notebook website](https://github.com/BUStools/BUS_notebooks_R).
