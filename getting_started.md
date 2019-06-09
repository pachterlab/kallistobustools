---
layout: page
title: "Getting Started"
group: navigation
---

{% include JB/setup %}

This page provides instructions for how to pre-process (the blah blah dataset) from (blah blah). Details for each of the steps are expanded on the [explanation page](www.google.com).

__Note:__ for the instructions, command line arguments are preceeded by`$`. For example, if you see `$ cd my_folder` then type `cd my_folder`. 

#### 0. Download and install software
Obtain ```kallisto``` from the [__kallisto__ installation page](https://pachterlab.github.io/kallisto/download), and ```bustools``` from the [bustools installation page](https://github.com/BUStools/bustools).

#### 1. Download materials
Prepare a folder:
```
$ mkdir kallisto_bustools_getting_started
$ cd kallisto_bustools_getting_started
```
Download the following files:

- Mouse transcriptome `Mus_musculus.GRCm38.cdna.all.fa.gz`
- 10x Chromium v2 chemistry barcode whitelist `10xv2_whitelist.txt`
- Transcripts to Genes map
- Read 1 fastq file `SRR8599150_S1_L001_R1_001.fastq.gz`
- Read 2 fastq file `SRR8599150_S1_L001_R2_001.fastq.gz`

```
$ wget ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
$ wget https://github.com/bustools/getting_started/releases/download/getting_started/10xv2_whitelist.txt
$ wget transcripts_to_genes.txt
$ wget https://github.com/bustools/getting_started/releases/download/getting_started/SRR8599150_S1_L001_R1_001.fastq.gz
$ wget https://github.com/pachterlab/bustools/getting_started/releases/download/getting_started/SRR8599150_S1_L001_R2_001.fastq.gz
```
#### 2. Build Index
Build the species index (alternatively download a pre-built index from the [kallisto transcriptome indices](https://github.com/pachterlab/kallisto-transcriptome-indices) page):
```
$ gunzip Mus_musculus.GRCm38.cdna.all.fa.gz
$ kallisto index -i Mus_musculus.GRCm38.cdna.all.idx -k 31 Mus_musculus.GRCm38.cdna.all.fa
```

#### 3. Run kallisto
Pseudoalign the reads:
```
$ kallisto bus -i Mus_musculus.GRCm38.cdna.all.idx -o bus_output/ -x 10xv2 -t 10 SRR8599150_S1_L001_R1_001.fastq.gz SRR8599150_S1_L001_R2_001.fastq.gz
```
#### 4. Run bustools
Correct, sort, and count the bus file. This creates the gene count matrix:
```
$ cd bus_output/
$ mkdir genecount
$ bustools correct -w ../10xv2_whitelist.txt -p output.bus | bustools sort -t 4 -p | bustools count -o genecount/gene -g ../transcripts_to_genes.txt -e matrix.ec -t transcripts.txt --genecounts output.correct.sort.bus
```

#### 5. Load count matrices into notebook
 

<!---
&nbsp;
&nbsp;
&nbsp;

## 0. Download and install the software
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

If you don't see this then you have either the programs are not installed the programs correctly, or you have not added the location of the binaries to your `$PATH` variable in your `~/.bashrc` file in order to execute from any location. See [here](http://pachterlab.github.io/kallisto/local_build.html) if you need to make a local build without root.  

&nbsp;
&nbsp;
&nbsp;

## 1. Download a reference, whitelist, gene map utility, and dataset
### 1a. Reference
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

```
$ wget ftp://ftp.ensembl.org/pub/release-96/gtf/mus_musculus/Mus_musculus.GRCm38.96.gtf.gz
```

### 1b. Barcode whitelist
Navigate to https://github.com/pachterlab/kallistobuspaper_2019/releases/tag/getting_started right-click on ```10xv2_whitelist.txt``` select```Copy Link Address``` and download this file on your terminal.

```
$ wget https://github.com/pachterlab/kallistobuspaper_2019/releases/download/getting_started/10xv2_whitelist.txt
```

### 1c. Gene map utility
Navigate to https://github.com/pachterlab/kallistobuspaper_2019/releases/tag/getting_started right-click on ```t2g.py``` select```Copy Link Address``` and download this file on your terminal.

```
$ wget https://github.com/pachterlab/kallistobuspaper_2019/releases/download/getting_started/t2g.py
```

and make the script executable
```
$ chmod +x t2g.py
```

### 1d. Dataset
Navigate to `https://github.com/pachterlab/kallistobuspaper_2019/releases/tag/getting_started`, right-click on ```SRR8599150_S1_L001_R1_001.fastq.gz``` select```Copy Link Address``` and download this file on your terminal, and do the same for ```SRR8599150_S1_L001_R2_001.fastq.gz```.

This is a single cell experiment on mouse retinal cells. You can find the original data (GSM3612831) deposited here: https://www.ncbi.nlm.nih.gov/sra/?term=GSM3612831

```
$ wget https://github.com/pachterlab/kallistobuspaper_2019/releases/download/getting_started/SRR8599150_S1_L001_R1_001.fastq.gz
$ wget https://github.com/pachterlab/kallistobuspaper_2019/releases/download/getting_started/SRR8599150_S1_L001_R2_001.fastq.gz
```

### 1e. Results
After downloading the transcriptome reference and GTF file from ensembl, the barcode whitelist and the data you should have the following file structure:

```
kallisto_bustools_getting_started/
├── Mus_musculus.GRCm38.96.gtf.gz
├── Mus_musculus.GRCm38.cdna.all.fa.gz
├── SRR8599150_S1_L001_I1_001.fastq.gz
├── SRR8599150_S1_L001_R1_001.fastq.gz
├── SRR8599150_S1_L001_R2_001.fastq.gz
├── t2g.py
└── 10xv2_whitelist.txt

0 directories, 6 files
```  

&nbsp;
&nbsp;
&nbsp;

## 2. Build the index and gene map
### 2a. Index
The index only needs to be built once for each species transcriptome for a given k-mer size (the default k-mer 31 is suggested). You have the option of downloading the index or building it yourself. 
Prebuilt indices constructed from [Ensembl reference transcriptomes](https://uswest.ensembl.org/info/data/ftp/index.html) can be download from the [kallisto transcriptome indices](https://github.com/pachterlab/kallisto-transcriptome-indices/releases) site. Building indices with __kallisto index__ will often be faster in practice than downloading index files. For example, the __kallisto__ index for the mouse transcriptome takes between 5--10 minutes to build on a standard desktop or laptop. Transcriptome fasta files for model organisms can be downloaded from the [Ensembl database](https://www.ensembl.org/info/data/ftp/index.html). We recommend using cDNA fasta, specifically the *.cdna.all.fa.gz files. __kallisto__ can build indices directly from gzipped files.

#### Downloading the index ( if downloading the index then skip to Step 2b.).
If you wish to download the index then navigate to https://github.com/pachterlab/kallistobuspaper_2019/releases/tag/getting_started right-click on ```Mus_musculus.GRCm38.cdna.all.idx.gz``` select```Copy Link Address``` and download this file on your terminal.

```
$ wget https://github.com/pachterlab/kallistobuspaper_2019/releases/download/getting_started/Mus_musculus.GRCm38.cdna.all.idx.gz
```
And the decompress (unzip) the index we just downloaded:

```
$ gunzip Mus_musculus.GRCm38.cdna.all.idx.gz
```

#### Building the index yourself
We first need to decompress (unzip) the reference fasta file we downloaded.

```
$ gunzip Mus_musculus.GRCm38.cdna.all.fa.gz
```

Now we can build the kallisto index. We recommend naming the index `Mus_musculus.GRCm38.cdna.all.idx` and using a kmer size of `31`. Note that a kmer size must always be odd and defaults to `31`.

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

### 2b. Gene map
As above, decompress (unzip) the genome GTF file we downloaded.
```
$ gunzip Mus_musculus.GRCm38.96.gtf.gz
```
Next make transcript to gene map using `t2g.py` script to parse the mouse GTF file
```
$ ./t2g.py --use_version < Mus_musculus.GRCm38.96.gtf > transcripts_to_genes.txt
```

### 2c. Results
After building the kallisto index from the reference fasta file and parting the GTF files to create the transcripts to genes map you should have the following files:

```
kallisto_bustools_getting_started/
├── Mus_musculus.GRCm38.96.gtf
├── Mus_musculus.GRCm38.cdna.all.fa
├── Mus_musculus.GRCm38.cdna.all.idx
├── SRR8599150_S1_L001_I1_001.fastq.gz
├── SRR8599150_S1_L001_R1_001.fastq.gz
├── SRR8599150_S1_L001_R2_001.fastq.gz
├── t2g.py
├── transcripts_to_genes.txt
└── 10xv2_whitelist.txt

0 directories, 8 files
```  

&nbsp;
&nbsp;
&nbsp;

## 3. Pseudoalign the reads with `kallisto bus`
The 10x Chromium v2 technology was used to generate the data we downloaded above. The technology dictates the Barcode/UMI structure and the whitelist used for barcode error correction. We have to specify the technology in the ```kallisto bus``` command and the whitelist in the ```bustools``` command. 

### 3a. Run `kallisto bus` to pseudoalign the reads
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
__Note:__ For running `kallisto bus` you need an even number of fastq files and the order of the ```.fastq``` files is important, ```R1``` comes first then ```R2``` goes second. Please see the __Tutorials__ page if you want to know how to process more than one set of fastq files in one go.

### 3b. Results
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
├── t2g.txt
├── transcripts_to_genes.txt
└── 10xv2_whitelist.txt

1 directories, 12 files
```  

&nbsp;
&nbsp;
&nbsp;

## 4. Process the BUS file with `bustools`
```bustools``` allows us to go from a __BUS__ file, to a equivalence-class-UMI count matrix or a gene-UMI count matrix that can be loaded directly into python for analysis. We will use __bustools__ to do the following: 

1. Correct the barcodes using `bustools correct`: fix the barcodes that are within one hamming distance of the barcodes in the whitelist using ```whitelist.txt```,
2. Sort the busfile using `bustools sort`: organize the busfile by barcode, UMI, set, and multiplicity, and
3. Count records in the BUS with `bustools count`: generate the UMI count matrix using ```transcripts_to_genes.txt```.

First navigate to your bus output directory. Your folder should contain the following items:
```
$ cd bus_output/
$ ls -1
matrix.ec
output.bus
run_info.json
transcripts.txt
```
### 4a. Correct the barcodes in the busfile with `bustools correct` and the `whitelist.txt`
This makes a corrected bus file ```output.correct.bus```
```
$ bustools correct -w ../10xv2_whitelist.txt -o output.correct.bus output.bus
Found 737280 barcodes in the whitelist
Number of hamming dist 1 barcodes = 20550336
Processed 3431849 bus records
In whitelist = 3281671
Corrected = 36927
Uncorrected = 113251
```
### 4b. Sort the busfile with `bustools sort`
This makes a sorted bus file ```output.correct.sort.bus```. This step __cannot__ be skipped. Sorting takes BUS records that are the same and rearranges them into one BUS record with multiplicity column indicating how many times those records were seen. Note that this is different than UMI collapsing and serves the purpose of making the bus file smaller and making UMI counting more efficient.
```
$ bustools sort -t 4 -o output.correct.sort.bus output.correct.bus
Read in 3318598 number of busrecords
```

### 4c. Count the UMIs in the busfile with `bustools count` and the `transcripts_to_genes.txt`
For organization first make the following two folders:
```
$ mkdir eqcount
$ mkdir genecount
```

To make the Transcript Compatibility Count (TCC) Matrix we want the default output of `bustools count`
```
$ bustools count -o eqcount/tcc -g ../transcripts_to_genes.txt -e matrix.ec -t transcripts.txt output.correct.sort.bus
bad counts = 0, rescued  =38627, compacted = 65899
```

To make the Gene Count Matrix we need to give it the `--genecounts ` option
```
$ bustools count -o genecount/gene -g ../transcripts_to_genes.txt -e matrix.ec -t transcripts.txt --genecounts output.correct.sort.bus
bad counts = 0, rescued  =0, compacted = 0
```
Bustools will output the matrices in `.mtx` format, and gene names in a file ending as `.genes.txt` and barcodes in a file ending with `.barcodes.txt`. These can then be loaded into python or R for further analysis.

### Results
After using ```bustools``` to correct, sort and count the entries in the BUS file, your files should be structured like this:
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
├── t2g.py
├── transcripts_to_genes.txt
└── 10xv2_whitelist.txt

3 directories, 20 files
```  

&nbsp;
&nbsp;
&nbsp;

## 5. Load and analyze the matrices with python/R
Analysis of matrices with python: python-notebook

Analysis of matrices with R: R-notebook

Other useful tutorial notebooks on the __BUStools__ repository include the [10x_hgmm_100 notebook](https://github.com/BUStools/BUS_notebooks_python/blob/master/dataset-notebooks/10x_hgmm_100_python/10x_hgmm_100.ipynb) which details the analysis of a small, and therefore easily downloadable dataset. Links to other tutorial notebooks are posted on the [__BUStools__ python notebook website](https://github.com/BUStools/BUS_notebooks_python) and the [__BUStools__ R notebook website](https://github.com/BUStools/BUS_notebooks_R).
-->
