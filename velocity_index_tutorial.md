---
layout: page
title: "RNA velocity tutorial"
---

{% include JB/setup %}

<p align="center">
  <a href="">
    <img src="assets/website_velocity.jpg" width="60%">
  </a>
</p>

This tutorial provides instructions for how to generate indicies to use with __kallisto &#124; bustools__ to perform an RNA velocity analysis. 

__Note:__ for the instructions, command line arguments are preceeded by`$`. For example, if you see `$ cd my_folder` then type `cd my_folder`. 

#### 0. Download and install software
Download and install ```bedtools``` from [here](https://bedtools.readthedocs.io/en/latest/content/installation.html).

#### 1. Determine your biological read length
Take your FASTQ file `R1.fastq.gz` find the length of the read
```
$ zcat R1.fastq.gz | head -2
@SRR8742283.1 NS500422:552:HJ5Y3BGX3:1:11101:21875:1038 length=61
CAGTCNTTTTTTTTAATTTAAAAAAAAAAAAAAGATTTATTAACAGTTTTAGAAGGCAGTT

$ echo -n CAGTCNTTTTTTTTAATTTAAAAAAAAAAAAAAGATTTATTAACAGTTTTAGAAGGCAGTT | wc -c
61
```
So `L = 61` as stated in the FASTQ header. `L` will vary so it is important that you check before proceeding.

#### 2. Download materials
Prepare a folder:
```
$ mkdir velocity_index/; cd velocity_index/
```
Download the **INTRONS BED file with L-1 flank**:

1. Go to the [UCSC table browser](http://genome.ucsc.edu/cgi-bin/hgTables?command=start).
2. Select desired species and assembly
3. Select group: `Genes and Gene Prediction Tracks`
4. Select track: UCSC Genes (or Refseq, Ensembl, etc.)
5. Select table: `knownGene`
6. Select region: `genome` (or you can test on a single chromosome or smaller region)
7. Select output format: `BED - browser extensible data`
8. Enter output file: `introns.tsv`
9. Select file type returned: `gzip compressed`
10. Select the 'get output' button
A second page of options relating to the BED file will appear.
11. Under 'create one BED record per:'. Select 'Introns plus'
12. Add flank `L - 1` flank
13. Hit the 'get BED' option

Download the **cDNA FASTA file**:

1. Go to the [UCSC table browser](http://genome.ucsc.edu/cgi-bin/hgTables?command=start).
2. Select desired species and assembly
3. Select group: `Genes and Gene Prediction Tracks`
4. Select track: UCSC Genes (or Refseq, Ensembl, etc.)
5. Select table: `knownGene`
6. Select region: `genome` (or you can test on a single chromosome or smaller region)
7. Select output format: `sequence`
8. Enter output file: `introns.tsv`
9. Select file type returned: `gzip compressed`
10. Hit the 'get output' button
11. Select `genomic` and click submit
A page of options relating to the FASTA file will appear.
12. Select `5' UTR Exons` & `CDS Exons` & `3' UTR Exons`
11. Select `One FASTA record per region (exon, intron, etc.) with  0 extra bases upstream (5') and  0 extra downstream (3')`
14. Select `All upper case`
13. Select `get sequence`

**Note**: You may ask why we don't just download the `sequence` of introns? The reason is because the FASTA file is large for complex organisms (you can do this for simple organisms) and the UCSC server times out after 20 minutes and results in a corrupted intron FASTA file.
