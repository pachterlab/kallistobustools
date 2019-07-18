---
layout: page
title: "Building a cDNA and intron index"
---

{% include JB/setup %}

This tutorial provides instructions for how to generate indicies to use with __kallisto &#124; bustools__ to perform an RNA velocity analysis. 

__Note:__ for the instructions, command line arguments are preceeded by`$`. For example, if you see `$ cd my_folder` then type `cd my_folder`. 

#### 0. Download and install software
Download and install ```bedtools``` from [here](https://bedtools.readthedocs.io/en/latest/content/installation.html). Download and install the transcripts to genes binary from [here](https://github.com/sbooeshaghi/tools/releases/tag/t2g_v0.24.0).

#### 1. Determine your biological read length
Take your FASTQ file `R1.fastq.gz` find the length of the read
```
$ zcat R1.fastq.gz | head -2 ## note on a mac you would do zcat < R1.fastq.gz | head
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
8. Enter output file: `introns.bed`
9. Select file type returned: `gzip compressed`
10. Select the 'get output' button
A second page of options relating to the BED file will appear.
11. Under 'create one BED record per:'. Select 'Introns plus'
12. Add flank `L - 1` flank
13. Select the 'get BED' option
14. Save as `introns.bed.gz` to `velocity_index/`

Download the **cDNA FASTA file**:

1. Go to the [UCSC table browser](http://genome.ucsc.edu/cgi-bin/hgTables?command=start).
2. Select desired species and assembly
3. Select group: `Genes and Gene Prediction Tracks`
4. Select track: UCSC Genes (or Refseq, Ensembl, etc.)
5. Select table: `knownGene`
6. Select region: `genome` (or you can test on a single chromosome or smaller region)
7. Select output format: `sequence`
8. Enter output file: `cDNA.fa.gz`
9. Select file type returned: `gzip compressed`
10. Hit the 'get output' button
11. Select `genomic` and click submit
A page of options relating to the FASTA file will appear.
12. Select `5' UTR Exons` & `CDS Exons` & `3' UTR Exons`
11. Select `One FASTA record per region (exon, intron, etc.) with  0 extra bases upstream (5') and  0 extra downstream (3')`
14. Select `All upper case`
13. Select `get sequence`
14. Save as `cDNA.fa.gz` to `velocity_index/`

**Note**: You may ask why we don't just download the `sequence` of introns? The reason is because the FASTA file is large for complex organisms (you can do this for simple organisms) and the UCSC server times out after 20 minutes and results in a corrupted intron FASTA file.

Download the Genome
1. Go to [the website specieid by the **track** in the UCSC table browser](https://www.gencodegenes.org/human/release_29.html) ## Example Gencode 29 GRCh38
2. Selected desired species
3. Right click `FASTA` next to `Genome sequence, primary assembly (GRCh38)`
4. Download species.dna.primary_assembly.fa.gz (where species will be your specific species) 
```
$ wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz
```

Download the GTF and make a transcripts to genes map
1. Go to [ensembl](http://uswest.ensembl.org/index.html) or the website for whichever **track** specified in the UCSC table browser
2. Selected desired species
3. Select `Download GTF`
5. Download species.gtf.gz (where species will be your specific species) 
```
$ wget  ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz ## Homo sapiens GRCh38 example
```

#### 3. Make transcripts to genes map
Gunzip (decompress) the GTF file
```
$ gunzip species.gtf
```
Use the `t2g` utility to make thr transcripts to genes map (takes a minute or two for larger GTF files). Pick **one** of the commands to run.
```
$ cat genes.gtf | t2g make -p - > tr2g.txt
$ t2g make -p - < genes.gtf > tr2g.txt
$ t2g make --version -p - < genes.gtf > tr2g.txt       # (with version number)
```
**Note:** for this tutorial we will use the transcript IDs _with a version number_.

#### 4a. Convert INTRONS BED file to INTRON FASTA file
Gunzip (decompress) files
```
$ gunzip species.primary_assembly.fa.gz
$ gunzip introns.bed.gz
$ head -4 introns.bed
chr1	12118	12721	ENST00000456328.2_intron_0_109_chr1_12228_f	0	+
chr1	12612	13329	ENST00000456328.2_intron_1_109_chr1_12722_f	0	+
chr1	11948	12287	ENST00000450305.2_intron_0_109_chr1_12058_f	0	+
chr1	12118	12721	ENST00000450305.2_intron_1_109_chr1_12228_f	0	+
```
Using `bedtools getfasta` we will slice up the primary assembly with the BED file to give us a FASTA file of introns
```
$ bedtools getfasta -name -fo introns.fa -fi species.primary_assmebly.fa -bed introns.bed
$ head -2 introns.fa ## Homo sapiens GRCh38
>ENST00000456328.2_intron_0_109_chr1_12228_f
CTGCATGTAACTTAATACCACAACCAGGCATAGGGGAAAGATTGGAGGAAAGATGAGTGAGAGCATCAACTTCTCTCACAACCTAGGCCAGTAAGTAGTGCTTGTGCTCATCTCCTTGGCTGTGATACGTGGCCGGCCCTCGCTCCAGCAGCTGGACCCCTACCTGCCGTCTGCTGCCATCGGAGCCCAAAGCCGGGCTGTGACTGCTCAGACCAGCCGGCTGGAGGGAGGGGCTCAGCAGGTCTGGCTTTGGCCCTGGGAGAGCAGGTGGAAGATCAGGCAGGCCATCGCTGCCACAGAACCCAGTGGATTGGCCTAGGTGGGATCTCTGAGCTCAACAAGCCCTCTCTGGGTGGTAGGTGCAGAGACGGGAGGGGCAGAGCCGCAGGCACAGCCAAGAGGGCTGAAGAAATGGTAGAACGGAGCAGCTGGTGATGTGTGGGCCCACCGGCCCCAGGCTCCTGTCTCCCCCCAGGTGTGTGGTGATGCCAGGCATGCCCTTCCCCAGCATCAGGTCTCCAGAGCTGCAGAAGACGACGGCCGACTTGGATCACACTCTTGTGAG
```

#### 4b. Get the transcripts to capture list and transcripts to genes for INTRONS
First we need to get a list of all of the intronic transcript IDs represented in our FASTA file, with or without version numbers. 
```
$ cat introns.fa | awk '/^>/ {print $0}' | tr "_" " " | awk '{print $3"."$4}' > introns_transcripts.txt
$ cat introns_transcripts.txt | tr "." " " | awk '{print $1}' > introns_transcripts_no_version.txt
```
**Note:** your FASTA header may not be structured exactly like the example. If that is the case, you can change the columns that are printed in the `awk` command. For example, `awk '{print $3}'` prints the 3rd column.

**Note:** its important that the transcript IDs be consistent with the `tr2g.txt` that you made earlier, i.e. keep both with or without version.

Now we add an identifier to the transcript IDs
```
$ cat introns_transcripts.txt | awk '{print $0"."NR"-I"}' > introns_transcripts.to_capture.txt
```
Next we have to map the transcripts to their respective genes.
```
$ awk 'NR==FNR{a[$1]=$2; b[$1]=$3;next} {$2=a[$1];$3=b[$1]} 1' tr2g.txt introns_transcripts.txt > introns_t2g.txt
```

#### 4c. Fix the INTRONS FASTA header
We need to fix all of the headers for the introns FASTA file so that they contain the transcript ID, an identifier specifying that the transcript is an "intronic" transcript, and a unique number to avoid duplicates. 
```
$ awk '{print ">"$1"."NR"-I"" gene_id:"$2" gene_name:"$3}' introns_t2g.txt > introns_fasta_header.txt
$ awk -v var=1 'FNR==NR{a[NR]=$0;next}{ if ($0~/^>/) {print a[var], var++} else {print $0}}' introns_fasta_header.txt introns.fa > introns.correct_header.fa
$ head -1 introns.correct_header.fa
>ENST00000456328.2.1-I gene_id:ENSG00000223972.5 gene_name:DDX11L1 1
```

#### 5a. Get the transcripts to capture list and transcripts to genes for cDNA
Gunzip the cDNA FASTA file
```
$ gunzip cDNA.fa.gz
$ head -1 cDNA.fa  ## Homo sapiens GRCh38
>hg38_knownGene_ENST00000456328.2_0 range=chr1:11869-12227 5'pad=0 3'pad=0 strand=+ repeatMasking=none ## GRCh38 example
```
Get list of transcripts
```
$ cat cDNA.fa | awk '/^>/ {print $0}' | tr "_" " " | awk '{print $3}' > cDNA_transcripts.txt
$ cat cDNA_transcripts.txt | tr "." " " | awk '{print $1}' > cDNA_transcripts_no_version.txt
```
**Note:** your FASTA header may not be structured exactly like the example. If that is the case, you can change the columns that are printed in the `awk` command. For example, `awk '{print $3}'` prints the 3rd column.

**Note:** its important that the transcript IDs be consistent with the `tr2g.txt` that you made earlier, i.e. keep both with or without version.

Add an identifier to thr transcript IDs
```
$ cat cDNA_transcripts.txt | awk '{print $0"."NR}' > cDNA_transcripts.to_capture.txt
```
Map the transcripts to genes
```
$ awk 'NR==FNR{a[$1]=$2; b[$1]=$3;next} {$2=a[$1];$3=b[$1]} 1' tr2g.txt cDNA_transcripts.txt > cDNA_t2g.txt
```

#### 5b. Fix the INTRONS FASTA header
```
$ awk '{print ">"$1"."NR" gene_id:"$2" gene_name:"$3}' cDNA_t2g.txt > cDNA_fasta_header.txt
$ awk -v var=1 'FNR==NR{a[NR]=$0;next}{ if ($0~/^>/) {print a[var], var++} else {print $0}}' cDNA_fasta_header.txt $cDNA_fa > cDNA.correct_header.fa
$ head -1 cDNA.correct_header.fa
>ENST00000456328.2.1 gene_id:ENSG00000223972.5 gene_name:DDX11L1 1
```

#### 6. Make the kallisto index
```
$ cat cDNA.correct_header.fa introns.correct_header.fa > cDNA_introns.fa
$ cat cDNA_t2g.txt introns_t2g.txt > cDNA_introns_t2g.txt
$ kallisto index -i cDNA_introns.idx -k 31 cDNA_introns.fa
```

#### 7. Align your reads
See [this](velocity_tutorial.html) tutorial for how to proceed.
