---
layout: page
title: "Data downloading"
---

{% include JB/setup %}

This tutorial provides information on where to find single-cell RNA-seq data, and how to download it for processing with the __kallisto &#124; bustools__ workflow.

__Note:__ for the instructions, command line arguments are preceeded by`$`. For example, if you see `$ cd my_folder` then type `cd my_folder`. 

#### Databases

There are four databases that are important repositories for sequencing data and metadata, and that are relevant for obtaining single-cell RNA-seq data. For each archive we provide an example of how the data is organized and how to download it.

- [__Biological Project Library__](https://bigd.big.ac.cn/bioproject/) (BioProject): The Biological Project Library organizes metadata for research projects involving genomic data types. This repository, which was started in 2016, is similar to the Gene Expression Omnibus. As an example, the data from the paper [Peng et al. 2019](https://www.nature.com/articles/s41422-019-0195-y) is organized under project accession [PRJCA001063](https://bigd.big.ac.cn/bioproject/browse/PRJCA001063). Each single-cell RNA-seq dataset has a "BioSample accession", e.g. [SAMC047103](https://bigd.big.ac.cn/biosample/browse/SAMC047103). A further link to the Genome Sequencing Archive provides access to FASTQ files.

- [__Genome Sequence Archive__](http://gsa.big.ac.cn/) (GSA): This repository contains reads for projects in FASTQ format. For example, reads for [SAMC047103](https://bigd.big.ac.cn/biosample/browse/SAMC047103) from the [PRJCA001063](https://bigd.big.ac.cn/bioproject/browse/PRJCA001063) in the BioProject repository are accessible under accession [CRA001160](https://bigd.big.ac.cn/gsa/browse/CRA001160). A specific run accession, e.g. [CRR034516](https://bigd.big.ac.cn/gsa/browse/CRA001160/CRR034516) provides direct access to FASTQ files.

- [__Gene Expression Omnibus__](https://www.ncbi.nlm.nih.gov/geo/) (GEO): The Gene Expression Omnibus is a repository for [MIAME (Minimum Infomration about a Microarray Experiment)](https://www.ncbi.nlm.nih.gov/geo/info/MIAME.html) compliant data. While the MIAME standards were established during a time when gene expression data was primarily collected with microarrays, the standards also apply to sequencing data and the GEO repository hosts project metadata for both types of research projects. As an example, the project link for the paper [Wolock et al. 2019](https://www.sciencedirect.com/science/article/pii/S2211124719307971) is [GSE132151](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132151). Most papers refer to their data via GEO accessions, so GEO is a useful repository for searching for data from projects.

- [__European Nucelotide Archive__](https://www.ebi.ac.uk/ena) (ENA): The ENA provides access to nucleotide sequences associated with genomic projects. In the case of [GSE132151](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132151) mentioned above, the nucleotide sequences are at [PRJNA546231](https://www.ebi.ac.uk/ena/data/view/PRJNA546231). The ENA provides direct access to FASTQ files from the project page. It also links to NCBI Sequence Read Archive format data.

- [__Sequence Read Archive__](https://www.ncbi.nlm.nih.gov/sra) (SRA):


#### Streaming
