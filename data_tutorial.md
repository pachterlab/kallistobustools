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

- [__Sequence Read Archive__](https://www.ncbi.nlm.nih.gov/sra) (SRA): The SRA is a sequence repository for genomic data. Files are stored in SRA format, which must be downloaded and converted to FASTQ format prior to pre-processing using the `fasterq-dump` program available as part of [SRA tools](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump). For example, the data in [Rossi et al., 2019](https://science.sciencemag.org/content/364/6447/1271) can be located in the SRA via [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130597), then to the [SRA](https://www.ncbi.nlm.nih.gov/sra?term=SRP194426), and finally a sequence data page for one of the runs, [SRX5779290](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR9000493) has information about the traces (reads). The SRA tools operate directly on SRA accessions.

#### Searching

The [__sra-explorer__](https://ewels.github.io/sra-explorer/) website is an effective and easy to use utility for searching the SRA and for downloading files. The utility finds SRA entires by keywords or accession numbers and produces links to the FASTQs and to commands for downloading them. 

#### Streaming

Single-cell RNA-seq data from sequence repositories can be streamed into __kallisto__ making possible a workflow that does not require saving files to disk prior to pre-processing. This can be done via process substitution or using `df`. For example, the following command can be used to stream data from the ENA for pre-processing:

`$ urlR1="https://github.com/bustools/getting_started/releases/download/getting_started/SRR8599150_S1_L001_R1_001.fastq.gz";`
`$ urlR2="https://github.com/bustools/getting_started/releases/download/getting_started/SRR8599150_S1_L001_R2_001.fastq.gz";`
`$ time kallisto bus -i Mus_musculus.GRCm38.cdna.all.idx -x 10xv2 -t 4 -o bus_out/ <(curl -Ls ${urlR1}) <(curl -Ls ${urlR2})`

The only required file that must be locally stored on disk prior to pre-processing is the transcriptome index [Mus_musculus.GRCm38.cdna.all.idx](https://github.com/pachterlab/kallisto-transcriptome-indices/releases). A complete tutorial for how to stream data, together with code based on `mkfifo` is available [here](https://sinabooeshaghi.com/2019/07/09/fasterq-to-count-matrices-for-single-cell-rna-seq/).
