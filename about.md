---
layout: page
title: "About"
group: navigation
---

{% include JB/setup %}

The analysis of single-cell RNA-Seq data involves a series of pre-processing steps that include: (1) association of reads with their cells of origin, (2) collapsing of reads according to unique molecular identifiers (UMIs), and (3) generation of gene or feature counts from the reads to generate a cell x gene matrix. We have developed two tools and one file format to acheive fast, efficient, and accurate single-cell preprocessing: [kallisto](https://www.nature.com/articles/nbt.3519) for pseudoalignment of single-cell reads, [bustools](https://www.biorxiv.org/content/10.1101/673285v2) for manipulating BUS files and the [BUS file format](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz279/5487510).

With kallisto | bustools you can 
* Generate a cell x gene or cell x transcript equivalence class matrix
* Perform RNA velocity and Single Nuclei RNA-seq
* Process numerous technologies such as 10x, inDrops, Dropseq 
* Process feature barcoding data such as Cite-seq, Reap-seq, Multi-seq, Perturb-seq, Clicktags
* Get rapid QC about your single-cell library

Click on the introduction tab to learn more or the tutorials tab for examples.

Páll Melsted, A. Sina Booeshaghi, Fan Gao, Eduardo Beltrame, Lambda Lu, Kristján Eldjárn Hjorleifsson, Jase Gehring and Lior Pachter, [Modular and efficient pre-processing of single-cell RNA-seq](https://www.biorxiv.org/content/10.1101/673285v1), bioRxiv, 2019.
