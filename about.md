---
layout: page
title: "About"
group: navigation
---

{% include JB/setup %}

<p align="center">
<img src="https://user-images.githubusercontent.com/10369156/58990086-06c81500-879a-11e9-886b-7e4a690c5862.png" width="70%" />
</p>


The analysis of single-cell RNA-Seq data involves a series of pre-processing steps that include: (1) association of reads with their cells of origin, (2) collapsing of reads according to unique molecular identifiers (UMIs), and (3) generation of feature counts from the reads to generate a feature-cell matrix.

We recently introduced the [__BUS__ file format](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz279/5487510) for single-cell RNA-seq data to facilitate the development of modular workflows for data pre-processing. It consists of a binary representation of barcode and UMI sequences from scRNA-seq reads, along with sets of equivalence classes obtained by pseudoalignment of reads to a reference transcriptome (hence the acronym Barcode, UMI, Set). We have implemented a command in __kallisto__ called `bus` that allows for the efficient generation of BUS format from any single-cell RNA-seq technology. Tools for manipulating BUS files are provided as part of the [__bustools__](https://bustools.github.io/) package.

This website provides tutorials and workflows to learn how to use the __kallisto__ and __bustools__ programs together to perform single-cell RNA-seq pre-processing. We suggest beginning with the [Getting Started](kb_getting_started.html) page. See the [__kallisto and applications__ Google group](https://groups.google.com/forum/#!forum/kallisto-and-applications) for answers to frequently asked questions. The kallisto &#124; bustools workflow is described in detail in

Páll Melsted, A. Sina Booeshaghi, Fan Gao, Eduardo Beltrame, Lambda Lu, Kristján Eldjárn Hjorleifsson, Jase Gehring and Lior Pachter, [Modular and efficient pre-processing of single-cell RNA-seq](https://www.biorxiv.org/content/10.1101/673285v2), bioRxiv, 2019.
