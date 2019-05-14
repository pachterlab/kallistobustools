---
layout: page
title: "About"
group: navigation
---

{% include JB/setup %}

__BUS__ format is a file format for single-cell RNA-seq data designed to facilitate the development of modular workflows for data processing. It consists of a binary representation of barcode and UMI sequences from scRNA-seq reads, along with sets of equivalence classes obtained by pseudoalignment of reads to a reference transcriptome (hence the acronym Barcode, UMI, Set). __BUS__ files are a convenient and useful checkpoint during single-cell RNA-seq processing. The format is described in detail in the [__BUStools__ __BUS__ format repository](https://github.com/BUStools/BUS-format) and in the preprint

[P. Melsted, V. Ntranos and L. Pachter, The Barcode, UMI, Set format and BUStools, bioRxiv, 2018](https://www.biorxiv.org/content/10.1101/472571v2).

We have implemented a new [feature in __kallisto__](https://pachterlab.github.io/kallisto/singlecell.html) that can be used to generate __BUS__ format files from a variety of different technologies. Once __BUS__ files have been produced, they can be manipulated with [__bustools__](https://github.com/BUStools/bustools). The __bustools__ programs can be used to perform tasks such as sorting __BUS__ files, which makes them easier to process into meaningful single-cell information in higher-level languages such as R or python. The end results are easy to develop, fast, and accurate workflows for a variety of single-cell RNA-seq and multimodal technologies.

To learn how to use the __BUStools__ suite of programs to perform single-cell RNA-seq analysis we suggested beginning with the [Getting Started](getting_started.html) page.