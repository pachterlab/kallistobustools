---
layout: page
title: "About"
group: navigation
---

{% include JB/setup %}


The analysis of single-cell RNA-Seq data involves a series of steps that include: (1) pre-processing of reads to associate them with their cells of origin, (2) possible collapsing of reads according to unique molecular identifiers (UMIs), (3) generation of feature counts from the reads to generate a feature-cell matrix and (4) analysis of the matrix to compare and contrast cells.

We have recently introduced the __BUS__ file format for single-cell RNA-seq data designed to facilitate the development of modular workflows for data processing. It consists of a binary representation of barcode and UMI sequences from scRNA-seq reads, along with sets of equivalence classes obtained by pseudoalignment of reads to a reference transcriptome (hence the acronym Barcode, UMI, Set). We have implemented a command in __kallisto__ version [0.45.0](http://pachterlab.github.io/kallisto//releases/2018/11/17/v0.45.0) called `bus` that allows for the efficient generation of BUS format from any single-cell RNA-seq technology. Tools for manipulating BUS files are provided as part of the [__bustools__](https://bustools.github.io/) package. 

This website provides tutorials and workflows to learn how to use the __kallisto__ and __BUStools__ programs to perform single-cell RNA-seq analysis. We suggested beginning with the [Getting Started](getting_started.html) page.
