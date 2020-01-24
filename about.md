---
layout: page
title: "About"
group: navigation
---

{% include JB/setup %}
kallisto bustools is a workflow for pre-processing single-cell RNA-seq data. 

Pre-processing single-cell rNA-seq involves: (1) association of reads with their cells of origin, (2) collapsing of reads according to unique molecular identifiers (UMIs), and (3) generation of gene or feature counts from the reads to generate a cell x gene matrix.

With kallisto | bustools you can 
* Generate a _cell x gene_ or _cell x transcript equivalence class_ matrix
* Perform RNA velocity and single-nuclei RNA-seq
* Quantify data from numerous technologies such as 10x, inDrops, and Dropseq.
* Customize the workflow for novel technologies and protocols.
* Process feature barcoding data such as CITE-seq, REAP-seq, MULTI-seq, Clicktags, and Perturb-seq.
* Obtain QC reports for single-cell libraries.

Click on the introduction tab to learn more or the tutorials tab for examples.

The kallisto bustools workflow is described in:

Páll Melsted, A. Sina Booeshaghi, Fan Gao, Eduardo Beltrame, Lambda Lu, Kristján Eldjárn Hjorleifsson, Jase Gehring and Lior Pachter, [Modular and efficient pre-processing of single-cell RNA-seq](https://www.biorxiv.org/content/10.1101/673285v1), bioRxiv, 2019.
