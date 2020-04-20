---
layout: page
title: "Computational analysis of sequencing based SARS-CoV-2 testing data"
---

## Overview

Several groups have been developing large scale testing for SARS-CoV-2 based on high-throughput sequencing. These methods have the potential to massively scale the extent of testing at low-cost and with assays that are highly sensitive.

The sequence data output of the proposed methods consists of two data types:

- (sequence) barcodes, each of which is associated to a sample.
- (biological) reads associated with sequences from genes. Genes may include viral genes, control genes, or spike-ins.

Pre-processing of the data involves several steps:

- Error correction of the barcode sequences. This is necessary in order to account for sequence errors that may have been introduced during sequencing.
- Association of the biological reads to the genes of origin.
- Collation of reads associated to a single sample in order to count the number of times different genes have been observed.

Following pre-processing of the data, analysis must be performed to determine which samples contained SARS-CoV-2, and which didn't. This involves setting thresholds based on the observed viral counts in the context of counts of spike-ins and controls.



{% include JB/setup %}


