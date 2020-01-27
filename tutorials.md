---
layout: page
title: "Tutorials"
group: navigation
---

{% include JB/setup %}

**Note:** All Google Colab notebooks can be run by selecting `Runtime > Run all > Run anyway` within the notebook.

---

#### Getting started [`python`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_standard.ipynb) 
- Download a mouse transcriptome index, pseudoalign reads, and perform downstream analysis

#### Build transcriptome index [`python`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_transcriptome_index.ipynb)
- Build a transcriptome index from a Genome and Genome annotation file

#### Downloading data [`command line`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/data_download.ipynb)
- Learn how to download data from different sequencing databases

#### Aggregating multiple count matrices [`python`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_aggregating_count_matrices.ipynb)
- Combine two count matrices for downstream analysis

#### Multiple FASTQs [`python`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_multiple_files.ipynb)
- Process multiple sets of FASTQs to create gene count matrices

#### Multi-species experiments [`python`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_species_mixing.ipynb), [`R`](https://bustools.github.io/BUS_notebooks_R/10xv2.html)
- Create a combine human mouse index and align single cell reads from a human mouse experiment

#### Feature barcoding [`python`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_kite.ipynb)
- Align feature barcoding reads from a 10x Feature Barcoding experiment

#### Parsing bus files [`R`](https://bustools.github.io/BUS_notebooks_R/10xv3.html)
- Parse a bus file 

#### Pre-processing single-nuclei RNA-seq [`python`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_single_nucleus.ipynb)
- Create a combine transciptome - intron index and process single-nuclei data

#### Pseudotime [`R - Monocole 2`](https://bustools.github.io/BUS_notebooks_R/monocle2.html), [`R - Monocole 3`](https://bustools.github.io/BUS_notebooks_R/monocle3.html), [`python - Monocole 3`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_monocle.ipynb), [`R - Slingshot`](https://bustools.github.io/BUS_notebooks_R/slingshot.html)
- Perform pseudotime to understand differentiation trajectories of cells

#### RNA velocity [`python - premade index`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_velocity.ipynb), [`R - custom index`](https://bustools.github.io/BUS_notebooks_R/velocity.html)
- Perform RNA velocity on the La Manno Human Forebrain data to create TSNE RNA velocity vectors

#### Custom Index for RNA velocity [`python`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_velocity_index.ipynb)
- Construct a custom transcriptome - intron index for RNA velocity
