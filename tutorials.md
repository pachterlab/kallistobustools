---
layout: page
title: "Tutorials"
group: navigation
---

{% include JB/setup %}

**Note:** All Google Colab notebooks can be run by selecting `Runtime > Run all > Run anyway` within the notebook.

---

| Tutorial | Description | Languages|
|:-----|:------------|---------:|
|Getting started | Download a mouse transcriptome index, pseudoalign reads, and perform downstream analysis | [`python`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_standard.ipynb)  |
|Build transcriptome index | Build a transcriptome index from a Genome and Genome annotation file | [`python`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_transcriptome_index.ipynb) |
| Downloading data | Learn how to download data from different sequencing databases | [`command line`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/data_download.ipynb)|
| Aggregating multiple count matrices | Combine two count matrices for downstream analysis | [`python`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_aggregating_count_matrices.ipynb)|
| Multiple FASTQs | Process multiple sets of FASTQs to create gene count matrices | [`python`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_multiple_files.ipynb) |
| Multi-species experiments | Create a combine human mouse index and align single cell reads from a human mouse experiment | [`python`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_species_mixing.ipynb), [`R`](https://bustools.github.io/BUS_notebooks_R/10xv2.html)|
| Feature barcoding | Align feature barcoding reads from a 10x Feature Barcoding experiment | [`python`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_kite.ipynb)|
| Parsing bus files | Parse a bus file | [`R`](https://bustools.github.io/BUS_notebooks_R/10xv3.html) |
| Pre-processing single-nuclei RNA-seq | Create a combine transciptome - intron index and process single-nuclei data | [`python`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_single_nucleus.ipynb)|
| Pseudotime | Perform pseudotime to understand differentiation trajectories of cells | [`R - Monocole 2`](https://bustools.github.io/BUS_notebooks_R/monocle2.html), [`R - Monocole 3`](https://bustools.github.io/BUS_notebooks_R/monocle3.html), [`python - Monocole 3`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_monocle.ipynb), [`R - Slingshot`](https://bustools.github.io/BUS_notebooks_R/slingshot.html)|
| RNA velocity | Perform RNA velocity on the La Manno Human Forebrain data to create TSNE RNA velocity vectors | [`python - premade index`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_velocity.ipynb), [`R - custom index`](https://bustools.github.io/BUS_notebooks_R/velocity.html)|
| Custom Index for RNA velocity | Construct a custom transcriptome - intron index for RNA velocity | [`python`](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_velocity_index.ipynb)|
