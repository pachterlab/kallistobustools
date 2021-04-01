<a href="https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/tutorials/docs/tutorials/kb_getting_started/python/kb_intro_2_python.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# Introduction to single-cell RNA-seq II: getting started with analysis

This notebook demonstrates pre-processing and basic analysis of the [mouse retinal cells GSE126783](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126783) dataset from [Koren et al., 2019](https://doi.org/10.1016/j.immuni.2019.02.007). Following pre-processing using kallisto and bustools and basic QC, the notebook demonstrates some initial analysis. The approximate running time of the notebook is about 13 minutes.

The notebook was written by Kyung Hoi (Joseph) Min, A. Sina Booeshaghi and Lior Pachter. If you use the methods in this notebook for your analysis please cite the following publications which describe the tools used in the notebook, as well as specific methods they run (these are cited inline in the notebook):

* Melsted, P., Booeshaghi, A.S. et al. Modular and efficient pre-processing of single-cell RNA-seq. bioRxiv (2019). doi:10.1101/673285
* Wolf, F. A., Angere, P. and Theis, F.J. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology (2018). doi:10.1186/s13059-017-1382-0

An R notebook implementing the same analysis is available [here](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_intro_2_R.ipynb). See the [kallistobus.tools tutorials](https://www.kallistobus.tools/tutorials) site for additional notebooks demonstrating other analyses.


## Setup


```
# This is  used to time the running of the notebook
import time
start_time = time.time()
```

### Install python packages


```
# These packages are pre-installed on Google Colab, but are included here to facilitate running this notebook locally
!pip install --quiet matplotlib
!pip install --quiet scikit-learn
!pip install --quiet numpy
!pip install --quiet scipy
```


```
%%time
# `kb` is a wrapper for the kallisto and bustools program, and the kb-python package contains the kallisto and bustools executables.
!pip install --quiet kb-python==0.24.1
```

    [K     |████████████████████████████████| 35.4MB 85kB/s 
    [K     |████████████████████████████████| 51kB 7.5MB/s 
    [K     |████████████████████████████████| 122kB 53.5MB/s 
    [K     |████████████████████████████████| 112kB 55.5MB/s 
    [?25h  Building wheel for loompy (setup.py) ... [?25l[?25hdone
      Building wheel for numpy-groupies (setup.py) ... [?25l[?25hdone
    CPU times: user 102 ms, sys: 14.3 ms, total: 117 ms
    Wall time: 9.38 s



```
%%time
# Install scanpy and other packages needed for single-cell RNA-seq analysis
!pip install --quiet scanpy python-igraph louvain MulticoreTSNE pybiomart
```

    [K     |████████████████████████████████| 10.3MB 17.3MB/s 
    [K     |████████████████████████████████| 3.2MB 48.3MB/s 
    [K     |████████████████████████████████| 2.2MB 52.2MB/s 
    [K     |████████████████████████████████| 81kB 11.6MB/s 
    [K     |████████████████████████████████| 1.2MB 33.4MB/s 
    [K     |████████████████████████████████| 51kB 7.4MB/s 
    [K     |████████████████████████████████| 71kB 10.0MB/s 
    [?25h  Building wheel for MulticoreTSNE (setup.py) ... [?25l[?25hdone
      Building wheel for umap-learn (setup.py) ... [?25l[?25hdone
      Building wheel for sinfo (setup.py) ... [?25l[?25hdone
      Building wheel for pynndescent (setup.py) ... [?25l[?25hdone
    CPU times: user 108 ms, sys: 17.5 ms, total: 125 ms
    Wall time: 15.1 s



```
# Import packages
import anndata
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import TruncatedSVD
from scipy import sparse, io

matplotlib.rcParams.update({'font.size': 12})
%config InlineBackend.figure_format = 'retina'
```

### Download the data

__Note:__ We use the `-O` option for `wget` to rename the files so they can be easily identified. The notebook requires reads in fastq format; the files can be processed in gzip compressed format.

In this example the reads are downloaded from a Box drive; see the [data download](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/data_download.ipynb) notebook for information on where to find publicly available single-cell RNA-seq data.


```
%%time
!wget -q https://caltech.box.com/shared/static/9u2qk1uqu6py03phooe1ti0kjd9v87pu.txt -O checksums.txt
!wget -q https://caltech.box.com/shared/static/w9ww8et5o029s2e3usjzpbq8lpot29rh.gz -O SRR8599150_S1_L001_R1_001.fastq.gz
!wget -q https://caltech.box.com/shared/static/ql00zyvqnpy7bf8ogdoe9zfy907guzy9.gz -O SRR8599150_S1_L001_R2_001.fastq.gz
```

    CPU times: user 163 ms, sys: 56.6 ms, total: 220 ms
    Wall time: 49.7 s




```
# This is formatted as code
```

Next, we verify the integrity of the files that were downloaded to confirm that they were not corrupted during the download.


```
!md5sum -c checksums.txt --ignore-missing
```

    SRR8599150_S1_L001_R1_001.fastq.gz: OK
    SRR8599150_S1_L001_R2_001.fastq.gz: OK


### Download an index

__Note:__ See [this notebook](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_transcriptome_index.ipynb) for a tutorial on how to build custom transcriptome indices or indices for RNA velocity.


```
%%time
# This downloads a pre-built index for kallisto to use when pseudoaligning the reads
!kb ref -d mouse -i index.idx -g t2g.txt
```

    [2021-03-31 21:50:23,293]    INFO Downloading files for mouse from https://caltech.box.com/shared/static/vcaz6cujop0xuapdmz0pplp3aoqc41si.gz to tmp/vcaz6cujop0xuapdmz0pplp3aoqc41si.gz
    [2021-03-31 21:52:03,259]    INFO Extracting files from tmp/vcaz6cujop0xuapdmz0pplp3aoqc41si.gz
    CPU times: user 449 ms, sys: 89.4 ms, total: 539 ms
    Wall time: 2min 16s


## Pseudoalignment and counting

### Run kallisto and bustools

The following command will generate an RNA count matrix of cells (rows) by genes (columns) in H5AD format, which is a binary format used to store [Anndata](https://anndata.readthedocs.io/en/stable/) objects. Notice that this requires providing the index and transcript-to-gene mapping downloaded in the previous step to the `-i` and `-g` arguments respectively. Also, since the reads were generated with the 10x Genomics Chromium Single Cell v2 Chemistry, the `-x 10xv2` argument is used. To view other supported technologies, run `kb --list`.

__Note:__ To output a [Loom](https://linnarssonlab.org/loompy/format/index.html) file instead, replace the `--h5ad` flag with `--loom`. To obtain the raw matrix output by `kb` instead of the H5AD or Loom converted files, omit these flags.


```
%%time
# This step runs `kb` to pseudoalign the reads, and then generate the cells x gene matrix in h5ad format.
!kb count -i index.idx -g t2g.txt -x 10xv2 --h5ad -t 2 \
SRR8599150_S1_L001_R1_001.fastq.gz SRR8599150_S1_L001_R2_001.fastq.gz
```

    [2021-03-31 21:52:45,896]    INFO Generating BUS file from
    [2021-03-31 21:52:45,896]    INFO         SRR8599150_S1_L001_R1_001.fastq.gz
    [2021-03-31 21:52:45,897]    INFO         SRR8599150_S1_L001_R2_001.fastq.gz
    [2021-03-31 21:55:04,845]    INFO Sorting BUS file ./output.bus to tmp/output.s.bus
    [2021-03-31 21:55:08,125]    INFO Whitelist not provided
    [2021-03-31 21:55:08,125]    INFO Copying pre-packaged 10XV2 whitelist to .
    [2021-03-31 21:55:08,226]    INFO Inspecting BUS file tmp/output.s.bus
    [2021-03-31 21:55:09,648]    INFO Correcting BUS records in tmp/output.s.bus to tmp/output.s.c.bus with whitelist ./10xv2_whitelist.txt
    [2021-03-31 21:55:24,007]    INFO Sorting BUS file tmp/output.s.c.bus to ./output.unfiltered.bus
    [2021-03-31 21:55:26,998]    INFO Generating count matrix ./counts_unfiltered/cells_x_genes from BUS file ./output.unfiltered.bus
    [2021-03-31 21:55:29,092]    INFO Converting matrix ./counts_unfiltered/cells_x_genes.mtx to h5ad ./counts_unfiltered/adata.h5ad
    CPU times: user 852 ms, sys: 103 ms, total: 955 ms
    Wall time: 2min 52s


### Exercise

- [Loom](https://linnarssonlab.org/loompy/format/index.html) is a alternative format for storing single-cell count matrices. Output a Loom file with kb by replacing the `--h5ad` flag with `--loom`, or obtain the raw matrix output by omitting the flags


```
%%time
# # This runs `kb` to pseudoalign the reads, and then generate the cells x gene matrix in Loom format.
# !kb count -i index.idx -g t2g.txt -x 10xv2 --h5ad -t 2 \
# SRR8599150_S1_L001_R1_001.fastq.gz SRR8599150_S1_L001_R2_001.fastq.gz
```

    CPU times: user 2 µs, sys: 1 µs, total: 3 µs
    Wall time: 4.53 µs


## Basic QC

First, the *cells x genes* matrix is imported into an H5AD-formatted [Anndata](https://icb-anndata.readthedocs-hosted.com/en/stable/anndata.AnnData.html) matrix. Anndata is a convenient format for storing single-cell count matrices in Python.


```
# import data
adata = anndata.read('counts_unfiltered/adata.h5ad')
adata
```




    AnnData object with n_obs × n_vars = 96775 × 55421



Represent the cells in 2D with PCA


```
# Perform SVD
tsvd = TruncatedSVD(n_components=2)
tsvd.fit(adata.X)
X = tsvd.transform(adata.X)

# Plot the cells in the 2D PCA projection
fig, ax = plt.subplots(figsize=(10, 7))

ax.scatter(X[:,0], X[:,1], alpha=0.5, c="green")

plt.axis('off')
plt.show()
```


![png](kb_intro_2_python_files/kb_intro_2_python_25_0.png)


### Test for library saturation


```
# Create a plot showing genes detected as a function of UMI counts.
fig, ax = plt.subplots(figsize=(10, 7))

x = np.asarray(adata.X.sum(axis=1))[:,0]
y = np.asarray(np.sum(adata.X>0, axis=1))[:,0]

ax.scatter(x, y, color="green", alpha=0.25)
ax.set_xlabel("UMI Counts")
ax.set_ylabel("Genes Detected")
ax.set_xscale('log')
ax.set_yscale('log', nonposy='clip')

ax.set_xlim((0.5, 4500))
ax.set_ylim((0.5,2000))


plt.show()
```


![png](kb_intro_2_python_files/kb_intro_2_python_27_0.png)


### Examine the knee plot

The "knee plot" was introduced in the Drop-seq paper: 
- Macosko et al., [Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets](https://www.cell.com/fulltext/S0092-8674(15)00549-8), 2015. DOI:10.1016/j.cell.2015.05.002

In this plot cells are ordered by the number of UMI counts associated to them (shown on the *x*-axis), and the fraction of droplets with at least that number of cells is shown on the *y*-axis:


```
#@title Threshold cells according to knee plot { run: "auto", vertical-output: true }
expected_num_cells =  3655#@param {type:"integer"}
knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]

fig, ax = plt.subplots(figsize=(10, 7))

ax.loglog(knee, range(len(knee)), linewidth=5, color="g")
ax.axvline(x=knee[expected_num_cells], linewidth=3, color="k")
ax.axhline(y=expected_num_cells, linewidth=3, color="k")

ax.set_xlabel("UMI Counts")
ax.set_ylabel("Set of Barcodes")

plt.grid(True, which="both")
plt.show()
```


![png](kb_intro_2_python_files/kb_intro_2_python_29_0.png)


The knee plot can be used to threshold cells based on the number of UMI counts they contain. The threshold is applied at the "knee", where there is a sharp dropoff in the number of UMIs per cell. In this example we use the nunber 3979 based on the publication describing the data.

### Filter empty droplets


```
adata
```




    AnnData object with n_obs × n_vars = 96775 × 55421




```
# Filter the cells according to the threshold determined from the knee plot
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, min_counts=knee[expected_num_cells])
```


```
adata
```




    AnnData object with n_obs × n_vars = 3536 × 55421
        obs: 'n_genes', 'n_counts'



### Filtering out by mitochondrial content


```
mito_ensembl_ids = sc.queries.mitochondrial_genes("mmusculus", attrname="ensembl_gene_id")
```


```
mito_genes = mito_ensembl_ids["ensembl_gene_id"].values
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1
```


```
sc.pl.scatter(adata, x='n_counts', y='percent_mito')
```


![png](kb_intro_2_python_files/kb_intro_2_python_38_0.png)



```
adata = adata[adata.obs.percent_mito < 0.03]
```

### Filter out genes that are not present in any cells


```
sc.pp.filter_genes(adata, min_cells=3)
```

    Trying to set attribute `.var` of view, copying.



```
adata
```




    AnnData object with n_obs × n_vars = 3507 × 15420
        obs: 'n_genes', 'n_counts', 'percent_mito'
        var: 'n_cells'



### Visualizing count distributions

Examination of the gene count and UMI count distributions is useful QC to evaluate the quality of the library and how deeply it was sequenced.


```
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True)
```


![png](kb_intro_2_python_files/kb_intro_2_python_45_0.png)


## Analysis

In this part of the tutorial, the cells are clustered by [Louvain community detection](https://en.wikipedia.org/wiki/Louvain_modularity).

### Normalize the counts


```
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
```

### Identify highly variable genes


```
# flavor="cell_ranger" is consistent with Seurat and flavor="suerat" is not consistent with Seurat
sc.pp.highly_variable_genes(adata, min_mean=0.01, max_mean=8, min_disp=1, n_top_genes=3000, n_bins=20, flavor="seurat")
sc.pl.highly_variable_genes(adata)
```


![png](kb_intro_2_python_files/kb_intro_2_python_50_0.png)



```
sc.pp.scale(adata, max_value=10)
```

### Clustering and visualization

There are many algorithms for clustering cells, and while they have been compared in detail in various benchmarks (see e.g., [Duo et al. 2018](https://f1000research.com/articles/7-1141/v2)), there is no univerally agreed upon method. Here we demonstrate clustering using [Louvain clustering](https://en.wikipedia.org/wiki/Louvain_modularity), which is a popular method for clustering single-cell RNA-seq data. The method was published in 

- Blondel, Vincent D; Guillaume, Jean-Loup; Lambiotte, Renaud; Lefebvre, Etienne (9 October 2008). "Fast unfolding of communities in large networks". Journal of Statistical Mechanics: Theory and Experiment. 2008 (10): P10008.


```
%%time
# Cluster the cells using Louvain clustering
sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True, n_comps=10)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=10, knn=True)
sc.tl.louvain(adata)
```

    CPU times: user 18 s, sys: 738 ms, total: 18.7 s
    Wall time: 18.5 s


It is useful to revisit the PCA project with points colored by cluster. Previously we computed the PCA projection directly; here we use a function in ScanPy which does the same.

### PCA


```
# Perform PCA and plot the projection to the first two dimensions, with points colored according to the Louvain clusters.
fig, ax = plt.subplots(figsize=(10, 7))
sc.pl.pca(adata, color='louvain', ax=ax)
```


![png](kb_intro_2_python_files/kb_intro_2_python_57_0.png)


The PCA representation is the result a *linear* map of the data from its ambient dimension, to low dimension (in the case above 2D). While such projections are useful, there are *non-linear* methods that can capture non-linear geometry in the data.

### tSNE

[t-SNE](https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding) is a non-linear dimensionality reduction technique described in:

- Maaten, Laurens van der, and Geoffrey Hinton. "Visualizing data using t-SNE." Journal of machine learning research 9.Nov (2008): 2579-2605.

Here it is applied to the 10-dimensional PCA projection of the cells.



```
%%time
# Visualize cells with t-SNE. The n_pcs parameter sets the number of principal components to project to prior to 
# performing t-SNE
sc.tl.tsne(adata, n_pcs=10)
```

    CPU times: user 19.7 s, sys: 7 ms, total: 19.7 s
    Wall time: 19.7 s



```
fig, ax = plt.subplots(figsize=(10, 7))
sc.pl.tsne(adata, color='louvain', ax=ax)
```


![png](kb_intro_2_python_files/kb_intro_2_python_61_0.png)


## UMAP

UMAP stands for Uniform Manifold Approximation and Projection is a non-linear dimensionality reduction techinque described in 

* Leland McInnes and John Healy and James Melville, "UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction" 2018 1802.03426 arXiv stat.ML



```
%%time
# Visualize cells with t-SNE. The n_pcs parameter sets the number of principal components to project to prior to 
# performing t-SNE
sc.tl.umap(adata)
```

    CPU times: user 11.6 s, sys: 185 ms, total: 11.8 s
    Wall time: 11.7 s



```
fig, ax = plt.subplots(figsize=(10, 7))
sc.pl.umap(adata, color='louvain', ax=ax)
```


![png](kb_intro_2_python_files/kb_intro_2_python_64_0.png)


### Exercises

- The variance explained by each principal component is a measure of how well a projection to that component represents the data. Compute the variance explained by each component.


```
# Compute and plot the variance explained by the PC subspaces.
# sc.pl.pca_variance_ratio(adata)
```

- In the notebook we visualized the data using t-SNE and UMAP. These two techniques can produce different two dimensional embeddings based with different parameters and even different random seeds. Change around these parameters to see how they affect the embedding.

For tSNE parameters see https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.tl.tsne.html

For UMAP parameters see https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.tl.umap.html


```
# %%time
# # Visualize cells with UMAP.
# sc.tl.tsne(adata, 
# n_pcs=None, 
# use_rep=None, 
# perplexity=30, 
# early_exaggeration=12, 
# learning_rate=1000, 
# random_state=0, 
# use_fast_tsne=True, 
# n_jobs=None, 
# copy=False)
# fig, ax = plt.subplots(figsize=(10, 7))
# sc.pl.tsne(adata, color='louvain', ax=ax)
```


```
# %%time
# # Visualize cells with UMAP.
# sc.tl.umap(adata, 
# min_dist=0.5, 
# spread=1.0, 
# n_components=2, 
# maxiter=None, 
# alpha=1.0, 
# gamma=1.0, 
# negative_sample_rate=5, 
# init_pos='spectral', 
# random_state=0, 
# a=None, 
# b=None, 
# copy=False, 
# method='umap')
# fig, ax = plt.subplots(figsize=(10, 7))
# sc.pl.umap(adata, color='louvain', ax=ax)
```

## Discussion

This notebook has demonstrated visualization of cells following pre-processing of single-cell RNA-seq data.


```
# Running time of the notebook
print("{:.2f} minutes".format((time.time()-start_time)/60))
```

    7.97 minutes


**Feedback**: please report any issues, or submit pull requests for improvements, in the [Github repository where this notebook is located](https://github.com/pachterlab/kallistobustools/blob/master/notebooks/kb_intro_2_python.ipynb).
