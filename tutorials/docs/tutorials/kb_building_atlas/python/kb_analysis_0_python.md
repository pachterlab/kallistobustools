<a href="https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/tutorials/docs/tutorials/kb_building_atlas/python/kb_analysis_0_python.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# Analysis of single-cell RNA-seq data: building and annotating an atlas
This Python notebook pre-processes the [pbmc_1k v3 dataset](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3) from 10X Genomics with kallisto and bustools using `kb`, and then performs an analysis of the cell types and their marker genes.

The notebook was written by A. Sina Booeshaghi and Lior Pachter and is based on three noteboks:
- The kallisto | bustools [Introduction to single-cell RNA-seq I](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_1_minute_intro.ipynb#scrollTo=wtwMjIjjCMcD) notebook.
- The kallisto | bustools [Introduction to single-cell RNA-seq II](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_standard.ipynb#scrollTo=ijU_u6uj3Sio) notebook.
- The Scanpy [Preprocessing and clustering 3k PBMCs" notebook](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html).

If you use the methods in this notebook for your analysis please cite the following publications which describe the tools used in the notebook:

* Melsted, P., Booeshaghi, A.S. et al. Modular and efficient pre-processing of single-cell RNA-seq. bioRxiv (2019). doi:10.1101/673285
* Wolf, F. A., Angere, P. and Theis, F.J. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology (2018). doi:10.1186/s13059-017-1382-0

An R notebook implementing the same analysis is available [here](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_analysis_0_R.ipynb). See the [kallistobus.tools tutorials](https://www.kallistobus.tools/tutorials) site for additional notebooks demonstrating other analyses.

## Setup


```python
# This is  used to time the running of the notebook
import time
start_time = time.time()
```

### Install python packages


```python
%%time
# These packages are pre-installed on Google Colab, but are included here to simplify running this notebook locally
!pip install --quiet matplotlib scikit-learn numpy scipy
!pip3 install --quiet leidenalg
!pip install --quiet louvain scanpy
```

    [K     |████████████████████████████████| 2.4MB 6.2MB/s 
    [K     |████████████████████████████████| 3.2MB 41.0MB/s 
    [K     |████████████████████████████████| 2.2MB 4.4MB/s 
    [K     |████████████████████████████████| 10.3MB 23.5MB/s 
    [K     |████████████████████████████████| 81kB 6.0MB/s 
    [K     |████████████████████████████████| 122kB 42.3MB/s 
    [K     |████████████████████████████████| 51kB 5.0MB/s 
    [K     |████████████████████████████████| 1.2MB 38.3MB/s 
    [K     |████████████████████████████████| 71kB 6.3MB/s 
    [?25h  Building wheel for umap-learn (setup.py) ... [?25l[?25hdone
      Building wheel for sinfo (setup.py) ... [?25l[?25hdone
      Building wheel for pynndescent (setup.py) ... [?25l[?25hdone
    CPU times: user 231 ms, sys: 57.1 ms, total: 288 ms
    Wall time: 19.6 s


### Install kb-python


```python
%%time
# install kb
!pip install --quiet kb-python 
```

    [K     |████████████████████████████████| 59.1MB 110kB/s 
    [K     |████████████████████████████████| 13.2MB 52.4MB/s 
    [K     |████████████████████████████████| 51kB 4.3MB/s 
    [K     |████████████████████████████████| 112kB 41.7MB/s 
    [?25h  Building wheel for loompy (setup.py) ... [?25l[?25hdone
      Building wheel for numpy-groupies (setup.py) ... [?25l[?25hdone
    CPU times: user 380 ms, sys: 61.1 ms, total: 441 ms
    Wall time: 31.2 s


### Download the data


```python
%%time
# Download the data from the 10x website
!wget -q http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar

# unpack the downloaded files
!tar -xf pbmc_1k_v3_fastqs.tar
```

    CPU times: user 1.41 s, sys: 199 ms, total: 1.6 s
    Wall time: 3min 3s


### Download an index

This data consists of peripheral blood mononuclear cells from a human, so we download the human index.


```python
!kb ref -d human -i index.idx -g t2g.txt -f1 transcriptome.fasta
```

    [2021-03-31 21:09:56,028]    INFO Downloading files for human from https://caltech.box.com/shared/static/v1nm7lpnqz5syh8dyzdk2zs8bglncfib.gz to tmp/v1nm7lpnqz5syh8dyzdk2zs8bglncfib.gz
    100% 2.23G/2.23G [01:44<00:00, 22.8MB/s]
    [2021-03-31 21:11:42,360]    INFO Extracting files from tmp/v1nm7lpnqz5syh8dyzdk2zs8bglncfib.gz


## Pseudoalignment and counting

### Run kallisto and bustools


```python
%%time
!kb count --h5ad -i index.idx -g t2g.txt -x 10xv3 -o output --filter bustools -t 2 \
pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz \
pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz \
pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz \
pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz
```

    [2021-03-31 21:12:22,466]    INFO Using index index.idx to generate BUS file to output from
    [2021-03-31 21:12:22,466]    INFO         pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz
    [2021-03-31 21:12:22,467]    INFO         pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz
    [2021-03-31 21:12:22,467]    INFO         pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz
    [2021-03-31 21:12:22,467]    INFO         pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz
    [2021-03-31 21:27:43,495]    INFO Sorting BUS file output/output.bus to output/tmp/output.s.bus
    [2021-03-31 21:28:15,485]    INFO Whitelist not provided
    [2021-03-31 21:28:15,485]    INFO Copying pre-packaged 10XV3 whitelist to output
    [2021-03-31 21:28:16,575]    INFO Inspecting BUS file output/tmp/output.s.bus
    [2021-03-31 21:28:39,940]    INFO Correcting BUS records in output/tmp/output.s.bus to output/tmp/output.s.c.bus with whitelist output/10xv3_whitelist.txt
    [2021-03-31 21:29:05,490]    INFO Sorting BUS file output/tmp/output.s.c.bus to output/output.unfiltered.bus
    [2021-03-31 21:29:34,684]    INFO Generating count matrix output/counts_unfiltered/cells_x_genes from BUS file output/output.unfiltered.bus
    [2021-03-31 21:29:51,979]    INFO Reading matrix output/counts_unfiltered/cells_x_genes.mtx
    [2021-03-31 21:30:00,389]    INFO Writing matrix to h5ad output/counts_unfiltered/adata.h5ad
    [2021-03-31 21:30:01,376]    INFO Filtering with bustools
    [2021-03-31 21:30:01,376]    INFO Generating whitelist output/filter_barcodes.txt from BUS file output/output.unfiltered.bus
    [2021-03-31 21:30:01,623]    INFO Correcting BUS records in output/output.unfiltered.bus to output/tmp/output.unfiltered.c.bus with whitelist output/filter_barcodes.txt
    [2021-03-31 21:30:15,206]    INFO Sorting BUS file output/tmp/output.unfiltered.c.bus to output/output.filtered.bus
    [2021-03-31 21:30:37,631]    INFO Generating count matrix output/counts_filtered/cells_x_genes from BUS file output/output.filtered.bus
    [2021-03-31 21:30:52,494]    INFO Reading matrix output/counts_filtered/cells_x_genes.mtx
    [2021-03-31 21:30:58,622]    INFO Writing matrix to h5ad output/counts_filtered/adata.h5ad
    CPU times: user 8.81 s, sys: 1.03 s, total: 9.84 s
    Wall time: 18min 38s


## Basic QC


```python
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from sklearn.decomposition import TruncatedSVD
import matplotlib
import matplotlib.pyplot as plt

```


```python
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80)
```


```python
# load the unfiltered matrix
results_file = 'pbmc1k.h5ad'  # the file that will store the analysis results
adata = anndata.read_h5ad("output/counts_unfiltered/adata.h5ad")
adata.var["gene_id"] = adata.var.index.values

t2g = pd.read_csv("t2g.txt", header=None, names=["tid", "gene_id", "gene_name"], sep="\t")
t2g.index = t2g.gene_id
t2g = t2g.loc[~t2g.index.duplicated(keep='first')]

adata.var["gene_name"] = adata.var.gene_id.map(t2g["gene_name"])
adata.var.index = adata.var["gene_name"]

adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
```

    /usr/local/lib/python3.7/dist-packages/anndata/utils.py:117: UserWarning: Suffix used (-[0-9]+) to deduplicate index values may make index values difficult to interpret. There values with a similar suffixes in the index. Consider using a different delimiter by passing `join={delimiter}`Example key collisions generated by the make_index_unique algorithm: ['SNORD116-1', 'SNORD116-2', 'SNORD116-3', 'SNORD116-4', 'SNORD116-5']
      + str(example_colliding_values)



```python
adata
```




    AnnData object with n_obs × n_vars = 259615 × 60623
        var: 'gene_name', 'gene_id'



### Test for library saturation


```python
# Create a plot showing genes detected as a function of UMI counts.
fig, ax = plt.subplots(figsize=(7, 7))

x = np.asarray(adata.X.sum(axis=1))[:,0]
y = np.asarray(np.sum(adata.X>0, axis=1))[:,0]

ax.scatter(x, y, color="green", alpha=0.25)
ax.set_xlabel("UMI Counts")
ax.set_ylabel("Genes Detected")
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim(1)
ax.set_ylim(1)


plt.show()
```


![png](kb_analysis_0_python_files/kb_analysis_0_python_22_0.png)


This plot is very misleading, as even the small alpha can't accurately show how many points are stacked at one location (This takes about a minute to run since there are a lot of points)


```python
fig, ax = plt.subplots(figsize=(7,7))

#histogram definition
bins = [1500, 1500] # number of bins

# histogram the data
hh, locx, locy = np.histogram2d(x, y, bins=bins)

# Sort the points by density, so that the densest points are plotted last
z = np.array([hh[np.argmax(a<=locx[1:]),np.argmax(b<=locy[1:])] for a,b in zip(x,y)])
idx = z.argsort()
x2, y2, z2 = x[idx], y[idx], z[idx]


s = ax.scatter(x2, y2, c=z2, cmap='Greens')  
fig.colorbar(s, ax=ax)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel("UMI Counts")
ax.set_ylabel("Genes Detected")

ax.set_xlim(1, 10**5)
ax.set_ylim(1, 10**4)

plt.show()
```


![png](kb_analysis_0_python_files/kb_analysis_0_python_24_0.png)


### Examine the knee plot

The "knee plot" was introduced in the Drop-seq paper: 
- Macosko et al., [Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets](https://www.cell.com/fulltext/S0092-8674(15)00549-8), 2015. DOI:10.1016/j.cell.2015.05.002

In this plot cells are ordered by the number of UMI counts associated to them (shown on the *x*-axis), and the fraction of droplets with at least that number of cells is shown on the *y*-axis:


```python
#@title Threshold cells according to knee plot { run: "auto", vertical-output: true }
expected_num_cells =  1178#@param {type:"integer"}
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


![png](kb_analysis_0_python_files/kb_analysis_0_python_26_0.png)


## Analysis

It is useful to examine mitochondrial genes, which are important for quality control. [(Lun, McCarthy & Marioni, 2017)](https://master.bioconductor.org/packages/release/workflows/html/simpleSingleCell.html#examining-gene-level-metrics) write that

> High proportions are indicative of poor-quality cells (Islam et al. 2014; Ilicic et al. 2016), possibly because of loss of cytoplasmic RNA from perforated cells. The reasoning is that mitochondria are larger than individual transcript molecules and less likely to escape through tears in the cell membrane.

Note you can also use the function `pp.calculate_qc_metrics` to compute the fraction of mitochondrial genes and additional measures.

Show those genes that yield the highest fraction of counts in each single cells, across all cells.


```python
sc.pl.highest_expr_genes(adata, n_top=20)
```

    normalizing counts per cell
    WARNING: Some cells have total count of genes equal to zero
        finished (0:00:00)



![png](kb_analysis_0_python_files/kb_analysis_0_python_31_1.png)


### Filter

Begin by filtering cells according to various criteria. First, a filter for genes and cells based on minimum thresholds:


```python
# Removes cells with less than 1070 umi counts
adata = adata[np.asarray(adata.X.sum(axis=1)).reshape(-1) > 1070]

# Removes genes with 0 umi counts
adata = adata[:, np.asarray(adata.X.sum(axis=0)).reshape(-1) > 0]
```


```python
adata
```




    View of AnnData object with n_obs × n_vars = 1180 × 31837
        var: 'gene_name', 'gene_id'




```python
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
```

    filtered out 5 cells that have less than 200 genes expressed
    Trying to set attribute `.obs` of view, copying.
    filtered out 5884 genes that are detected in less than 3 cells



```python
adata
```




    AnnData object with n_obs × n_vars = 1175 × 25953
        obs: 'n_genes'
        var: 'gene_name', 'gene_id', 'n_cells'



Next, filter by mitochondrial gene content


```python
mito_genes = adata.var_names.str.startswith('MT-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1
```

Perform a QC check of the counts post-filtering


```python
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True)
```

    ... storing 'gene_name' as categorical



![png](kb_analysis_0_python_files/kb_analysis_0_python_41_1.png)



```python
#examine mitochondrial content 
sc.pl.scatter(adata, x='n_counts', y='percent_mito')
sc.pl.scatter(adata, x='n_counts', y='n_genes')
```


![png](kb_analysis_0_python_files/kb_analysis_0_python_42_0.png)



![png](kb_analysis_0_python_files/kb_analysis_0_python_42_1.png)



```python
# Create a mask to filter out cells with more than 6500 genes, less than 200 genes or less than 0.2 mitochondrial umi counts
mask = np.logical_or((adata.obs.n_genes < 6500).values, (adata.obs.n_genes > 200).values, (adata.obs.percent_mito < 0.2).values)
```


```python
#filter
adata = adata[mask, :]
```


```python
adata
```




    View of AnnData object with n_obs × n_vars = 1175 × 25953
        obs: 'n_genes', 'percent_mito', 'n_counts'
        var: 'gene_name', 'gene_id', 'n_cells'



### Normalize counts

Total-count normalize (library-size correct) the data matrix $\mathbf{X}$ to 10,000 reads per cell, so that counts become comparable among cells.


```python
# normalize counts in each cell to be equal
sc.pp.normalize_total(adata, target_sum=10**4)
```

    /usr/local/lib/python3.7/dist-packages/scanpy/preprocessing/_normalization.py:138: UserWarning: Revieved a view of an AnnData. Making a copy.
      view_to_actual(adata)
    normalizing counts per cell
        finished (0:00:00)


Log the counts


```python
# Replace raw counts with their logarithm
sc.pp.log1p(adata)
```

Lets now look at the highest expressed genes after filtering, normalization, and log


```python
sc.pl.highest_expr_genes(adata, n_top=20)
```

    normalizing counts per cell
        finished (0:00:00)



![png](kb_analysis_0_python_files/kb_analysis_0_python_52_1.png)


Set the `.raw` attribute of AnnData object to the normalized and logarithmized raw gene expression for later use in differential testing and visualizations of gene expression. This simply freezes the state of the AnnData object.

The unnormalized data is stored in `.raw`.


```python
adata.raw = adata
```

<div class="alert alert-info">

**Note**
The result of the following highly-variable-genes detection is stored as an annotation in `.var.highly_variable` and auto-detected by PCA and hence, `sc.pp.neighbors` and subsequent manifold/graph tools.

</div>

### Identify highly-variable genes.


```python
# flavor="cell_ranger" is consistent with Seurat and flavor="suerat" is not consistent with Seurat
sc.pp.highly_variable_genes(adata, min_mean=0.01, max_mean=8, min_disp=1, n_top_genes=2000, flavor="cell_ranger", n_bins=20)
```

    If you pass `n_top_genes`, all cutoffs are ignored.
    extracting highly variable genes
    /usr/local/lib/python3.7/dist-packages/statsmodels/tools/_testing.py:19: FutureWarning: pandas.util.testing is deprecated. Use the functions in the public API at pandas.testing instead.
      import pandas.util.testing as tm
        finished (0:00:01)
    --> added
        'highly_variable', boolean vector (adata.var)
        'means', float vector (adata.var)
        'dispersions', float vector (adata.var)
        'dispersions_norm', float vector (adata.var)



```python
sc.pl.highly_variable_genes(adata)
```


![png](kb_analysis_0_python_files/kb_analysis_0_python_58_0.png)


Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.

We do not regress out as per https://github.com/theislab/scanpy/issues/526


```python
# sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
```

### Scaling the data
Scale each gene to unit variance. Clip values exceeding standard deviation 10. 


```python
sc.pp.scale(adata, max_value=10)
```

    ... as `zero_center=True`, sparse input is densified and may lead to large memory consumption


### Principal component analysis

Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.


```python
# We perform PCA on just the highly variable genes
sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True)
```

    computing PCA
        on highly variable genes
        with n_comps=50
        finished (0:00:00)



```python
sc.pl.pca(adata, color='CST3')
```


![png](kb_analysis_0_python_files/kb_analysis_0_python_67_0.png)


Let us inspect the contribution of single PCs to the total variance in the data. This gives us information about how many PCs we should consider in order to compute the neighborhood relations of cells, e.g. used in the clustering function  `sc.tl.leiden()` or tSNE `sc.tl.tsne()`. In our experience, often, a rough estimate of the number of PCs does fine.


```python
sc.pl.pca_variance_ratio(adata, log=True)
```


![png](kb_analysis_0_python_files/kb_analysis_0_python_69_0.png)



```python
adata
```




    AnnData object with n_obs × n_vars = 1175 × 25953
        obs: 'n_genes', 'percent_mito', 'n_counts'
        var: 'gene_name', 'gene_id', 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'mean', 'std'
        uns: 'log1p', 'hvg', 'pca'
        obsm: 'X_pca'
        varm: 'PCs'



### Compute the neighborhood graph

Next we compute the neighborhood graph of cells using the PCA representation of the data matrix. You might simply use default values here. In order to be consistent with Seurat's results, we use the following values.


```python
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=10)
```

    computing neighbors
        using 'X_pca' with n_pcs = 10
        finished: added to `.uns['neighbors']`
        `.obsp['distances']`, distances for each pair of neighbors
        `.obsp['connectivities']`, weighted adjacency matrix (0:00:29)


### Embed the neighborhood graph

### UMAP

UMAP (UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction) is a manifold learning technique that can also be used to visualize cells. It was published in:

- McInnes, Leland, John Healy, and James Melville. "Umap: Uniform manifold approximation and projection for dimension reduction." arXiv preprint arXiv:1802.03426 (2018).

We run that to visualize the results:

```
tl.paga(adata)
pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
tl.umap(adata, init_pos='paga')
```


```python
sc.tl.umap(adata)
```

    computing UMAP
        finished: added
        'X_umap', UMAP coordinates (adata.obsm) (0:00:04)



```python
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])
```


![png](kb_analysis_0_python_files/kb_analysis_0_python_77_0.png)


### Cluster the neighborhood graph

### Clustering

There are many algorithms for clustering cells, and while they have been compared in detail in various benchmarks (see e.g., [Duo et al. 2018](https://f1000research.com/articles/7-1141/v2)), there is no univerally agreed upon method. Here we demonstrate clustering using [Louvain clustering](https://en.wikipedia.org/wiki/Louvain_modularity), which is a popular method for clustering single-cell RNA-seq data. The method was published in 

- Blondel, Vincent D; Guillaume, Jean-Loup; Lambiotte, Renaud; Lefebvre, Etienne (9 October 2008). "Fast unfolding of communities in large networks". Journal of Statistical Mechanics: Theory and Experiment. 2008 (10): P10008.

Note that Louvain clustering directly clusters the neighborhood graph of cells, which we already computed in the previous section.


```python
sc.tl.louvain(adata,resolution=0.5, random_state=42)
```

    running Louvain clustering
        using the "louvain" package of Traag (2017)
        finished: found 8 clusters and added
        'louvain', the cluster labels (adata.obs, categorical) (0:00:00)


A plot of the clusters is shown below:


```python
sc.pl.umap(adata, color=['louvain', 'CST3', 'NKG7'])
```


![png](kb_analysis_0_python_files/kb_analysis_0_python_82_0.png)


### Find marker genes

A key aspect of annotating a cell atlas is identifying "marker genes". These are genes specific to individual clusters that "mark" them, and are important both for assigning functions to cell clusters, and for designing downstream experiments to probe activity of clusters. 

A gene marker analysis begins with ranking genes in each cluster according to how different they are relative to other clusters. Typically the t-test is used for this purpose.


```python
sc.tl.rank_genes_groups(adata, 'louvain', method='t-test', corr_method="bonferroni")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
```

    ranking genes
        finished: added to `.uns['rank_genes_groups']`
        'names', sorted np.recarray to be indexed by group ids
        'scores', sorted np.recarray to be indexed by group ids
        'logfoldchanges', sorted np.recarray to be indexed by group ids
        'pvals', sorted np.recarray to be indexed by group ids
        'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:00)



![png](kb_analysis_0_python_files/kb_analysis_0_python_85_1.png)



```python
sc.settings.verbosity = 2  # reduce the verbosity
```

An alternative to the parametric t-test is the non-parametric [Wilcoxon rank-sum (Mann-Whitney-U)](https://de.wikipedia.org/wiki/Wilcoxon-Mann-Whitney-Test) test.


```python
sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon', corr_method="bonferroni")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
```

    ranking genes
        finished (0:00:04)



![png](kb_analysis_0_python_files/kb_analysis_0_python_88_1.png)


As an alternative, genes can be ranked using logistic regression (see [Natranos et al. (2018)](https://doi.org/10.1101/258566)).


```python
tmp = adata.copy()
```


```python
sc.tl.rank_genes_groups(tmp, 'louvain', method='logreg')
sc.pl.rank_genes_groups(tmp, n_genes=25, sharey=False)
```

    ranking genes
    /usr/local/lib/python3.7/dist-packages/sklearn/linear_model/_logistic.py:940: ConvergenceWarning: lbfgs failed to converge (status=1):
    STOP: TOTAL NO. of ITERATIONS REACHED LIMIT.
    
    Increase the number of iterations (max_iter) or scale the data as shown in:
        https://scikit-learn.org/stable/modules/preprocessing.html
    Please also refer to the documentation for alternative solver options:
        https://scikit-learn.org/stable/modules/linear_model.html#logistic-regression
      extra_warning_msg=_LOGISTIC_SOLVER_CONVERGENCE_MSG)
        finished (0:00:20)



![png](kb_analysis_0_python_files/kb_analysis_0_python_91_1.png)


With the exceptions of *IL7R*, which is only found by the t-test and *FCER1A*, which is only found by the other two appraoches, all marker genes are recovered in all approaches.

Louvain Group | Markers | Cell Type
---|---|---
0 | IL7R | CD4 T cells
1 | CD14, LYZ | CD14+ Monocytes
2 | MS4A1 |	B cells
3 | GNLY, NKG7 | 	NK cells
4 | FCGR3A, MS4A7 |	FCGR3A+ Monocytes
5 | CD8A |	CD8 T cells
6 | MS4A1 |	B cells

Let us also define a list of marker genes for later reference.


```python
marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',  
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']
```

Reload the object that has been save with the Wilcoxon Rank-Sum test result.

Show the 10 top ranked genes per cluster 0, 1, ..., 7 in a dataframe.


```python
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
      <th>1</th>
      <th>2</th>
      <th>3</th>
      <th>4</th>
      <th>5</th>
      <th>6</th>
      <th>7</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>S100A8</td>
      <td>IL7R</td>
      <td>RPL32</td>
      <td>IGHD</td>
      <td>NKG7</td>
      <td>CST3</td>
      <td>KLRB1</td>
      <td>CD79A</td>
    </tr>
    <tr>
      <th>1</th>
      <td>S100A9</td>
      <td>IL32</td>
      <td>RPS12</td>
      <td>IGHM</td>
      <td>CTSW</td>
      <td>HLA-DPA1</td>
      <td>KLRG1</td>
      <td>IGKC</td>
    </tr>
    <tr>
      <th>2</th>
      <td>VCAN</td>
      <td>TRAC</td>
      <td>RPL30</td>
      <td>CD37</td>
      <td>GZMA</td>
      <td>NPC2</td>
      <td>NKG7</td>
      <td>CD79B</td>
    </tr>
    <tr>
      <th>3</th>
      <td>S100A12</td>
      <td>LDHB</td>
      <td>RPL35A</td>
      <td>CD79A</td>
      <td>CST7</td>
      <td>HLA-DPB1</td>
      <td>IL32</td>
      <td>BANK1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>FOS</td>
      <td>BCL11B</td>
      <td>RPS15A</td>
      <td>MS4A1</td>
      <td>GZMM</td>
      <td>FCER1G</td>
      <td>GZMA</td>
      <td>TNFRSF13C</td>
    </tr>
  </tbody>
</table>
</div>



Get a table with the scores and groups.


```python
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(5)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0_n</th>
      <th>0_p</th>
      <th>1_n</th>
      <th>1_p</th>
      <th>2_n</th>
      <th>2_p</th>
      <th>3_n</th>
      <th>3_p</th>
      <th>4_n</th>
      <th>4_p</th>
      <th>5_n</th>
      <th>5_p</th>
      <th>6_n</th>
      <th>6_p</th>
      <th>7_n</th>
      <th>7_p</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>S100A8</td>
      <td>4.325218e-141</td>
      <td>IL7R</td>
      <td>3.851653e-57</td>
      <td>RPL32</td>
      <td>7.151780e-62</td>
      <td>IGHD</td>
      <td>1.303809e-75</td>
      <td>NKG7</td>
      <td>9.495682e-64</td>
      <td>CST3</td>
      <td>1.595268e-42</td>
      <td>KLRB1</td>
      <td>8.274802e-41</td>
      <td>CD79A</td>
      <td>3.916006e-28</td>
    </tr>
    <tr>
      <th>1</th>
      <td>S100A9</td>
      <td>1.357507e-140</td>
      <td>IL32</td>
      <td>3.390931e-49</td>
      <td>RPS12</td>
      <td>4.320062e-61</td>
      <td>IGHM</td>
      <td>1.624935e-74</td>
      <td>CTSW</td>
      <td>1.281002e-58</td>
      <td>HLA-DPA1</td>
      <td>9.721044e-42</td>
      <td>KLRG1</td>
      <td>1.245774e-34</td>
      <td>IGKC</td>
      <td>7.115329e-27</td>
    </tr>
    <tr>
      <th>2</th>
      <td>VCAN</td>
      <td>1.698042e-136</td>
      <td>TRAC</td>
      <td>5.952130e-49</td>
      <td>RPL30</td>
      <td>2.823714e-59</td>
      <td>CD37</td>
      <td>1.704719e-73</td>
      <td>GZMA</td>
      <td>5.700633e-57</td>
      <td>NPC2</td>
      <td>5.443673e-33</td>
      <td>NKG7</td>
      <td>7.556990e-32</td>
      <td>CD79B</td>
      <td>6.053017e-26</td>
    </tr>
    <tr>
      <th>3</th>
      <td>S100A12</td>
      <td>2.736086e-136</td>
      <td>LDHB</td>
      <td>6.487053e-46</td>
      <td>RPL35A</td>
      <td>6.098855e-58</td>
      <td>CD79A</td>
      <td>2.184785e-70</td>
      <td>CST7</td>
      <td>1.401230e-56</td>
      <td>HLA-DPB1</td>
      <td>7.522829e-33</td>
      <td>IL32</td>
      <td>1.135243e-30</td>
      <td>BANK1</td>
      <td>6.905314e-26</td>
    </tr>
    <tr>
      <th>4</th>
      <td>FOS</td>
      <td>3.264834e-132</td>
      <td>BCL11B</td>
      <td>4.040622e-42</td>
      <td>RPS15A</td>
      <td>6.949561e-58</td>
      <td>MS4A1</td>
      <td>6.806009e-70</td>
      <td>GZMM</td>
      <td>5.222455e-51</td>
      <td>FCER1G</td>
      <td>3.106843e-32</td>
      <td>GZMA</td>
      <td>1.885113e-30</td>
      <td>TNFRSF13C</td>
      <td>7.946763e-24</td>
    </tr>
  </tbody>
</table>
</div>



Compare to a single cluster. 


```python
sc.tl.rank_genes_groups(adata, 'louvain', groups=['0'], reference='1', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['0'], n_genes=20)
```

    ranking genes
        finished (0:00:01)



![png](kb_analysis_0_python_files/kb_analysis_0_python_101_1.png)


If we want a more detailed view for a certain group, use `sc.pl.rank_genes_groups_violin`.


```python
sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)
```


![png](kb_analysis_0_python_files/kb_analysis_0_python_103_0.png)


If you want to compare a certain gene across groups, use the following.


```python
sc.pl.violin(adata, marker_genes, groupby='louvain')
```


![png](kb_analysis_0_python_files/kb_analysis_0_python_105_0.png)


Actually mark the cell types.


```python
pd.DataFrame(result["names"]).head(10)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
      <th>1</th>
      <th>2</th>
      <th>3</th>
      <th>4</th>
      <th>5</th>
      <th>6</th>
      <th>7</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>S100A8</td>
      <td>IL7R</td>
      <td>RPL32</td>
      <td>IGHD</td>
      <td>NKG7</td>
      <td>CST3</td>
      <td>KLRB1</td>
      <td>CD79A</td>
    </tr>
    <tr>
      <th>1</th>
      <td>S100A9</td>
      <td>IL32</td>
      <td>RPS12</td>
      <td>IGHM</td>
      <td>CTSW</td>
      <td>HLA-DPA1</td>
      <td>KLRG1</td>
      <td>IGKC</td>
    </tr>
    <tr>
      <th>2</th>
      <td>VCAN</td>
      <td>TRAC</td>
      <td>RPL30</td>
      <td>CD37</td>
      <td>GZMA</td>
      <td>NPC2</td>
      <td>NKG7</td>
      <td>CD79B</td>
    </tr>
    <tr>
      <th>3</th>
      <td>S100A12</td>
      <td>LDHB</td>
      <td>RPL35A</td>
      <td>CD79A</td>
      <td>CST7</td>
      <td>HLA-DPB1</td>
      <td>IL32</td>
      <td>BANK1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>FOS</td>
      <td>BCL11B</td>
      <td>RPS15A</td>
      <td>MS4A1</td>
      <td>GZMM</td>
      <td>FCER1G</td>
      <td>GZMA</td>
      <td>TNFRSF13C</td>
    </tr>
    <tr>
      <th>5</th>
      <td>MNDA</td>
      <td>CD3D</td>
      <td>RPS6</td>
      <td>CD74</td>
      <td>PRF1</td>
      <td>COTL1</td>
      <td>CTSW</td>
      <td>MS4A1</td>
    </tr>
    <tr>
      <th>6</th>
      <td>LYZ</td>
      <td>CD3E</td>
      <td>RPS14</td>
      <td>TCL1A</td>
      <td>HCST</td>
      <td>AIF1</td>
      <td>CCL5</td>
      <td>TCF4</td>
    </tr>
    <tr>
      <th>7</th>
      <td>SRGN</td>
      <td>CD3G</td>
      <td>RPS3A</td>
      <td>CD79B</td>
      <td>KLRD1</td>
      <td>CTSZ</td>
      <td>IL7R</td>
      <td>CD74</td>
    </tr>
    <tr>
      <th>8</th>
      <td>FCN1</td>
      <td>TCF7</td>
      <td>RPS28</td>
      <td>LINC00926</td>
      <td>CCL5</td>
      <td>HLA-DRA</td>
      <td>AQP3</td>
      <td>HLA-DQA1</td>
    </tr>
    <tr>
      <th>9</th>
      <td>DUSP1</td>
      <td>CD2</td>
      <td>RPS27</td>
      <td>BCL11A</td>
      <td>B2M</td>
      <td>HLA-DRB1</td>
      <td>SLC4A10</td>
      <td>LINC00926</td>
    </tr>
  </tbody>
</table>
</div>



Recall

Louvain Group | Markers | Cell Type
---|---|---
0 | IL7R | CD4 T cells
1 | CD14, LYZ | CD14+ Monocytes
2 | MS4A1 |	B cells
3 | GNLY, NKG7 | 	NK cells
4 | FCGR3A, MS4A7 |	FCGR3A+ Monocytes
5 | CD8A |	CD8 T cells
6 | MS4A1 |	B cells


```python
## Map the gene ids to a cluster label
# Each cluster (index in the list) corresponds to a cell type
new_cluster_names = [
                     "CD4 T",
                     "B Cells",
                     "NK",
                     "CD14 Monocytes",
                     "FCGR3A Monocytes",
                     "CD8 T",
                     "B-2",
]
```


```python
# Relabel the clusters
# adata.rename_categories('louvain', new_cluster_names)
```


```python
sc.pl.umap(adata, color='louvain', legend_loc='on data', title='', frameon=False)
```


![png](kb_analysis_0_python_files/kb_analysis_0_python_111_0.png)


Now that we annotated the cell types, let us visualize the marker genes.


```python
ax = sc.pl.dotplot(adata, marker_genes, groupby='louvain')
```


![png](kb_analysis_0_python_files/kb_analysis_0_python_113_0.png)


There is also a very compact violin plot.


```python
ax = sc.pl.stacked_violin(adata, marker_genes, groupby='louvain', rotation=90)
```


![png](kb_analysis_0_python_files/kb_analysis_0_python_115_0.png)


Note that as a result of the analysis the adata object has accumulated several annotations:


```python
adata
```




    AnnData object with n_obs × n_vars = 1175 × 25953
        obs: 'n_genes', 'percent_mito', 'n_counts', 'louvain'
        var: 'gene_name', 'gene_id', 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'mean', 'std'
        uns: 'log1p', 'hvg', 'pca', 'neighbors', 'umap', 'louvain', 'louvain_colors', 'rank_genes_groups'
        obsm: 'X_pca', 'X_umap'
        varm: 'PCs'
        obsp: 'distances', 'connectivities'



If you want to export to "csv", you have the following options:


```python
# Export single fields of the annotation of observations
# adata.obs[['n_counts', 'louvain_groups']].to_csv(
#     './write/pbmc3k_corrected_louvain_groups.csv')

# Export single columns of the multidimensional annotation
# adata.obsm.to_df()[['X_pca1', 'X_pca2']].to_csv(
#     './write/pbmc3k_corrected_X_pca.csv')

# Or export everything except the data using `.write_csvs`.
# Set `skip_data=False` if you also want to export the data.
# adata.write_csvs(results_file[:-5], )
```


```python
# Running time of the notebook
print("{:.2f} minutes".format((time.time()-start_time)/60))
```

    26.77 minutes


**Feedback**: please report any issues, or submit pull requests for improvements, in the [Github repository where this notebook is located](https://github.com/pachterlab/kallistobustools/).


```python

```
