<a href="https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_analysis_0_R.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# Analysis of single-cell RNA-seq data: building and annotating an atlas
This R notebook pre-processes the [pbmc_1k v3 dataset](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3) from 10X Genomics with kallisto and bustools using `kb`, and then performs an analysis of the cell types and their marker genes.

The notebook was written by A. Sina Booeshaghi, Lambda Lu and Lior Pachter and is based on three noteboks:
- The kallisto | bustools [Introduction to single-cell RNA-seq I](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_1_minute_intro.ipynb#scrollTo=wtwMjIjjCMcD) notebook.
- The kallisto | bustools [Introduction to single-cell RNA-seq II](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_standard.ipynb#scrollTo=ijU_u6uj3Sio) notebook.
- The Seurat [Guided Clustering Tutorial](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html).

If you use the methods in this notebook for your analysis please cite the following publications which describe the tools used in the notebook:

* Melsted, P., Booeshaghi, A.S. et al. Modular and efficient pre-processing of single-cell RNA-seq. bioRxiv (2019). doi:10.1101/673285
* Stuart, Butler et al. Comprehensive Integration of Single-cell Data. Cell (2019). doi:10.1016/j.cell.2019.05.031

A Python notebook implementing the same analysis is available [here](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_analysis_0_python.ipynb#scrollTo=3kZG9UyUPE_B). See the [kallistobus.tools tutorials](https://www.kallistobus.tools/tutorials) site for additional notebooks demonstrating other analyses.


## Setup


```R
# This is used to time the running of the notebook
start_time <- Sys.time()
```

### Install R packages
A large fraction of the running time of this notebook is in installing the Seurat R package, since it has lots of dependencies and many of them use Rcpp which results in the need to compile lots of C++ code.


```R
system.time({
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(c('multtest', "DropletUtils"), Ncpus = 2)
  install.packages(c("Seurat", "scico"), Ncpus = 2)
})
```

    Installing package into ‘/usr/local/lib/R/site-library’
    (as ‘lib’ is unspecified)
    
    Bioconductor version 3.10 (BiocManager 1.30.10), R 3.6.2 (2019-12-12)
    
    Installing package(s) 'BiocVersion', 'multtest', 'DropletUtils'
    
    also installing the dependencies ‘zlibbioc’, ‘bitops’, ‘XVector’, ‘RCurl’, ‘GenomeInfoDbData’, ‘formatR’, ‘GenomicRanges’, ‘GenomeInfoDb’, ‘lambda.r’, ‘futile.options’, ‘matrixStats’, ‘SummarizedExperiment’, ‘futile.logger’, ‘snow’, ‘limma’, ‘locfit’, ‘DelayedArray’, ‘IRanges’, ‘R.oo’, ‘R.methodsS3’, ‘sitmo’, ‘BiocGenerics’, ‘Biobase’, ‘SingleCellExperiment’, ‘S4Vectors’, ‘BiocParallel’, ‘edgeR’, ‘rhdf5’, ‘HDF5Array’, ‘R.utils’, ‘dqrng’, ‘beachmat’, ‘Rhdf5lib’
    
    
    Old packages: 'curl', 'DT', 'farver', 'jsonlite', 'knitr', 'mime', 'rprojroot',
      'rstudioapi', 'svglite', 'xfun', 'xtable', 'nlme'
    
    Installing packages into ‘/usr/local/lib/R/site-library’
    (as ‘lib’ is unspecified)
    
    also installing the dependencies ‘mnormt’, ‘numDeriv’, ‘TH.data’, ‘sandwich’, ‘lsei’, ‘bibtex’, ‘gbRd’, ‘sn’, ‘mvtnorm’, ‘plotrix’, ‘multcomp’, ‘gtools’, ‘gdata’, ‘caTools’, ‘npsurv’, ‘globals’, ‘listenv’, ‘zoo’, ‘Rdpack’, ‘TFisher’, ‘mutoss’, ‘hexbin’, ‘data.table’, ‘rappdirs’, ‘gplots’, ‘gridExtra’, ‘RcppEigen’, ‘FNN’, ‘RSpectra’, ‘RcppParallel’, ‘RcppProgress’
    
    
    Downloading GitHub repo satijalab/Seurat@master
    


    curl     (4.2   -> 4.3  ) [CRAN]
    jsonlite (1.6   -> 1.6.1) [CRAN]
    mime     (0.8   -> 0.9  ) [CRAN]
    farver   (2.0.1 -> 2.0.3) [CRAN]
    xtable   (1.8-3 -> 1.8-4) [CRAN]


    Skipping 3 packages ahead of CRAN: multtest, BiocGenerics, Biobase
    
    Installing 5 packages: curl, jsonlite, mime, farver, xtable
    
    Installing packages into ‘/usr/local/lib/R/site-library’
    (as ‘lib’ is unspecified)
    


    [32m✔[39m  [90mchecking for file ‘/tmp/RtmpxYnpP2/remotes766bb59517/satijalab-seurat-49a1be0/DESCRIPTION’[39m[36m[39m
    [90m─[39m[90m  [39m[90mpreparing ‘Seurat’:[39m[36m[39m
    [32m✔[39m  [90mchecking DESCRIPTION meta-information[39m[36m[39m
    [90m─[39m[90m  [39m[90mcleaning src[39m[36m[39m
    [90m─[39m[90m  [39m[90mchecking for LF line-endings in source and make files and shell scripts[39m[36m[39m
    [90m─[39m[90m  [39m[90mchecking for empty or unneeded directories[39m[36m[39m
    [90m─[39m[90m  [39m[90mlooking to see if a ‘data/datalist’ file should be added[39m[36m[39m
    [90m─[39m[90m  [39m[90mbuilding ‘Seurat_3.1.2.tar.gz’[39m[36m[39m
       


    Installing package into ‘/usr/local/lib/R/site-library’
    (as ‘lib’ is unspecified)
    



        user   system  elapsed 
    2467.462  219.185 1542.399 


The package installation took 26 minutes which is almost half the running time of the notebook.

### Install kb-python


```R
# Install kb (includes installing kallisto and bustools)
system("pip3 install kb-python", intern=TRUE)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Collecting kb-python'</li><li><span style=white-space:pre-wrap>'  Downloading https://files.pythonhosted.org/packages/62/c9/2e5b8fa2cd873a23ae1aeb128b33165d6a9387a2f56ea1fafec1d6d32477/kb_python-0.24.4-py3-none-any.whl (35.4MB)'</span></li><li>'Collecting loompy&gt;=3.0.6'</li><li><span style=white-space:pre-wrap>'  Downloading https://files.pythonhosted.org/packages/36/52/74ed37ae5988522fbf87b856c67c4f80700e6452410b4cd80498c5f416f9/loompy-3.0.6.tar.gz (41kB)'</span></li><li>'Collecting anndata&gt;=0.6.22.post1'</li><li><span style=white-space:pre-wrap>'  Downloading https://files.pythonhosted.org/packages/5b/c8/5c594a95ba293433dfe1cf188075ccbabe495bf2d291be744974aca85ffc/anndata-0.7.1-py3-none-any.whl (97kB)'</span></li><li>'Requirement already satisfied: h5py in /usr/local/lib/python3.6/dist-packages (from loompy&gt;=3.0.6-&gt;kb-python) (2.8.0)'</li><li>'Requirement already satisfied: numpy in /usr/local/lib/python3.6/dist-packages (from loompy&gt;=3.0.6-&gt;kb-python) (1.17.5)'</li><li>'Requirement already satisfied: scipy in /usr/local/lib/python3.6/dist-packages (from loompy&gt;=3.0.6-&gt;kb-python) (1.4.1)'</li><li>'Requirement already satisfied: setuptools in /usr/local/lib/python3.6/dist-packages (from loompy&gt;=3.0.6-&gt;kb-python) (45.1.0)'</li><li>'Requirement already satisfied: numba in /usr/local/lib/python3.6/dist-packages (from loompy&gt;=3.0.6-&gt;kb-python) (0.47.0)'</li><li>'Requirement already satisfied: click in /usr/local/lib/python3.6/dist-packages (from loompy&gt;=3.0.6-&gt;kb-python) (7.0)'</li><li>'Collecting numpy-groupies'</li><li><span style=white-space:pre-wrap>'  Downloading https://files.pythonhosted.org/packages/57/ae/18217b57ba3e4bb8a44ecbfc161ed065f6d1b90c75d404bd6ba8d6f024e2/numpy_groupies-0.9.10.tar.gz (43kB)'</span></li><li>'Requirement already satisfied: pandas&gt;=0.23.0 in /usr/local/lib/python3.6/dist-packages (from anndata&gt;=0.6.22.post1-&gt;kb-python) (0.25.3)'</li><li>'Requirement already satisfied: packaging in /usr/local/lib/python3.6/dist-packages (from anndata&gt;=0.6.22.post1-&gt;kb-python) (20.1)'</li><li>'Requirement already satisfied: importlib-metadata&gt;=0.7; python_version &lt; "3.8" in /usr/local/lib/python3.6/dist-packages (from anndata&gt;=0.6.22.post1-&gt;kb-python) (1.5.0)'</li><li>'Requirement already satisfied: natsort in /usr/local/lib/python3.6/dist-packages (from anndata&gt;=0.6.22.post1-&gt;kb-python) (5.5.0)'</li><li>'Requirement already satisfied: six in /usr/local/lib/python3.6/dist-packages (from h5py-&gt;loompy&gt;=3.0.6-&gt;kb-python) (1.12.0)'</li><li>'Requirement already satisfied: llvmlite&gt;=0.31.0dev0 in /usr/local/lib/python3.6/dist-packages (from numba-&gt;loompy&gt;=3.0.6-&gt;kb-python) (0.31.0)'</li><li>'Requirement already satisfied: pytz&gt;=2017.2 in /usr/local/lib/python3.6/dist-packages (from pandas&gt;=0.23.0-&gt;anndata&gt;=0.6.22.post1-&gt;kb-python) (2018.9)'</li><li>'Requirement already satisfied: python-dateutil&gt;=2.6.1 in /usr/local/lib/python3.6/dist-packages (from pandas&gt;=0.23.0-&gt;anndata&gt;=0.6.22.post1-&gt;kb-python) (2.6.1)'</li><li>'Requirement already satisfied: pyparsing&gt;=2.0.2 in /usr/local/lib/python3.6/dist-packages (from packaging-&gt;anndata&gt;=0.6.22.post1-&gt;kb-python) (2.4.6)'</li><li>'Requirement already satisfied: zipp&gt;=0.5 in /usr/local/lib/python3.6/dist-packages (from importlib-metadata&gt;=0.7; python_version &lt; "3.8"-&gt;anndata&gt;=0.6.22.post1-&gt;kb-python) (2.1.0)'</li><li>'Building wheels for collected packages: loompy, numpy-groupies'</li><li><span style=white-space:pre-wrap>'  Building wheel for loompy (setup.py): started'</span></li><li><span style=white-space:pre-wrap>'  Building wheel for loompy (setup.py): finished with status \'done\''</span></li><li><span style=white-space:pre-wrap>'  Created wheel for loompy: filename=loompy-3.0.6-cp36-none-any.whl size=47895 sha256=cb353db9cb57063d9db53d8ed434b4f4fdb15a77787b3eee9e60797c60c5cb57'</span></li><li><span style=white-space:pre-wrap>'  Stored in directory: /root/.cache/pip/wheels/f9/a4/90/5a98ad83419732b0fba533b81a2a52ba3dbe230a936ca4cdc9'</span></li><li><span style=white-space:pre-wrap>'  Building wheel for numpy-groupies (setup.py): started'</span></li><li><span style=white-space:pre-wrap>'  Building wheel for numpy-groupies (setup.py): finished with status \'done\''</span></li><li><span style=white-space:pre-wrap>'  Created wheel for numpy-groupies: filename=numpy_groupies-0+unknown-cp36-none-any.whl size=28042 sha256=db9efb43f40b5edece6fc869edb49644efb79cb201a2ae423a38223d51d8124a'</span></li><li><span style=white-space:pre-wrap>'  Stored in directory: /root/.cache/pip/wheels/30/ac/83/64d5f9293aeaec63f9539142fc629a41af064cae1b3d8d94aa'</span></li><li>'Successfully built loompy numpy-groupies'</li><li>'Installing collected packages: numpy-groupies, loompy, anndata, kb-python'</li><li>'Successfully installed anndata-0.7.1 kb-python-0.24.4 loompy-3.0.6 numpy-groupies-0+unknown'</li></ol>



### Download the data


```R
# Download the data from the 10x website
system("wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar", intern=TRUE)
system("tar -xvf pbmc_1k_v3_fastqs.tar", intern=TRUE)
```






<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'pbmc_1k_v3_fastqs/'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_I1_001.fastq.gz'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_I1_001.fastq.gz'</li></ol>



### Download an index


```R
system("kb ref -d human -i index.idx -g t2g.txt -f1 transcriptome.fasta",intern=TRUE)
```





## Pseudoalignment and counting

### Run kallisto and bustools


```R
system("kb count -i index.idx -g t2g.txt -x 10xv3 -o output --filter bustools -t 2 pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz",intern=TRUE)
```





## Basic QC


```R
library(DropletUtils)
library(Seurat)
library(Matrix)
library(tidyverse)
library(scico)
theme_set(theme_bw())
```

    Loading required package: SingleCellExperiment
    
    Loading required package: SummarizedExperiment
    
    Loading required package: GenomicRanges
    
    Loading required package: stats4
    
    Loading required package: BiocGenerics
    
    Loading required package: parallel
    
    
    Attaching package: ‘BiocGenerics’
    
    
    The following objects are masked from ‘package:parallel’:
    
        clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
        clusterExport, clusterMap, parApply, parCapply, parLapply,
        parLapplyLB, parRapply, parSapply, parSapplyLB
    
    
    The following objects are masked from ‘package:stats’:
    
        IQR, mad, sd, var, xtabs
    
    
    The following objects are masked from ‘package:base’:
    
        anyDuplicated, append, as.data.frame, basename, cbind, colnames,
        dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
        grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
        order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
        rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
        union, unique, unsplit, which, which.max, which.min
    
    
    Loading required package: S4Vectors
    
    
    Attaching package: ‘S4Vectors’
    
    
    The following object is masked from ‘package:base’:
    
        expand.grid
    
    
    Loading required package: IRanges
    
    Loading required package: GenomeInfoDb
    
    Loading required package: Biobase
    
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    
    Loading required package: DelayedArray
    
    Loading required package: matrixStats
    
    
    Attaching package: ‘matrixStats’
    
    
    The following objects are masked from ‘package:Biobase’:
    
        anyMissing, rowMedians
    
    
    Loading required package: BiocParallel
    
    
    Attaching package: ‘DelayedArray’
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges
    
    
    The following objects are masked from ‘package:base’:
    
        aperm, apply, rowsum
    
    
    
    Attaching package: ‘Seurat’
    
    
    The following object is masked from ‘package:SummarizedExperiment’:
    
        Assays
    
    
    
    Attaching package: ‘Matrix’
    
    
    The following object is masked from ‘package:S4Vectors’:
    
        expand
    
    
    ── [1mAttaching packages[22m ─────────────────────────────────────── tidyverse 1.3.0 ──
    
    [32m✔[39m [34mggplot2[39m 3.2.1     [32m✔[39m [34mpurrr  [39m 0.3.3
    [32m✔[39m [34mtibble [39m 2.1.3     [32m✔[39m [34mdplyr  [39m 0.8.4
    [32m✔[39m [34mtidyr  [39m 1.0.2     [32m✔[39m [34mstringr[39m 1.4.0
    [32m✔[39m [34mreadr  [39m 1.3.1     [32m✔[39m [34mforcats[39m 0.4.0
    
    ── [1mConflicts[22m ────────────────────────────────────────── tidyverse_conflicts() ──
    [31m✖[39m [34mdplyr[39m::[32mcollapse()[39m   masks [34mIRanges[39m::collapse()
    [31m✖[39m [34mdplyr[39m::[32mcombine()[39m    masks [34mBiobase[39m::combine(), [34mBiocGenerics[39m::combine()
    [31m✖[39m [34mdplyr[39m::[32mcount()[39m      masks [34mmatrixStats[39m::count()
    [31m✖[39m [34mdplyr[39m::[32mdesc()[39m       masks [34mIRanges[39m::desc()
    [31m✖[39m [34mtidyr[39m::[32mexpand()[39m     masks [34mMatrix[39m::expand(), [34mS4Vectors[39m::expand()
    [31m✖[39m [34mdplyr[39m::[32mfilter()[39m     masks [34mstats[39m::filter()
    [31m✖[39m [34mdplyr[39m::[32mfirst()[39m      masks [34mS4Vectors[39m::first()
    [31m✖[39m [34mdplyr[39m::[32mlag()[39m        masks [34mstats[39m::lag()
    [31m✖[39m [34mtidyr[39m::[32mpack()[39m       masks [34mMatrix[39m::pack()
    [31m✖[39m [34mggplot2[39m::[32mPosition()[39m masks [34mBiocGenerics[39m::Position(), [34mbase[39m::Position()
    [31m✖[39m [34mpurrr[39m::[32mreduce()[39m     masks [34mGenomicRanges[39m::reduce(), [34mIRanges[39m::reduce()
    [31m✖[39m [34mdplyr[39m::[32mrename()[39m     masks [34mS4Vectors[39m::rename()
    [31m✖[39m [34mpurrr[39m::[32msimplify()[39m   masks [34mDelayedArray[39m::simplify()
    [31m✖[39m [34mdplyr[39m::[32mslice()[39m      masks [34mIRanges[39m::slice()
    [31m✖[39m [34mtidyr[39m::[32munpack()[39m     masks [34mMatrix[39m::unpack()
    



```R
list.files(".", recursive = TRUE)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'index.idx'</li><li>'output/10xv3_whitelist.txt'</li><li>'output/counts_filtered/cells_x_genes.barcodes.txt'</li><li>'output/counts_filtered/cells_x_genes.genes.txt'</li><li>'output/counts_filtered/cells_x_genes.mtx'</li><li>'output/counts_unfiltered/cells_x_genes.barcodes.txt'</li><li>'output/counts_unfiltered/cells_x_genes.genes.txt'</li><li>'output/counts_unfiltered/cells_x_genes.mtx'</li><li>'output/filter_barcodes.txt'</li><li>'output/inspect.json'</li><li>'output/matrix.ec'</li><li>'output/output.bus'</li><li>'output/output.filtered.bus'</li><li>'output/output.unfiltered.bus'</li><li>'output/run_info.json'</li><li>'output/transcripts.txt'</li><li>'pbmc_1k_v3_fastqs.tar'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_I1_001.fastq.gz'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_I1_001.fastq.gz'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz'</li><li>'sample_data/anscombe.json'</li><li>'sample_data/california_housing_test.csv'</li><li>'sample_data/california_housing_train.csv'</li><li>'sample_data/mnist_test.csv'</li><li>'sample_data/mnist_train_small.csv'</li><li>'sample_data/README.md'</li><li>'t2g.txt'</li></ol>




```R
# Slightly modified from BUSpaRse, just to avoid installing a few dependencies not used here
read_count_output <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cells in rows
  ge <- ".genes.txt"
  genes <- readLines(file(paste0(dir, "/", name, ge)))
  barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}
```


```R
res_mat <- read_count_output("./output/counts_unfiltered", name = "cells_x_genes")
```


```R
dim(res_mat)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>60623</li><li>259615</li></ol>



### Test for library saturation


```R
tot_counts <- colSums(res_mat)
lib_sat <- tibble(nCount = tot_counts,
                  nGene = colSums(res_mat > 0))
```


```R
options(repr.plot.width=9, repr.plot.height=6)
ggplot(lib_sat, aes(nCount, nGene)) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_x_log10() + scale_y_log10() + annotation_logticks()
```

    Warning message:
    “Transformation introduced infinite values in continuous x-axis”
    Warning message:
    “Transformation introduced infinite values in continuous y-axis”



![png](kb_analysis_0_R_files/kb_analysis_0_R_24_1.png)


This plot is very misleading, as even the small alpha can't accurately show how many points are stacked at one location.


```R
ggplot(lib_sat, aes(nCount, nGene)) +
  geom_bin2d(bins = 50) +
  scale_fill_scico(palette = "devon", direction = -1, end = 0.95) +
  scale_x_log10() + scale_y_log10() + annotation_logticks()
```

    Warning message:
    “Transformation introduced infinite values in continuous x-axis”
    Warning message:
    “Transformation introduced infinite values in continuous y-axis”
    Warning message:
    “Removed 19583 rows containing non-finite values (stat_bin2d).”



![png](kb_analysis_0_R_files/kb_analysis_0_R_26_1.png)


Lots of points are piled at around 1 gene and 1 count, and correspond to empty or near empty droplets.

### Examine the knee plot

The "knee plot" was introduced in the Drop-seq paper: 
- Macosko et al., [Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets](https://www.cell.com/fulltext/S0092-8674(15)00549-8), 2015. DOI:10.1016/j.cell.2015.05.002

In this plot cells are ordered by the number of UMI counts associated to them (shown on the *x*-axis), and the fraction of droplets with at least that number of cells is shown on the *y*-axis:


```R
summary(tot_counts)
```


        Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
        0.00     1.00     1.00    43.64     6.00 60120.00 



```R
bc_rank <- barcodeRanks(res_mat, lower = 1000)
```


```R
#' Knee plot for filtering empty droplets
#' 
#' Visualizes the inflection point to filter empty droplets. This function plots 
#' different datasets with a different color. Facets can be added after calling
#' this function with `facet_*` functions. Will be added to the next release
#' version of BUSpaRse.
#' 
#' @param bc_rank A `DataFrame` output from `DropletUtil::barcodeRanks`.
#' @return A ggplot2 object.
knee_plot <- function(bc_rank) {
  knee_plt <- tibble(rank = bc_rank[["rank"]],
                     total = bc_rank[["total"]]) %>% 
    distinct() %>% 
    dplyr::filter(total > 0)
  annot <- tibble(inflection = metadata(bc_rank)[["inflection"]],
                  rank_cutoff = max(bc_rank$rank[bc_rank$total > metadata(bc_rank)[["inflection"]]]))
  p <- ggplot(knee_plt, aes(total, rank)) +
    geom_line() +
    geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2) +
    geom_vline(aes(xintercept = inflection), data = annot, linetype = 2) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    labs(y = "Rank", x = "Total UMIs")
  return(p)
}
```


```R
options(repr.plot.width=9, repr.plot.height=6)
knee_plot(bc_rank)
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_32_0.png)


## Analysis
We begin by asking for genes with the highest proportions in droplets (prior to filtering out empty droplets).


```R
tr2g <- read_tsv("t2g.txt", col_names = c("transcript", "gene", "gene_name"))
tr2g <- distinct(tr2g[, c("gene", "gene_name")])
```

    Parsed with column specification:
    cols(
      transcript = [31mcol_character()[39m,
      gene = [31mcol_character()[39m,
      gene_name = [31mcol_character()[39m
    )
    



```R
plot_pct_genes <- function(mat, tr2g, top_n = 20, symbol = "ensembl") {
  pct_tx <- rowSums(mat)
  gs <- rownames(mat)[order(-pct_tx)]
  df <- as.data.frame(t(mat[gs[1:20],]))
  df <- df %>%
    mutate_all(function(x) x/colSums(mat)) %>%
    pivot_longer(everything(), names_to = "gene")
  if (symbol == "ensembl") {
    df <- left_join(df, tr2g, by = "gene")
  } else {
    df <- rename(df, gene_name = gene)
  }
    df %>%
    mutate(gene = fct_reorder(gene_name, value, .fun = median)) %>%
    ggplot(aes(gene, value)) +
    geom_boxplot() +
    labs(x = "", y = "Proportion of total counts") +
    coord_flip()
}
```


```R
options(repr.plot.width=6, repr.plot.height=10)
plot_pct_genes(res_mat, tr2g)
```

    Warning message:
    “Removed 391660 rows containing non-finite values (stat_boxplot).”



![png](kb_analysis_0_R_files/kb_analysis_0_R_36_1.png)


For many barcodes, the top genes by proportion of all counts are ribosomal or mitochondrial genes. Also, the proportions plotted above seem to have some discrete values; this effect is a result of computing fractions with small denominator, which happens when droplets produce very few UMI counts.

### Filter


```R
res_mat <- res_mat[, tot_counts > metadata(bc_rank)$inflection]
res_mat <- res_mat[Matrix::rowSums(res_mat) > 0,]
dim(res_mat)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>31832</li><li>1178</li></ol>




```R
# Convert from Ensembl gene ID to gene symbol
rownames(res_mat) <- tr2g$gene_name[match(rownames(res_mat), tr2g$gene)]
```


```R
(pbmc <- CreateSeuratObject(counts = res_mat, project = "pbmc1k", min.cells = 3, min.features = 200))
```

    Warning message:
    “Non-unique features (rownames) present in the input matrix, making unique”
    Warning message:
    “Feature names cannot have underscores ('_'), replacing with dashes ('-')”



    An object of class Seurat 
    25950 features across 1173 samples within 1 assay 
    Active assay: RNA (25950 features)


The steps below constitute a standard analysis worklow for single-cell RNA-seq data. 


```R
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```

The number of unique genes and total molecules are automatically calculated when running the `CreateSeuratObject` command.
The associated data is stored in the object metadata.



```R
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)
```


<table>
<caption>A data.frame: 5 × 4</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_RNA</th><th scope=col>nFeature_RNA</th><th scope=col>percent.mt</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACCCAAGGAGAGTA</th><td>pbmc1k</td><td>9289</td><td>3198</td><td>11.271396</td></tr>
	<tr><th scope=row>AAACGCTTCAGCCCAG</th><td>pbmc1k</td><td>6483</td><td>2513</td><td> 8.252352</td></tr>
	<tr><th scope=row>AAAGAACAGACGACTG</th><td>pbmc1k</td><td>5011</td><td>2082</td><td> 6.166434</td></tr>
	<tr><th scope=row>AAAGAACCAATGGCAG</th><td>pbmc1k</td><td>3264</td><td>1555</td><td> 6.893382</td></tr>
	<tr><th scope=row>AAAGAACGTCTGCAAT</th><td>pbmc1k</td><td>7488</td><td>2508</td><td> 6.610577</td></tr>
</tbody>
</table>



Next, we visualize some QC metrics and use the results to set filtering criteria.


```R
# Visualize QC metrics as a violin plot
options(repr.plot.width=12, repr.plot.height=6)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_47_0.png)



```R
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_48_0.png)



```R
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 20)
```

### Normalize

After removing unwanted cells from the dataset, the next step is to normalize the data. A standard choice is `LogNormalize` which normalizes the UMI counts for each cell by the total counts, multiplies this by a scale factor (10,000 by default), and finally log-transforms the result. Normalized values are stored in pbmc[["RNA"]]@data. 

We recommend the preprint 
- Breda, J., Zavolan, M. and van Nimwegen, E. Bayesian inference of the gene expression states of single cells from scRNA-seq data. bioRxiv (2019). doi.org/10.1101/2019.12.28.889956 

for a thorough discussion of normalization.


```R
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```

For clarity, in this previous line of code (and in future commands), we provide the default values for certain parameters in the function call. However, this isn’t required and the same behavior can be achieved with:


```R
# pbmc <- NormalizeData(pbmc)
```

### Highly expressed genes


To identify a subset of genes that exhibit high cell-to-cell variation in the dataset we apply a procedure implemented in the `FindVariableFeatures` function. By default, it returns 2,000 genes per dataset. These will be used in downstream analysis.

Seurat documentation describes the method used to find highly variable genes here as such:

> First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess). Then standardizes the feature values using the observed mean and expected variance (given by the fitted line). Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).


```R
options(repr.plot.width=6, repr.plot.height=10)
plot_pct_genes(GetAssayData(pbmc, slot = "counts"), tr2g, symbol = "symbol")
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_57_0.png)



```R
options(repr.plot.width=9, repr.plot.height=6)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc, log = FALSE)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
```

    Warning message:
    “Using `as.character()` on a quosure is deprecated as of rlang 0.3.0.
    Please use `as_label()` or `as_name()` instead.
    [90mThis warning is displayed once per session.[39m”
    When using repel, set xnudge and ynudge to 0 for optimal results
    



![png](kb_analysis_0_R_files/kb_analysis_0_R_58_1.png)


### Scaling the data
Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData function shifts the expression of each gene, so that the mean expression across cells is 0 and the variance across cells is 1
This step gives equal weight to genes in downstream analyses, so that highly-expressed genes do not dominate. The results of this are stored in pbmc[["RNA"]]@scale.data


```R
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

    Centering and scaling data matrix
    


We apply this only to the genes identified as highly variable:


```R
# pbmc <- ScaleData(pbmc)

```

The scaling does not affect PCA or clustering results. However, Seurat heatmaps (produced as shown below with DoHeatmap) require genes in the heatmap to be scaled so that highly-expressed genes don’t dominate. To make sure we don’t leave any genes out of the heatmap later, we are scaling all genes in this tutorial.



In Seurat v2 we also use the ScaleData function to remove unwanted sources of variation from a single-cell dataset. For example, we could ‘regress out’ heterogeneity associated with (for example) cell cycle stage, or mitochondrial contamination. These features are still supported in ScaleData in Seurat v3, i.e.:


```R
# pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
```

### Principal component analysis

Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input.


```R
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
```

    PC_ 1 
    Positive:  S100A9, FCN1, MNDA, CST3, FGL2, LYZ, S100A8, CTSS, NCF2, SERPINA1 
    	   PSAP, AIF1, TYMP, VCAN, CSTA, KLF4, GRN, MPEG1, CPVL, MS4A6A 
    	   CLEC7A, LST1, TNFAIP2, FCER1G, CD14, CYBB, LGALS1, S100A12, TYROBP, CD36 
    Negative:  LTB, TRAC, CD3E, TRBC2, CD3D, IL32, BCL11B, CD3G, IL7R, TCF7 
    	   CD69, ISG20, CD247, CD27, SPOCK2, ARL4C, CD7, CD2, GZMM, TRBC1 
    	   CD6, PRKCQ-AS1, NOSIP, AC058791.1, RORA, CTSW, CCR7, AQP3, PEBP1, ITM2A 
    PC_ 2 
    Positive:  CD3E, IL32, CD247, GZMM, CD7, CTSW, CD3D, S100A4, GZMA, NKG7 
    	   ANXA1, TRAC, BCL11B, CD3G, IL7R, PRF1, CST7, KLRB1, ARL4C, SAMD3 
    	   CCL5, CD2, TRBC1, KLRG1, A2M, ITGB2, MT2A, RORA, ACTG1, TCF7 
    Negative:  CD79A, MS4A1, IGHM, BANK1, BCL11A, LINC00926, CD79B, TNFRSF13C, IGHD, CD74 
    	   HLA-DQB1, CD22, HLA-DQA1, HLA-DRB1, HLA-DRA, HLA-DPA1, HLA-DPB1, TCL1A, FCER2, AFF3 
    	   PAX5, IGKC, VPREB3, SPIB, MEF2C, RALGPS2, HVCN1, FCRL1, CD40, HLA-DOB 
    PC_ 3 
    Positive:  PRKAR2B, CAVIN2, GNG11, PPBP, PF4, PTCRA, TUBB1, GP9, LY6G6F, ITGA2B 
    	   CLU, TREML1, CD9, C2orf88, MED12L, LINC00989, CTTN, GMPR, MTURN, CMTM5 
    	   ESAM, TMEM158, CLDN5, HIST1H2AC, TNNC2, MEIS1, SPARC, NAT8B, GP1BA, CXCR2P1 
    Negative:  FOS, DUSP1, VIM, NEAT1, EVI2B, ZFP36L1, CYBA, MCL1, AC020916.1, KLF6 
    	   NCF1C, LCP1, LTB, HNRNPU, S100A10, LTA4H, RHOB, ITGB2, RBP7, APLP2 
    	   ATP2B1-AS1, SCPEP1, TKT, NFKBIA, S100A6, CTSS, VMP1, TSPO, MEGF9, CASP1 
    PC_ 4 
    Positive:  LEF1, TCF7, IL7R, MAL, CCR7, BCL11B, CD3D, NOSIP, LTB, TRAC 
    	   CD3G, RGS10, CAMK4, PASK, NELL2, CD27, RGCC, EGR1, SLC2A3, BEX3 
    	   CD40LG, FHIT, CD6, INPP4B, TRAT1, ADTRP, NOG, VIM, TSHZ2, PRKCQ-AS1 
    Negative:  GZMB, GNLY, NKG7, CLIC3, KLRF1, PRF1, CST7, SPON2, FGFBP2, KLRD1 
    	   GZMA, ADGRG1, CCL4, TRDC, HOPX, MATK, IL2RB, APOBEC3G, CTSW, TTC38 
    	   TBX21, RHOC, PTGDR, S1PR5, FCGR3A, SH2D1B, C12orf75, MYOM2, CMC1, GZMH 
    PC_ 5 
    Positive:  GNLY, FGFBP2, KLRF1, PRF1, NKG7, CST7, CCL4, KLRD1, MS4A1, IGHD 
    	   LINC00926, ADGRG1, CD79B, CD79A, TRDC, CD22, GZMA, TBX21, SPON2, TNFRSF13C 
    	   MATK, FCER2, IL2RB, MYOM2, PAX5, HOPX, S1PR5, SH2D1B, TTC38, BANK1 
    Negative:  LILRA4, SCT, PACSIN1, SMPD3, LRRC26, SERPINF1, TPM2, AL096865.1, IL3RA, DNASE1L3 
    	   TNFRSF21, CLEC4C, PLD4, CUX2, ITM2C, GAS6, MYBL2, EPHA2, PPP1R14B, UGCG 
    	   PPP1R14B-AS1, LAMP5, CUEDC1, RUNX2, PPM1J, SERPINF2, NRP1, DERL3, CIB2, LINC02812 
    


Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction, DimPlot, and DimHeatmap




```R
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
```

    PC_ 1 
    Positive:  S100A9, FCN1, MNDA, CST3, FGL2 
    Negative:  LTB, TRAC, CD3E, TRBC2, CD3D 
    PC_ 2 
    Positive:  CD3E, IL32, CD247, GZMM, CD7 
    Negative:  CD79A, MS4A1, IGHM, BANK1, BCL11A 
    PC_ 3 
    Positive:  PRKAR2B, CAVIN2, GNG11, PPBP, PF4 
    Negative:  FOS, DUSP1, VIM, NEAT1, EVI2B 
    PC_ 4 
    Positive:  LEF1, TCF7, IL7R, MAL, CCR7 
    Negative:  GZMB, GNLY, NKG7, CLIC3, KLRF1 
    PC_ 5 
    Positive:  GNLY, FGFBP2, KLRF1, PRF1, NKG7 
    Negative:  LILRA4, SCT, PACSIN1, SMPD3, LRRC26 


Which genes are contributing the most to the first 2 PCs?


```R
options(repr.plot.width=6, repr.plot.height=8)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_71_0.png)



```R
options(repr.plot.width=7, repr.plot.height=6)
FeaturePlot(pbmc, reduction = "pca", feature = "CST3")
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_72_0.png)


### Determining dimensionality
To overcome the extensive technical noise in any single feature for scRNA-seq data, one can cluster cells based on their PCA projections, with each PC essentially representing a ‘metafeature’ that combines information across a correlated feature set.

A common heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one (ElbowPlot function). In this example, we can observe an ‘elbow’ around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs.


```R
options(repr.plot.width=9, repr.plot.height=6)
ElbowPlot(pbmc)
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_75_0.png)


### The neighborhood graph

We cluster cells using the Louvain algorithm (a default in Seurat), which iteratively group cells together, with the goal of optimizing the standard modularity function. The FindClusters function implements this procedure, and contains a resolution parameter that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters. We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. The clusters can be found using the Idents function.


```R
pbmc <- FindNeighbors(pbmc, dims = 1:10, k.param = 20)
pbmc <- FindClusters(pbmc, resolution = 0.6)
```

    Computing nearest neighbor graph
    
    Computing SNN
    


    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 1113
    Number of edges: 36080
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8763
    Number of communities: 9
    Elapsed time: 0 seconds



```R
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>AAACCCAAGGAGAGTA</dt><dd>1</dd><dt>AAACGCTTCAGCCCAG</dt><dd>4</dd><dt>AAAGAACAGACGACTG</dt><dd>5</dd><dt>AAAGAACCAATGGCAG</dt><dd>5</dd><dt>AAAGAACGTCTGCAAT</dt><dd>0</dd></dl>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'0'</li><li>'1'</li><li>'2'</li><li>'3'</li><li>'4'</li><li>'5'</li><li>'6'</li><li>'7'</li><li>'8'</li></ol>
</details>


### UMAP and t-SNE
tSNE and UMAP can be used to visualize and explore non-linear aspects of high-dimensional data. Here we apply these methods to the PC projection of the data (with same dimension as used for clustering).

[UMAP](https://en.wikipedia.org/wiki/Nonlinear_dimensionality_reduction) (UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction) is a manifold learning technique that can also be used to visualize cells. It was published in:

- McInnes, Leland, John Healy, and James Melville. "Umap: Uniform manifold approximation and projection for dimension reduction." arXiv preprint arXiv:1802.03426 (2018).

[t-SNE](https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding) is a non-linear dimensionality reduction technique described in:

- Maaten, Laurens van der, and Geoffrey Hinton. "Visualizing data using t-SNE." Journal of machine learning research 9.Nov (2008): 2579-2605.




```R
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
```

    Warning message:
    “The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    This message will be shown once per session”



```R
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_82_0.png)



```R
options(repr.plot.width=16, repr.plot.height=5)
FeaturePlot(pbmc, reduction = "umap", features = c("CST3", "NKG7", "PPBP"),
ncol = 3)
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_83_0.png)


### Finding differentially expressed features (cluster biomarkers)
A key follow-up step to clustering cells is to find gene markers that are associated with them. We used Seurat's FindAllMarkers function which automates the process for all clusters.

The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups. You can set both of these to 0, but with a dramatic increase in time - since this will test a large number of features that are unlikely to be highly discriminatory. As another option to speed up these computations, max.cells.per.ident can be set. This will downsample each identity class to have no more cells than whatever this is set to. While there is generally going to be a loss in power, the speed increases can be significiant and the most highly differentially expressed features will likely still rise to the top.


```R
# Scanpy style gene rank plot
plot_gene_rank <- function(markers, n) {
  df_plot <- markers %>%
    group_by(cluster) %>%
    top_n(25, avg_logFC) %>%
    mutate(rank = factor(row_number(desc(avg_logFC))))
  ggplot(df_plot, aes(rank, avg_logFC)) +
    geom_text(aes(label = gene), angle = -90, hjust = 1) +
    facet_wrap(~ cluster) +
    scale_y_continuous(expand = expand_scale(mult = c(0.05, 0.25)))
}
```

Several methods for differential expression are supported by Seurat. The default is Wilcoxon rank sum test.


```R
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, test.use = "wilcox", only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)
```

    Calculating cluster 0
    
    Calculating cluster 1
    
    Calculating cluster 2
    
    Calculating cluster 3
    
    Calculating cluster 4
    
    Calculating cluster 5
    
    Calculating cluster 6
    
    Calculating cluster 7
    
    Calculating cluster 8
    



```R
options(repr.plot.width=12, repr.plot.height=8)
plot_gene_rank(pbmc.markers, 25)
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_88_0.png)


Student's t test is also supported


```R
pbmc.markers.t <- FindAllMarkers(pbmc, test.use = "t", only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)
```

    Calculating cluster 0
    
    Calculating cluster 1
    
    Calculating cluster 2
    
    Calculating cluster 3
    
    Calculating cluster 4
    
    Calculating cluster 5
    
    Calculating cluster 6
    
    Calculating cluster 7
    
    Calculating cluster 8
    



```R
plot_gene_rank(pbmc.markers.t, 25)
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_91_0.png)


Also logistic regression to test how good each gene is for deciding whether a cell is in a cluster.


```R
pbmc.markers.lr <- FindAllMarkers(pbmc, test.use = "LR", only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)
```

    Calculating cluster 0
    
    Calculating cluster 1
    
    Calculating cluster 2
    
    Calculating cluster 3
    
    Calculating cluster 4
    
    Calculating cluster 5
    
    Calculating cluster 6
    
    Calculating cluster 7
    
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Calculating cluster 8
    
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”



```R
plot_gene_rank(pbmc.markers.lr, 25)
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_94_0.png)


Seurat includes several tools for visualizing marker expression. VlnPlot (shows expression probability distributions across clusters), and FeaturePlot (visualizes feature expression on a tSNE or PCA plot) are our most commonly used visualizations. We also suggest exploring RidgePlot, CellScatter, and DotPlot as additional methods to view your dataset.




```R
marker_genes <- sort(c('IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',  
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP', 'CCR7',
                'S100A4'))
```


```R
options(repr.plot.width=16, repr.plot.height=20)
VlnPlot(pbmc, features = marker_genes, ncol = 4)
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_97_0.png)



```R
options(repr.plot.width=16, repr.plot.height=20)
FeaturePlot(pbmc, features = marker_genes, ncol = 4)
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_98_0.png)


### Assigning cell type identity to clusters
In this dataset, we can use canonical markers to easily match the *de novo* clustering to known cell types:

Cluster ID | Markers | Cell Type
-----------|---------|-------------
0	| IL7R, S100A4 |	Memory CD4+ T
1 |	CD14, LYZ |	CD14+ Mono
2 |	IL7R, CCR7 |	Naive CD4+
3 |	FCGR3A, MS4A7 |	FCGR3A+ Mono
4 |	MS4A1, CD79A |	B
5 |	GNLY, NKG7 | NK
6 |	CD8A | CD8+ T
7 | MS4A1, CD79A | B
8 |	PPBP | Platelet


```R
options(repr.plot.width=6, repr.plot.height=7)
DotPlot(pbmc, assay = "RNA", features = marker_genes, scale.by = "size") +
  coord_flip()
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_101_0.png)



```R
options(repr.plot.width=9, repr.plot.height=6)
new.cluster.ids <- c("Memory CD4 T", "CD14+ Mono", "Naive CD4 T", "FCGR3A+ Mono", 
    "B1", "NK", "CD8 T", "B2", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6) + NoLegend()
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_102_0.png)



```R
Sys.time() - start_time
```


    Time difference of 54.48928 mins



```R
sessionInfo()
```


    R version 3.6.2 (2019-12-12)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: Ubuntu 18.04.3 LTS
    
    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
    LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so
    
    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
     [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
     [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
     [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    
    attached base packages:
    [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    [8] methods   base     
    
    other attached packages:
     [1] scico_1.1.0                 forcats_0.4.0              
     [3] stringr_1.4.0               dplyr_0.8.4                
     [5] purrr_0.3.3                 readr_1.3.1                
     [7] tidyr_1.0.2                 tibble_2.1.3               
     [9] ggplot2_3.2.1               tidyverse_1.3.0            
    [11] Matrix_1.2-18               Seurat_3.1.2               
    [13] DropletUtils_1.6.1          SingleCellExperiment_1.8.0 
    [15] SummarizedExperiment_1.16.1 DelayedArray_0.12.2        
    [17] BiocParallel_1.20.1         matrixStats_0.55.0         
    [19] Biobase_2.46.0              GenomicRanges_1.38.0       
    [21] GenomeInfoDb_1.22.0         IRanges_2.20.2             
    [23] S4Vectors_0.24.3            BiocGenerics_0.32.0        
    
    loaded via a namespace (and not attached):
      [1] reticulate_1.14        R.utils_2.9.2          tidyselect_1.0.0      
      [4] htmlwidgets_1.5.1      grid_3.6.2             Rtsne_0.15            
      [7] devtools_2.2.1         munsell_0.5.0          codetools_0.2-16      
     [10] mutoss_0.1-12          ica_1.0-2              pbdZMQ_0.3-3          
     [13] future_1.16.0          withr_2.1.2            colorspace_1.4-1      
     [16] uuid_0.1-2             rstudioapi_0.10        ROCR_1.0-7            
     [19] gbRd_0.4-11            listenv_0.8.0          labeling_0.3          
     [22] Rdpack_0.11-1          repr_1.1.0             GenomeInfoDbData_1.2.2
     [25] mnormt_1.5-6           farver_2.0.3           rhdf5_2.30.1          
     [28] rprojroot_1.2          vctrs_0.2.2            generics_0.0.2        
     [31] TH.data_1.0-10         R6_2.4.1               rsvd_1.0.2            
     [34] locfit_1.5-9.1         bitops_1.0-6           assertthat_0.2.1      
     [37] SDMTools_1.1-221.2     scales_1.1.0           multcomp_1.4-12       
     [40] gtable_0.3.0           npsurv_0.4-0           globals_0.12.5        
     [43] processx_3.4.1         sandwich_2.5-1         rlang_0.4.4           
     [46] splines_3.6.2          lazyeval_0.2.2         broom_0.5.4           
     [49] BiocManager_1.30.10    reshape2_1.4.3         modelr_0.1.5          
     [52] backports_1.1.5        tools_3.6.2            usethis_1.5.1         
     [55] ellipsis_0.3.0         gplots_3.0.1.2         RColorBrewer_1.1-2    
     [58] sessioninfo_1.1.1      ggridges_0.5.2         TFisher_0.2.0         
     [61] Rcpp_1.0.3             plyr_1.8.5             base64enc_0.1-3       
     [64] zlibbioc_1.32.0        RCurl_1.98-1.1         ps_1.3.0              
     [67] prettyunits_1.1.1      pbapply_1.4-2          cowplot_1.0.0         
     [70] zoo_1.8-7              haven_2.2.0            ggrepel_0.8.1         
     [73] cluster_2.1.0          fs_1.3.1               magrittr_1.5          
     [76] RSpectra_0.16-0        data.table_1.12.8      lmtest_0.9-37         
     [79] reprex_0.3.0           RANN_2.6.1             mvtnorm_1.0-12        
     [82] fitdistrplus_1.0-14    pkgload_1.0.2          hms_0.5.3             
     [85] lsei_1.2-0             evaluate_0.14          readxl_1.3.1          
     [88] gridExtra_2.3          testthat_2.3.1         compiler_3.6.2        
     [91] KernSmooth_2.23-16     crayon_1.3.4           R.oo_1.23.0           
     [94] htmltools_0.4.0        RcppParallel_4.4.4     lubridate_1.7.4       
     [97] DBI_1.1.0              dbplyr_1.4.2           MASS_7.3-51.5         
    [100] rappdirs_0.3.1         cli_2.0.1              R.methodsS3_1.7.1     
    [103] gdata_2.18.0           metap_1.3              igraph_1.2.4.2        
    [106] pkgconfig_2.0.3        sn_1.5-5               numDeriv_2016.8-1.1   
    [109] IRdisplay_0.7.0        plotly_4.9.1           xml2_1.2.2            
    [112] dqrng_0.2.1            multtest_2.42.0        XVector_0.26.0        
    [115] rvest_0.3.5            bibtex_0.4.2.2         callr_3.4.1           
    [118] digest_0.6.23          sctransform_0.2.1      RcppAnnoy_0.0.14      
    [121] tsne_0.1-3             cellranger_1.1.0       leiden_0.3.3          
    [124] uwot_0.1.5             edgeR_3.28.0           curl_4.2              
    [127] gtools_3.8.1           lifecycle_0.1.0        nlme_3.1-143          
    [130] jsonlite_1.6           Rhdf5lib_1.8.0         desc_1.2.0            
    [133] viridisLite_0.3.0      limma_3.42.2           fansi_0.4.1           
    [136] pillar_1.4.3           lattice_0.20-38        httr_1.4.1            
    [139] plotrix_3.7-7          pkgbuild_1.0.6         survival_3.1-8        
    [142] glue_1.3.1             remotes_2.1.0          png_0.1-7             
    [145] stringi_1.4.5          HDF5Array_1.14.2       caTools_1.18.0        
    [148] memoise_1.1.0          IRkernel_1.1           irlba_2.3.3           
    [151] future.apply_1.4.0     ape_5.3               


**Feedback**: please report any issues, or submit pull requests for improvements, in the [Github repository where this notebook is located](https://github.com/pachterlab/kallistobustools/blob/master/notebooks/kb_analysis_0_R.ipynb).
