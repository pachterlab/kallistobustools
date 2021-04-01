<a href="https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/tutorials/docs/tutorials/kb_building_atlas/R/kb_analysis_0_R.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

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
A large fraction of the running time of this notebook is in installing the Seurat R package, since it has lots of dependencies and many of them use Rcpp which results in the need to compile lots of C++ code. Compilation is required because CRAN does not distribute binaries for Linux, which is the operating system here.


```R
system.time({
  install.packages("Seurat", Ncpus = 2)
})
```

    Installing package into ‘/usr/local/lib/R/site-library’
    (as ‘lib’ is unspecified)
    
    also installing the dependencies ‘bitops’, ‘gtools’, ‘caTools’, ‘sass’, ‘jquerylib’, ‘sitmo’, ‘globals’, ‘listenv’, ‘parallelly’, ‘plyr’, ‘zoo’, ‘data.table’, ‘gplots’, ‘reshape2’, ‘gridExtra’, ‘RcppArmadillo’, ‘httpuv’, ‘xtable’, ‘sourcetools’, ‘bslib’, ‘spatstat.data’, ‘spatstat.utils’, ‘spatstat.sparse’, ‘abind’, ‘tensor’, ‘goftest’, ‘deldir’, ‘polyclip’, ‘FNN’, ‘RSpectra’, ‘dqrng’, ‘cowplot’, ‘fitdistrplus’, ‘future’, ‘future.apply’, ‘ggrepel’, ‘ggridges’, ‘ica’, ‘igraph’, ‘irlba’, ‘leiden’, ‘lmtest’, ‘matrixStats’, ‘miniUI’, ‘patchwork’, ‘pbapply’, ‘plotly’, ‘png’, ‘RANN’, ‘RcppAnnoy’, ‘reticulate’, ‘ROCR’, ‘Rtsne’, ‘scattermore’, ‘sctransform’, ‘SeuratObject’, ‘shiny’, ‘spatstat.core’, ‘spatstat.geom’, ‘uwot’, ‘RcppEigen’, ‘RcppProgress’
    
    



        user   system  elapsed 
    1923.767  192.391 1168.018 


The package installation took 20 minutes (elapsed) which is 40% of the running time of the entire notebook. The user time is nearly twice the elapsed time here because 2 cores were used to install all those packages.

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
<ol class=list-inline><li>'Collecting kb-python'</li><li><span style=white-space:pre-wrap>'  Downloading https://files.pythonhosted.org/packages/29/94/855ed1c11110a65a466cd95a6fef64958bad055f2678270b80a32e42cdb1/kb_python-0.25.1-py3-none-any.whl (59.1MB)'</span></li><li>'Requirement already satisfied: tqdm&gt;=4.39.0 in /usr/local/lib/python3.7/dist-packages (from kb-python) (4.41.1)'</li><li>'Collecting loompy&gt;=3.0.6'</li><li><span style=white-space:pre-wrap>'  Downloading https://files.pythonhosted.org/packages/36/52/74ed37ae5988522fbf87b856c67c4f80700e6452410b4cd80498c5f416f9/loompy-3.0.6.tar.gz (41kB)'</span></li><li>'Requirement already satisfied: nbformat&gt;=4.4.0 in /usr/local/lib/python3.7/dist-packages (from kb-python) (5.1.2)'</li><li>'Requirement already satisfied: h5py&gt;=2.10.0 in /usr/local/lib/python3.7/dist-packages (from kb-python) (2.10.0)'</li><li>'Requirement already satisfied: requests&gt;=2.19.0 in /usr/local/lib/python3.7/dist-packages (from kb-python) (2.23.0)'</li><li>'Requirement already satisfied: nbconvert&gt;=5.6.0 in /usr/local/lib/python3.7/dist-packages (from kb-python) (5.6.1)'</li><li>'Requirement already satisfied: numpy&gt;=1.17.2 in /usr/local/lib/python3.7/dist-packages (from kb-python) (1.19.5)'</li><li>'Collecting anndata&gt;=0.6.22.post1'</li><li><span style=white-space:pre-wrap>'  Downloading https://files.pythonhosted.org/packages/81/b1/743cc79f89d9db6dccbfb7e6000795acb218a6c6320b7a2337cad99bd047/anndata-0.7.5-py3-none-any.whl (119kB)'</span></li><li>'Collecting scanpy&gt;=1.4.4.post1'</li><li><span style=white-space:pre-wrap>'  Downloading https://files.pythonhosted.org/packages/c1/f4/a7848e6f990cf5bcbedade93702baf3e99ae704714563fe9bdceb3d597c7/scanpy-1.7.1-py3-none-any.whl (10.3MB)'</span></li><li>'Requirement already satisfied: Jinja2&gt;2.10.1 in /usr/local/lib/python3.7/dist-packages (from kb-python) (2.11.3)'</li><li>'Collecting plotly&gt;=4.5.0'</li><li><span style=white-space:pre-wrap>'  Downloading https://files.pythonhosted.org/packages/1f/f6/bd3c17c8003b6641df1228e80e1acac97ed8402635e46c2571f8e1ef63af/plotly-4.14.3-py2.py3-none-any.whl (13.2MB)'</span></li><li>'Requirement already satisfied: scikit-learn&gt;=0.21.3 in /usr/local/lib/python3.7/dist-packages (from kb-python) (0.22.2.post1)'</li><li>'Requirement already satisfied: scipy in /usr/local/lib/python3.7/dist-packages (from loompy&gt;=3.0.6-&gt;kb-python) (1.4.1)'</li><li>'Requirement already satisfied: setuptools in /usr/local/lib/python3.7/dist-packages (from loompy&gt;=3.0.6-&gt;kb-python) (54.2.0)'</li><li>'Requirement already satisfied: numba in /usr/local/lib/python3.7/dist-packages (from loompy&gt;=3.0.6-&gt;kb-python) (0.51.2)'</li><li>'Requirement already satisfied: click in /usr/local/lib/python3.7/dist-packages (from loompy&gt;=3.0.6-&gt;kb-python) (7.1.2)'</li><li>'Collecting numpy-groupies'</li><li><span style=white-space:pre-wrap>'  Downloading https://files.pythonhosted.org/packages/99/eb/fc72b507219957cffdf2c5952e396cc04a30c2223e2fd789f4a744ffc52f/numpy_groupies-0.9.13.tar.gz (109kB)'</span></li><li>'Requirement already satisfied: jsonschema!=2.5.0,&gt;=2.4 in /usr/local/lib/python3.7/dist-packages (from nbformat&gt;=4.4.0-&gt;kb-python) (2.6.0)'</li><li>'Requirement already satisfied: traitlets&gt;=4.1 in /usr/local/lib/python3.7/dist-packages (from nbformat&gt;=4.4.0-&gt;kb-python) (5.0.5)'</li><li>'Requirement already satisfied: jupyter-core in /usr/local/lib/python3.7/dist-packages (from nbformat&gt;=4.4.0-&gt;kb-python) (4.7.1)'</li><li>'Requirement already satisfied: ipython-genutils in /usr/local/lib/python3.7/dist-packages (from nbformat&gt;=4.4.0-&gt;kb-python) (0.2.0)'</li><li>'Requirement already satisfied: six in /usr/local/lib/python3.7/dist-packages (from h5py&gt;=2.10.0-&gt;kb-python) (1.15.0)'</li><li>'Requirement already satisfied: urllib3!=1.25.0,!=1.25.1,&lt;1.26,&gt;=1.21.1 in /usr/local/lib/python3.7/dist-packages (from requests&gt;=2.19.0-&gt;kb-python) (1.24.3)'</li><li>'Requirement already satisfied: idna&lt;3,&gt;=2.5 in /usr/local/lib/python3.7/dist-packages (from requests&gt;=2.19.0-&gt;kb-python) (2.10)'</li><li>'Requirement already satisfied: certifi&gt;=2017.4.17 in /usr/local/lib/python3.7/dist-packages (from requests&gt;=2.19.0-&gt;kb-python) (2020.12.5)'</li><li>'Requirement already satisfied: chardet&lt;4,&gt;=3.0.2 in /usr/local/lib/python3.7/dist-packages (from requests&gt;=2.19.0-&gt;kb-python) (3.0.4)'</li><li>'Requirement already satisfied: entrypoints&gt;=0.2.2 in /usr/local/lib/python3.7/dist-packages (from nbconvert&gt;=5.6.0-&gt;kb-python) (0.3)'</li><li>'Requirement already satisfied: pygments in /usr/local/lib/python3.7/dist-packages (from nbconvert&gt;=5.6.0-&gt;kb-python) (2.6.1)'</li><li>'Requirement already satisfied: bleach in /usr/local/lib/python3.7/dist-packages (from nbconvert&gt;=5.6.0-&gt;kb-python) (3.3.0)'</li><li>'Requirement already satisfied: pandocfilters&gt;=1.4.1 in /usr/local/lib/python3.7/dist-packages (from nbconvert&gt;=5.6.0-&gt;kb-python) (1.4.3)'</li><li>'Requirement already satisfied: testpath in /usr/local/lib/python3.7/dist-packages (from nbconvert&gt;=5.6.0-&gt;kb-python) (0.4.4)'</li><li>'Requirement already satisfied: mistune&lt;2,&gt;=0.8.1 in /usr/local/lib/python3.7/dist-packages (from nbconvert&gt;=5.6.0-&gt;kb-python) (0.8.4)'</li><li>'Requirement already satisfied: defusedxml in /usr/local/lib/python3.7/dist-packages (from nbconvert&gt;=5.6.0-&gt;kb-python) (0.7.1)'</li><li>'Requirement already satisfied: packaging in /usr/local/lib/python3.7/dist-packages (from anndata&gt;=0.6.22.post1-&gt;kb-python) (20.9)'</li><li>'Requirement already satisfied: importlib-metadata&gt;=0.7; python_version &lt; "3.8" in /usr/local/lib/python3.7/dist-packages (from anndata&gt;=0.6.22.post1-&gt;kb-python) (3.8.1)'</li><li>'Requirement already satisfied: pandas!=1.1,&gt;=1.0 in /usr/local/lib/python3.7/dist-packages (from anndata&gt;=0.6.22.post1-&gt;kb-python) (1.1.5)'</li><li>'Requirement already satisfied: natsort in /usr/local/lib/python3.7/dist-packages (from anndata&gt;=0.6.22.post1-&gt;kb-python) (5.5.0)'</li><li>'Collecting umap-learn&gt;=0.3.10'</li><li><span style=white-space:pre-wrap>'  Downloading https://files.pythonhosted.org/packages/75/69/85e7f950bb75792ad5d666d86c5f3e62eedbb942848e7e3126513af9999c/umap-learn-0.5.1.tar.gz (80kB)'</span></li><li>'Requirement already satisfied: joblib in /usr/local/lib/python3.7/dist-packages (from scanpy&gt;=1.4.4.post1-&gt;kb-python) (1.0.1)'</li><li>'Requirement already satisfied: seaborn in /usr/local/lib/python3.7/dist-packages (from scanpy&gt;=1.4.4.post1-&gt;kb-python) (0.11.1)'</li><li>'Requirement already satisfied: statsmodels&gt;=0.10.0rc2 in /usr/local/lib/python3.7/dist-packages (from scanpy&gt;=1.4.4.post1-&gt;kb-python) (0.10.2)'</li><li>'Requirement already satisfied: patsy in /usr/local/lib/python3.7/dist-packages (from scanpy&gt;=1.4.4.post1-&gt;kb-python) (0.5.1)'</li><li>'Requirement already satisfied: matplotlib&gt;=3.1.2 in /usr/local/lib/python3.7/dist-packages (from scanpy&gt;=1.4.4.post1-&gt;kb-python) (3.2.2)'</li><li>'Collecting sinfo'</li><li><span style=white-space:pre-wrap>'  Downloading https://files.pythonhosted.org/packages/e1/4c/aef8456284f1a1c3645b938d9ca72388c9c4878e6e67b8a349c7d22fac78/sinfo-0.3.1.tar.gz'</span></li><li>'Collecting legacy-api-wrap'</li><li><span style=white-space:pre-wrap>'  Downloading https://files.pythonhosted.org/packages/a4/68/da997bc56bb69dcdcee4054f0bc42266909307b905389fbc54c9158f42da/legacy_api_wrap-1.2-py3-none-any.whl'</span></li><li>'Requirement already satisfied: tables in /usr/local/lib/python3.7/dist-packages (from scanpy&gt;=1.4.4.post1-&gt;kb-python) (3.4.4)'</li><li>'Requirement already satisfied: networkx&gt;=2.3 in /usr/local/lib/python3.7/dist-packages (from scanpy&gt;=1.4.4.post1-&gt;kb-python) (2.5)'</li><li>'Requirement already satisfied: MarkupSafe&gt;=0.23 in /usr/local/lib/python3.7/dist-packages (from Jinja2&gt;2.10.1-&gt;kb-python) (1.1.1)'</li><li>'Requirement already satisfied: retrying&gt;=1.3.3 in /usr/local/lib/python3.7/dist-packages (from plotly&gt;=4.5.0-&gt;kb-python) (1.3.3)'</li><li>'Requirement already satisfied: llvmlite&lt;0.35,&gt;=0.34.0.dev0 in /usr/local/lib/python3.7/dist-packages (from numba-&gt;loompy&gt;=3.0.6-&gt;kb-python) (0.34.0)'</li><li>'Requirement already satisfied: webencodings in /usr/local/lib/python3.7/dist-packages (from bleach-&gt;nbconvert&gt;=5.6.0-&gt;kb-python) (0.5.1)'</li><li>'Requirement already satisfied: pyparsing&gt;=2.0.2 in /usr/local/lib/python3.7/dist-packages (from packaging-&gt;anndata&gt;=0.6.22.post1-&gt;kb-python) (2.4.7)'</li><li>'Requirement already satisfied: typing-extensions&gt;=3.6.4; python_version &lt; "3.8" in /usr/local/lib/python3.7/dist-packages (from importlib-metadata&gt;=0.7; python_version &lt; "3.8"-&gt;anndata&gt;=0.6.22.post1-&gt;kb-python) (3.7.4.3)'</li><li>'Requirement already satisfied: zipp&gt;=0.5 in /usr/local/lib/python3.7/dist-packages (from importlib-metadata&gt;=0.7; python_version &lt; "3.8"-&gt;anndata&gt;=0.6.22.post1-&gt;kb-python) (3.4.1)'</li><li>'Requirement already satisfied: python-dateutil&gt;=2.7.3 in /usr/local/lib/python3.7/dist-packages (from pandas!=1.1,&gt;=1.0-&gt;anndata&gt;=0.6.22.post1-&gt;kb-python) (2.8.1)'</li><li>'Requirement already satisfied: pytz&gt;=2017.2 in /usr/local/lib/python3.7/dist-packages (from pandas!=1.1,&gt;=1.0-&gt;anndata&gt;=0.6.22.post1-&gt;kb-python) (2018.9)'</li><li>'Collecting pynndescent&gt;=0.5'</li><li><span style=white-space:pre-wrap>'  Downloading https://files.pythonhosted.org/packages/af/65/8189298dd3a05bbad716ee8e249764ff8800e365d8dc652ad2192ca01b4a/pynndescent-0.5.2.tar.gz (1.1MB)'</span></li><li>'Requirement already satisfied: cycler&gt;=0.10 in /usr/local/lib/python3.7/dist-packages (from matplotlib&gt;=3.1.2-&gt;scanpy&gt;=1.4.4.post1-&gt;kb-python) (0.10.0)'</li><li>'Requirement already satisfied: kiwisolver&gt;=1.0.1 in /usr/local/lib/python3.7/dist-packages (from matplotlib&gt;=3.1.2-&gt;scanpy&gt;=1.4.4.post1-&gt;kb-python) (1.3.1)'</li><li>'Collecting stdlib_list'</li><li><span style=white-space:pre-wrap>'  Downloading https://files.pythonhosted.org/packages/7a/b1/52f59dcf31ead2f0ceff8976288449608d912972b911f55dff712cef5719/stdlib_list-0.8.0-py3-none-any.whl (63kB)'</span></li><li>'Collecting get-version&gt;=2.0.4'</li><li><span style=white-space:pre-wrap>'  Downloading https://files.pythonhosted.org/packages/23/48/7610e884e62fff2183e7bc8592397c39a020267fb5147905fcd3f9cc820c/get_version-2.1-py3-none-any.whl (43kB)'</span></li><li>'Requirement already satisfied: numexpr&gt;=2.5.2 in /usr/local/lib/python3.7/dist-packages (from tables-&gt;scanpy&gt;=1.4.4.post1-&gt;kb-python) (2.7.3)'</li><li>'Requirement already satisfied: decorator&gt;=4.3.0 in /usr/local/lib/python3.7/dist-packages (from networkx&gt;=2.3-&gt;scanpy&gt;=1.4.4.post1-&gt;kb-python) (4.4.2)'</li><li>'Building wheels for collected packages: loompy, numpy-groupies, umap-learn, sinfo, pynndescent'</li><li><span style=white-space:pre-wrap>'  Building wheel for loompy (setup.py): started'</span></li><li><span style=white-space:pre-wrap>'  Building wheel for loompy (setup.py): finished with status \'done\''</span></li><li><span style=white-space:pre-wrap>'  Created wheel for loompy: filename=loompy-3.0.6-cp37-none-any.whl size=47896 sha256=9a9ab62e41414913f9bbfb2eef13596a38bc9dd39316cabab892f51143ed4d10'</span></li><li><span style=white-space:pre-wrap>'  Stored in directory: /root/.cache/pip/wheels/f9/a4/90/5a98ad83419732b0fba533b81a2a52ba3dbe230a936ca4cdc9'</span></li><li><span style=white-space:pre-wrap>'  Building wheel for numpy-groupies (setup.py): started'</span></li><li><span style=white-space:pre-wrap>'  Building wheel for numpy-groupies (setup.py): finished with status \'done\''</span></li><li><span style=white-space:pre-wrap>'  Created wheel for numpy-groupies: filename=numpy_groupies-0.9.13-cp37-none-any.whl size=24068 sha256=bb4d7f7f96af13b58c827bf94f338dd2559b83522f750c17546667d33088ca2e'</span></li><li><span style=white-space:pre-wrap>'  Stored in directory: /root/.cache/pip/wheels/ef/97/d7/270bc85eb8b1b84629caac97a3900bff23edb7f834f6ed729e'</span></li><li><span style=white-space:pre-wrap>'  Building wheel for umap-learn (setup.py): started'</span></li><li><span style=white-space:pre-wrap>'  Building wheel for umap-learn (setup.py): finished with status \'done\''</span></li><li><span style=white-space:pre-wrap>'  Created wheel for umap-learn: filename=umap_learn-0.5.1-cp37-none-any.whl size=76569 sha256=627ba759b3e6d1787390e59be46a8ee88ea2f8e3961907069df003122d8e052c'</span></li><li><span style=white-space:pre-wrap>'  Stored in directory: /root/.cache/pip/wheels/ad/df/d5/a3691296ff779f25cd1cf415a3af954b987fb53111e3392cf4'</span></li><li><span style=white-space:pre-wrap>'  Building wheel for sinfo (setup.py): started'</span></li><li><span style=white-space:pre-wrap>'  Building wheel for sinfo (setup.py): finished with status \'done\''</span></li><li><span style=white-space:pre-wrap>'  Created wheel for sinfo: filename=sinfo-0.3.1-cp37-none-any.whl size=7012 sha256=0a7ca4b4b981f9ff8c0ebc0f677119dfa8e32cb5746fca18b261ae01f0f18e7e'</span></li><li><span style=white-space:pre-wrap>'  Stored in directory: /root/.cache/pip/wheels/11/f0/23/347d6d8e59787c2bc272162d18223dc3b45bd6dc40aceee6af'</span></li><li><span style=white-space:pre-wrap>'  Building wheel for pynndescent (setup.py): started'</span></li><li><span style=white-space:pre-wrap>'  Building wheel for pynndescent (setup.py): finished with status \'done\''</span></li><li><span style=white-space:pre-wrap>'  Created wheel for pynndescent: filename=pynndescent-0.5.2-cp37-none-any.whl size=51351 sha256=46d4d6df8cbd2dabd69c0f051626ff9b8978fb9141ff68278e7c7d83ec086277'</span></li><li><span style=white-space:pre-wrap>'  Stored in directory: /root/.cache/pip/wheels/ba/52/4e/4c28d04d144a28f89e2575fb63628df6e6d49b56c5ddd0c74e'</span></li><li>'Successfully built loompy numpy-groupies umap-learn sinfo pynndescent'</li><li>'Installing collected packages: numpy-groupies, loompy, anndata, pynndescent, umap-learn, stdlib-list, sinfo, get-version, legacy-api-wrap, scanpy, plotly, kb-python'</li><li><span style=white-space:pre-wrap>'  Found existing installation: plotly 4.4.1'</span></li><li><span style=white-space:pre-wrap>'    Uninstalling plotly-4.4.1:'</span></li><li><span style=white-space:pre-wrap>'      Successfully uninstalled plotly-4.4.1'</span></li><li>'Successfully installed anndata-0.7.5 get-version-2.1 kb-python-0.25.1 legacy-api-wrap-1.2 loompy-3.0.6 numpy-groupies-0.9.13 plotly-4.14.3 pynndescent-0.5.2 scanpy-1.7.1 sinfo-0.3.1 stdlib-list-0.8.0 umap-learn-0.5.1'</li></ol>



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
library(Seurat)
library(Matrix)
library(tidyverse)
library(patchwork)
theme_set(theme_bw())
```

    Attaching SeuratObject
    
    Warning message in system("timedatectl", intern = TRUE):
    “running command 'timedatectl' had status 1”
    Registered S3 method overwritten by 'cli':
      method     from         
      print.boxx spatstat.geom
    
    ── [1mAttaching packages[22m ─────────────────────────────────────── tidyverse 1.3.0 ──
    
    [32m✔[39m [34mggplot2[39m 3.3.3     [32m✔[39m [34mpurrr  [39m 0.3.4
    [32m✔[39m [34mtibble [39m 3.1.0     [32m✔[39m [34mdplyr  [39m 1.0.5
    [32m✔[39m [34mtidyr  [39m 1.1.3     [32m✔[39m [34mstringr[39m 1.4.0
    [32m✔[39m [34mreadr  [39m 1.4.0     [32m✔[39m [34mforcats[39m 0.5.1
    
    ── [1mConflicts[22m ────────────────────────────────────────── tidyverse_conflicts() ──
    [31m✖[39m [34mtidyr[39m::[32mexpand()[39m masks [34mMatrix[39m::expand()
    [31m✖[39m [34mdplyr[39m::[32mfilter()[39m masks [34mstats[39m::filter()
    [31m✖[39m [34mdplyr[39m::[32mlag()[39m    masks [34mstats[39m::lag()
    [31m✖[39m [34mtidyr[39m::[32mpack()[39m   masks [34mMatrix[39m::pack()
    [31m✖[39m [34mtidyr[39m::[32munpack()[39m masks [34mMatrix[39m::unpack()
    



```R
list.files(".", recursive = TRUE)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'index.idx'</li><li>'output/10xv3_whitelist.txt'</li><li>'output/counts_filtered/cells_x_genes.barcodes.txt'</li><li>'output/counts_filtered/cells_x_genes.genes.txt'</li><li>'output/counts_filtered/cells_x_genes.mtx'</li><li>'output/counts_unfiltered/cells_x_genes.barcodes.txt'</li><li>'output/counts_unfiltered/cells_x_genes.genes.txt'</li><li>'output/counts_unfiltered/cells_x_genes.mtx'</li><li>'output/filter_barcodes.txt'</li><li>'output/inspect.json'</li><li>'output/kb_info.json'</li><li>'output/matrix.ec'</li><li>'output/output.bus'</li><li>'output/output.filtered.bus'</li><li>'output/output.unfiltered.bus'</li><li>'output/run_info.json'</li><li>'output/transcripts.txt'</li><li>'pbmc_1k_v3_fastqs.tar'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_I1_001.fastq.gz'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_I1_001.fastq.gz'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz'</li><li>'pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz'</li><li>'sample_data/anscombe.json'</li><li>'sample_data/california_housing_test.csv'</li><li>'sample_data/california_housing_train.csv'</li><li>'sample_data/mnist_test.csv'</li><li>'sample_data/mnist_train_small.csv'</li><li>'sample_data/README.md'</li><li>'t2g.txt'</li></ol>




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
  scale_fill_distiller(palette = "Blues", direction = 1) +
  scale_x_log10() + scale_y_log10() + annotation_logticks()
```

    Warning message:
    “Transformation introduced infinite values in continuous x-axis”
    Warning message:
    “Transformation introduced infinite values in continuous y-axis”
    Warning message:
    “Removed 19583 rows containing non-finite values (stat_bin2d).”



![png](kb_analysis_0_R_files/kb_analysis_0_R_26_1.png)


Lots of points are piled at around 1 gene and 1 count, and those with 0 gene or count were removed for introducing -Inf in log transform. These correspond to empty or near empty droplets.

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
#' @rdname knee_plot
#' @param mat Gene count matrix, a dgCMatrix.
#' @return `get_knee_df` returns a tibble with two columns: \code{total} for 
#' total UMI counts for each barcode, and \code{rank} for rank of the total 
#' counts, with number 1 for the barcode with the most counts.
#' @export
#' @importFrom dplyr row_number desc arrange
#' @importFrom Matrix colSums
get_knee_df <- function(mat) {
  total <- rank <- NULL
  tibble(total = Matrix::colSums(mat),
         rank = row_number(desc(total))) %>%
    distinct() %>%
    dplyr::filter(total > 0) %>% 
    arrange(rank)
}

#' @rdname knee_plot
#' @param df The data frame from \code{\link{get_knee_df}}.
#' @param lower Minimum total UMI counts for barcode for it to be considered
#' when calculating the inflection point; this helps to avoid the noisy part of
#' the curve for barcodes with very few counts.
#' @return `get_inflection` returns a \code{numeric(1)} for the total UMI count 
#' at the inflection point.
#' @note Code in part adapted from \code{barcodeRanks} from \code{DropetUtils}.
#' @export
#' @importFrom dplyr transmute
#' 
get_inflection <- function(df, lower = 100) {
  log_total <- log_rank <- total <-  NULL
  df_fit <- df %>% 
    dplyr::filter(total > lower) %>% 
    transmute(log_total = log10(total),
              log_rank = log10(rank))
  d1n <- diff(df_fit$log_total)/diff(df_fit$log_rank)
  right.edge <- which.min(d1n)
  10^(df_fit$log_total[right.edge])
}

#' Plot the transposed knee plot and inflection point
#' 
#' Plot a transposed knee plot, showing the inflection point and
#' the number of remaining cells after inflection point filtering. It's
#' transposed since it's more generalizable to multi-modal data. Taken from the 
#' BUSpaRse package.
#' 
#' @param df The data frame from \code{\link{get_knee_df}}.
#' @param inflection Output of \code{\link{get_inflection}}.
#' @return `knee_plot` returns a \code{ggplot2} object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_path geom_vline geom_hline 
#' scale_x_log10 scale_y_log10 labs annotation_logticks geom_text
knee_plot <- function(df, inflection) {
  total <- rank_cutoff <- NULL
  annot <- tibble(inflection = inflection,
                  rank_cutoff = max(df$rank[df$total > inflection]))
  ggplot(df, aes(total, rank)) +
    geom_path() +
    geom_vline(aes(xintercept = inflection), data = annot, linetype = 2, 
               color = "gray40") +
    geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2, 
               color = "gray40") +
    geom_text(aes(inflection, rank_cutoff, 
                  label = paste(rank_cutoff, "'cells'")),
              data = annot, vjust = 1) +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = "Rank", x = "Total UMIs") +
    annotation_logticks()
}
```


```R
options(repr.plot.width=9, repr.plot.height=6)
knee_df <- get_knee_df(res_mat)
inflection <- get_inflection(knee_df)
knee_plot(knee_df, inflection)
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_31_0.png)


## Analysis
We begin by asking for genes with the highest proportions in droplets (prior to filtering out empty droplets).


```R
tr2g <- read_tsv("t2g.txt", col_names = c("transcript", "gene", "gene_name"))
tr2g <- distinct(tr2g[, c("gene", "gene_name")])
```

    
    [36m──[39m [1m[1mColumn specification[1m[22m [36m────────────────────────────────────────────────────────[39m
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



![png](kb_analysis_0_R_files/kb_analysis_0_R_35_1.png)


For many barcodes, the top genes by proportion of all counts are ribosomal or mitochondrial genes. Also, the proportions plotted above seem to have some discrete values; this effect is a result of computing fractions with small denominator, which happens when droplets produce very few UMI counts.

### Filter


```R
res_mat <- res_mat[, tot_counts > inflection]
res_mat <- res_mat[Matrix::rowSums(res_mat) > 0,]
dim(res_mat)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>31861</li><li>1322</li></ol>




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
    25966 features across 1208 samples within 1 assay 
    Active assay: RNA (25966 features, 0 variable features)


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


<table class="dataframe">
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


![png](kb_analysis_0_R_files/kb_analysis_0_R_46_0.png)



```R
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_47_0.png)



```R
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 20)
```

### Normalize

After removing unwanted cells from the dataset, the next step is to normalize the data. A standard choice is `LogNormalize` which normalizes the UMI counts for each cell by the total counts, multiplies this by a scale factor (10,000 by default), and finally log-transforms the result. 

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


![png](kb_analysis_0_R_files/kb_analysis_0_R_56_0.png)



```R
options(repr.plot.width=9, repr.plot.height=6)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc, log = FALSE)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
```

    When using repel, set xnudge and ynudge to 0 for optimal results
    



![png](kb_analysis_0_R_files/kb_analysis_0_R_57_1.png)


### Scaling the data
Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData function shifts the expression of each gene, so that the mean expression across cells is 0 and the variance across cells is 1
This step gives equal weight to genes in downstream analyses, so that highly-expressed genes do not dominate.


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
    Positive:  S100A9, FCN1, MNDA, FGL2, S100A8, CTSS, CST3, SERPINA1, PSAP, NCF2 
    	   LYZ, AIF1, TYMP, VCAN, KLF4, GRN, CSTA, MPEG1, CPVL, CLEC7A 
    	   LST1, MS4A6A, CD14, LGALS1, S100A12, TYROBP, TNFAIP2, FCER1G, CD36, CSF3R 
    Negative:  LTB, TRAC, TRBC2, CD3D, IL32, BCL11B, CD3G, IL7R, TCF7, CD69 
    	   ISG20, CD247, CD27, SPOCK2, ARL4C, CD7, CD2, GZMM, TRBC1, CD6 
    	   PRKCQ-AS1, NOSIP, AC058791.1, RORA, CTSW, CCR7, AQP3, ITM2A, PEBP1, SAMD3 
    PC_ 2 
    Positive:  CD79A, MS4A1, IGHM, BANK1, BCL11A, LINC00926, CD79B, TNFRSF13C, IGHD, CD74 
    	   HLA-DQB1, CD22, HLA-DQA1, HLA-DRB1, HLA-DRA, HLA-DPA1, HLA-DPB1, TCL1A, FCER2, AFF3 
    	   PAX5, IGKC, VPREB3, SPIB, MEF2C, RALGPS2, HVCN1, FCRL1, HLA-DOB, HLA-DMA 
    Negative:  IL32, CD247, GZMM, CD7, CTSW, CD3D, GZMA, NKG7, S100A4, TRAC 
    	   ANXA1, BCL11B, PRF1, CST7, KLRB1, CD3G, IL7R, ARL4C, SAMD3, CD2 
    	   TRBC1, CCL5, KLRG1, A2M, MT2A, RORA, ITGB2, GNLY, TCF7, MATK 
    PC_ 3 
    Positive:  CAVIN2, GP9, PF4, GNG11, PPBP, CD9, TREML1, CMTM5, TUBB1, SPARC 
    	   CLU, HIST1H2AC, ACRBP, PTCRA, PRKAR2B, NRGN, ITGA2B, CTTN, TMEM40, TSC22D1 
    	   AC147651.1, GMPR, PF4V1, CLDN5, CA2, MAP3K7CL, PGRMC1, CXCR2P1, HIST1H3H, MMD 
    Negative:  CYBA, VIM, FOS, ITGB2, NEAT1, HNRNPU, CALR, LSP1, LCP1, DUSP1 
    	   S100A10, S100A6, KLF6, CD74, PLAC8, LTB, ZFP36L1, S100A4, IFITM2, ISG20 
    	   SPCS1, SEC61B, ANXA1, MCL1, EVI2B, HSPA5, APOBEC3G, HSP90B1, PEBP1, AC020916.1 
    PC_ 4 
    Positive:  LEF1, TCF7, IL7R, MAL, CCR7, BCL11B, CD3D, NOSIP, LTB, CD3G 
    	   TRAC, CAMK4, NELL2, PASK, CD27, EGR1, SLC2A3, RGCC, FHIT, RGS10 
    	   CD6, CD40LG, VIM, INPP4B, ADTRP, TRAT1, NOG, TSHZ2, PRKCQ-AS1, TESPA1 
    Negative:  GZMB, GNLY, CLIC3, NKG7, KLRF1, PRF1, CST7, SPON2, FGFBP2, KLRD1 
    	   GZMA, ADGRG1, CCL4, TRDC, HOPX, MATK, IL2RB, TTC38, APOBEC3G, CTSW 
    	   TBX21, RHOC, C12orf75, S1PR5, FCGR3A, SH2D1B, PTGDR, MYOM2, CMC1, GZMH 
    PC_ 5 
    Positive:  LILRA4, SCT, PACSIN1, SMPD3, LRRC26, SERPINF1, TPM2, AL096865.1, IL3RA, DNASE1L3 
    	   TNFRSF21, CUX2, PLD4, ITM2C, GAS6, MYBL2, CLEC4C, PPP1R14B, EPHA2, UGCG 
    	   PPP1R14B-AS1, CUEDC1, LAMP5, RUNX2, PPM1J, SERPINF2, NRP1, DERL3, LINC02812, CIB2 
    Negative:  GNLY, FGFBP2, KLRF1, PRF1, NKG7, CCL4, CST7, KLRD1, MS4A1, IGHD 
    	   LINC00926, ADGRG1, CD79B, CD79A, TRDC, CD22, GZMA, TBX21, SPON2, TNFRSF13C 
    	   MATK, FCER2, IL2RB, PAX5, MYOM2, HOPX, S1PR5, TTC38, SH2D1B, BANK1 
    


Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction, DimPlot, and DimHeatmap




```R
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
```

    PC_ 1 
    Positive:  S100A9, FCN1, MNDA, FGL2, S100A8 
    Negative:  LTB, TRAC, TRBC2, CD3D, IL32 
    PC_ 2 
    Positive:  CD79A, MS4A1, IGHM, BANK1, BCL11A 
    Negative:  IL32, CD247, GZMM, CD7, CTSW 
    PC_ 3 
    Positive:  CAVIN2, GP9, PF4, GNG11, PPBP 
    Negative:  CYBA, VIM, FOS, ITGB2, NEAT1 
    PC_ 4 
    Positive:  LEF1, TCF7, IL7R, MAL, CCR7 
    Negative:  GZMB, GNLY, CLIC3, NKG7, KLRF1 
    PC_ 5 
    Positive:  LILRA4, SCT, PACSIN1, SMPD3, LRRC26 
    Negative:  GNLY, FGFBP2, KLRF1, PRF1, NKG7 


Which genes are contributing the most to the first 2 PCs?


```R
options(repr.plot.width=6, repr.plot.height=8)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_70_0.png)



```R
options(repr.plot.width=7, repr.plot.height=6)
FeaturePlot(pbmc, reduction = "pca", feature = "CST3")
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_71_0.png)


### Determining dimensionality
To overcome the extensive technical noise in any single feature for scRNA-seq data, one can cluster cells based on their PCA projections, with each PC essentially representing a ‘metafeature’ that combines information across a correlated feature set.

A common heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one (ElbowPlot function). In this example, we can observe an ‘elbow’ around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs.


```R
options(repr.plot.width=9, repr.plot.height=6)
ElbowPlot(pbmc)
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_74_0.png)


### The neighborhood graph

We cluster cells using the Louvain algorithm (a default in Seurat), which iteratively group cells together, with the goal of optimizing the standard modularity function. The FindClusters function implements this procedure, and contains a resolution parameter that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters. We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. The clusters can be found using the Idents function.


```R
pbmc <- FindNeighbors(pbmc, dims = 1:10, k.param = 20)
pbmc <- FindClusters(pbmc, resolution = 0.6)
```

    Computing nearest neighbor graph
    
    Computing SNN
    


    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 1138
    Number of edges: 35679
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8720
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
</style><dl class=dl-inline><dt>AAACCCAAGGAGAGTA</dt><dd>0</dd><dt>AAACGCTTCAGCCCAG</dt><dd>3</dd><dt>AAAGAACAGACGACTG</dt><dd>5</dd><dt>AAAGAACCAATGGCAG</dt><dd>5</dd><dt>AAAGAACGTCTGCAAT</dt><dd>1</dd></dl>

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


![png](kb_analysis_0_R_files/kb_analysis_0_R_81_0.png)



```R
options(repr.plot.width=16, repr.plot.height=5)
FeaturePlot(pbmc, reduction = "umap", features = c("CST3", "NKG7", "PPBP"),
ncol = 3)
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_82_0.png)


### Finding differentially expressed features (cluster biomarkers)
A key follow-up step to clustering cells is to find gene markers that are associated with them. We used Seurat's FindAllMarkers function which automates the process for all clusters.

The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups. You can set both of these to 0, but with a dramatic increase in time - since this will test a large number of features that are unlikely to be highly discriminatory. As another option to speed up these computations, max.cells.per.ident can be set. This will downsample each identity class to have no more cells than whatever this is set to. While there is generally going to be a loss in power, the speed increases can be significiant and the most highly differentially expressed features will likely still rise to the top.


```R
# Scanpy style gene rank plot
plot_gene_rank <- function(markers, n) {
  df_plot <- markers %>%
    group_by(cluster) %>%
    top_n(25, avg_log2FC) %>%
    mutate(rank = factor(row_number(desc(avg_log2FC))))
  ggplot(df_plot, aes(rank, avg_log2FC)) +
    geom_text(aes(label = gene), angle = -90, hjust = 1) +
    facet_wrap(~ cluster) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
}
```

Several methods for differential expression are supported by Seurat. The default is Wilcoxon rank sum test.


```R
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, test.use = "wilcox", only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)
```

    Calculating cluster 0
    
    For a more efficient implementation of the Wilcoxon Rank Sum Test,
    (default method for FindMarkers) please install the limma package
    --------------------------------------------
    install.packages('BiocManager')
    BiocManager::install('limma')
    --------------------------------------------
    After installation of limma, Seurat will automatically use the more 
    efficient implementation (no further action necessary).
    This message will be shown once per session
    
    Calculating cluster 1
    
    Calculating cluster 2
    
    Calculating cluster 3
    
    Calculating cluster 4
    
    Calculating cluster 5
    
    Calculating cluster 6
    
    Calculating cluster 7
    
    Calculating cluster 8
    



```R
head(pbmc.markers)
```


<table class="dataframe">
<caption>A data.frame: 6 × 7</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>cluster</th><th scope=col>gene</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>S100A12</th><td>4.064804e-197</td><td>4.213263</td><td>0.978</td><td>0.076</td><td>1.055467e-192</td><td>0</td><td>S100A12</td></tr>
	<tr><th scope=row>VCAN</th><td>3.352042e-189</td><td>3.484980</td><td>0.996</td><td>0.104</td><td>8.703913e-185</td><td>0</td><td>VCAN   </td></tr>
	<tr><th scope=row>S100A8</th><td>4.318649e-186</td><td>5.520490</td><td>1.000</td><td>0.143</td><td>1.121380e-181</td><td>0</td><td>S100A8 </td></tr>
	<tr><th scope=row>CD14</th><td>7.936720e-173</td><td>2.574749</td><td>0.944</td><td>0.088</td><td>2.060849e-168</td><td>0</td><td>CD14   </td></tr>
	<tr><th scope=row>S100A9</th><td>7.548587e-171</td><td>4.954452</td><td>0.996</td><td>0.197</td><td>1.960066e-166</td><td>0</td><td>S100A9 </td></tr>
	<tr><th scope=row>MNDA</th><td>4.262766e-170</td><td>3.039318</td><td>0.981</td><td>0.130</td><td>1.106870e-165</td><td>0</td><td>MNDA   </td></tr>
</tbody>
</table>




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
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
    Warning message:
    “glm.fit: algorithm did not converge”
    Warning message:
    “glm.fit: fitted probabilities numerically 0 or 1 occurred”
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
0 |	CD14, LYZ |	CD14+ Mono
1	| IL7R, S100A4 |	Memory CD4+ T
2 |	IL7R, CCR7 |	Naive CD4+
3 |	MS4A1, CD79A |	B
4 |	FCGR3A, MS4A7 |	FCGR3A+ Mono
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
new.cluster.ids <- c("CD14+ Mono", "Memory CD4 T", "Naive CD4 T", "B1", "FCGR3A+ Mono", 
    "NK", "CD8+ T", "B2", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6) + NoLegend()
```


![png](kb_analysis_0_R_files/kb_analysis_0_R_102_0.png)



```R
Sys.time() - start_time
```


    Time difference of 50.73962 mins



```R
sessionInfo()
```


    R version 4.0.4 (2021-02-15)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: Ubuntu 18.04.5 LTS
    
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
    [1] stats     graphics  grDevices utils     datasets  methods   base     
    
    other attached packages:
     [1] patchwork_1.1.1    forcats_0.5.1      stringr_1.4.0      dplyr_1.0.5       
     [5] purrr_0.3.4        readr_1.4.0        tidyr_1.1.3        tibble_3.1.0      
     [9] ggplot2_3.3.3      tidyverse_1.3.0    Matrix_1.3-2       SeuratObject_4.0.0
    [13] Seurat_4.0.1      
    
    loaded via a namespace (and not attached):
      [1] Rtsne_0.15            colorspace_2.0-0      deldir_0.2-10        
      [4] ellipsis_0.3.1        ggridges_0.5.3        IRdisplay_1.0        
      [7] fs_1.5.0              base64enc_0.1-3       rstudioapi_0.13      
     [10] spatstat.data_2.1-0   farver_2.1.0          leiden_0.3.7         
     [13] listenv_0.8.0         ggrepel_0.9.1         RSpectra_0.16-0      
     [16] lubridate_1.7.10      fansi_0.4.2           xml2_1.3.2           
     [19] codetools_0.2-18      splines_4.0.4         polyclip_1.10-0      
     [22] IRkernel_1.1.1        jsonlite_1.7.2        broom_0.7.5          
     [25] ica_1.0-2             cluster_2.1.1         dbplyr_2.1.0         
     [28] png_0.1-7             uwot_0.1.10           shiny_1.6.0          
     [31] sctransform_0.3.2     spatstat.sparse_2.0-0 compiler_4.0.4       
     [34] httr_1.4.2            backports_1.2.1       assertthat_0.2.1     
     [37] fastmap_1.1.0         lazyeval_0.2.2        cli_2.3.1            
     [40] later_1.1.0.1         htmltools_0.5.1.1     tools_4.0.4          
     [43] igraph_1.2.6          gtable_0.3.0          glue_1.4.2           
     [46] RANN_2.6.1            reshape2_1.4.4        Rcpp_1.0.6           
     [49] scattermore_0.7       cellranger_1.1.0      vctrs_0.3.6          
     [52] nlme_3.1-152          lmtest_0.9-38         ps_1.6.0             
     [55] globals_0.14.0        rvest_1.0.0           mime_0.10            
     [58] miniUI_0.1.1.1        lifecycle_1.0.0       irlba_2.3.3          
     [61] goftest_1.2-2         future_1.21.0         MASS_7.3-53.1        
     [64] zoo_1.8-9             scales_1.1.1          spatstat.core_2.0-0  
     [67] hms_1.0.0             promises_1.2.0.1      spatstat.utils_2.1-0 
     [70] parallel_4.0.4        RColorBrewer_1.1-2    reticulate_1.18      
     [73] pbapply_1.4-3         gridExtra_2.3         rpart_4.1-15         
     [76] stringi_1.5.3         repr_1.1.3            rlang_0.4.10         
     [79] pkgconfig_2.0.3       matrixStats_0.58.0    evaluate_0.14        
     [82] lattice_0.20-41       ROCR_1.0-11           tensor_1.5           
     [85] labeling_0.4.2        htmlwidgets_1.5.3     cowplot_1.1.1        
     [88] tidyselect_1.1.0      parallelly_1.24.0     RcppAnnoy_0.0.18     
     [91] plyr_1.8.6            magrittr_2.0.1        R6_2.5.0             
     [94] generics_0.1.0        pbdZMQ_0.3-5          DBI_1.1.1            
     [97] withr_2.4.1           haven_2.3.1           pillar_1.5.1         
    [100] mgcv_1.8-34           fitdistrplus_1.1-3    survival_3.2-10      
    [103] abind_1.4-5           future.apply_1.7.0    modelr_0.1.8         
    [106] crayon_1.4.1          uuid_0.1-4            KernSmooth_2.23-18   
    [109] utf8_1.2.1            spatstat.geom_2.0-1   plotly_4.9.3         
    [112] readxl_1.3.1          grid_4.0.4            data.table_1.14.0    
    [115] reprex_1.0.0          digest_0.6.27         xtable_1.8-4         
    [118] httpuv_1.5.5          munsell_0.5.0         viridisLite_0.3.0    


**Feedback**: please report any issues, or submit pull requests for improvements, in the [Github repository where this notebook is located](https://github.com/pachterlab/kallistobustools/blob/master/notebooks/kb_analysis_0_R.ipynb).
