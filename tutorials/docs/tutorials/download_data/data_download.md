# Data downloading

This tutorial provides information on where to find single-cell RNA-seq data, and how to download it for processing with the **kallisto | bustools** workflow.

## Databases

There are four databases that are important repositories for sequencing data and metadata, and that are relevant for obtaining single-cell RNA-seq data. For each archive we provide an example of how the data is organized and how to download it.

* **[Biological Project Library](https://bigd.big.ac.cn/bioproject/)** (BioProject): The Biological Project Library organizes metadata for research projects involving genomic data types. This repository, which was started in 2016, is similar to the Gene Expression Omnibus. As an example, the data from the paper [Peng et al. 2019](https://www.nature.com/articles/s41422-019-0195-y) is organized under project accession [PRJCA001063](https://bigd.big.ac.cn/bioproject/browse/PRJCA001063). Each single-cell RNA-seq dataset has a “BioSample accession”, e.g. [SAMC047103](https://bigd.big.ac.cn/biosample/browse/SAMC047103). A further link to the Genome Sequencing Archive provides access to FASTQ files.

* **[Genome Sequence Archive](http://gsa.big.ac.cn/)** (GSA): This repository contains reads for projects in FASTQ format. For example, reads for [SAMC047103](https://bigd.big.ac.cn/biosample/browse/SAMC047103) from the [PRJCA001063](https://bigd.big.ac.cn/bioproject/browse/PRJCA001063) in the BioProject repository are accessible under accession [CRA001160](https://bigd.big.ac.cn/gsa/browse/CRA001160). A specific run accession, e.g. [CRR034516](https://bigd.big.ac.cn/gsa/browse/CRA001160/CRR034516) provides direct access to FASTQ files.

* **[Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/)** (GEO): The Gene Expression Omnibus is a repository for [MIAME (Minimum Infomration about a Microarray Experiment)](https://www.ncbi.nlm.nih.gov/geo/info/MIAME.html) compliant data. While the MIAME standards were established during a time when gene expression data was primarily collected with microarrays, the standards also apply to sequencing data and the GEO repository hosts project metadata for both types of research projects. As an example, the project link for the paper [Wolock et al. 2019](https://www.sciencedirect.com/science/article/pii/S2211124719307971) is [GSE132151](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132151). Most papers refer to their data via GEO accessions, so GEO is a useful repository for searching for data from projects.

* **[European Nucelotide Archive](https://www.ebi.ac.uk/ena)** (ENA): The ENA provides access to nucleotide sequences associated with genomic projects. In the case of [GSE132151](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132151) mentioned above, the nucleotide sequences are at [PRJNA546231](https://www.ebi.ac.uk/ena/data/view/PRJNA546231). The ENA provides direct access to FASTQ files from the project page. It also links to NCBI Sequence Read Archive format data.

* **[Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra)** (SRA): The SRA is a sequence repository for genomic data. Files are stored in SRA format, which must be downloaded and converted to FASTQ format prior to pre-processing using the `fasterq-dump` program available as part of [SRA tools](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump). For example, the data in [Rossi et al., 2019](https://science.sciencemag.org/content/364/6447/1271) can be located in the SRA via [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130597), then to the [SRA](https://www.ncbi.nlm.nih.gov/sra?term=SRP194426), and finally a sequence data page for one of the runs, [SRX5779290](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR9000493) has information about the traces (reads). The SRA tools operate directly on SRA accessions.

## Searching

The [sra-explorer](https://ewels.github.io/sra-explorer/) website is an effective and easy to use utility for searching the SRA and for downloading files. The utility finds SRA entires by keywords or accession numbers and produces links to the FASTQs and to commands for downloading them.

## Streaming

Single-cell RNA-seq data from sequence repositories can be streamed into `kb` making possible a workflow that does not require saving files to disk prior to pre-processing. For example, the following command can be used to stream data from the a URL:

__Note__: Streaming is not supported on Windows.

### Install `kb`


```
!pip install kb-python
```

    Collecting kb-python
    [?25l  Downloading https://files.pythonhosted.org/packages/62/c9/2e5b8fa2cd873a23ae1aeb128b33165d6a9387a2f56ea1fafec1d6d32477/kb_python-0.24.4-py3-none-any.whl (35.4MB)
    [K     |████████████████████████████████| 35.4MB 124kB/s 
    [?25hCollecting loompy>=3.0.6
    [?25l  Downloading https://files.pythonhosted.org/packages/36/52/74ed37ae5988522fbf87b856c67c4f80700e6452410b4cd80498c5f416f9/loompy-3.0.6.tar.gz (41kB)
    [K     |████████████████████████████████| 51kB 8.2MB/s 
    [?25hCollecting anndata>=0.6.22.post1
    [?25l  Downloading https://files.pythonhosted.org/packages/2b/72/87196c15f68d9865c31a43a10cf7c50bcbcedd5607d09f9aada0b3963103/anndata-0.6.22.post1-py3-none-any.whl (47kB)
    [K     |████████████████████████████████| 51kB 8.5MB/s 
    [?25hRequirement already satisfied: h5py in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (2.8.0)
    Requirement already satisfied: numpy in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (1.17.5)
    Requirement already satisfied: scipy in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (1.4.1)
    Requirement already satisfied: setuptools in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (42.0.2)
    Requirement already satisfied: numba in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (0.47.0)
    Requirement already satisfied: click in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (7.0)
    Collecting numpy-groupies
    [?25l  Downloading https://files.pythonhosted.org/packages/57/ae/18217b57ba3e4bb8a44ecbfc161ed065f6d1b90c75d404bd6ba8d6f024e2/numpy_groupies-0.9.10.tar.gz (43kB)
    [K     |████████████████████████████████| 51kB 8.9MB/s 
    [?25hRequirement already satisfied: pandas>=0.23.0 in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python) (0.25.3)
    Requirement already satisfied: natsort in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python) (5.5.0)
    Requirement already satisfied: six in /usr/local/lib/python3.6/dist-packages (from h5py->loompy>=3.0.6->kb-python) (1.12.0)
    Requirement already satisfied: llvmlite>=0.31.0dev0 in /usr/local/lib/python3.6/dist-packages (from numba->loompy>=3.0.6->kb-python) (0.31.0)
    Requirement already satisfied: pytz>=2017.2 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata>=0.6.22.post1->kb-python) (2018.9)
    Requirement already satisfied: python-dateutil>=2.6.1 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata>=0.6.22.post1->kb-python) (2.6.1)
    Building wheels for collected packages: loompy, numpy-groupies
      Building wheel for loompy (setup.py) ... [?25l[?25hdone
      Created wheel for loompy: filename=loompy-3.0.6-cp36-none-any.whl size=47896 sha256=b4583c69463e8d433e513784273fd2c82a98415cca57c9b732877b018d860909
      Stored in directory: /root/.cache/pip/wheels/f9/a4/90/5a98ad83419732b0fba533b81a2a52ba3dbe230a936ca4cdc9
      Building wheel for numpy-groupies (setup.py) ... [?25l[?25hdone
      Created wheel for numpy-groupies: filename=numpy_groupies-0+unknown-cp36-none-any.whl size=28044 sha256=01158d0fd3640057dc39fdbbbeaebc47c415aef3d55b09394bc19205410b8f03
      Stored in directory: /root/.cache/pip/wheels/30/ac/83/64d5f9293aeaec63f9539142fc629a41af064cae1b3d8d94aa
    Successfully built loompy numpy-groupies
    Installing collected packages: numpy-groupies, loompy, anndata, kb-python
    Successfully installed anndata-0.6.22.post1 kb-python-0.24.4 loompy-3.0.6 numpy-groupies-0+unknown


### Download a pre-built mouse index

The only required file that must be locally stored on disk prior to pre-processing is the index, which is why we download it here.


```
%%time
!kb ref -d mouse -i index.idx -g t2g.txt
```

    [2020-01-13 22:10:41,612]    INFO Downloading files for mouse from https://caltech.box.com/shared/static/vcaz6cujop0xuapdmz0pplp3aoqc41si.gz to tmp/vcaz6cujop0xuapdmz0pplp3aoqc41si.gz
    [2020-01-13 22:13:24,857]    INFO Extracting files from tmp/vcaz6cujop0xuapdmz0pplp3aoqc41si.gz
    CPU times: user 664 ms, sys: 123 ms, total: 787 ms
    Wall time: 3min 16s



```
%%time
!kb count -i index.idx -g t2g.txt -x 10xv2 --h5ad -t 2 \
https://caltech.box.com/shared/static/w9ww8et5o029s2e3usjzpbq8lpot29rh.gz \
https://caltech.box.com/shared/static/ql00zyvqnpy7bf8ogdoe9zfy907guzy9.gz
```

    [2020-01-13 22:14:01,325]    INFO Piping https://caltech.box.com/shared/static/w9ww8et5o029s2e3usjzpbq8lpot29rh.gz to tmp/w9ww8et5o029s2e3usjzpbq8lpot29rh.gz
    [2020-01-13 22:14:01,327]    INFO Piping https://caltech.box.com/shared/static/ql00zyvqnpy7bf8ogdoe9zfy907guzy9.gz to tmp/ql00zyvqnpy7bf8ogdoe9zfy907guzy9.gz
    [2020-01-13 22:14:01,327]    INFO Generating BUS file from
    [2020-01-13 22:14:01,327]    INFO         tmp/w9ww8et5o029s2e3usjzpbq8lpot29rh.gz
    [2020-01-13 22:14:01,327]    INFO         tmp/ql00zyvqnpy7bf8ogdoe9zfy907guzy9.gz
    [2020-01-13 22:16:39,231]    INFO Sorting BUS file ./output.bus to tmp/output.s.bus
    [2020-01-13 22:16:42,864]    INFO Whitelist not provided
    [2020-01-13 22:16:42,864]    INFO Copying pre-packaged 10XV2 whitelist to .
    [2020-01-13 22:16:42,983]    INFO Inspecting BUS file tmp/output.s.bus
    [2020-01-13 22:16:44,697]    INFO Correcting BUS records in tmp/output.s.bus to tmp/output.s.c.bus with whitelist ./10xv2_whitelist.txt
    [2020-01-13 22:17:01,096]    INFO Sorting BUS file tmp/output.s.c.bus to ./output.unfiltered.bus
    [2020-01-13 22:17:04,468]    INFO Generating count matrix ./counts_unfiltered/cells_x_genes from BUS file ./output.unfiltered.bus
    [2020-01-13 22:17:06,930]    INFO Converting matrix ./counts_unfiltered/cells_x_genes.mtx to h5ad ./counts_unfiltered/adata.h5ad
    CPU times: user 828 ms, sys: 85.2 ms, total: 913 ms
    Wall time: 3min 10s



```

```
