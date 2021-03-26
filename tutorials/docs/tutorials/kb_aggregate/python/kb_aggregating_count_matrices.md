# Aggregating multiple count matrices tutorial

This tutorial describes how to aggregate multiple count matrices by concatenating them into a single [AnnData](https://anndata.readthedocs.io/en/latest/anndata.AnnData.html) object with batch labels for different samples.

This is similar to the Cell Ranger aggr function, however no normalization is performed. cellranger aggr is described at https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate

For this tutorial we use dataset E-MTAB-6108.

## Download the raw data

The raw data for E-MTAB-6108 is available at https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6108/


```
%%time
!wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample1_L005_R1.fastq.gz
!wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample1_L005_R2.fastq.gz
!wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample2_L005_R1.fastq.gz
!wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample2_L005_R2.fastq.gz
```

    --2020-01-14 22:03:03--  https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample1_L005_R1.fastq.gz
    Resolving www.ebi.ac.uk (www.ebi.ac.uk)... 193.62.193.80
    Connecting to www.ebi.ac.uk (www.ebi.ac.uk)|193.62.193.80|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 2812327595 (2.6G) [application/x-gzip]
    Saving to: ‘iPSC_RGCscRNAseq_Sample1_L005_R1.fastq.gz’
    
    iPSC_RGCscRNAseq_Sa 100%[===================>]   2.62G  25.9MB/s    in 1m 46s  
    
    2020-01-14 22:04:55 (25.4 MB/s) - ‘iPSC_RGCscRNAseq_Sample1_L005_R1.fastq.gz’ saved [2812327595/2812327595]
    
    --2020-01-14 22:04:56--  https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample1_L005_R2.fastq.gz
    Resolving www.ebi.ac.uk (www.ebi.ac.uk)... 193.62.193.80
    Connecting to www.ebi.ac.uk (www.ebi.ac.uk)|193.62.193.80|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 10626337032 (9.9G) [application/x-gzip]
    Saving to: ‘iPSC_RGCscRNAseq_Sample1_L005_R2.fastq.gz’
    
    iPSC_RGCscRNAseq_Sa 100%[===================>]   9.90G  26.7MB/s    in 6m 33s  
    
    2020-01-14 22:11:29 (25.8 MB/s) - ‘iPSC_RGCscRNAseq_Sample1_L005_R2.fastq.gz’ saved [10626337032/10626337032]
    
    --2020-01-14 22:11:30--  https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample2_L005_R1.fastq.gz
    Resolving www.ebi.ac.uk (www.ebi.ac.uk)... 193.62.193.80
    Connecting to www.ebi.ac.uk (www.ebi.ac.uk)|193.62.193.80|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 1938548485 (1.8G) [application/x-gzip]
    Saving to: ‘iPSC_RGCscRNAseq_Sample2_L005_R1.fastq.gz’
    
    iPSC_RGCscRNAseq_Sa 100%[===================>]   1.80G  27.6MB/s    in 72s     
    
    2020-01-14 22:12:43 (25.5 MB/s) - ‘iPSC_RGCscRNAseq_Sample2_L005_R1.fastq.gz’ saved [1938548485/1938548485]
    
    --2020-01-14 22:12:44--  https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample2_L005_R2.fastq.gz
    Resolving www.ebi.ac.uk (www.ebi.ac.uk)... 193.62.193.80
    Connecting to www.ebi.ac.uk (www.ebi.ac.uk)|193.62.193.80|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 7563578426 (7.0G) [application/x-gzip]
    Saving to: ‘iPSC_RGCscRNAseq_Sample2_L005_R2.fastq.gz’
    
    iPSC_RGCscRNAseq_Sa 100%[===================>]   7.04G  27.6MB/s    in 4m 42s  
    
    2020-01-14 22:17:27 (25.6 MB/s) - ‘iPSC_RGCscRNAseq_Sample2_L005_R2.fastq.gz’ saved [7563578426/7563578426]
    
    CPU times: user 6.92 s, sys: 1.37 s, total: 8.29 s
    Wall time: 14min 24s


## Install `kb`

Install `kb` for running the kallisto|bustools workflow.


```
!pip install kb-python
```

    Collecting kb-python
    [?25l  Downloading https://files.pythonhosted.org/packages/62/c9/2e5b8fa2cd873a23ae1aeb128b33165d6a9387a2f56ea1fafec1d6d32477/kb_python-0.24.4-py3-none-any.whl (35.4MB)
    [K     |████████████████████████████████| 35.4MB 75kB/s 
    [?25hCollecting loompy>=3.0.6
    [?25l  Downloading https://files.pythonhosted.org/packages/36/52/74ed37ae5988522fbf87b856c67c4f80700e6452410b4cd80498c5f416f9/loompy-3.0.6.tar.gz (41kB)
    [K     |████████████████████████████████| 51kB 5.6MB/s 
    [?25hCollecting anndata>=0.6.22.post1
    [?25l  Downloading https://files.pythonhosted.org/packages/2b/72/87196c15f68d9865c31a43a10cf7c50bcbcedd5607d09f9aada0b3963103/anndata-0.6.22.post1-py3-none-any.whl (47kB)
    [K     |████████████████████████████████| 51kB 5.8MB/s 
    [?25hRequirement already satisfied: h5py in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (2.8.0)
    Requirement already satisfied: numpy in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (1.17.5)
    Requirement already satisfied: scipy in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (1.4.1)
    Requirement already satisfied: setuptools in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (42.0.2)
    Requirement already satisfied: numba in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (0.47.0)
    Requirement already satisfied: click in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (7.0)
    Collecting numpy-groupies
    [?25l  Downloading https://files.pythonhosted.org/packages/57/ae/18217b57ba3e4bb8a44ecbfc161ed065f6d1b90c75d404bd6ba8d6f024e2/numpy_groupies-0.9.10.tar.gz (43kB)
    [K     |████████████████████████████████| 51kB 6.1MB/s 
    [?25hRequirement already satisfied: pandas>=0.23.0 in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python) (0.25.3)
    Requirement already satisfied: natsort in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python) (5.5.0)
    Requirement already satisfied: six in /usr/local/lib/python3.6/dist-packages (from h5py->loompy>=3.0.6->kb-python) (1.12.0)
    Requirement already satisfied: llvmlite>=0.31.0dev0 in /usr/local/lib/python3.6/dist-packages (from numba->loompy>=3.0.6->kb-python) (0.31.0)
    Requirement already satisfied: python-dateutil>=2.6.1 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata>=0.6.22.post1->kb-python) (2.6.1)
    Requirement already satisfied: pytz>=2017.2 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata>=0.6.22.post1->kb-python) (2018.9)
    Building wheels for collected packages: loompy, numpy-groupies
      Building wheel for loompy (setup.py) ... [?25l[?25hdone
      Created wheel for loompy: filename=loompy-3.0.6-cp36-none-any.whl size=47896 sha256=47ebc61e06159e73ae729b6d3f65f4e13422a67ac13e5d7aab3c15af4fe93961
      Stored in directory: /root/.cache/pip/wheels/f9/a4/90/5a98ad83419732b0fba533b81a2a52ba3dbe230a936ca4cdc9
      Building wheel for numpy-groupies (setup.py) ... [?25l[?25hdone
      Created wheel for numpy-groupies: filename=numpy_groupies-0+unknown-cp36-none-any.whl size=28044 sha256=50bf49ee83e33e8c5312562dd519758710d97415c67a9e857aee497922da19d8
      Stored in directory: /root/.cache/pip/wheels/30/ac/83/64d5f9293aeaec63f9539142fc629a41af064cae1b3d8d94aa
    Successfully built loompy numpy-groupies
    Installing collected packages: numpy-groupies, loompy, anndata, kb-python
    Successfully installed anndata-0.6.22.post1 kb-python-0.24.4 loompy-3.0.6 numpy-groupies-0+unknown


## Download a pre-built human index

__Note:__ See [this notebook]() for a tutorial on how to build custom transcriptome or RNA velocity indices.


```
%%time
!kb ref -d human -i index.idx -g t2g.txt
```

    [2020-01-14 22:17:40,464]    INFO Downloading files for human from https://caltech.box.com/shared/static/v1nm7lpnqz5syh8dyzdk2zs8bglncfib.gz to tmp/v1nm7lpnqz5syh8dyzdk2zs8bglncfib.gz
    [2020-01-14 22:19:31,668]    INFO Extracting files from tmp/v1nm7lpnqz5syh8dyzdk2zs8bglncfib.gz
    CPU times: user 578 ms, sys: 77.4 ms, total: 655 ms
    Wall time: 2min 32s


## Generate an RNA count matrices in H5AD format

The following command will generate an RNA count matrix of cells (rows) by genes (columns) in H5AD format, which is a binary format used to store [Anndata](https://anndata.readthedocs.io/en/stable/) objects. Notice we are providing the index and transcript-to-gene mapping we downloaded in the previous step to the `-i` and `-g` arguments respectively. Also, these reads were generated with the 10x Genomics Chromium Single Cell v2 Chemistry, hence the `-x 10xv2` argument. To view other supported technologies, run `kb --list`.

The `--filter` flag is used to filter out barcodes with low UMI counts. This will generate two matrices, one in the `counts_unfiltered` directory and another in the `counts_filtered` directory.

__Note:__ If you would like a Loom file instead, replace the `--h5ad` flag with `--loom`. If you want to use the raw matrix output by `kb` instead of their H5AD or Loom converted files, omit these flags.

### Sample 1


```
%%time
!kb count -i index.idx -g t2g.txt -x 10xv2 -o sample1 --h5ad -t 2 --filter bustools \
iPSC_RGCscRNAseq_Sample1_L005_R1.fastq.gz \
iPSC_RGCscRNAseq_Sample1_L005_R2.fastq.gz
```

    [2020-01-14 22:55:56,693]    INFO Skipping kallisto bus because output files already exist. Use the --overwrite flag to overwrite.
    [2020-01-14 22:55:56,693]    INFO Sorting BUS file sample1/output.bus to tmp/output.s.bus
    [2020-01-14 22:57:31,354]    INFO Whitelist not provided
    [2020-01-14 22:57:31,354]    INFO Copying pre-packaged 10XV2 whitelist to sample1
    [2020-01-14 22:57:35,347]    INFO Inspecting BUS file tmp/output.s.bus
    [2020-01-14 22:57:47,155]    INFO Correcting BUS records in tmp/output.s.bus to tmp/output.s.c.bus with whitelist sample1/10xv2_whitelist.txt
    [2020-01-14 22:58:13,462]    INFO Sorting BUS file tmp/output.s.c.bus to sample1/output.unfiltered.bus
    [2020-01-14 22:58:48,480]    INFO Generating count matrix sample1/counts_unfiltered/cells_x_genes from BUS file sample1/output.unfiltered.bus
    [2020-01-14 22:59:01,129]    INFO Converting matrix sample1/counts_unfiltered/cells_x_genes.mtx to h5ad sample1/counts_unfiltered/adata.h5ad
    [2020-01-14 22:59:11,951]    INFO Filtering with bustools
    [2020-01-14 22:59:11,952]    INFO Generating whitelist sample1/filter_barcodes.txt from BUS file sample1/output.unfiltered.bus
    [2020-01-14 22:59:12,274]    INFO Capturing records from BUS file sample1/output.unfiltered.bus to tmp/output.filtered.bus with capture list sample1/filter_barcodes.txt
    [2020-01-14 22:59:15,831]    INFO Sorting BUS file tmp/output.filtered.bus to sample1/output.filtered.bus
    [2020-01-14 22:59:52,828]    INFO Generating count matrix sample1/counts_filtered/cells_x_genes from BUS file sample1/output.filtered.bus
    [2020-01-14 23:00:03,942]    INFO Converting matrix sample1/counts_filtered/cells_x_genes.mtx to h5ad sample1/counts_filtered/adata.h5ad
    CPU times: user 1.24 s, sys: 161 ms, total: 1.4 s
    Wall time: 4min 16s


### Sample 2


```
%%time
!kb count -i index.idx -g t2g.txt -x 10xv2 -o sample2 --h5ad -t 2 --filter bustools \
iPSC_RGCscRNAseq_Sample2_L005_R1.fastq.gz \
iPSC_RGCscRNAseq_Sample2_L005_R2.fastq.gz
```

    [2020-01-14 23:00:13,871]    INFO Skipping kallisto bus because output files already exist. Use the --overwrite flag to overwrite.
    [2020-01-14 23:00:13,871]    INFO Sorting BUS file sample2/output.bus to tmp/output.s.bus
    [2020-01-14 23:01:14,475]    INFO Whitelist not provided
    [2020-01-14 23:01:14,475]    INFO Copying pre-packaged 10XV2 whitelist to sample2
    [2020-01-14 23:01:14,681]    INFO Inspecting BUS file tmp/output.s.bus
    [2020-01-14 23:01:21,144]    INFO Correcting BUS records in tmp/output.s.bus to tmp/output.s.c.bus with whitelist sample2/10xv2_whitelist.txt
    [2020-01-14 23:01:45,237]    INFO Sorting BUS file tmp/output.s.c.bus to sample2/output.unfiltered.bus
    [2020-01-14 23:01:54,119]    INFO Generating count matrix sample2/counts_unfiltered/cells_x_genes from BUS file sample2/output.unfiltered.bus
    [2020-01-14 23:01:59,981]    INFO Converting matrix sample2/counts_unfiltered/cells_x_genes.mtx to h5ad sample2/counts_unfiltered/adata.h5ad
    [2020-01-14 23:02:03,635]    INFO Filtering with bustools
    [2020-01-14 23:02:03,635]    INFO Generating whitelist sample2/filter_barcodes.txt from BUS file sample2/output.unfiltered.bus
    [2020-01-14 23:02:03,803]    INFO Capturing records from BUS file sample2/output.unfiltered.bus to tmp/output.filtered.bus with capture list sample2/filter_barcodes.txt
    [2020-01-14 23:02:05,366]    INFO Sorting BUS file tmp/output.filtered.bus to sample2/output.filtered.bus
    [2020-01-14 23:02:12,500]    INFO Generating count matrix sample2/counts_filtered/cells_x_genes from BUS file sample2/output.filtered.bus
    [2020-01-14 23:02:17,853]    INFO Converting matrix sample2/counts_filtered/cells_x_genes.mtx to h5ad sample2/counts_filtered/adata.h5ad
    CPU times: user 626 ms, sys: 82.7 ms, total: 709 ms
    Wall time: 2min 7s


# Install `anndata`


```
!pip install anndata
```

    Requirement already satisfied: anndata in /usr/local/lib/python3.6/dist-packages (0.6.22.post1)
    Requirement already satisfied: numpy~=1.14 in /usr/local/lib/python3.6/dist-packages (from anndata) (1.17.5)
    Requirement already satisfied: scipy~=1.0 in /usr/local/lib/python3.6/dist-packages (from anndata) (1.4.1)
    Requirement already satisfied: pandas>=0.23.0 in /usr/local/lib/python3.6/dist-packages (from anndata) (0.25.3)
    Requirement already satisfied: natsort in /usr/local/lib/python3.6/dist-packages (from anndata) (5.5.0)
    Requirement already satisfied: h5py in /usr/local/lib/python3.6/dist-packages (from anndata) (2.8.0)
    Requirement already satisfied: pytz>=2017.2 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata) (2018.9)
    Requirement already satisfied: python-dateutil>=2.6.1 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata) (2.6.1)
    Requirement already satisfied: six in /usr/local/lib/python3.6/dist-packages (from h5py->anndata) (1.12.0)


# Read sample1 and sample2 gene counts into anndata


```
import anndata
sample1 = anndata.read_h5ad('sample1/counts_filtered/adata.h5ad')
sample2 = anndata.read_h5ad('sample2/counts_filtered/adata.h5ad')
```


```
sample1
```




    AnnData object with n_obs × n_vars = 1396 × 60623 




```
sample1.X
```




    <1396x60623 sparse matrix of type '<class 'numpy.float32'>'
    	with 4828398 stored elements in Compressed Sparse Row format>




```
sample1.obs.head()
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
    </tr>
    <tr>
      <th>index</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>AAACCTGAGCTGTTCA</th>
    </tr>
    <tr>
      <th>AAACCTGCAATTCCTT</th>
    </tr>
    <tr>
      <th>AAACCTGGTCTACCTC</th>
    </tr>
    <tr>
      <th>AAACCTGGTTTCCACC</th>
    </tr>
    <tr>
      <th>AAACCTGTCGGAGCAA</th>
    </tr>
  </tbody>
</table>
</div>




```
sample1.var.head()
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
    </tr>
    <tr>
      <th>index</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ENSG00000223972.5</th>
    </tr>
    <tr>
      <th>ENSG00000227232.5</th>
    </tr>
    <tr>
      <th>ENSG00000278267.1</th>
    </tr>
    <tr>
      <th>ENSG00000243485.5</th>
    </tr>
    <tr>
      <th>ENSG00000284332.1</th>
    </tr>
  </tbody>
</table>
</div>




```
sample2
```




    AnnData object with n_obs × n_vars = 279 × 60623 




```
sample2.X
```




    <279x60623 sparse matrix of type '<class 'numpy.float32'>'
    	with 1282741 stored elements in Compressed Sparse Row format>




```
sample2.obs.head()
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
    </tr>
    <tr>
      <th>index</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>AAACCTGAGACCACGA</th>
    </tr>
    <tr>
      <th>AAACCTGTCTGATACG</th>
    </tr>
    <tr>
      <th>AAACGGGAGTGTTGAA</th>
    </tr>
    <tr>
      <th>AAAGATGTCCGAACGC</th>
    </tr>
    <tr>
      <th>AAAGTAGGTTAGTGGG</th>
    </tr>
  </tbody>
</table>
</div>




```
sample2.var.head()
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
    </tr>
    <tr>
      <th>index</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ENSG00000223972.5</th>
    </tr>
    <tr>
      <th>ENSG00000227232.5</th>
    </tr>
    <tr>
      <th>ENSG00000278267.1</th>
    </tr>
    <tr>
      <th>ENSG00000243485.5</th>
    </tr>
    <tr>
      <th>ENSG00000284332.1</th>
    </tr>
  </tbody>
</table>
</div>



## Concatenate the anndatas


```
concat_samples = sample1.concatenate(
    sample2, join='outer', batch_categories=['sample1', 'sample2'], index_unique='-'
)
```


```
concat_samples
```




    AnnData object with n_obs × n_vars = 1675 × 60623 
        obs: 'batch'




```
concat_samples.var.head()
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
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ENSG00000223972.5</th>
    </tr>
    <tr>
      <th>ENSG00000227232.5</th>
    </tr>
    <tr>
      <th>ENSG00000278267.1</th>
    </tr>
    <tr>
      <th>ENSG00000243485.5</th>
    </tr>
    <tr>
      <th>ENSG00000284332.1</th>
    </tr>
  </tbody>
</table>
</div>




```
concat_samples.obs
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
      <th>batch</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>AAACCTGAGCTGTTCA-sample1</th>
      <td>sample1</td>
    </tr>
    <tr>
      <th>AAACCTGCAATTCCTT-sample1</th>
      <td>sample1</td>
    </tr>
    <tr>
      <th>AAACCTGGTCTACCTC-sample1</th>
      <td>sample1</td>
    </tr>
    <tr>
      <th>AAACCTGGTTTCCACC-sample1</th>
      <td>sample1</td>
    </tr>
    <tr>
      <th>AAACCTGTCGGAGCAA-sample1</th>
      <td>sample1</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
    </tr>
    <tr>
      <th>TTTGCGCGTTGACGTT-sample2</th>
      <td>sample2</td>
    </tr>
    <tr>
      <th>TTTGCGCGTTGTCTTT-sample2</th>
      <td>sample2</td>
    </tr>
    <tr>
      <th>TTTGGTTGTCATGCAT-sample2</th>
      <td>sample2</td>
    </tr>
    <tr>
      <th>TTTGTCAGTGAGTGAC-sample2</th>
      <td>sample2</td>
    </tr>
    <tr>
      <th>TTTGTCATCTTCATGT-sample2</th>
      <td>sample2</td>
    </tr>
  </tbody>
</table>
<p>1675 rows × 1 columns</p>
</div>




```

```
