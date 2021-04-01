<a href="https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/tutorials/docs/tutorials/kb_aggregate/python/kb_aggregating_count_matrices.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# Aggregating multiple count matrices tutorial

This tutorial describes how to aggregate multiple count matrices by concatenating them into a single [AnnData](https://anndata.readthedocs.io/en/latest/anndata.AnnData.html) object with batch labels for different samples.

This is similar to the Cell Ranger aggr function, however no normalization is performed. cellranger aggr is described at https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate

For this tutorial we use dataset E-MTAB-6108.

The notebook will take some time to run. To ensure that Google Colab does not shut down because of inactivity paste the following code into the console of this tab (*Cntrl [Mac: Cmd]  + Option + i  -> Console tab -> paste code -> press Enter*).

```javascript
function ConnectButton(){
    console.log("Connect pushed"); 
    document.querySelector("#top-toolbar > colab-connect-button").shadowRoot.querySelector("#connect").click() 
}
setInterval(ConnectButton,60000);
```

## Download the raw data

The raw data for E-MTAB-6108 is available at https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6108/


```
%%time
!wget -q https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample1_L005_R1.fastq.gz
!wget -q https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample1_L005_R2.fastq.gz
!wget -q https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample2_L005_R1.fastq.gz
!wget -q https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample2_L005_R2.fastq.gz
```

    CPU times: user 5.42 s, sys: 839 ms, total: 6.25 s
    Wall time: 15min 30s


## Install `kb`

Install `kb` for running the kallisto|bustools workflow.


```
!pip install --quiet kb-python
```

    [K     |████████████████████████████████| 59.1MB 77kB/s 
    [K     |████████████████████████████████| 10.3MB 34.3MB/s 
    [K     |████████████████████████████████| 13.2MB 50.1MB/s 
    [K     |████████████████████████████████| 51kB 5.6MB/s 
    [K     |████████████████████████████████| 81kB 6.8MB/s 
    [K     |████████████████████████████████| 112kB 56.7MB/s 
    [K     |████████████████████████████████| 71kB 6.7MB/s 
    [K     |████████████████████████████████| 1.2MB 50.2MB/s 
    [K     |████████████████████████████████| 51kB 5.0MB/s 
    [?25h  Building wheel for loompy (setup.py) ... [?25l[?25hdone
      Building wheel for sinfo (setup.py) ... [?25l[?25hdone
      Building wheel for umap-learn (setup.py) ... [?25l[?25hdone
      Building wheel for numpy-groupies (setup.py) ... [?25l[?25hdone
      Building wheel for pynndescent (setup.py) ... [?25l[?25hdone


## Download a pre-built human index

__Note:__ See [this notebook]() for a tutorial on how to build custom transcriptome or RNA velocity indices.


```
%%time
!kb ref -d human -i index.idx -g t2g.txt
```

    [2021-03-31 20:50:10,750]    INFO Downloading files for human from https://caltech.box.com/shared/static/v1nm7lpnqz5syh8dyzdk2zs8bglncfib.gz to tmp/v1nm7lpnqz5syh8dyzdk2zs8bglncfib.gz
    100% 2.23G/2.23G [01:35<00:00, 25.0MB/s]
    [2021-03-31 20:51:47,840]    INFO Extracting files from tmp/v1nm7lpnqz5syh8dyzdk2zs8bglncfib.gz
    CPU times: user 1.51 s, sys: 288 ms, total: 1.8 s
    Wall time: 2min 15s


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

    [2021-03-31 20:52:24,861]    INFO Using index index.idx to generate BUS file to sample1 from
    [2021-03-31 20:52:24,861]    INFO         iPSC_RGCscRNAseq_Sample1_L005_R1.fastq.gz
    [2021-03-31 20:52:24,861]    INFO         iPSC_RGCscRNAseq_Sample1_L005_R2.fastq.gz
    [2021-03-31 21:15:29,824]    INFO Sorting BUS file sample1/output.bus to sample1/tmp/output.s.bus
    [2021-03-31 21:16:21,667]    INFO Whitelist not provided
    [2021-03-31 21:16:21,668]    INFO Copying pre-packaged 10XV2 whitelist to sample1
    [2021-03-31 21:16:22,484]    INFO Inspecting BUS file sample1/tmp/output.s.bus
    [2021-03-31 21:16:36,519]    INFO Correcting BUS records in sample1/tmp/output.s.bus to sample1/tmp/output.s.c.bus with whitelist sample1/10xv2_whitelist.txt
    [2021-03-31 21:16:46,573]    INFO Sorting BUS file sample1/tmp/output.s.c.bus to sample1/output.unfiltered.bus
    [2021-03-31 21:17:25,239]    INFO Generating count matrix sample1/counts_unfiltered/cells_x_genes from BUS file sample1/output.unfiltered.bus
    [2021-03-31 21:17:48,992]    INFO Reading matrix sample1/counts_unfiltered/cells_x_genes.mtx
    [2021-03-31 21:18:00,111]    INFO Writing matrix to h5ad sample1/counts_unfiltered/adata.h5ad
    [2021-03-31 21:18:00,914]    INFO Filtering with bustools
    [2021-03-31 21:18:00,915]    INFO Generating whitelist sample1/filter_barcodes.txt from BUS file sample1/output.unfiltered.bus
    [2021-03-31 21:18:01,259]    INFO Correcting BUS records in sample1/output.unfiltered.bus to sample1/tmp/output.unfiltered.c.bus with whitelist sample1/filter_barcodes.txt
    [2021-03-31 21:18:09,292]    INFO Sorting BUS file sample1/tmp/output.unfiltered.c.bus to sample1/output.filtered.bus
    [2021-03-31 21:18:46,457]    INFO Generating count matrix sample1/counts_filtered/cells_x_genes from BUS file sample1/output.filtered.bus
    [2021-03-31 21:19:09,724]    INFO Reading matrix sample1/counts_filtered/cells_x_genes.mtx
    [2021-03-31 21:19:18,183]    INFO Writing matrix to h5ad sample1/counts_filtered/adata.h5ad
    CPU times: user 9.95 s, sys: 1.39 s, total: 11.3 s
    Wall time: 26min 56s


### Sample 2


```
%%time
!kb count -i index.idx -g t2g.txt -x 10xv2 -o sample2 --h5ad -t 2 --filter bustools \
iPSC_RGCscRNAseq_Sample2_L005_R1.fastq.gz \
iPSC_RGCscRNAseq_Sample2_L005_R2.fastq.gz
```

    [2021-03-31 21:19:22,185]    INFO Using index index.idx to generate BUS file to sample2 from
    [2021-03-31 21:19:22,185]    INFO         iPSC_RGCscRNAseq_Sample2_L005_R1.fastq.gz
    [2021-03-31 21:19:22,185]    INFO         iPSC_RGCscRNAseq_Sample2_L005_R2.fastq.gz
    [2021-03-31 21:37:11,095]    INFO Sorting BUS file sample2/output.bus to sample2/tmp/output.s.bus
    [2021-03-31 21:37:35,255]    INFO Whitelist not provided
    [2021-03-31 21:37:35,255]    INFO Copying pre-packaged 10XV2 whitelist to sample2
    [2021-03-31 21:37:35,379]    INFO Inspecting BUS file sample2/tmp/output.s.bus
    [2021-03-31 21:37:43,363]    INFO Correcting BUS records in sample2/tmp/output.s.bus to sample2/tmp/output.s.c.bus with whitelist sample2/10xv2_whitelist.txt
    [2021-03-31 21:37:47,960]    INFO Sorting BUS file sample2/tmp/output.s.c.bus to sample2/output.unfiltered.bus
    [2021-03-31 21:37:58,445]    INFO Generating count matrix sample2/counts_unfiltered/cells_x_genes from BUS file sample2/output.unfiltered.bus
    [2021-03-31 21:38:08,901]    INFO Reading matrix sample2/counts_unfiltered/cells_x_genes.mtx
    [2021-03-31 21:38:13,045]    INFO Writing matrix to h5ad sample2/counts_unfiltered/adata.h5ad
    [2021-03-31 21:38:13,797]    INFO Filtering with bustools
    [2021-03-31 21:38:13,798]    INFO Generating whitelist sample2/filter_barcodes.txt from BUS file sample2/output.unfiltered.bus
    [2021-03-31 21:38:13,965]    INFO Correcting BUS records in sample2/output.unfiltered.bus to sample2/tmp/output.unfiltered.c.bus with whitelist sample2/filter_barcodes.txt
    [2021-03-31 21:38:16,943]    INFO Sorting BUS file sample2/tmp/output.unfiltered.c.bus to sample2/output.filtered.bus
    [2021-03-31 21:38:25,772]    INFO Generating count matrix sample2/counts_filtered/cells_x_genes from BUS file sample2/output.filtered.bus
    [2021-03-31 21:38:33,900]    INFO Reading matrix sample2/counts_filtered/cells_x_genes.mtx
    [2021-03-31 21:38:36,553]    INFO Writing matrix to h5ad sample2/counts_filtered/adata.h5ad
    CPU times: user 7.29 s, sys: 1.04 s, total: 8.33 s
    Wall time: 19min 17s


# Install `anndata`


```
!pip install --quiet anndata
```

# Read sample1 and sample2 gene counts into anndata


```
import anndata
sample1 = anndata.read_h5ad('sample1/counts_filtered/adata.h5ad')
sample2 = anndata.read_h5ad('sample2/counts_filtered/adata.h5ad')
```


```
sample1
```




    AnnData object with n_obs × n_vars = 1424 × 60623
        var: 'gene_name'




```
sample1.X
```




    <1424x60623 sparse matrix of type '<class 'numpy.float32'>'
    	with 4829530 stored elements in Compressed Sparse Row format>




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
      <th>barcode</th>
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
      <th>AAACCTGTCCTATGTT</th>
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
      <th>gene_name</th>
    </tr>
    <tr>
      <th>gene_id</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ENSG00000223972.5</th>
      <td>DDX11L1</td>
    </tr>
    <tr>
      <th>ENSG00000227232.5</th>
      <td>WASH7P</td>
    </tr>
    <tr>
      <th>ENSG00000278267.1</th>
      <td>MIR6859-1</td>
    </tr>
    <tr>
      <th>ENSG00000243485.5</th>
      <td>MIR1302-2HG</td>
    </tr>
    <tr>
      <th>ENSG00000284332.1</th>
      <td>MIR1302-2</td>
    </tr>
  </tbody>
</table>
</div>




```
sample2
```




    AnnData object with n_obs × n_vars = 281 × 60623
        var: 'gene_name'




```
sample2.X
```




    <281x60623 sparse matrix of type '<class 'numpy.float32'>'
    	with 1282359 stored elements in Compressed Sparse Row format>




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
      <th>barcode</th>
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
      <th>AAAGATGTCAGAGACG</th>
    </tr>
    <tr>
      <th>AAAGATGTCCGAACGC</th>
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
      <th>gene_name</th>
    </tr>
    <tr>
      <th>gene_id</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ENSG00000223972.5</th>
      <td>DDX11L1</td>
    </tr>
    <tr>
      <th>ENSG00000227232.5</th>
      <td>WASH7P</td>
    </tr>
    <tr>
      <th>ENSG00000278267.1</th>
      <td>MIR6859-1</td>
    </tr>
    <tr>
      <th>ENSG00000243485.5</th>
      <td>MIR1302-2HG</td>
    </tr>
    <tr>
      <th>ENSG00000284332.1</th>
      <td>MIR1302-2</td>
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




    AnnData object with n_obs × n_vars = 1705 × 60623
        obs: 'batch'
        var: 'gene_name'




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
      <th>gene_name</th>
    </tr>
    <tr>
      <th>gene_id</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ENSG00000223972.5</th>
      <td>DDX11L1</td>
    </tr>
    <tr>
      <th>ENSG00000227232.5</th>
      <td>WASH7P</td>
    </tr>
    <tr>
      <th>ENSG00000278267.1</th>
      <td>MIR6859-1</td>
    </tr>
    <tr>
      <th>ENSG00000243485.5</th>
      <td>MIR1302-2HG</td>
    </tr>
    <tr>
      <th>ENSG00000284332.1</th>
      <td>MIR1302-2</td>
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
    <tr>
      <th>barcode</th>
      <th></th>
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
      <th>AAACCTGTCCTATGTT-sample1</th>
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
<p>1705 rows × 1 columns</p>
</div>




```

```
