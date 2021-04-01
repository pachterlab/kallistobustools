<a href="https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/tutorials/docs/tutorials/kb_kite/python/kb_kite.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# Pre-processing and analysis of feature barcode single-cell RNA-seq data with KITE.

In this notebook, we will perform pre-processing and analysis of [10x Genomics pbmc_1k_protein_v3](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_protein_v3) feature barcoding dataset using the **Kallisto Indexing and Tag Extraction (KITE)** workflow, implemented with a wrapper called `kb`. It was developed by Kyung Hoi (Joseph) Min and A. Sina Booeshaghi.

In Feature Barcoding assays, cellular data are recorded as short DNA sequences using procedures adapted from single-cell RNA-seq. The **KITE** workflow generates a "Mismatch Map" containing the sequences of all Feature Barcodes used in the experiment as well as all of their single-base mismatches. The Mismatch Map is used to produce transcipt-to-gene (.t2g) and fasta (.fa) files to be used as inputs for kallisto. An index is made with kallisto index, then kallisto | bustools effectively searches the sequencing data for the sequences in the Mismatch Map.

## Pre-processing

### Download the data

__Note:__ We use the `-O` option for `wget` to rename the files to easily identify them.


```
%%time
!wget -q https://caltech.box.com/shared/static/asmj4nu90ydhsrk3pm7aaxu00cnnfige.txt -O checksums.txt
!wget -q https://caltech.box.com/shared/static/mp2vr3p6dztdyatuag8ir3cektmrztg8.gz -O pbmc_1k_protein_v3_antibody_S2_L001_R1_001.fastq.gz
!wget -q https://caltech.box.com/shared/static/f3payi1za7mn0jfai7vm10sy3yqwgpqh.gz -O pbmc_1k_protein_v3_antibody_S2_L001_R2_001.fastq.gz
!wget -q https://caltech.box.com/shared/static/e112bbczh9o1rl6gfin36bqp0ga7uvdy.gz -O pbmc_1k_protein_v3_antibody_S2_L002_R1_001.fastq.gz
!wget -q https://caltech.box.com/shared/static/3ve2axc8dr8v5nnrhmynrdgpqj6xg42k.gz -O pbmc_1k_protein_v3_antibody_S2_L002_R2_001.fastq.gz
```

    CPU times: user 477 ms, sys: 86.2 ms, total: 563 ms
    Wall time: 1min


Then, we verify the integrity of the files we downloaded to make sure they were not corrupted during the download.


```
!md5sum -c checksums.txt --ignore-missing
```

    pbmc_1k_protein_v3_antibody_S2_L001_R1_001.fastq.gz: OK
    pbmc_1k_protein_v3_antibody_S2_L001_R2_001.fastq.gz: OK
    pbmc_1k_protein_v3_antibody_S2_L002_R1_001.fastq.gz: OK
    pbmc_1k_protein_v3_antibody_S2_L002_R2_001.fastq.gz: OK


### Install `kb`

Install `kb` for running the kallisto|bustools workflow.


```
!pip install --quiet git+https://github.com/pachterlab/kb_python@devel
```

    [K     |████████████████████████████████| 51kB 2.6MB/s 
    [K     |████████████████████████████████| 13.2MB 314kB/s 
    [K     |████████████████████████████████| 112kB 43.4MB/s 
    [?25h  Building wheel for kb-python (setup.py) ... [?25l[?25hdone
      Building wheel for loompy (setup.py) ... [?25l[?25hdone
      Building wheel for numpy-groupies (setup.py) ... [?25l[?25hdone


### Build the feature barcode mismatch index

`kb` is able to generate a FASTA file containing all hamming distance < 2 variants of the feature barcodes and create a kallisto index of these sequences. But it in order to do so, we first need to prepare a TSV containing feature barcode sequences in the first column and the feature barcode names in the second.

First, we download the feature reference file provided by 10x Genomics.


```
!wget -q http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_feature_ref.csv
```

Let's load it in as a Pandas DataFrame.


```
import pandas as pd

df = pd.read_csv('pbmc_1k_protein_v3_feature_ref.csv')
df
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
      <th>id</th>
      <th>name</th>
      <th>read</th>
      <th>pattern</th>
      <th>sequence</th>
      <th>feature_type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>CD3</td>
      <td>CD3_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>AACAAGACCCTTGAG</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>1</th>
      <td>CD4</td>
      <td>CD4_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>TACCCGTAATAGCGT</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>2</th>
      <td>CD8a</td>
      <td>CD8a_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>ATTGGCACTCAGATG</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>3</th>
      <td>CD14</td>
      <td>CD14_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>GAAAGTCAAAGCACT</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>4</th>
      <td>CD15</td>
      <td>CD15_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>ACGAATCAATCTGTG</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>5</th>
      <td>CD16</td>
      <td>CD16_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>GTCTTTGTCAGTGCA</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>6</th>
      <td>CD56</td>
      <td>CD56_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>GTTGTCCGACAATAC</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>7</th>
      <td>CD19</td>
      <td>CD19_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>TCAACGCTTGGCTAG</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>8</th>
      <td>CD25</td>
      <td>CD25_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>GTGCATTCAACAGTA</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>9</th>
      <td>CD45RA</td>
      <td>CD45RA_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>GATGAGAACAGGTTT</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>10</th>
      <td>CD45RO</td>
      <td>CD45RO_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>TGCATGTCATCGGTG</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>11</th>
      <td>PD-1</td>
      <td>PD-1_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>AAGTCGTGAGGCATG</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>12</th>
      <td>TIGIT</td>
      <td>TIGIT_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>TGAAGGCTCATTTGT</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>13</th>
      <td>CD127</td>
      <td>CD127_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>ACATTGACGCAACTA</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>14</th>
      <td>IgG2a</td>
      <td>IgG2a_control_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>CTCTATTCAGACCAG</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>15</th>
      <td>IgG1</td>
      <td>IgG1_control_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>ACTCACTGGAGTCTC</td>
      <td>Antibody Capture</td>
    </tr>
    <tr>
      <th>16</th>
      <td>IgG2b</td>
      <td>IgG2b_control_TotalSeqB</td>
      <td>R2</td>
      <td>5PNNNNNNNNNN(BC)NNNNNNNNN</td>
      <td>ATCACATCGTTGCCA</td>
      <td>Antibody Capture</td>
    </tr>
  </tbody>
</table>
</div>



We'll convert this dataframe into a TSV format that `kb` requires.


```
df[['sequence', 'id']].to_csv('features.tsv', index=None, header=None, sep='\t')
!cat features.tsv
```

    AACAAGACCCTTGAG	CD3
    TACCCGTAATAGCGT	CD4
    ATTGGCACTCAGATG	CD8a
    GAAAGTCAAAGCACT	CD14
    ACGAATCAATCTGTG	CD15
    GTCTTTGTCAGTGCA	CD16
    GTTGTCCGACAATAC	CD56
    TCAACGCTTGGCTAG	CD19
    GTGCATTCAACAGTA	CD25
    GATGAGAACAGGTTT	CD45RA
    TGCATGTCATCGGTG	CD45RO
    AAGTCGTGAGGCATG	PD-1
    TGAAGGCTCATTTGT	TIGIT
    ACATTGACGCAACTA	CD127
    CTCTATTCAGACCAG	IgG2a
    ACTCACTGGAGTCTC	IgG1
    ATCACATCGTTGCCA	IgG2b


Finally, we use `kb` to generate the mismatch kallisto index.


```
!kb ref -i mismatch.idx -f1 mismatch.fa -g t2g.txt --workflow kite features.tsv
```

    [2021-03-31 23:30:15,970]    INFO Generating mismatch FASTA at mismatch.fa
    [2021-03-31 23:30:15,982]    INFO Creating transcript-to-gene mapping at t2g.txt
    [2021-03-31 23:30:15,986]    INFO Indexing mismatch.fa to mismatch.idx


### Generate a feature count matrix in H5AD format

The following command will generate an RNA count matrix of cells (rows) by genes (columns) in H5AD format, which is a binary format used to store [Anndata](https://anndata.readthedocs.io/en/stable/) objects. Notice we are providing the index and transcript-to-gene mapping we generated in the previous step to the `-i` and `-g` arguments respectively. Also, these reads were generated with the 10x Genomics Chromium Single Cell v3 Chemistry, hence the `-x 10xv3` argument. To view other supported technologies, run `kb --list`.

__Note:__ If you would like a Loom file instead, replace the `--h5ad` flag with `--loom`. If you want to use the raw matrix output by `kb` instead of their H5AD or Loom converted files, omit these flags.


```
%%time
!kb count --h5ad -i mismatch.idx -g t2g.txt -x 10xv3 --workflow kite -t 2 \
pbmc_1k_protein_v3_antibody_S2_L001_R1_001.fastq.gz \
pbmc_1k_protein_v3_antibody_S2_L001_R2_001.fastq.gz \
pbmc_1k_protein_v3_antibody_S2_L002_R1_001.fastq.gz \
pbmc_1k_protein_v3_antibody_S2_L002_R2_001.fastq.gz
```

    [2021-03-31 23:30:17,474]    INFO Using index mismatch.idx to generate BUS file to . from
    [2021-03-31 23:30:17,474]    INFO         pbmc_1k_protein_v3_antibody_S2_L001_R1_001.fastq.gz
    [2021-03-31 23:30:17,474]    INFO         pbmc_1k_protein_v3_antibody_S2_L001_R2_001.fastq.gz
    [2021-03-31 23:30:17,474]    INFO         pbmc_1k_protein_v3_antibody_S2_L002_R1_001.fastq.gz
    [2021-03-31 23:30:17,474]    INFO         pbmc_1k_protein_v3_antibody_S2_L002_R2_001.fastq.gz
    [2021-03-31 23:32:02,567]    INFO Sorting BUS file ./output.bus to ./tmp/output.s.bus
    [2021-03-31 23:32:18,183]    INFO Whitelist not provided
    [2021-03-31 23:32:18,183]    INFO Copying pre-packaged 10XV3 whitelist to .
    [2021-03-31 23:32:19,157]    INFO Inspecting BUS file ./tmp/output.s.bus
    [2021-03-31 23:32:31,284]    INFO Correcting BUS records in ./tmp/output.s.bus to ./tmp/output.s.c.bus with whitelist ./10xv3_whitelist.txt
    [2021-03-31 23:32:49,746]    INFO Sorting BUS file ./tmp/output.s.c.bus to ./output.unfiltered.bus
    [2021-03-31 23:33:01,416]    INFO Generating count matrix ./counts_unfiltered/cells_x_features from BUS file ./output.unfiltered.bus
    [2021-03-31 23:33:03,924]    INFO Reading matrix ./counts_unfiltered/cells_x_features.mtx
    [2021-03-31 23:33:05,698]    INFO Writing matrix to h5ad ./counts_unfiltered/adata.h5ad
    CPU times: user 1.3 s, sys: 180 ms, total: 1.48 s
    Wall time: 2min 49s


## Analysis

In this part of the tutorial, we will load the RNA count matrix generated by `kb count` into Python and cluster the cells with Leiden.

### Install packages

Google Colab does not come with `Scanpy`, `python-igraph`, or `louvain` (but comes with `matplotlib`, `numpy`, `pandas`, and `scipy`).


```
!pip --quiet install leidenalg scanpy MulticoreTSNE
```

### Import packages


```
import anndata
import numpy as np
import scanpy as sc
```


```
adata = anndata.read_h5ad('counts_unfiltered/adata.h5ad')
```


```
adata
```




    AnnData object with n_obs × n_vars = 124716 × 17
        var: 'feature_name'




```
adata.obs.head()
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
      <th>AAACCCAAGAAACCCA</th>
    </tr>
    <tr>
      <th>AAACCCAAGACGAGGA</th>
    </tr>
    <tr>
      <th>AAACCCAAGAGTGTGT</th>
    </tr>
    <tr>
      <th>AAACCCAAGAGTGTTG</th>
    </tr>
    <tr>
      <th>AAACCCAAGATAGCAC</th>
    </tr>
  </tbody>
</table>
</div>




```
adata.var
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
      <th>feature_name</th>
    </tr>
    <tr>
      <th>feature_id</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>CD3</th>
      <td>CD3</td>
    </tr>
    <tr>
      <th>CD4</th>
      <td>CD4</td>
    </tr>
    <tr>
      <th>CD8a</th>
      <td>CD8a</td>
    </tr>
    <tr>
      <th>CD14</th>
      <td>CD14</td>
    </tr>
    <tr>
      <th>CD15</th>
      <td>CD15</td>
    </tr>
    <tr>
      <th>CD16</th>
      <td>CD16</td>
    </tr>
    <tr>
      <th>CD56</th>
      <td>CD56</td>
    </tr>
    <tr>
      <th>CD19</th>
      <td>CD19</td>
    </tr>
    <tr>
      <th>CD25</th>
      <td>CD25</td>
    </tr>
    <tr>
      <th>CD45RA</th>
      <td>CD45RA</td>
    </tr>
    <tr>
      <th>CD45RO</th>
      <td>CD45RO</td>
    </tr>
    <tr>
      <th>PD-1</th>
      <td>PD-1</td>
    </tr>
    <tr>
      <th>TIGIT</th>
      <td>TIGIT</td>
    </tr>
    <tr>
      <th>CD127</th>
      <td>CD127</td>
    </tr>
    <tr>
      <th>IgG2a</th>
      <td>IgG2a</td>
    </tr>
    <tr>
      <th>IgG1</th>
      <td>IgG1</td>
    </tr>
    <tr>
      <th>IgG2b</th>
      <td>IgG2b</td>
    </tr>
  </tbody>
</table>
</div>



### Plot counts


```
sc.pp.filter_cells(adata, min_counts=0)
```


```
sc.pp.filter_genes(adata, min_counts=0)
```


```
sc.pl.violin(adata, keys='n_counts')
```


![png](kb_kite_files/kb_kite_32_0.png)



```
adata.obs['n_countslog'] = np.log1p(adata.obs['n_counts'])
```


```
sc.pl.violin(adata, keys='n_countslog')
```


![png](kb_kite_files/kb_kite_34_0.png)



```
adata.obs.index
```




    Index(['AAACCCAAGAAACCCA', 'AAACCCAAGACGAGGA', 'AAACCCAAGAGTGTGT',
           'AAACCCAAGAGTGTTG', 'AAACCCAAGATAGCAC', 'AAACCCAAGATGAGTC',
           'AAACCCAAGATGGTAC', 'AAACCCAAGATTCGTT', 'AAACCCAAGATTTGGG',
           'AAACCCAAGCAAGCAT',
           ...
           'TTTGTTGTCGTCGCCT', 'TTTGTTGTCGTTGACG', 'TTTGTTGTCTAACCGG',
           'TTTGTTGTCTATGTAG', 'TTTGTTGTCTCAACAA', 'TTTGTTGTCTCACTCA',
           'TTTGTTGTCTCTTCGA', 'TTTGTTGTCTCTTGGT', 'TTTGTTGTCTGCACTT',
           'TTTGTTGTCTGCGACA'],
          dtype='object', name='barcode', length=124716)




```
sc.pp.filter_cells(adata, min_counts=1000)
sc.pl.violin(adata, keys='n_countslog', title="kallisto UMI counts")
adata
```


![png](kb_kite_files/kb_kite_36_0.png)





    AnnData object with n_obs × n_vars = 725 × 17
        obs: 'n_counts', 'n_countslog'
        var: 'feature_name', 'n_counts'



Here are violin plots for each Feature Barcode (antibody-oligo conjugates, x-axis) across all cells.


```
sc.pl.violin(adata, keys=list(adata.var.index)[-17:], xlabel='kallisto')
```


![png](kb_kite_files/kb_kite_38_0.png)


### Cluster with Leiden


```
sc.pp.normalize_per_cell(adata, counts_per_cell_after=10000)
```


```
sc.pp.neighbors(adata)
sc.tl.umap(adata)
```


```
sc.tl.leiden(adata, resolution=0.05)
```


```
sc.pl.umap(adata, color='leiden', palette='tab10')
```


![png](kb_kite_files/kb_kite_43_0.png)


### Embedding and Antibody Quantification


```
sc.pl.umap(adata, color=adata.var.index)
```


![png](kb_kite_files/kb_kite_45_0.png)



```
sc.pl.violin(adata, keys=list(adata.var.index[:2]), groupby='leiden')
```


![png](kb_kite_files/kb_kite_46_0.png)



```

```
