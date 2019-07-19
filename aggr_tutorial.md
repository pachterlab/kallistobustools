# Aggregating multiple  count matrices

This tutorial describes how to aggregate multiple count matrices by concatenating them into a single [AnnData](https://anndata.readthedocs.io/en/latest/anndata.AnnData.html) object with batch labels for different samples. A notebook showing the entire workflow (including running kallisto and bsutools) is available [here](https://github.com/BUStools/getting_started/blob/master/aggr_tutorial.ipynb).

This is similar to the Cell Ranger `aggr` function, however no normalization is performed. `cellranger aggr` is described at https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate

For this tutorial we use dataset E-MTAB-6108. We provide the count matrices as an 80mb zip file at https://github.com/BUStools/getting_started/releases/download/aggr/E-MTAB-6108_sample1_sample2_genecounts.zip

The raw data for E-MTAB-6108 is available at https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6108/

**We assume that you downloaded the zip file and have your work directory structure with the files as shown below**
```
.
├── sample1
│   ├── genecounts
│   │   ├── genes.barcodes.txt
│   │   ├── genes.genes.txt
│   │   └── genes.mtx
│   ├── matrix.ec
│   ├── run_info.json
│   └── transcripts.txt
└── sample2
    ├── genecounts
    │   ├── genes.barcodes.txt
    │   ├── genes.genes.txt
    │   └── genes.mtx
    ├── matrix.ec
    ├── run_info.json
    └── transcripts.txt
```



## Imports
```python
from anndata import AnnData
import anndata
from scipy import sparse
import scipy
import anndata
import scipy.io
import os
```

### Read sample1 and sample2 gene count matrices into anndata and concatenate them

```python
## load sample1 on anndata as sparse crs matrix
sample1 = anndata.AnnData(scipy.io.mmread('./sample1/genecounts/genes.mtx').tocsr())
sample1.obs= pd.read_csv('./sample1/genecounts/genes.barcodes.txt', index_col = 0, header = None, names = ['barcode'])
sample1.var = pd.read_csv('./sample1/genecounts/genes.genes.txt', header = None, index_col = 0, names =['ensembl_id'], sep = '\t')
print('Loaded sample1 mtx:',sample1.X.shape)
```
```
>>> Loaded sample1 mtx: (226612, 35606)

```

```python
## load sample2 on anndata as sparse crs matrix
sample2 = anndata.AnnData(scipy.io.mmread('./sample2/genecounts/genes.mtx').tocsr())
sample2.obs= pd.read_csv('./sample2/genecounts/genes.barcodes.txt', index_col = 0, header = None, names = ['barcode'])
sample2.var = pd.read_csv('./sample2/genecounts/genes.genes.txt', header = None, index_col = 0, names =['ensembl_id'], sep = '\t')
print('Loaded sample2 mtx:',sample2.X.shape)
```

```
>>>  Loaded sample2 mtx: (135582, 35606)
```

```python
concat_samples = AnnData.concatenate(sample1, sample2, join='outer', batch_categories=['sample1','sample2'],index_unique='-')
```
The resuting AnnData `concat_samples`:
```
AnnData object with n_obs × n_vars = 362194 × 35606 
    obs: 'batch'
```
With the ensembl gene names under `concat_samples.var`
```
ENSG00000223972.5
ENSG00000227232.5
ENSG00000268020.3
ENSG00000240361.2
ENSG00000186092.6
```

And it will have the following observation labels under `concat_samples.obs`
```
                            batch
AAACCTGAGAAACCAT-sample1  sample1
AAACCTGAGAAACCGC-sample1  sample1
AAACCTGAGAAACCTA-sample1  sample1
AAACCTGAGAAACGAG-sample1  sample1
AAACCTGAGAAAGTGG-sample1  sample1
...                           ...
TTTGTCATCTTACCTA-sample2  sample2
TTTGTCATCTTCATGT-sample2  sample2
TTTGTCATCTTGCCGT-sample2  sample2
TTTGTCATCTTTACAC-sample2  sample2
TTTGTCATCTTTCCTC-sample2  sample2
```

