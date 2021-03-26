# Constructing a transcriptome index with `kb`

This tutorial provides instructions for how to generate a transcriptome index to use with **kallisto | bustools** using `kb`.

## Download reference files

Download the genomic (DNA) FASTA and GTF annotations for your desired organism from the database of your choice. This tutorial uses mouse reference files downloaded from [Ensembl](https://uswest.ensembl.org/info/data/ftp/index.html).


```
%%time
!wget ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
!wget ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz
```

    --2020-01-13 22:23:27--  ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
               => ‘Mus_musculus.GRCm38.dna.primary_assembly.fa.gz’
    Resolving ftp.ensembl.org (ftp.ensembl.org)... 193.62.193.8
    Connecting to ftp.ensembl.org (ftp.ensembl.org)|193.62.193.8|:21... connected.
    Logging in as anonymous ... Logged in!
    ==> SYST ... done.    ==> PWD ... done.
    ==> TYPE I ... done.  ==> CWD (1) /pub/release-98/fasta/mus_musculus/dna ... done.
    ==> SIZE Mus_musculus.GRCm38.dna.primary_assembly.fa.gz ... 805984352
    ==> PASV ... done.    ==> RETR Mus_musculus.GRCm38.dna.primary_assembly.fa.gz ... done.
    Length: 805984352 (769M) (unauthoritative)
    
    Mus_musculus.GRCm38 100%[===================>] 768.65M  95.3MB/s    in 8.0s    
    
    2020-01-13 22:23:35 (95.9 MB/s) - ‘Mus_musculus.GRCm38.dna.primary_assembly.fa.gz’ saved [805984352]
    
    --2020-01-13 22:23:37--  ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz
               => ‘Mus_musculus.GRCm38.98.gtf.gz’
    Resolving ftp.ensembl.org (ftp.ensembl.org)... 193.62.193.8
    Connecting to ftp.ensembl.org (ftp.ensembl.org)|193.62.193.8|:21... connected.
    Logging in as anonymous ... Logged in!
    ==> SYST ... done.    ==> PWD ... done.
    ==> TYPE I ... done.  ==> CWD (1) /pub/release-98/gtf/mus_musculus ... done.
    ==> SIZE Mus_musculus.GRCm38.98.gtf.gz ... 30256597
    ==> PASV ... done.    ==> RETR Mus_musculus.GRCm38.98.gtf.gz ... done.
    Length: 30256597 (29M) (unauthoritative)
    
    Mus_musculus.GRCm38 100%[===================>]  28.85M  76.4MB/s    in 0.4s    
    
    2020-01-13 22:23:37 (76.4 MB/s) - ‘Mus_musculus.GRCm38.98.gtf.gz’ saved [30256597]
    
    CPU times: user 130 ms, sys: 21.8 ms, total: 152 ms
    Wall time: 11.2 s


## Install `kb`


```
!pip install kb-python
```

    Collecting kb-python
    [?25l  Downloading https://files.pythonhosted.org/packages/62/c9/2e5b8fa2cd873a23ae1aeb128b33165d6a9387a2f56ea1fafec1d6d32477/kb_python-0.24.4-py3-none-any.whl (35.4MB)
    [K     |████████████████████████████████| 35.4MB 119kB/s 
    [?25hCollecting loompy>=3.0.6
    [?25l  Downloading https://files.pythonhosted.org/packages/36/52/74ed37ae5988522fbf87b856c67c4f80700e6452410b4cd80498c5f416f9/loompy-3.0.6.tar.gz (41kB)
    [K     |████████████████████████████████| 51kB 202kB/s 
    [?25hCollecting anndata>=0.6.22.post1
    [?25l  Downloading https://files.pythonhosted.org/packages/2b/72/87196c15f68d9865c31a43a10cf7c50bcbcedd5607d09f9aada0b3963103/anndata-0.6.22.post1-py3-none-any.whl (47kB)
    [K     |████████████████████████████████| 51kB 6.8MB/s 
    [?25hRequirement already satisfied: h5py in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (2.8.0)
    Requirement already satisfied: numpy in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (1.17.5)
    Requirement already satisfied: scipy in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (1.4.1)
    Requirement already satisfied: setuptools in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (42.0.2)
    Requirement already satisfied: numba in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (0.47.0)
    Requirement already satisfied: click in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (7.0)
    Collecting numpy-groupies
    [?25l  Downloading https://files.pythonhosted.org/packages/57/ae/18217b57ba3e4bb8a44ecbfc161ed065f6d1b90c75d404bd6ba8d6f024e2/numpy_groupies-0.9.10.tar.gz (43kB)
    [K     |████████████████████████████████| 51kB 7.5MB/s 
    [?25hRequirement already satisfied: natsort in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python) (5.5.0)
    Requirement already satisfied: pandas>=0.23.0 in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python) (0.25.3)
    Requirement already satisfied: six in /usr/local/lib/python3.6/dist-packages (from h5py->loompy>=3.0.6->kb-python) (1.12.0)
    Requirement already satisfied: llvmlite>=0.31.0dev0 in /usr/local/lib/python3.6/dist-packages (from numba->loompy>=3.0.6->kb-python) (0.31.0)
    Requirement already satisfied: pytz>=2017.2 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata>=0.6.22.post1->kb-python) (2018.9)
    Requirement already satisfied: python-dateutil>=2.6.1 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata>=0.6.22.post1->kb-python) (2.6.1)
    Building wheels for collected packages: loompy, numpy-groupies
      Building wheel for loompy (setup.py) ... [?25l[?25hdone
      Created wheel for loompy: filename=loompy-3.0.6-cp36-none-any.whl size=47896 sha256=dfb1f98246b6c01b60795a6ac413732febd33adc6e87b0a96573ae9e4a300484
      Stored in directory: /root/.cache/pip/wheels/f9/a4/90/5a98ad83419732b0fba533b81a2a52ba3dbe230a936ca4cdc9
      Building wheel for numpy-groupies (setup.py) ... [?25l[?25hdone
      Created wheel for numpy-groupies: filename=numpy_groupies-0+unknown-cp36-none-any.whl size=28044 sha256=5f21c80ba168876f1f6deb4bfe4fdd037dc7c0a8c447417663ebe14ed33a0c53
      Stored in directory: /root/.cache/pip/wheels/30/ac/83/64d5f9293aeaec63f9539142fc629a41af064cae1b3d8d94aa
    Successfully built loompy numpy-groupies
    Installing collected packages: numpy-groupies, loompy, anndata, kb-python
    Successfully installed anndata-0.6.22.post1 kb-python-0.24.4 loompy-3.0.6 numpy-groupies-0+unknown


## Build the index

`kb` automatically splits the genome into a cDNA FASTA file and uses that to build a kallisto index.


```
%%time
!kb ref -i transcriptome.idx -g transcripts_to_genes.txt -f1 cdna.fa \
Mus_musculus.GRCm38.dna.primary_assembly.fa.gz \
Mus_musculus.GRCm38.98.gtf.gz
```

    [2020-01-13 22:26:03,328]    INFO Decompressing Mus_musculus.GRCm38.98.gtf.gz to tmp
    [2020-01-13 22:26:06,666]    INFO Creating transcript-to-gene mapping at transcripts_to_genes.txt
    [2020-01-13 22:26:36,616]    INFO Decompressing Mus_musculus.GRCm38.dna.primary_assembly.fa.gz to tmp
    [2020-01-13 22:27:01,812]    INFO Sorting tmp/Mus_musculus.GRCm38.dna.primary_assembly.fa
    [2020-01-13 22:34:09,011]    INFO Sorting tmp/Mus_musculus.GRCm38.98.gtf
    [2020-01-13 22:34:57,702]    INFO Splitting genome into cDNA at cdna.fa
    [2020-01-13 22:36:02,242]    INFO Indexing to transcriptome.idx
    CPU times: user 5.08 s, sys: 886 ms, total: 5.97 s
    Wall time: 19min 38s



```

```
