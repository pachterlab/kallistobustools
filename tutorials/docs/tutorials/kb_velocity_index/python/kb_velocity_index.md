# Constructing a velocity index with `kb`

This tutorial provides instructions for how to generate a velocity index to use with **kallisto | bustools** using `kb`.

## Download reference files

Download the genomic (DNA) FASTA and GTF annotations for your desired organism from the database of your choice. This tutorial uses mouse reference files downloaded from [Ensembl](https://uswest.ensembl.org/info/data/ftp/index.html).


```
%%time
!wget ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
!wget ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz
```

    --2020-01-16 00:33:28--  ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
               => ‘Mus_musculus.GRCm38.dna.primary_assembly.fa.gz’
    Resolving ftp.ensembl.org (ftp.ensembl.org)... 193.62.193.8
    Connecting to ftp.ensembl.org (ftp.ensembl.org)|193.62.193.8|:21... connected.
    Logging in as anonymous ... Logged in!
    ==> SYST ... done.    ==> PWD ... done.
    ==> TYPE I ... done.  ==> CWD (1) /pub/release-98/fasta/mus_musculus/dna ... done.
    ==> SIZE Mus_musculus.GRCm38.dna.primary_assembly.fa.gz ... 805984352
    ==> PASV ... done.    ==> RETR Mus_musculus.GRCm38.dna.primary_assembly.fa.gz ... done.
    Length: 805984352 (769M) (unauthoritative)
    
    Mus_musculus.GRCm38 100%[===================>] 768.65M  3.78MB/s    in 4m 49s  
    
    2020-01-16 00:38:18 (2.66 MB/s) - ‘Mus_musculus.GRCm38.dna.primary_assembly.fa.gz’ saved [805984352]
    
    --2020-01-16 00:38:18--  ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz
               => ‘Mus_musculus.GRCm38.98.gtf.gz’
    Resolving ftp.ensembl.org (ftp.ensembl.org)... 193.62.193.8
    Connecting to ftp.ensembl.org (ftp.ensembl.org)|193.62.193.8|:21... connected.
    Logging in as anonymous ... Logged in!
    ==> SYST ... done.    ==> PWD ... done.
    ==> TYPE I ... done.  ==> CWD (1) /pub/release-98/gtf/mus_musculus ... done.
    ==> SIZE Mus_musculus.GRCm38.98.gtf.gz ... 30256597
    ==> PASV ... done.    ==> RETR Mus_musculus.GRCm38.98.gtf.gz ... done.
    Length: 30256597 (29M) (unauthoritative)
    
    Mus_musculus.GRCm38 100%[===================>]  28.85M  17.6MB/s    in 1.6s    
    
    2020-01-16 00:38:21 (17.6 MB/s) - ‘Mus_musculus.GRCm38.98.gtf.gz’ saved [30256597]
    
    CPU times: user 2.91 s, sys: 500 ms, total: 3.41 s
    Wall time: 4min 53s


## Install `kb`


```
!pip install git+https://github.com/pachterlab/kb_python@count-kite
```

    Collecting git+https://github.com/pachterlab/kb_python@count-kite
      Cloning https://github.com/pachterlab/kb_python (to revision count-kite) to /tmp/pip-req-build-0qd_9um1
      Running command git clone -q https://github.com/pachterlab/kb_python /tmp/pip-req-build-0qd_9um1
      Running command git checkout -b count-kite --track origin/count-kite
      Switched to a new branch 'count-kite'
      Branch 'count-kite' set up to track remote branch 'count-kite' from 'origin'.
    Requirement already satisfied: anndata>=0.6.22.post1 in /usr/local/lib/python3.6/dist-packages (from kb-python==0.24.4) (0.6.22.post1)
    Requirement already satisfied: loompy>=3.0.6 in /usr/local/lib/python3.6/dist-packages (from kb-python==0.24.4) (3.0.6)
    Requirement already satisfied: requests>=2.19.0 in /usr/local/lib/python3.6/dist-packages (from kb-python==0.24.4) (2.21.0)
    Requirement already satisfied: tqdm>=4.39.0 in /usr/local/lib/python3.6/dist-packages (from kb-python==0.24.4) (4.41.1)
    Requirement already satisfied: h5py in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python==0.24.4) (2.8.0)
    Requirement already satisfied: pandas>=0.23.0 in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python==0.24.4) (0.25.3)
    Requirement already satisfied: scipy~=1.0 in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python==0.24.4) (1.4.1)
    Requirement already satisfied: natsort in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python==0.24.4) (5.5.0)
    Requirement already satisfied: numpy~=1.14 in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python==0.24.4) (1.17.5)
    Requirement already satisfied: setuptools in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python==0.24.4) (42.0.2)
    Requirement already satisfied: click in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python==0.24.4) (7.0)
    Requirement already satisfied: numpy-groupies in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python==0.24.4) (0+unknown)
    Requirement already satisfied: numba in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python==0.24.4) (0.47.0)
    Requirement already satisfied: certifi>=2017.4.17 in /usr/local/lib/python3.6/dist-packages (from requests>=2.19.0->kb-python==0.24.4) (2019.11.28)
    Requirement already satisfied: idna<2.9,>=2.5 in /usr/local/lib/python3.6/dist-packages (from requests>=2.19.0->kb-python==0.24.4) (2.8)
    Requirement already satisfied: chardet<3.1.0,>=3.0.2 in /usr/local/lib/python3.6/dist-packages (from requests>=2.19.0->kb-python==0.24.4) (3.0.4)
    Requirement already satisfied: urllib3<1.25,>=1.21.1 in /usr/local/lib/python3.6/dist-packages (from requests>=2.19.0->kb-python==0.24.4) (1.24.3)
    Requirement already satisfied: six in /usr/local/lib/python3.6/dist-packages (from h5py->anndata>=0.6.22.post1->kb-python==0.24.4) (1.12.0)
    Requirement already satisfied: pytz>=2017.2 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata>=0.6.22.post1->kb-python==0.24.4) (2018.9)
    Requirement already satisfied: python-dateutil>=2.6.1 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata>=0.6.22.post1->kb-python==0.24.4) (2.6.1)
    Requirement already satisfied: llvmlite>=0.31.0dev0 in /usr/local/lib/python3.6/dist-packages (from numba->loompy>=3.0.6->kb-python==0.24.4) (0.31.0)
    Building wheels for collected packages: kb-python
      Building wheel for kb-python (setup.py) ... [?25l[?25hdone
      Created wheel for kb-python: filename=kb_python-0.24.4-cp36-none-any.whl size=80991434 sha256=4dc5ecd507f54a3aa99e9f4862f2ea29dd9a0ce442af305bf0eff1177f35e22f
      Stored in directory: /tmp/pip-ephem-wheel-cache-7b9entel/wheels/8e/56/56/c89223de74af26792675e82f4bb5223e7cf0d653a33038e34c
    Successfully built kb-python
    Installing collected packages: kb-python
    Successfully installed kb-python-0.24.4


## Build the index

`kb` automatically splits the genome into cDNA and intron FASTA files. Because Google Colab has limited memory, we need to split the index into parts (here, we use `-n 4`). This will reduce the maximum memory `kb` uses, but the runtime of `kb count` will increase, which is a fair tradeoff in favor of less memory.


```
%%time
!kb ref -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa -c1 cdna_t2c.txt -c2 intron_t2c.txt --workflow lamanno -n 8 \
Mus_musculus.GRCm38.dna.primary_assembly.fa.gz \
Mus_musculus.GRCm38.98.gtf.gz
```

    [2020-01-16 03:31:03,222]    INFO Preparing Mus_musculus.GRCm38.dna.primary_assembly.fa.gz, Mus_musculus.GRCm38.98.gtf.gz
    [2020-01-16 03:31:03,222]    INFO Decompressing Mus_musculus.GRCm38.dna.primary_assembly.fa.gz to tmp
    [2020-01-16 03:31:30,853]    INFO Sorting tmp/Mus_musculus.GRCm38.dna.primary_assembly.fa to /content/tmp/tmpl3plby1k
    [2020-01-16 03:38:59,002]    INFO Decompressing Mus_musculus.GRCm38.98.gtf.gz to tmp
    [2020-01-16 03:39:03,235]    INFO Sorting tmp/Mus_musculus.GRCm38.98.gtf to /content/tmp/tmp7ebqamug
    [2020-01-16 03:40:00,940]    INFO Splitting genome tmp/Mus_musculus.GRCm38.dna.primary_assembly.fa into cDNA at /content/tmp/tmp1zp3oo9w
    [2020-01-16 03:40:00,940] WARNING The following chromosomes were found in the FASTA but doens't have any "transcript" features in the GTF: GL456213.1, JH584302.1, GL456392.1, JH584300.1, GL456396.1, GL456383.1, GL456389.1, GL456379.1, GL456378.1, GL456359.1, JH584301.1, GL456366.1, GL456360.1, GL456370.1, GL456368.1, GL456390.1, GL456393.1, GL456387.1, GL456382.1, GL456394.1, GL456367.1. No sequences will be generated for these chromosomes.
    [2020-01-16 03:41:14,163]    INFO Wrote 142446 cDNA transcripts
    [2020-01-16 03:41:14,168]    INFO Creating cDNA transcripts-to-capture at /content/tmp/tmpxjcopm_m
    [2020-01-16 03:41:15,248]    INFO Splitting genome into introns at /content/tmp/tmpmo19l2ry
    [2020-01-16 03:45:44,829]    INFO Wrote 647972 intron sequences
    [2020-01-16 03:45:44,836]    INFO Creating intron transcripts-to-capture at /content/tmp/tmprztg9wry
    [2020-01-16 03:46:51,138]    INFO Concatenating 1 cDNA FASTAs to cdna.fa
    [2020-01-16 03:46:55,954]    INFO Concatenating 1 cDNA transcripts-to-captures to cdna_t2c.txt
    [2020-01-16 03:46:56,028]    INFO Concatenating 1 intron FASTAs to intron.fa
    [2020-01-16 03:47:46,054]    INFO Concatenating 1 intron transcripts-to-captures to intron_t2c.txt
    [2020-01-16 03:47:46,396]    INFO Concatenating cDNA and intron FASTAs to /content/tmp/tmpn9ihglyn
    [2020-01-16 03:49:08,305]    INFO Creating transcript-to-gene mapping at t2g.txt
    [2020-01-16 03:50:17,445]    INFO Splitting /content/tmp/tmpn9ihglyn into 8 parts
    [2020-01-16 03:51:11,012]    INFO Indexing /content/tmp/tmphffimvvn to index.idx.0
    [2020-01-16 04:09:50,394]    INFO Indexing /content/tmp/tmpc8mhcmkw to index.idx.1
    [2020-01-16 04:25:30,897]    INFO Indexing /content/tmp/tmpcpy9f0cj to index.idx.2
    [2020-01-16 04:41:43,273]    INFO Indexing /content/tmp/tmp6bmc3x9s to index.idx.3
    [2020-01-16 04:57:03,671]    INFO Indexing /content/tmp/tmpj3nyxnh1 to index.idx.4
    [2020-01-16 05:11:56,815]    INFO Indexing /content/tmp/tmpwg6y882o to index.idx.5
    [2020-01-16 05:27:03,548]    INFO Indexing /content/tmp/tmp00l853b4 to index.idx.6
    [2020-01-16 05:42:58,221]    INFO Indexing /content/tmp/tmpijfn8tph to index.idx.7
    CPU times: user 50.3 s, sys: 6.14 s, total: 56.4 s
    Wall time: 2h 27min 11s



```

```
