# Processing Multiple Lanes at Once

This tutorial provides instructions for how to pre-process the [mouse T cells SRR8206317](https://www.ncbi.nlm.nih.gov/sra/?term=SRR8206317) dataset from [Miller & Sen et al., 2019](https://doi.org/10.1038/s41590-019-0312-6) using the **kallisto | bustools** workflow.

## Download the data


```
!wget ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR820/SRR8206317/d10_Tet_possorted_genome_bam.bam
```

    --2020-01-14 22:38:35--  ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR820/SRR8206317/d10_Tet_possorted_genome_bam.bam
               => ‘d10_Tet_possorted_genome_bam.bam’
    Resolving ftp.sra.ebi.ac.uk (ftp.sra.ebi.ac.uk)... 193.62.192.7
    Connecting to ftp.sra.ebi.ac.uk (ftp.sra.ebi.ac.uk)|193.62.192.7|:21... connected.
    Logging in as anonymous ... Logged in!
    ==> SYST ... done.    ==> PWD ... done.
    ==> TYPE I ... done.  ==> CWD (1) /vol1/run/SRR820/SRR8206317 ... 
    Error in server response, closing control connection.
    Retrying.
    
    --2020-01-14 22:43:41--  ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR820/SRR8206317/d10_Tet_possorted_genome_bam.bam
      (try: 2) => ‘d10_Tet_possorted_genome_bam.bam’
    Connecting to ftp.sra.ebi.ac.uk (ftp.sra.ebi.ac.uk)|193.62.192.7|:21... connected.
    Error in server response. Closing.
    Retrying.
    
    --2020-01-14 22:43:43--  ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR820/SRR8206317/d10_Tet_possorted_genome_bam.bam
      (try: 3) => ‘d10_Tet_possorted_genome_bam.bam’
    Connecting to ftp.sra.ebi.ac.uk (ftp.sra.ebi.ac.uk)|193.62.192.7|:21... connected.
    Logging in as anonymous ... Logged in!
    ==> SYST ... done.    ==> PWD ... done.
    ==> TYPE I ... done.  ==> CWD (1) /vol1/run/SRR820/SRR8206317 ... done.
    ==> SIZE d10_Tet_possorted_genome_bam.bam ... 9998262104
    ==> PASV ... done.    ==> RETR d10_Tet_possorted_genome_bam.bam ... done.
    Length: 9998262104 (9.3G) (unauthoritative)
    
    d10_Tet_possorted_g 100%[===================>]   9.31G  17.2MB/s    in 52m 55s 
    
    2020-01-14 23:37:02 (3.00 MB/s) - Control connection closed.
    Retrying.
    
    --2020-01-14 23:52:05--  ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR820/SRR8206317/d10_Tet_possorted_genome_bam.bam
      (try: 4) => ‘d10_Tet_possorted_genome_bam.bam’
    Connecting to ftp.sra.ebi.ac.uk (ftp.sra.ebi.ac.uk)|193.62.192.7|:21... connected.
    Logging in as anonymous ... Logged in!
    ==> SYST ... done.    ==> PWD ... done.
    ==> TYPE I ... done.  ==> CWD (1) /vol1/run/SRR820/SRR8206317 ... done.
    ==> SIZE d10_Tet_possorted_genome_bam.bam ... 9998262104
    File has already been retrieved.
    2020-01-14 23:52:07 (0.00 B/s) - ‘d10_Tet_possorted_genome_bam.bam’ saved [9998262104]
    


## Install `kb` and `bamtofastq`

We will be using `bamtofastq` to generate the original FASTQ files from the BAM files provided by the authors.


```
!pip install kb-python
```

    Collecting kb-python
    [?25l  Downloading https://files.pythonhosted.org/packages/62/c9/2e5b8fa2cd873a23ae1aeb128b33165d6a9387a2f56ea1fafec1d6d32477/kb_python-0.24.4-py3-none-any.whl (35.4MB)
    [K     |████████████████████████████████| 35.4MB 118kB/s 
    [?25hCollecting anndata>=0.6.22.post1
    [?25l  Downloading https://files.pythonhosted.org/packages/2b/72/87196c15f68d9865c31a43a10cf7c50bcbcedd5607d09f9aada0b3963103/anndata-0.6.22.post1-py3-none-any.whl (47kB)
    [K     |████████████████████████████████| 51kB 6.4MB/s 
    [?25hCollecting loompy>=3.0.6
    [?25l  Downloading https://files.pythonhosted.org/packages/36/52/74ed37ae5988522fbf87b856c67c4f80700e6452410b4cd80498c5f416f9/loompy-3.0.6.tar.gz (41kB)
    [K     |████████████████████████████████| 51kB 6.6MB/s 
    [?25hRequirement already satisfied: pandas>=0.23.0 in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python) (0.25.3)
    Requirement already satisfied: natsort in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python) (5.5.0)
    Requirement already satisfied: h5py in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python) (2.8.0)
    Requirement already satisfied: scipy~=1.0 in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python) (1.4.1)
    Requirement already satisfied: numpy~=1.14 in /usr/local/lib/python3.6/dist-packages (from anndata>=0.6.22.post1->kb-python) (1.17.5)
    Requirement already satisfied: setuptools in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (42.0.2)
    Requirement already satisfied: numba in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (0.47.0)
    Requirement already satisfied: click in /usr/local/lib/python3.6/dist-packages (from loompy>=3.0.6->kb-python) (7.0)
    Collecting numpy-groupies
    [?25l  Downloading https://files.pythonhosted.org/packages/57/ae/18217b57ba3e4bb8a44ecbfc161ed065f6d1b90c75d404bd6ba8d6f024e2/numpy_groupies-0.9.10.tar.gz (43kB)
    [K     |████████████████████████████████| 51kB 6.7MB/s 
    [?25hRequirement already satisfied: python-dateutil>=2.6.1 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata>=0.6.22.post1->kb-python) (2.6.1)
    Requirement already satisfied: pytz>=2017.2 in /usr/local/lib/python3.6/dist-packages (from pandas>=0.23.0->anndata>=0.6.22.post1->kb-python) (2018.9)
    Requirement already satisfied: six in /usr/local/lib/python3.6/dist-packages (from h5py->anndata>=0.6.22.post1->kb-python) (1.12.0)
    Requirement already satisfied: llvmlite>=0.31.0dev0 in /usr/local/lib/python3.6/dist-packages (from numba->loompy>=3.0.6->kb-python) (0.31.0)
    Building wheels for collected packages: loompy, numpy-groupies
      Building wheel for loompy (setup.py) ... [?25l[?25hdone
      Created wheel for loompy: filename=loompy-3.0.6-cp36-none-any.whl size=47896 sha256=fdf56d2f34f5776d33b9aea58f602517bfa0c29689bdd3de852d912e2d3444fd
      Stored in directory: /root/.cache/pip/wheels/f9/a4/90/5a98ad83419732b0fba533b81a2a52ba3dbe230a936ca4cdc9
      Building wheel for numpy-groupies (setup.py) ... [?25l[?25hdone
      Created wheel for numpy-groupies: filename=numpy_groupies-0+unknown-cp36-none-any.whl size=28044 sha256=7e582790055188f3b8d5132cbae763ebe112e210588928a8ed09dfa07140f063
      Stored in directory: /root/.cache/pip/wheels/30/ac/83/64d5f9293aeaec63f9539142fc629a41af064cae1b3d8d94aa
    Successfully built loompy numpy-groupies
    Installing collected packages: anndata, numpy-groupies, loompy, kb-python
    Successfully installed anndata-0.6.22.post1 kb-python-0.24.4 loompy-3.0.6 numpy-groupies-0+unknown



```
!wget http://cf.10xgenomics.com/misc/bamtofastq-1.2.0
!chmod +x bamtofastq-1.2.0
```

    --2020-01-14 22:29:02--  http://cf.10xgenomics.com/misc/bamtofastq-1.2.0
    Resolving cf.10xgenomics.com (cf.10xgenomics.com)... 13.224.29.56, 13.224.29.52, 13.224.29.102, ...
    Connecting to cf.10xgenomics.com (cf.10xgenomics.com)|13.224.29.56|:80... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 13288280 (13M) [binary/octet-stream]
    Saving to: ‘bamtofastq-1.2.0.1’
    
    bamtofastq-1.2.0.1  100%[===================>]  12.67M  52.9MB/s    in 0.2s    
    
    2020-01-14 22:29:02 (52.9 MB/s) - ‘bamtofastq-1.2.0.1’ saved [13288280/13288280]
    



```
!./bamtofastq-1.2.0
```

    bamtofastq v1.2.0
    Invalid arguments.
    
    Usage:
      bamtofastq [options] <bam> <output-path>
      bamtofastq (-h | --help)


## Download a pre-built mouse index


```
%%time
!kb ref -d mouse -i index.idx -g t2g.txt
```

    [2020-01-15 00:45:51,708]    INFO Downloading files for mouse from https://caltech.box.com/shared/static/vcaz6cujop0xuapdmz0pplp3aoqc41si.gz to tmp/vcaz6cujop0xuapdmz0pplp3aoqc41si.gz
    [2020-01-15 00:47:16,807]    INFO Extracting files from tmp/vcaz6cujop0xuapdmz0pplp3aoqc41si.gz
    CPU times: user 649 ms, sys: 118 ms, total: 767 ms
    Wall time: 2min 3s


## Generate the FASTQs from the BAM file

Use the `bamtofastq` utility to generate the FASTQs.


```
%%time
!./bamtofastq-1.2.0 --reads-per-fastq=500000000 d10_Tet_possorted_genome_bam.bam ./fastqs
```

    bamtofastq v1.2.0
    Args { arg_bam: "d10_Tet_possorted_genome_bam.bam", arg_output_path: "./fastqs", flag_nthreads: 4, flag_locus: None, flag_bx_list: None, flag_reads_per_fastq: 500000000, flag_gemcode: false, flag_lr20: false, flag_cr11: false }
    Writing finished.  Observed 85992089 read pairs. Wrote 85992089 read pairs
    CPU times: user 3.56 s, sys: 359 ms, total: 3.92 s
    Wall time: 13min 3s


## Generate an RNA count matrix in H5AD Format

The following command will generate an RNA count matrix of cells (rows) by genes (columns) in H5AD format, which is a binary format used to store [Anndata](https://anndata.readthedocs.io/en/stable/) objects. Notice we are providing the index and transcript-to-gene mapping we downloaded in the previous step to the `-i` and `-g` arguments respectively. Also, these reads were generated with the 10x Genomics Chromium Single Cell v2 Chemistry, hence the `-x 10xv2` argument. To view other supported technologies, run `kb --list`.

__Note:__ If you would like a Loom file instead, replace the `--h5ad` flag with `--loom`. If you want to use the raw matrix output by `kb` instead of their H5AD or Loom converted files, omit these flags.


```
!kb count -i index.idx -g t2g.txt -x 10xv2 -o output -t 2 \
fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L001_R1_001.fastq.gz \
fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L001_R2_001.fastq.gz \
fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L002_R1_001.fastq.gz \
fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L002_R2_001.fastq.gz \
fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L003_R1_001.fastq.gz \
fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L003_R2_001.fastq.gz \
fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L004_R1_001.fastq.gz \
fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L004_R2_001.fastq.gz
```

    [2020-01-15 01:31:02,446]    INFO Generating BUS file from
    [2020-01-15 01:31:02,446]    INFO         fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L001_R1_001.fastq.gz
    [2020-01-15 01:31:02,447]    INFO         fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L001_R2_001.fastq.gz
    [2020-01-15 01:31:02,447]    INFO         fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L002_R1_001.fastq.gz
    [2020-01-15 01:31:02,447]    INFO         fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L002_R2_001.fastq.gz
    [2020-01-15 01:31:02,447]    INFO         fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L003_R1_001.fastq.gz
    [2020-01-15 01:31:02,447]    INFO         fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L003_R2_001.fastq.gz
    [2020-01-15 01:31:02,447]    INFO         fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L004_R1_001.fastq.gz
    [2020-01-15 01:31:02,447]    INFO         fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L004_R2_001.fastq.gz
    [2020-01-15 01:37:42,296]    INFO Sorting BUS file output/output.bus to tmp/output.s.bus
    [2020-01-15 01:39:09,358]    INFO Whitelist not provided
    [2020-01-15 01:39:09,362]    INFO Copying pre-packaged 10XV2 whitelist to output
    [2020-01-15 01:39:09,479]    INFO Inspecting BUS file tmp/output.s.bus
    [2020-01-15 01:39:27,669]    INFO Correcting BUS records in tmp/output.s.bus to tmp/output.s.c.bus with whitelist output/10xv2_whitelist.txt
    [2020-01-15 01:40:00,727]    INFO Sorting BUS file tmp/output.s.c.bus to output/output.unfiltered.bus
    [2020-01-15 01:41:04,474]    INFO Generating count matrix output/counts_unfiltered/cells_x_genes from BUS file output/output.unfiltered.bus


## Load the count matrices into a notebook

See the getting started tutorial for how to load the count matrices into [ScanPy](https://scanpy.readthedocs.io/en/latest/index.html) for analysis.


```

```
