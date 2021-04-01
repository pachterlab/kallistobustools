<a href="https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/tutorials/docs/tutorials/kb_multiple_fastqs/python/kb_multiple_files.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# Processing Multiple Lanes at Once

This tutorial provides instructions for how to pre-process the [mouse T cells SRR8206317](https://www.ncbi.nlm.nih.gov/sra/?term=SRR8206317) dataset from [Miller & Sen et al., 2019](https://doi.org/10.1038/s41590-019-0312-6) using the **kallisto | bustools** workflow.

## Download the data


```
!wget -q ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR820/SRR8206317/d10_Tet_possorted_genome_bam.bam
```

## Install `kb` and `bamtofastq`

We will be using `bamtofastq` to generate the original FASTQ files from the BAM files provided by the authors.


```
!pip install --quiet kb-python
```

    [K     |████████████████████████████████| 59.1MB 131kB/s 
    [K     |████████████████████████████████| 51kB 6.3MB/s 
    [K     |████████████████████████████████| 13.2MB 37.1MB/s 
    [K     |████████████████████████████████| 10.3MB 19.2MB/s 
    [K     |████████████████████████████████| 122kB 45.0MB/s 
    [K     |████████████████████████████████| 112kB 44.8MB/s 
    [K     |████████████████████████████████| 81kB 8.0MB/s 
    [K     |████████████████████████████████| 51kB 6.0MB/s 
    [K     |████████████████████████████████| 71kB 7.7MB/s 
    [K     |████████████████████████████████| 1.2MB 39.6MB/s 
    [?25h  Building wheel for loompy (setup.py) ... [?25l[?25hdone
      Building wheel for numpy-groupies (setup.py) ... [?25l[?25hdone
      Building wheel for sinfo (setup.py) ... [?25l[?25hdone
      Building wheel for umap-learn (setup.py) ... [?25l[?25hdone
      Building wheel for pynndescent (setup.py) ... [?25l[?25hdone



```
!wget -q http://cf.10xgenomics.com/misc/bamtofastq-1.2.0
!chmod +x bamtofastq-1.2.0
```


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

    [2021-03-31 23:49:33,545]    INFO Downloading files for mouse from https://caltech.box.com/shared/static/vcaz6cujop0xuapdmz0pplp3aoqc41si.gz to tmp/vcaz6cujop0xuapdmz0pplp3aoqc41si.gz
    100% 1.89G/1.89G [01:26<00:00, 23.4MB/s]
    [2021-03-31 23:51:01,788]    INFO Extracting files from tmp/vcaz6cujop0xuapdmz0pplp3aoqc41si.gz
    CPU times: user 1.29 s, sys: 298 ms, total: 1.58 s
    Wall time: 2min 1s


## Generate the FASTQs from the BAM file

Use the `bamtofastq` utility to generate the FASTQs.


```
%%time
!./bamtofastq-1.2.0 --reads-per-fastq=500000000 d10_Tet_possorted_genome_bam.bam ./fastqs
```

    bamtofastq v1.2.0
    Args { arg_bam: "d10_Tet_possorted_genome_bam.bam", arg_output_path: "./fastqs", flag_nthreads: 4, flag_locus: None, flag_bx_list: None, flag_reads_per_fastq: 500000000, flag_gemcode: false, flag_lr20: false, flag_cr11: false }
    Writing finished.  Observed 85992089 read pairs. Wrote 85992089 read pairs
    CPU times: user 4.46 s, sys: 491 ms, total: 4.95 s
    Wall time: 12min 12s


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

    [2021-04-01 00:03:48,050]    INFO Using index index.idx to generate BUS file to output from
    [2021-04-01 00:03:48,051]    INFO         fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L001_R1_001.fastq.gz
    [2021-04-01 00:03:48,051]    INFO         fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L001_R2_001.fastq.gz
    [2021-04-01 00:03:48,051]    INFO         fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L002_R1_001.fastq.gz
    [2021-04-01 00:03:48,051]    INFO         fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L002_R2_001.fastq.gz
    [2021-04-01 00:03:48,051]    INFO         fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L003_R1_001.fastq.gz
    [2021-04-01 00:03:48,051]    INFO         fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L003_R2_001.fastq.gz
    [2021-04-01 00:03:48,051]    INFO         fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L004_R1_001.fastq.gz
    [2021-04-01 00:03:48,051]    INFO         fastqs/MAR1_POOL_2_10_tet_SI-GA-C8_MissingLibrary_1_HJNKJBGX5/bamtofastq_S1_L004_R2_001.fastq.gz
    [2021-04-01 00:10:48,408]    INFO Sorting BUS file output/output.bus to output/tmp/output.s.bus
    [2021-04-01 00:12:03,976]    INFO Whitelist not provided
    [2021-04-01 00:12:03,976]    INFO Copying pre-packaged 10XV2 whitelist to output
    [2021-04-01 00:12:04,105]    INFO Inspecting BUS file output/tmp/output.s.bus
    [2021-04-01 00:12:22,834]    INFO Correcting BUS records in output/tmp/output.s.bus to output/tmp/output.s.c.bus with whitelist output/10xv2_whitelist.txt
    [2021-04-01 00:12:39,443]    INFO Sorting BUS file output/tmp/output.s.c.bus to output/output.unfiltered.bus
    [2021-04-01 00:13:57,946]    INFO Generating count matrix output/counts_unfiltered/cells_x_genes from BUS file output/output.unfiltered.bus


## Load the count matrices into a notebook

See the getting started tutorial for how to load the count matrices into [ScanPy](https://scanpy.readthedocs.io/en/latest/index.html) for analysis.


```

```
