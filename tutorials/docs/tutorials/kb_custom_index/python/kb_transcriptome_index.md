<a href="https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/tutorials/docs/tutorials/kb_custom_index/python/kb_transcriptome_index.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# Constructing a transcriptome index with `kb`

This tutorial provides instructions for how to generate a transcriptome index to use with **kallisto | bustools** using `kb`.

## Download reference files

Download the genomic (DNA) FASTA and GTF annotations for your desired organism from the database of your choice. This tutorial uses mouse reference files downloaded from [Ensembl](https://uswest.ensembl.org/info/data/ftp/index.html).


```
%%time
!wget -q ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
!wget -q ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz
```

    CPU times: user 173 ms, sys: 37.3 ms, total: 211 ms
    Wall time: 28.8 s


## Install `kb`


```
!pip install --quiet kb-python
```

    [K     |████████████████████████████████| 59.1MB 71kB/s 
    [K     |████████████████████████████████| 51kB 4.9MB/s 
    [K     |████████████████████████████████| 122kB 46.5MB/s 
    [K     |████████████████████████████████| 13.2MB 42.9MB/s 
    [K     |████████████████████████████████| 10.3MB 30.2MB/s 
    [K     |████████████████████████████████| 112kB 61.3MB/s 
    [K     |████████████████████████████████| 81kB 7.9MB/s 
    [K     |████████████████████████████████| 1.2MB 42.2MB/s 
    [K     |████████████████████████████████| 51kB 3.8MB/s 
    [K     |████████████████████████████████| 71kB 5.2MB/s 
    [?25h  Building wheel for loompy (setup.py) ... [?25l[?25hdone
      Building wheel for numpy-groupies (setup.py) ... [?25l[?25hdone
      Building wheel for umap-learn (setup.py) ... [?25l[?25hdone
      Building wheel for sinfo (setup.py) ... [?25l[?25hdone
      Building wheel for pynndescent (setup.py) ... [?25l[?25hdone


## Build the index

`kb` automatically splits the genome into a cDNA FASTA file and uses that to build a kallisto index.


```
%%time
!kb ref -i transcriptome.idx -g transcripts_to_genes.txt -f1 cdna.fa \
Mus_musculus.GRCm38.dna.primary_assembly.fa.gz \
Mus_musculus.GRCm38.98.gtf.gz
```

    [2021-03-31 19:42:16,748]    INFO Preparing Mus_musculus.GRCm38.dna.primary_assembly.fa.gz, Mus_musculus.GRCm38.98.gtf.gz
    [2021-03-31 19:42:16,748]    INFO Decompressing Mus_musculus.GRCm38.98.gtf.gz to tmp
    [2021-03-31 19:42:20,207]    INFO Creating transcript-to-gene mapping at /content/tmp/tmp8orgc74k
    [2021-03-31 19:43:01,135]    INFO Decompressing Mus_musculus.GRCm38.dna.primary_assembly.fa.gz to tmp
    [2021-03-31 19:43:25,710]    INFO Sorting tmp/Mus_musculus.GRCm38.dna.primary_assembly.fa to /content/tmp/tmpees_4ry7
    [2021-03-31 19:50:50,364]    INFO Sorting tmp/Mus_musculus.GRCm38.98.gtf to /content/tmp/tmpssw7nu7e
    [2021-03-31 19:51:51,290]    INFO Splitting genome tmp/Mus_musculus.GRCm38.dna.primary_assembly.fa into cDNA at /content/tmp/tmpg65jndu0
    [2021-03-31 19:51:51,290] WARNING The following chromosomes were found in the FASTA but does not have any "transcript" features in the GTF: JH584302.1, GL456387.1, GL456396.1, GL456367.1, GL456366.1, GL456394.1, GL456383.1, GL456382.1, GL456393.1, GL456368.1, GL456379.1, GL456390.1, GL456378.1, GL456360.1, GL456389.1, JH584301.1, JH584300.1, GL456392.1, GL456370.1, GL456359.1, GL456213.1. No sequences will be generated for these chromosomes.
    [2021-03-31 19:53:04,043]    INFO Wrote 142446 cDNA transcripts
    [2021-03-31 19:53:04,047]    INFO Concatenating 1 transcript-to-gene mappings to transcripts_to_genes.txt
    [2021-03-31 19:53:04,264]    INFO Concatenating 1 cDNAs to cdna.fa
    [2021-03-31 19:53:05,204]    INFO Indexing cdna.fa to transcriptome.idx
    CPU times: user 7.92 s, sys: 1.04 s, total: 8.96 s
    Wall time: 21min 57s



```

```
