**kallisto** and **bustools** are wrapped in an easy-to-use program called `kb` which is part of the `kb-python` package ([developer documentation](https://kb-python.readthedocs.io/en/latest/index.html)). `kb-python` can be installed on any machine by typing `pip install kb-python` on the command line; this installs everything needed to process single-cell RNA-seq reads. 

`kb` has three commands:

```
$ kb
usage: kb [-h] [--list] <CMD> ...

kb_python 0.26.4

positional arguments:
  <CMD>
    info      Display package and citation information
    ref       Build a kallisto index and transcript-to-gene mapping
    count     Generate count matrices from a set of single-cell FASTQ files

optional arguments:
  -h, --help  Show this help message and exit
  --list      Display list of supported single-cell technologies
```

With two simple commands, `kb ref` and `kb count` you can process single-cell RNA-seq reads:

```
$ kb ref -i transcriptome.idx -g transcripts_to_genes.txt -f1 cdna.fa dna.primary_assembly.fa.gz gtf.gz
$ kb count -i index.idx -g t2g.txt -x 10xv2 --h5ad -t 2 read_1.fastq.gz read_2.fastq.gz
```

The last command, `kb info`, print out a description of the program:

```
$ kb info
kb_python 0.26.4
    kallisto: 0.46.2
    bustools: 0.40.0

kb is a python package for rapidly pre-processing single-cell RNA-seq data. It
is a wrapper for the methods described in [2].

The goal of the wrapper is to simplify downloading and running of the kallisto
[1] and bustools [2] programs. It was inspired by Sten Linnarsson’s loompy
fromfq command (http://linnarssonlab.org/loompy/kallisto/index.html)

The kb program consists of two parts:

The `kb ref` command builds or downloads a species-specific index for
pseudoalignment of reads. This command must be run prior to `kb count`, and it
runs the `kallisto index` [1].

The `kb count` command runs the kallisto [1] and bustools [2] programs. It can
be used for pre-processing of data from a variety of single-cell RNA-seq
technologies, and for a number of different workflows (e.g. production of gene
count matrices, RNA velocity analyses, etc.). The output can be saved in a
variety of formats including mix and loom. Examples are provided below.

Examples are available at: https://www.kallistobus.tools/tutorials

References
==========
[1] Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal
probabilistic RNA-seq quantification. Nature biotechnology, 34(5), 525.
[2] Melsted, P., Booeshaghi, A. S., Liu, L., Gao, F., Lu, L., Min, K. H., da
Veiga Beltrame, E., Hjorleifsson, K. E., Gehring, J., & Pachter, L. (2021).
Modular and efficient pre-processing of single-cell RNA-seq. Nature
Biotechnology.
```

#### kallisto and bustools

`kb-python` utilizes the `kallisto` and `bustools` programs. Details on the use of these tools can be found below:

The `kallisto` manual is available at: [https://pachterlab.github.io/kallisto/download](https://pachterlab.github.io/kallisto/download)

The `bustools` manual is available at: [https://bustools.github.io/manual](https://bustools.github.io/manual)
