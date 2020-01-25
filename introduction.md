---
layout: page
title: "Introduction"
group: navigation
---

{% include JB/setup %}

__kallisto__ and __bustools__ are wrapped in an easy-to-use program called `kb` which is part of the `kb-python` package, and that can be installed on any machine by typing `pip install kb-python` on the command line. This installs everything needed to process single-cell RNA-seq reads with two simple commands. The first is `kb ref` and the second `kb count`:

```
$ kb ref -i transcriptome.idx -g transcripts_to_genes.txt -f1 cdna.fa dna.primary_assembly.fa.gz gtf.gz
$ kb count -i index.idx -g t2g.txt -x 10xv2 --h5ad -t 2 read_1.fastq.gz read_2.fastq.gz
```

To learn how to explore a dataset using `kb-python` begin with the [getting started tutorial](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_species_mixing.ipynb). For an in-depth overview of `kb` see the [docs](https://kb-python.readthedocs.io/en/latest/index.html).

## Details
### 1. Building a reference
Starting with a genome and a genome annotation a transcriptome index  can be built with `kallisto` via `kb ref`. This allows flexibility in building a transcriptomes from genomes and associated genome annotations. In addition `kb ref` can be used to download pre-built indices:

```
$ kb ref
usage: kb ref [-h] [--tmp TMP] [--keep-tmp] [--verbose] -i INDEX -g T2G -f1
              FASTA [-f2 FASTA] [-c1 T2C] [-c2 T2C] [-n N]
              [-d {human,mouse,linnarsson}] [-k K]
              [--workflow {standard,lamanno,nucleus,kite}] [--lamanno]
              [--overwrite]
              fasta gtf [feature]

Build a kallisto index and transcript-to-gene mapping

positional arguments:
  fasta                 Genomic FASTA file(s), comma-delimited
  gtf                   Reference GTF file(s), comma-delimited
  feature               [`kite` workflow only] Path to TSV containing barcodes
                        and feature names.

optional arguments:
  -h, --help            Show this help message and exit
  --tmp TMP             Override default temporary directory
  --keep-tmp            Do not delete the tmp directory
  --verbose             Print debugging information
  -n N                  Number of files to split the index into. If this
                        option is specified, the FASTA that is normally used
                        to create an index is split into `N` approximately-
                        equal parts. Each of these FASTAs are indexed
                        separately.
  -d {human,mouse,linnarsson}
                        Download a pre-built kallisto index (along with all
                        necessary files) instead of building it locally
  -k K                  Use this option to override the k-mer length of the
                        index. Usually, the k-mer length automatically
                        calculated by `kb` provides the best results.
  --workflow {standard,lamanno,nucleus,kite}
                        Type of workflow to prepare files for. Use `lamanno`
                        for RNA velocity based on La Manno et al. 2018 logic.
                        Use `nucleus` for RNA velocity on single-nucleus RNA-
                        seq reads. Use `kite` for feature barcoding. (default:
                        standard)
  --lamanno             Deprecated. Use `--workflow lamanno` instead.
  --overwrite           Overwrite existing kallisto index

required arguments:
  -i INDEX              Path to the kallisto index to be constructed. If `-n`
                        is also specified, this is the prefix for the n
                        indices to construct.
  -g T2G                Path to transcript-to-gene mapping to be generated
  -f1 FASTA             [Optional with -d] Path to the cDNA FASTA (lamanno,
                        nucleus) or mismatch FASTA (kite) to be generated

required arguments for `lamanno` and `nucleus` workflows:
  -f2 FASTA             Path to the intron FASTA to be generated
  -c1 T2C               Path to generate cDNA transcripts-to-capture
  -c2 T2C               Path to generate intron transcripts-to-capture
```

### 2. Quantifying a dataset
Once an index has been generated or downloaded, `kb count` uses `kallisto` to pseudoalign reads and `bustools` to quantify the data: 
```
usage: kb count [-h] [--tmp TMP] [--keep-tmp] [--verbose] -i INDEX -g T2G -x
                TECHNOLOGY [-o OUT] [-w WHITELIST] [-t THREADS] [-m MEMORY]
                [--workflow {standard,lamanno,nucleus,kite,kite:10xFB}]
                [--mm | --tcc] [--filter [{bustools}]] [-c1 T2C] [-c2 T2C]
                [--overwrite] [--dry-run] [--lamanno | --nucleus]
                [--loom | --h5ad]
                fastqs [fastqs ...]

Generate count matrices from a set of single-cell FASTQ files. Run `kb --list`
to view single-cell technology information.

positional arguments:
  fastqs                FASTQ files

optional arguments:
  -h, --help            Show this help message and exit
  --tmp TMP             Override default temporary directory
  --keep-tmp            Do not delete the tmp directory
  --verbose             Print debugging information
  -o OUT                Path to output directory (default: current directory)
  -w WHITELIST          Path to file of whitelisted barcodes to correct to. If
                        not provided and bustools supports the technology, a
                        pre-packaged whitelist is used. If not, the bustools
                        whitelist command is used. (`kb --list` to view
                        whitelists)
  -t THREADS            Number of threads to use (default: 8)
  -m MEMORY             Maximum memory used (default: 4G)
  --workflow {standard,lamanno,nucleus,kite,kite:10xFB}
                        Type of workflow. Use `lamanno` for RNA velocity based
                        on La Manno et al. 2018 logic. Use `nucleus` for RNA
                        velocity on single-nucleus RNA-seq reads. Use `kite`
                        for feature barcoding. Use `kite:10xFB` for 10x
                        Genomics Feature Barcoding technology. (default:
                        standard)
  --mm                  Include reads that pseudoalign to multiple genes.
  --tcc                 Generate a TCC matrix instead of a gene count matrix.
  --filter [{bustools}]
                        Produce a filtered gene count matrix (default:
                        bustools)
  --overwrite           Overwrite existing output.bus file
  --dry-run             Dry run
  --lamanno             Deprecated. Use `--workflow lamanno` instead.
  --nucleus             Deprecated. Use `--workflow nucleus` instead.
  --loom                Generate loom file from count matrix
  --h5ad                Generate h5ad file from count matrix

required arguments:
  -i INDEX              Path to kallisto index/indices, comma-delimited
  -g T2G                Path to transcript-to-gene mapping
  -x TECHNOLOGY         Single-cell technology used (`kb --list` to view)

required arguments for `lamanno` and `nucleus` workflows:
  -c1 T2C               Path to cDNA transcripts-to-capture
  -c2 T2C               Path to intron transcripts-to-captured
```

### 3. kallisto and bustools
`kb-python` utilizes the `kallisto` and `bustools` programs. Details on the use of these tools can be found below:

The `kallisto` manual is available at: [https://pachterlab.github.io/kallisto/download](https://pachterlab.github.io/kallisto/download)

The `bustools` manual is available at: [https://bustools.github.io/manual](https://bustools.github.io/manual)
