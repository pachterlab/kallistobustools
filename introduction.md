---
layout: page
title: "Introduction"
group: navigation
---

{% include JB/setup %}

We have wrapped kallisto and bustools into an easy-to-use program called `kb-python` which can be installed on any machine by typing `pip install kb-python` on the command line. This installs everything you need to process your single-cell RNA-seq reads with two simple commands. The first command is `kb ref` and the second command `kb count`.

```
$ kb ref -i transcriptome.idx -g transcripts_to_genes.txt -f1 cdna.fa dna.primary_assembly.fa.gz gtf.gz
$ kb count -i index.idx -g t2g.txt -x 10xv2 --h5ad -t 2 read_1.fastq.gz read_2.fastq.gz
```

To play around with a real data set using `kb-python`, check out this [tutorial](https://colab.research.google.com/github/pachterlab/kallistobustools/blob/master/notebooks/kb_species_mixing.ipynb). For an in-depth overview of `kb` checkout the [docs](https://kb-python.readthedocs.io/en/latest/index.html).

## `kb-python` explained
### 1. Building a reference
takes a genome and a genome annotation and builds a transcriptome index to use with `kallisto` under the hood of `kb-python`. `kb ref` allows flexibility in making a transcriptome from a custom genome and genome annotation, making custom references for RNA velocity/single-nuclei RNA-seq, and it allows one to download premade references.

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

### 2. Running a workflow
Once a reference is made, FASTQ reads must be aligned to it. `kb count` is responsible for running a variety of workflows that align reads to a reference for standard count matrix generation to feature barcoding.
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

## `kb-python` is built on top of `kallisto` and `bustools`

### kallisto and bustools manuals

The kallisto &#124; bustools single-cell RNA-seq workflow requires two programs: `kallisto` and `bustools`

The `kallisto` manual is available at: [https://pachterlab.github.io/kallisto/download](https://pachterlab.github.io/kallisto/download)

The `bustools` manual is available at: [https://bustools.github.io/manual](https://bustools.github.io/manual)
