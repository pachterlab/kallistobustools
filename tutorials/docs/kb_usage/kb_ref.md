#### Building a reference

Starting with a genome and a genome annotation a transcriptome index can be built with `kallisto` via `kb ref`. This allows flexibility in building a transcriptomes from genomes and associated genome annotations. In addition `kb ref` can be used to download pre-built indices:

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
