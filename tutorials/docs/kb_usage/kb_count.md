#### Quantifying a dataset

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
