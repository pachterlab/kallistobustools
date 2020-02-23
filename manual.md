---
layout: page
title: "Manual"
---

{% include JB/setup %}

Manual of all single-cell RNA-seq workflows supported by `kb`. This document is organized by workflow and then by the commands that are run, followed by descriptions of each required and optional argument.

### Table of Contents

1. [Standard workflow](#standard)
2. [RNA velocity workflow](#velocity)
3. [Single-nucleus workflow](#nucleus)
4. [KITE feature barcoding workflow](#kite)

-------------------------------------------

### 1. Standard workflow<a name='standard'></a>
### ref
**Options that apply to all commands**

The following options apply to all commands when running `kb ref` standard workflow.

|:---|:---|
| `--tmp TMP`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Set the temporary directory, which defaults to `./tmp`, to `TMP`. |
| `--keep-tmp` | Do not delete the temporary directory once `kb` finishes running. |
| `--verbose` | Output debugging information. |
| `--overwrite` | Overwrite any existing files. |

&nbsp;

#### Download a pre-built mouse or human index
This will download a `tar.gz` archive containing a kallisto index and transcript-to-gene mapping.
```
kb ref -i INDEX -g T2G -d ORGANISM [options]
```
**Required arguments**

|:---|:---|
| `INDEX`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Path to the kallisto index. The extracted kallisto index will be placed here. |
| `T2G` | Path to transcript-to-gene mapping. The extracted mapping will be placed here. |
| `ORGANISM` | Specifies which organism's index to download. Either `human` or `mouse`. |

&nbsp;

#### Build an index with local or remote files
This will build a kallisto index and a transcript-to-gene mapping using local or remote files. Note that Windows only supports local files.
```
kb ref -i INDEX -g T2G -f1 CDNA [options] FASTA GTF
```
**Required arguments**

|:---|:---|
| `INDEX`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Path to the kallisto index to be constructed. |
| `T2G` | Path to transcript-to-gene mapping to be generated. |
| `CDNA` | Path to the cDNA FASTA to be generated. |
| `FASTA` | Path or URL to the input genomic FASTA file. |
| `GTF` | Path or URL to the input GTF file. |

&nbsp;

**Additional options**

|:---|:---|
| `-n N`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Split the kallisto index into `N` files. Specifying this option will internally split the cDNA FASTA into `N` approximately-equal parts, and each of these FASTAs are indexed separately. |
| `-k K` | Use this option to override the k-mer length of the index to `K`. The default k-mer length is 31. |

### count
**Options that apply to all commands**

The following options apply to all commands when running `kb count` standard workflow.

|:---|:---|
| `-t THREADS`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Number of threads to use (default: `8`). |
| `-m MEMORY` | Maximum memory used (default: `4G`). This option may not be honored when using the `--report` flag. |
| `-w WHITELIST` | Path to plaintext file of whitelisted barcodes to correct to. If not provided and bustools supports the technology, a pre-packaged whitelist is used. If not, the `bustools whitelist` command is used. Use `kb --list` to view whitelists. |
| `--tmp TMP` | Set the temporary directory, which defaults to `./tmp`, to `TMP`. |
| `--keep-tmp` | Do not delete the temporary directory once `kb` finishes running. |
| `--verbose` | Output debugging information. |
| `--overwrite` | Overwrite any existing files. |
| `--dry-run` | Perform a dry run. This prints out the commands that would be run, such that the output of `kb` can be directly used as a bash (MacOS and Linux) or batch (Windows) script. This option does not support conversion of sparse matrices to other formats and report generation. |
| `--report` | Generate a Jupyter notebook and corresponding HTML report containing run statistics and basic plots. Using this option may cause `kb` to use more memory than specified with the `-m` option. It may also cause it to crash due to memory. |

&nbsp;

#### Generate an unfiltered (+ filtered) sparse gene count matrix
This command generates an unfiltered (+ filtered) gene count matrix in sparse matrix format.
```
kb count -i INDEX -g T2G -x TECHNOLOGY [options] FASTQS
```

|:---|:---|
| `INDEX`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Path to the kallisto index generated with `kb ref`. |
| `T2G` | Path to transcript-to-gene mapping generated with `kb ref`. |
| `TECHNOLOGY` | Sequencing/barcoding technology used. By default, `kb` supports some of the popular sequencing methods, such as 10X and inDrops. The full list can be viewed with `kb --list`. It is also possible to specify a custom technology in the format supported by kallisto. See [here](https://pachterlab.github.io/kallisto/manual#bus) for specifications. |
| `FASTQS` | Space-delimited list of input FASTQs. |

&nbsp;

**Additional options**

|:---|:---|
| `--filter FILTER` | Produce an additional, filtered, count matrix (default: `bustools`). Currently, only `bustools` is supported, which runs `bustools whitelist` to detect and filter empty droplets. |


#### Generate an unfiltered (+ filtered) sparse TCC matrix
This command generates an unfiltered (+ filtered) TCC matrix in sparse matrix format.
```
kb count -i INDEX -g T2G -x TECHNOLOGY --tcc [options] FASTQS
```

|:---|:---|
| `INDEX`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Path to the kallisto index generated with `kb ref`. |
| `T2G` | Path to transcript-to-gene mapping generated with `kb ref`. |
| `TECHNOLOGY` | Sequencing/barcoding technology used. By default, `kb` supports some of the popular sequencing methods, such as 10X and inDrops. The full list can be viewed with `kb --list`. It is also possible to specify a custom technology in the format supported by kallisto. See [here](https://pachterlab.github.io/kallisto/manual#bus) for specifications. |
| `FASTQS` | Space-delimited list of input FASTQs. |

&nbsp;

**Additional options**

|:---|:---|
| `--filter FILTER` | Produce an additional, filtered, count matrix (default: `bustools`). Currently, only `bustools` is supported, which runs `bustools whitelist` to detect and filter empty droplets. |

&nbsp;
<!-- #### Generate an unfiltered (+ filtered) sparse isoform count matrix -->

#### Convert sparse matrices to other formats
The following options convert the final unfiltered (+ filtered) sparse matrices to other popular formats.

|:---|:---|
| `--h5ad`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Convert the final matrices to H5AD format. The count matrices are loaded as an [Anndata](https://icb-anndata.readthedocs-hosted.com/en/stable/anndata.AnnData.html) object, which is then saved as `adata.h5ad`. |
| `--loom` | Convert the final matrices to Loom format. The count matrices are loaded as an [Anndata](https://icb-anndata.readthedocs-hosted.com/en/stable/anndata.AnnData.html) object, which is then saved as `adata.loom`. |
| `--cellranger` | Convert count matrices to cellranger-compatible format. This means the matrix is transposed so that the cells (barcodes) are the columns and observations (genes or transcripts) are the rows. The files are written to a separate `cellranger` subdirectory. |


### 2. RNA velocity workflow<a name='velocity'></a>
### ref

### count

### 3. Single-nucleus workflow<a name='nucleus'></a>
### ref

### count

### 4. KITE feature barcoding workflow<a name='kite'></a>
### ref

### count
