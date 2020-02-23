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
|:--------------|:------------------------------------------------------------------|
| `--tmp TMP`   | Set the temporary directory, which defaults to `./tmp`, to `TMP`. |
| `--keep-tmp`  | Do not delete the temporary directory once `kb` finishes running. |
| `--verbose`   | Output debugging information.                                     |
| `--overwrite` | Overwrite any existing files.                                     |

#### Download a pre-built mouse or human index
This will download a `tar.gz` archive containing a kallisto index and transcript-to-gene mapping.
```
kb ref -i INDEX -g T2G -d ORGANISM [options]
```
**Required arguments**

|:-----------|:-------------------------------------------------------------------------------|
| `INDEX`    | Path to the kallisto index. The extracted kallisto index will be placed here.  |
| `T2G`      | Path to transcript-to-gene mapping. The extracted mapping will be placed here. |
| `ORGANISM` | Specifies which organism's index to download. Either `human` or `mouse`.       |

#### Build an index with local or remote files
This will build a kallisto index and a transcript-to-gene mapping using local or remote files. Note that Windows only supports local files.
```
kb ref -i INDEX -g T2G -f1 CDNA [options] FASTA GTF
```
**Required arguments**

|:-----------|:-------------------------------------------------------------------------------|
| `INDEX`     | Path to the kallisto index to be constructed.  |
| `T2G`       | Path to transcript-to-gene mapping to be generated. |
| `CDNA`      | Path to the cDNA FASTA to be generated.  |
| `FASTA`     | Path or URL to the input genomic FASTA file.  |
| `GTF`       | Path or URL to the input GTF file.  |

&nbsp;

**Options**

|:-------|:-------------------------------------------------------------------------------|
| `-n N` | Split the kallisto index into `N` files. Specifying this option will internally split the cDNA FASTA into `N` approximately-equal parts, and each of these FASTAs are indexed separately. |
| `-k K` | Use this option to override the k-mer length of the index to `K`. The default k-mer length is 31. |

### count

### 2. RNA velocity workflow<a name='velocity'></a>
### ref

### count

### 3. Single-nucleus workflow<a name='nucleus'></a>
### ref

### count

### 4. KITE feature barcoding workflow<a name='kite'></a>
### ref

### count
