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

## 1. Standard workflow<a name='standard'></a>
### ref
#### Options that apply to all commands

|:--------------|:------------------------------------------------------------------|
| `--tmp TMP`   | Set the temporary directory, which defaults to `./tmp`, to `TMP`. |
| `--keep-tmp`  | Do not delete the temporary directory once `kb` finishes running. |
| `--verbose`   | Output debugging information.                                     |
| `--overwrite` | Overwrite any existing files.                                     |

#### Download a pre-built mouse or human index
```
kb ref -i INDEX -g T2G -d ORGANISM [options]
```
* `INDEX`
* `T2G`
* `ORGANISM`

#### Build an index with local files or files streamed from the internet
```
kb ref -i INDEX -g T2G -f1 CDNA [options] FASTA GTF
```
* `INDEX`
* `T2G`
* `CDNA`
* `FASTA`
* `GTF`
* `[options]`
  * `-n N`
  * `-k K`

### count

## 2. RNA velocity workflow<a name='velocity'></a>
### ref

### count

## 3. Single-nucleus workflow<a name='nucleus'></a>
### ref

### count

## 4. KITE feature barcoding workflow<a name='kite'></a>
### ref

### count
