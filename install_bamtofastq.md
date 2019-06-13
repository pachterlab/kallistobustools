---
layout: page
title: "Install bamtofastq utility"
group: navigation
---

{% include JB/setup %}

Since most authors upload only BAM files to the SRA, we need a utility to convert BAM files back into FASTQ files for processing with kallisto and bustools.

### Download bamtofastq
Navigate to the [10x website](https://support.10xgenomics.com/docs/bamtofastq). Click on the link `Download bamtofastq 1.1.2` and download the file to your computer. `bamtofastq` is a single executable that can be run directly and does not need to be compiled or installed.

If you want to run the bamtofastq utility from anywhere on your computer (in your terminal) Google "adding binaries to your path [insert Windows/Mac/Linux here]". Otherwise you can run the binary by specifiying where it is i.e.
```
$ /home/sina/tools/bamtofastq [options here]
```

### Using bamtofastq
Given a bam file `my_file.bam` we can convert it into FASTQs using this utility. Open up your terminal and navigate to your BAM file. The run the following command:

```
$ bamtofastq --reads-per-fastq=500000000 my_file.bam ./fastqs
```
And your FASTQs will be in the folder `./fastqs`.
