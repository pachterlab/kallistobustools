**kallisto** and **bustools** are wrapped in an easy-to-use program called `kb` which is part of the `kb-python` package, and that can be installed on any machine by typing `pip install kb-python` on the command line. This installs everything needed to process single-cell RNA-seq reads with two simple commands. The first is `kb ref` and the second `kb count`:

```
$ kb ref -i transcriptome.idx -g transcripts_to_genes.txt -f1 cdna.fa dna.primary_assembly.fa.gz gtf.gz
$ kb count -i index.idx -g t2g.txt -x 10xv2 --h5ad -t 2 read_1.fastq.gz read_2.fastq.gz
```

For an in-depth overview of `kb` see the [docs](https://kb-python.readthedocs.io/en/latest/index.html).

#### kallisto and bustools

`kb-python` utilizes the `kallisto` and `bustools` programs. Details on the use of these tools can be found below:

The `kallisto` manual is available at: [https://pachterlab.github.io/kallisto/download](https://pachterlab.github.io/kallisto/download)

The `bustools` manual is available at: [https://bustools.github.io/manual](https://bustools.github.io/manual)
