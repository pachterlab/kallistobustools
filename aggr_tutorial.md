# Aggregating multiple  count matrices

This tutorial describes how to aggregate multiple count matrices by concatenating them into a single AnnData object with batch labels for different samples. A notebook showing the entire workflow (including running kallisto and bsutools) is available [here](here).

This is similar to the Cell Ranger `aggr` function, however no normalization is performed. `cellranger aggr` is described at https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate

For this tutorial we use dataset E-MTAB-6108. We provide the count matrices as an 80mb zip file at https://github.com/BUStools/getting_started/releases/download/aggr/E-MTAB-6108_sample1_sample2_genecounts.zip

If you download the zip file is has the following structure: 
```
E-MTAB-6108_sample1_sample2_genecounts.zip
├── sample1
│   ├── genecounts
│   │   ├── genes.barcodes.txt
│   │   ├── genes.genes.txt
│   │   └── genes.mtx
│   ├── matrix.ec
│   ├── run_info.json
│   └── transcripts.txt
└── sample2
    ├── genecounts
    │   ├── genes.barcodes.txt
    │   ├── genes.genes.txt
    │   └── genes.mtx
    ├── matrix.ec
    ├── run_info.json
    └── transcripts.txt
```

The raw data for E-MTAB-6108 is available at:
https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6108/
