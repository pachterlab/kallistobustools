# kallisto | bustools for single-cell RNA-seq pre-processing

---

**kallisto &#124; bustools** is a workflow for pre-processing single-cell RNA-seq data. Pre-processing single-cell RNA-seq involves: (1) association of reads with their cells of origin, (2) collapsing of reads according to unique molecular identifiers (UMIs), and (3) generation of gene or feature counts from the reads to generate a _cell x gene_ matrix.

With **kallisto &#124; bustools** you can

- Generate a _cell x gene_ or _cell x transcript equivalence class_ count matrix
- Perform RNA velocity and single-nuclei RNA-seq analsis
- Quantify data from numerous technologies such as 10x, inDrops, and Dropseq.
- Customize workflows for new technologies and protocols.
- Process feature barcoding data such as CITE-seq, REAP-seq, MULTI-seq, Clicktags, and Perturb-seq.
- Obtain QC reports from single-cell RNA-seq data

The **kallisto &#124; bustools** workflow is described in:

Páll Melsted, A. Sina Booeshaghi, Fan Gao, Eduardo Beltrame, Lambda Lu, Kristján Eldjárn Hjorleifsson, Jase Gehring and Lior Pachter, <a href="https://www.biorxiv.org/content/10.1101/673285v2" class="external-link" target="_blank">Modular and efficient pre-processing of single-cell RNA-seq</a>, bioRxiv, 2019.
