---
layout: slim_banner
title: GLD-1 Binding Motif (GBM) finder
---

## Description
In the article listed below, we characterized the detailed binding specificity
of the RNA binding protein GLD-1 which is specifically expressed in the *C.elegans*
germline. We incorporated these binding determinants into a small R script which
can be used to scan any given sequence for GLD-1 binding motifs (GBM's). To make
use of the software, install R (http://www.r-project.org) as well as the R package
`Biostrings` (http://www.bioconductor.org). Download and extract the **GBM finder**
(zip file listed below) into an empty directory and start R. Execute the example
script `example.r` that scans the test sequences provided in `seqs.fa` by typing
`source('example.r')`.

The function which performs the scan (scanForGBMs) is called as follows:

```
hits <- scanForGBMs("seqs.fa","GBM_scores.tab")
```

where `seqs.fa` is the fasta file with the sequences that should be scanned and
`GBM_scores.tab` is the table which contains the search parameters (provided in
the zip file). The result of the scan is returned in the table hits.

## Download GBM finder
- [GBM finder](GBM_finder.zip)  

## Citation
**A quantitative RNA code for mRNA target selection by the germline fate determinant GLD-1.**  
Jane E Wright, Dimos Gaidatzis, Mathias Senften, Brian M Farley, Eric Westhof, Sean P Ryder and Rafal Ciosk  
*EMBO Journal* 2010; doi:10.1038/emboj.2010.334  

[PubMed](https://www.ncbi.nlm.nih.gov/pubmed/21169991) |
[free HTML](http://emboj.embopress.org/content/30/3/533.long)  
