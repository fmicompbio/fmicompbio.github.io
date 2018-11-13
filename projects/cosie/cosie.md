---
layout: slim_banner
title: COSIE - COrrected Splicing Indices for Exon arrays
---

## Description
Affymetrix exon arrays were designed as a tool for monitoring the
relative expression levels of hundreds of thousands of known and predicted
exons with a view to detecting alternative splicing events. In the article
referenced below, we have characterised a systematic bias of the exon array
platform that leads to an overestimation of alternative splicing events
in genes that are differentially expressed.

`COSIE` is an [R](http://www.r-project.org) function that for a given set of
exon arrays corrects for the observed bias and improves the detection of
alternative splicing. It adjusts splicing indices for exons, especially for
those that belong to differentially expressed genes. For this adjustment,
`COSIE` uses parameters that are specific for each probeset (see download
section below) which were trained from a large number of published exon arrays.
The downside of this approach is that such parameters cannot be estimated
for all probesets on the microarray. Based on our training set, `COSIE`
corrects 95.1% of the probesets. Separate parameter files are provided for both
the *full* and *core* sets, including all probesets that are linked to transcripts.
We recommend the use of the *core* set that was also used in the cited study
below. The *full* set is not as well characterized.

`SI_*` returned by COSIE contains the net probeset expressions after factoring
out gene expression and exon array bias (pre splicing index). To obtain the
differential probeset inclusion rates (final splicing indices), two columns
of `SI_*` simply need to be subtracted from one another. `tclevel_*` contains
the transcript levels used internally by COSIE (in case someone needs them)
but they are not required by the user when detecting alternatively included probesets.

A typical exon array data analysis workflow may look as
follows:

1. normalize and condense exon arrays  
2. adjust pre-splicing indices (`COSIE`, see also usage example below)  
3. compare pre-splicing indices in different samples to identify alternative exons  

## Download `COSIE` Software
- [COSIE software:](cosie.R) implementation of COSIE in [R](http://www.r-project.org)  
- [usage example:](cosie_example.R) how to calculate corrected splicing indices in R using `COSIE`  

## Download COSIE Parameter Files
- Human Exon Array  
  * [Core probesets v1](cosie_human_Core_v1.txt.gz) (based on HuEx-1_0-st-v2.r2.dt1.hg18, HuEx-1_0-st-v2.na24.hg18.transcript)  
  * [Full probesets v1](cosie_human_Full_v1.txt.gz) (based on HuEx-1_0-st-v2.r2.dt1.hg18, HuEx-1_0-st-v2.na24.hg18.transcript)  
- Mouse Exon Array  
  * [Core probesets v1](cosie_mouse_Core_v1.txt.gz) (based on MoEx-1_0-st-v1.r2.dt1.mm8, MoEx-1_0-st-v1.na24.mm8)  
  * [Full probesets v1](cosie_mouse_Full_v1.txt.gz) (based on MoEx-1_0-st-v1.r2.dt1.mm8, MoEx-1_0-st-v1.na24.mm8)  

## Citation
**Overestimation of alternative splicing caused by variable probe characteristics in exon arrays.**  
Dimos Gaidatzis; Kirsten Jacobeit; Edward J. Oakeley; Michael B. Stadler  
*Nucleic Acids Research* 2009; doi: 10.1093/nar/gkp508  
[Abstract](https://www.ncbi.nlm.nih.gov/pubmed/19528075) |
[Full text (free)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2760813/) |
[PDF (free)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2760813/pdf/gkp508.pdf) |
[PDF (local)](gkp508v1_screen.pdf)  
