---
layout: slim_banner
title: MethylSeekR supplementary data
---

## Description
`MethylSeekR` is a Bionductor/R package for the identification of unmethylated
(UMRs) and low-methylated regions (LMRs) from whole-genome bisulfite-sequencing
data. It can be downloaded from Bioconductor (http://www.bioconductor.org).

MethylSeekR segmentation results for 18 published human methylomes can be downloaded here:
[Identified segments (27.2MB)](IdentifiedSegments.tar.gz)  

Examples profiles showing the segmentation results for a representative subset
of the data can be downloaded here:  
[Segmentation Examples (5.1MB)](SegmentationExamples.tar.gz)  

The segmentation examples are multi-page pdf files for a representative set of
methylomes. Each page shows the segmentation for a randomly chosen region.
Each region displayed is broken up into 3 pairs of panels, where in each pair
the same region is shown twice, once with raw methylation levels (top) and once
with methylation levels smoothed over 3 consecutive CpGs (bottom). In both cases
only CpGs with a coverage of at least 5 reads are shown. The raw data best illustrates the
disordered nature of methylation levels in partially methylated domains(PMDs), whereas the smoothed
methylation levels more clearly show UMRs and LMRs. In all figures, UMRs are
shown as blue squares (placed at the middle of the identified segment), LMRs
as red triangles (placed at the middle of the identified segment) and PMDs as
green bars (extending over the entire PMD). The cut-off on methylation to
determine UMRs and LMRs is shown as a grey dashed line.
