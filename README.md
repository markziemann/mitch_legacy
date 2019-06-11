# Mitch
Mitch is a tool for multi-dimensional enrichment analysis. At it's heart, it uses a rank-MANOVA based statistical approach to detect sets of genes that exhibit enrichment in the multidimensional space as compared to the background. Mitch is useful for pathway analysis of profiling studies with two to or more contrasts, or in studies with multiple omics profiling, for example proteomic, transcriptomic, epigenomic analysis of the same samples.

## Installation
```
install.packages("devtools")
library("devtools")
devtools::install_github("markziemann/Mitch")
library("mitch")
```

## Workflow overview
### Importing gene sets
Mitch has a function to import GMT files to R lists, which was adapted from work by [Yu et al, 2012](https://dx.doi.org/10.1089%2Fomi.2011.0118) in the clusterProfiler package.
### Importing profiling data
Mitch accepts pre-ranked data supplied by the user, but also has functions for importing tables generated by Limma, edgeR, DESeq2, ABSSeq and Sleuth. By default, only the genes that are detected in all contrasts are included, but this behaviour can be modified. The below example imports two edgeR tables called "dge1" and "dge2".
```
x<-list("dge1"=dge1,"dge2"=dge2)
y<-mitch_import(x,"edger")
```
By default, differential gene activiy is scored using the directional nominal p-value.

S=-log10(p-value) * sgn(logFC)

If this is not desired, then users can perform their own custom scoring procedure.



