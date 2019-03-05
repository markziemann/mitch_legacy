library("biomaRt")
library("plyr")

source("mitch.R")

#get gene identifier information from ensembl
ensembl=useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
gt<-getBM(attributes = c('ensembl_gene_id', 'external_gene_name') , mart=ensembl)

#read in the seurat files
file_list<-list.files(".",pattern="WCvWD.tsv")
x = lapply(file_list, read.delim)
names(x)<-file_list

# run the mitch import function for seurat
y<-mitch_import(x,DEtype="seurat",geneTable=gt)

# import the reactome gene sets for mouse
genesets<-gmt_import("../reactome.v5.2.symbols_mouse.gmt")

# remove genes with >1 NA value
#x5<-x4[which(apply(x4, 1, function(x) sum(is.na(x)))<2),]

# run the enrichment calculation
res<-mitch_calc(y,genesets,resrows=10,bootstraps=100,priority="confidence",minsetsize=5)

# generate the pdf report
mitch_plots(res,outfile="seurat_mitch_plots.pdf")

# generate the html report
mitch_report(res,"seurat_mitch_report.html")


