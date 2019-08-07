library("getDEE2")
library(DESeq2)
#library(mitch)
library(fgsea)
source("mitch.dev.R")
# contiinue from mid fig1

##################################################
# get the gene sets and import
##################################################
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip", destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets<-gmt_import("ReactomePathways.gmt")




##################################################
# Obtain gene expression counts and run DESeq2
##################################################
mdat<-getDee2Metadata("hsapiens")
samplesheet<-mdat[which(mdat$GSE_accession=="GSE109140"),]
samplesheet<-samplesheet[order(samplesheet$SRR_accession),]
samplesheet$HG<-as.factor(as.numeric(grepl("high",samplesheet$experiment_title)))
samplesheet$VPA<-as.factor(as.numeric(grepl("VPA",samplesheet$experiment_title)))

# samplesheets for the 2 contrasts
s1<-samplesheet[grep("VPA",samplesheet$experiment_title,invert=T),]
s2<-samplesheet[grep("high",samplesheet$experiment_title),]

w<-getDEE2("hsapiens",samplesheet$SRR_accession)
x<-Tx2Gene(w)
x<-x$Tx2Gene

# save the genetable for later
gt<-w$GeneInfo[,1,drop=FALSE]
gt$accession<-rownames(gt)

# filter out lowly expressed genes
x<-x[which(rowSums(x)/ncol(x)>=(10)),]

# counts for the 2 contrasts
x1<-x[,which(colnames(x) %in% s1$SRR_accession)]
x2<-x[,which(colnames(x) %in% s2$SRR_accession)]

#run DESeq2 for LG vs HG
y <- DESeqDataSetFromMatrix(countData = round(x1), colData = s1, design = ~ HG)
y <- DESeq(y)
LGvHG <- as.data.frame(results(y))
rownames(LGvHG)<-sapply(strsplit(rownames(LGvHG),"\\."),"[[",1)

##################################################
# FGSEA
##################################################
z<-merge(LGvHG,gt,by=0)

z<-z[,which(colnames(z) %in% c("stat","GeneSymbol"))]

z<-aggregate(. ~ GeneSymbol,z,mean)

zz<-as.vector(z$stat)

names(zz)<-z$GeneSymbol

fgseaRes<-fgsea(pathways=genesets,stats=zz,nperm=10000,minSize=10)

fgseaRes<-fgseaRes[order(fgseaRes$pval),]

##################################################
# MITCH
##################################################
a<-list(LGvHG)

names(a)="LGvHG"

d<-mitch_import(a,DEtype="deseq2",geneTable=gt)

res<-mitch_calc(d,genesets,resrows=50,bootstraps=0,priority="effect")

res<-mitch_calc(d,genesets,resrows=50,bootstraps=0,priority="significance")

mitch_plots(res,outfile="1dcharts.pdf")

save.image("1d_test.Rdata")
