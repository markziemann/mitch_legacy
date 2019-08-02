#install.packages("devtools")
#devtools::install_github("markziemann/dee2/getDEE2")
#devtools::install_github("markziemann/Mitch")
#devtools::install_github("hrbrmstr/taucharts")

library("devtools")
library("getDEE2")
library("mitch")
library("taucharts")
library("DESeq2")


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
LGvHG <- results(y)
LGvHG<-as.data.frame(LGvHG[order(LGvHG$pvalue),])
rownames(LGvHG)<-sapply(strsplit(rownames(LGvHG),"\\."),"[[",1)

# run DESeq2 for HG vs HGVPA
y <- DESeqDataSetFromMatrix(countData = round(x2), colData = s2, design = ~ VPA)
y <- DESeq(y)
HGvHGVPA <- results(y)
HGvHGVPA<-as.data.frame(HGvHGVPA[order(HGvHGVPA$pvalue),])
rownames(HGvHGVPA)<-sapply(strsplit(rownames(HGvHGVPA),"\\."),"[[",1)

# import
d<-list("LGvHG"=LGvHG,"HGvHGVPA"=HGvHGVPA)
d<-mitch_import(d,DEtype="deseq2",geneTable=gt)

##################################################
# get the gene sets and import
##################################################
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip", destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets<-gmt_import("ReactomePathways.gmt")

##################################################
# run the analysis significance vs effect
##################################################
res<-mitch_calc(d,genesets,resrows=50,bootstraps=500,priority="effect")
mitch_plots(res,outfile="HGVPA_eff.pdf")
mitch_report(res,"HGVPA_eff.html")

res<-mitch_calc(d,genesets,resrows=50,bootstraps=500,priority="significance")
mitch_plots(res,outfile="HGVPA_sig.pdf")
mitch_report(res,"HGVPA_sig.html")

##################################################
# bootstrap saturation analysis reactome
##################################################

bb=NULL
for (i in 1:10) {
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=5,priority="confidence")
  b0005<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=10,priority="confidence")
  b0010<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=20,priority="confidence")
  b0020<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=50,priority="confidence")
  b0050<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=100,priority="confidence")
  b0100<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=200,priority="confidence")
  b0200<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=500,priority="confidence")
  b0500<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=1000,priority="confidence")
  b1000<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=2000,priority="confidence")
  b2000<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=5000,priority="confidence")
  b5000<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=10000,priority="confidence")
  b10000<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=100000,priority="confidence")
  b100000<-res$manova_result

  b0005r<-b0005[,c(1,9)] 
  b0010r<-b0010[,c(1,9)]
  b0020r<-b0020[,c(1,9)]
  b0050r<-b0050[,c(1,9)]
  b0100r<-b0100[,c(1,9)]
  b0200r<-b0200[,c(1,9)]
  b0500r<-b0500[,c(1,9)]
  b1000r<-b1000[,c(1,9)]
  b2000r<-b2000[,c(1,9)]
  b5000r<-b5000[,c(1,9)]
  b10000r<-b10000[,c(1,9)]
  b100000r<-b100000[,c(1,9)]

  b<-join_all(list(b0005r,b0010r,b0020r,b0050r,b0100r,b0200r,b0500r,b1000r,b2000r,b5000r,b10000r,b100000r),by="set")
  colnames(b)<-c("set","b0005","b0010","b0020","b0050","b0100","b0200","b0500","b1000","b2000","b5000","b10000","b100000")
  rownames(b)<-b$set
  b$set=NULL

  bb<-c(bb,cor(b)[,ncol(b)])
}

save.image("fig1.RData")

b0005m<-bb[which(names(bb)=="b0005")]
b0010m<-bb[which(names(bb)=="b0010")]
b0020m<-bb[which(names(bb)=="b0020")]
b0050m<-bb[which(names(bb)=="b0050")]
b0100m<-bb[which(names(bb)=="b0100")]
b0200m<-bb[which(names(bb)=="b0200")]
b0500m<-bb[which(names(bb)=="b0500")]
b1000m<-bb[which(names(bb)=="b1000")]
b2000m<-bb[which(names(bb)=="b2000")]
b5000m<-bb[which(names(bb)=="b5000")]
b10000m<-bb[which(names(bb)=="b10000")]
b100000m<-bb[which(names(bb)=="b100000m")]

pdf("bootstrap_saturation.pdf",width=10,height=10)
boxplot(b0005m,b0010m,b0020m,b0050m,b0100m,b0200m,b0500m,b1000m,b2000m,b5000m,b10000m, 
 names = c("5", "10", "20", "50", "100", "200", "500", "1000", "2000", "5000", "10000") ,
 ylab="Pearson correlation with 50k bootstraps",
 xlab="no. bootstraps")
dev.off()


q()
##################################################
# get the msigdb gene sets and import
##################################################
library(msigdbr)
m_df = msigdbr(species = "Homo sapiens")
genesets = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

##################################################
# run the analysis significance vs effect
##################################################
res<-mitch_calc(d,genesets,resrows=50,bootstraps=500,priority="effect")
mitch_plots(res,outfile="HGVPA_m_eff.pdf")
mitch_report(res,"HGVPA_m_eff.html")

res<-mitch_calc(d,genesets,resrows=50,bootstraps=500,priority="significance")
mitch_plots(res,outfile="HGVPA_m_sig.pdf")
mitch_report(res,"HGVPA_m_sig.html")

##################################################
# bootstrap saturation analysis reactome
##################################################

bb=NULL
for (i in 1:10) {
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=5,priority="confidence")
  b0005<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=10,priority="confidence")
  b0010<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=20,priority="confidence")
  b0020<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=50,priority="confidence")
  b0050<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=100,priority="confidence")
  b0100<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=200,priority="confidence")
  b0200<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=500,priority="confidence")
  b0500<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=1000,priority="confidence")
  b1000<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=2000,priority="confidence")
  b2000<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=5000,priority="confidence")
  b5000<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=10000,priority="confidence")
  b10000<-res$manova_result
  res<-mitch_calc(d,genesets,resrows=50,bootstraps=100000,priority="confidence")
  b100000<-res$manova_result

  b0005r<-b0005[,c(1,9)]
  b0010r<-b0010[,c(1,9)]
  b0020r<-b0020[,c(1,9)]
  b0050r<-b0050[,c(1,9)]
  b0100r<-b0100[,c(1,9)]
  b0200r<-b0200[,c(1,9)]
  b0500r<-b0500[,c(1,9)]
  b1000r<-b1000[,c(1,9)]
  b2000r<-b2000[,c(1,9)]
  b5000r<-b5000[,c(1,9)]
  b10000r<-b10000[,c(1,9)]
  b100000r<-b100000[,c(1,9)]

  b<-join_all(list(b0005r,b0010r,b0020r,b0050r,b0100r,b0200r,b0500r,b1000r,b2000r,b5000r,b10000r,b100000r),by="set")
  colnames(b)<-c("set","b0005","b0010","b0020","b0050","b0100","b0200","b0500","b1000","b2000","b5000","b10000","b100000")
  rownames(b)<-b$set
  b$set=NULL

  bb<-c(bb,cor(b)[,ncol(b)])
}

save.image("fig1.RData")

b0005m<-bb[which(names(bb)=="b0005")]
b0010m<-bb[which(names(bb)=="b0010")]
b0020m<-bb[which(names(bb)=="b0020")]
b0050m<-bb[which(names(bb)=="b0050")]
b0100m<-bb[which(names(bb)=="b0100")]
b0200m<-bb[which(names(bb)=="b0200")]
b0500m<-bb[which(names(bb)=="b0500")]
b1000m<-bb[which(names(bb)=="b1000")]
b2000m<-bb[which(names(bb)=="b2000")]
b5000m<-bb[which(names(bb)=="b5000")]
b10000m<-bb[which(names(bb)=="b10000")]
b100000m<-bb[which(names(bb)=="b100000m")]

pdf("bootstrap_saturation.pdf",width=10,height=10)
boxplot(b0005m,b0010m,b0020m,b0050m,b0100m,b0200m,b0500m,b1000m,b2000m,b5000m,b10000m,
 names = c("5", "10", "20", "50", "100", "200", "500", "1000", "2000", "5000", "10000") ,
 ylab="Pearson correlation with 100k bootstraps",
 xlab="no. bootstraps")
dev.off()


