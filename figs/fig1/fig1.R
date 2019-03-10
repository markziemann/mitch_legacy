install.packages("devtools")
library("devtools")

devtools::install_github("markziemann/dee2/getDEE2")
library("getDEE2")

devtools::install_github("markziemann/Mitch")
library("mitch")

devtools::install_github("hrbrmstr/taucharts")
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
download.file("https://raw.githubusercontent.com/markziemann/Mitch/master/getReactomeGMT.sh", destfile="getReactomeGMT.sh")
system("chmod +x getReactomeGMT.sh")
system("./getReactomeGMT.sh")
genesets<-gmt_import("ReactomePathways.gmt")

##################################################
# run the analysis significance vs effect
##################################################
res<-mitch_calc(d,genesets,resrows=50,bootstraps=100,priority="effect")
mitch_plots(res,outfile="HGVPA_eff.pdf")
mitch_report(res,"HGVPA_eff.html")

res<-mitch_calc(d,genesets,resrows=50,bootstraps=100,priority="significance")
mitch_plots(res,outfile="HGVPA_sig.pdf")
mitch_report(res,"HGVPA_sig.html")

##################################################
# bootstrap saturation analysis
##################################################
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


b0005r<-as.data.frame( cbind( b0005$set , rank(-b0005$confESp) ) ) 
b0010r<-as.data.frame( cbind( b0010$set , rank(-b0010$confESp) ) )
b0020r<-as.data.frame( cbind( b0020$set , rank(-b0020$confESp) ) )
b0050r<-as.data.frame( cbind( b0050$set , rank(-b0050$confESp) ) )
b0100r<-as.data.frame( cbind( b0100$set , rank(-b0100$confESp) ) )
b0200r<-as.data.frame( cbind( b0200$set , rank(-b0200$confESp) ) )
b0500r<-as.data.frame( cbind( b0500$set , rank(-b0500$confESp) ) )
b1000r<-as.data.frame( cbind( b1000$set , rank(-b1000$confESp) ) )
b2000r<-as.data.frame( cbind( b2000$set , rank(-b2000$confESp) ) )
b5000r<-as.data.frame( cbind( b5000$set , rank(-b5000$confESp) ) )

b<-join_all(list(b0005r,b0010r,b0020r,b0050r,b0100r,b0200r,b0500r,b1000r,b2000r,b5000r),by="V1")
colnames(b)<-c("Setname","b0005","b0010","b0020","b0050","b0100","b0200","b0500","b1000","b2000")
rownames(b)<-b$Setname
b$Setname=NULL

bb<-apply(b,2,function(x) as.numeric(as.character(x)) )
plot(cor(bb)[,ncol(bb)])
