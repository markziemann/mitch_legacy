library("getDEE2")
library("edgeR")
library("DESeq2")
source("nDrich.R")

mdat<-getDee2Metadata("hsapiens")
samples<-mdat[which(mdat$SRP_accession=="SRP096178"),]
samples<-samples[order(samples$SRR_accession),]

samples$group<-factor(c("Ctrl","Ctrl","Ctrl","Ami1","Ami1","Ami1","Ami5","Ami5","Ami5"),levels=c("Ctrl","Ami1","Ami5"))
samples$label<-c("Ctrl_1","Ctrl_2","Ctrl_3","Ami1_1","Ami1_2","Ami1_3","Ami5_1","Ami5_2","Ami5_3")
SRRlist<-as.vector(samples$SRR_accession)
x<-getDEE2("hsapiens",SRRlist)
x<-Tx2Gene(x)
y<-x$Tx2Gene

#here attach genenames
genenames<-unique(x$TxInfo[,1:2])
rownames(genenames)=genenames$GeneID
genenames$GeneID=NULL
y<-merge(y,genenames,by=0)
rownames(y)<-paste(y$Row.names,y$GeneSymbol,sep="|")
y$Row.names=y$GeneSymbol=NULL

#filter very low genes
y<-y[which(rowSums(y)/ncol(y)>=(10)),]

#EdgeR analysis
design<-model.matrix(~samples$group)
rownames(design)<-samples$SRR_accession
z<-DGEList(counts=y)
z<-calcNormFactors(y)
z <- estimateDisp(y, design,robust=TRUE)

#DE for AMI1
dge=NULL
fit <- glmFit(z, design[,1:2])
lrt <- glmLRT(fit)
dge<-as.data.frame(topTags(lrt,n=Inf))
dge$dispersion<-lrt$dispersion
dge<-merge(dge,lrt$fitted.values,by='row.names')
rownames(dge)=dge$Row.names
dge$Row.names=NULL
dge<-merge(dge,z$counts,by='row.names')
ami1_edger<-dge[order(dge$PValue),]
ami1_edger$GeneSymbol<-sapply((strsplit(ami1_edger$Row.names,"\\|")),"[[",2)

#DE for AMI5
dge=NULL
fit <- glmFit(z, design[,c(1,3)])
lrt <- glmLRT(fit)
dge<-as.data.frame(topTags(lrt,n=Inf))
dge$dispersion<-lrt$dispersion
dge<-merge(dge,lrt$fitted.values,by='row.names')
rownames(dge)=dge$Row.names
dge$Row.names=NULL
dge<-merge(dge,z$counts,by='row.names')
ami5_edger<-dge[order(dge$PValue),]
ami5_edger$GeneSymbol<-sapply((strsplit(ami5_edger$Row.names,"\\|")),"[[",2)

save.image("ami1ami5_dge.RData")

#now do DESeq2
y<-round(y)

#ami1
des<-samples[1:6,]
des$grp<-as.numeric(grepl("Ami",des$group))
dds <- DESeqDataSetFromMatrix(countData =y[,1:6], colData = des, design = ~ grp)
res <- DESeq(dds)
z<- results(res)
vsd <- vst(dds, blind=FALSE)
#stick on the normalised expression values to the table
zz<-cbind(z,assay(vsd))
#sort by adjusted p-value
zz<-as.data.frame(zz[order(zz$pvalue),])
ami1_deseq2<-zz
ami1_deseq2$GeneSymbol<-sapply((strsplit(rownames(ami1_deseq2),"\\|")),"[[",2)


#AMI5
des<-samples[c(1:3,7:9),]
des$grp<-as.numeric(grepl("Ami",des$group))
dds <- DESeqDataSetFromMatrix(countData =y[,c(1:3,7:9)], colData = des, design = ~ grp)
res <- DESeq(dds)
z<- results(res)
vsd <- vst(dds, blind=FALSE)
#stick on the normalised expression values to the table
zz<-cbind(z,assay(vsd))
#sort by adjusted p-value
zz<-as.data.frame(zz[order(zz$pvalue),])
ami5_deseq2<-zz
ami5_deseq2$GeneSymbol<-sapply((strsplit(rownames(ami5_deseq2),"\\|")),"[[",2)

save.image("ami1ami5_dge.RData")

#read in GMT file
genesets<-gmt_import("ReactomePathways.gmt")

#edger analysis
x1<-list("ami1_edger"=ami1_edger,"ami5_edger"=ami5_edger)
y1<-ndrich_import(x1,"edger",geneIDcol="GeneSymbol")
res<-endrich(y1,genesets,resrows=25)
plotSets(res,outfile="ami1ami5_edger.pdf")
render_report(res,"ami1ami5_edger.html")

#deseq analysis
x2<-list("ami1_deseq2"=ami1_deseq2,"ami5_deseq2"=ami5_deseq2)
y2<-ndrich_import(x2,"deseq2",geneIDcol="GeneSymbol")
res<-endrich(y2,genesets,resrows=25)
plotSets(res,outfile="ami1ami5_deseq.pdf")
render_report(res,"ami1ami5_deseq.html")


