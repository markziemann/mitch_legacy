library("getDEE2")
library("edgeR")
library("DESeq2")

mdat<-getDee2Metadata("hsapiens")
samples<-mdat[which(mdat$SRP_accession=="SRP096178"),]
samples<-samples[order(samples$SRR_accession),]

samples$group<-factor(c("Ctrl","Ctrl","Ctrl","Ami1","Ami1","Ami1","Ami5","Ami5","Ami5"),levels=c("Ctrl","Ami1","Ami5"))
samples$label<-c("Ctrl_1","Ctrl_2","Ctrl_3","Ami1_1","Ami1_2","Ami1_3","Ami5_1","Ami5_2","Ami5_3")
SRRlist<-as.vector(samples$SRR_accession)
x<-getDEE2("hsapiens",SRRlist)
x<-Tx2Gene(x)
y<-x$Tx2Gene
y<-y[which(rowSums(y)/ncol(y)>=(10)),]

#EdgeR analysis
design<-model.matrix(~samples$group)
rownames(design)<-samples$SRR_accession
y<-DGEList(counts=y)
y<-calcNormFactors(y)
y <- estimateDisp(y, design,robust=TRUE)

#DE for AMI1
dge=NULL
fit <- glmFit(y, design[,1:2])
lrt <- glmLRT(fit)
dge<-as.data.frame(topTags(lrt,n=Inf))
dge$dispersion<-lrt$dispersion
dge<-merge(dge,lrt$fitted.values,by='row.names')
rownames(dge)=dge$Row.names
dge$Row.names=NULL
dge<-merge(dge,y$counts,by='row.names')
ami1_edger<-dge[order(dge$PValue),]


#DE for AMI5
dge=NULL
fit <- glmFit(y, design[,c(1,3)])
lrt <- glmLRT(fit)
dge<-as.data.frame(topTags(lrt,n=Inf))
dge$dispersion<-lrt$dispersion
dge<-merge(dge,lrt$fitted.values,by='row.names')
rownames(dge)=dge$Row.names
dge$Row.names=NULL
dge<-merge(dge,y$counts,by='row.names')
ami5_edger<-dge[order(dge$PValue),]

save.image("ami1ami5_dge.RData")

#now do DESeq2
y<-x$Tx2Gene
y<-round(y[which(rowSums(y)/ncol(y)>=(10)),])

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

save.image("ami1ami5_dge.RData")

