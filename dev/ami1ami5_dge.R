library("getDEE2")
mdat<-getDee2Metadata("hsapiens")
samples<-mdat[which(mdat$SRP_accession=="SRP096178"),]
samples<-samples[order(samples$SRR_accession),]

samples$group<-factor(c("Ctrl","Ctrl","Ctrl","Ami1","Ami1","Ami1","Ami5","Ami5","Ami5"),levels=c("Ctrl","Ami1","Ami5"))


samples$label<-c("Ctrl_1","Ctrl_2","Ctrl_3","Ami1_1","Ami1_2","Ami1_3","Ami5_1","Ami5_2","Ami5_3")
SRRlist<-as.vector(samples$SRR_accession)
x<-getDEE2("hsapiens",SRRlist)
x<-Tx2Gene(x)

library("edgeR")
y<-x$Tx2Gene
y<-y[which(rowSums(y)/ncol(y)>=(10)),]
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

