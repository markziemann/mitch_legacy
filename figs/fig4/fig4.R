#SRP004777  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3063777/
#comparison of Evaluating Gene Expression in C57BL/6J and DBA/2J Mouse Striatum Using RNA-Seq and Microarrays

library("getDEE2")
library(DESeq2)

mdat<-getDee2Metadata("mmusculus")
m1<-mdat[which(mdat$SRP_accession %in% "SRP004777"),]
SRRlist<-as.vector(m1$SRR_accession)
SRRlist<-SRRlist[order(SRRlist)]
x<-getDEE2("mmusculus",SRRlist)
x<-Tx2Gene(x)
tx<-x$Tx2Gene

#ctrl=c57 trt=DBA
s<-c(1,1,1,1,0,0,0,0,0,0,0,1,1,1,0,0,0,1,1,1,1)
s<-as.data.frame(s)
rownames(s)<-colnames(tx)


tx1<-x$GeneCounts[which(rowSums(tx)/ncol(tx)>=(10)),]

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(tx1), colData = s, design = ~ s)
res <- DESeq(dds)
z<- results(res)

