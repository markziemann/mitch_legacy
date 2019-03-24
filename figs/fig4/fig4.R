#SRP004777  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3063777/
#comparison of Evaluating Gene Expression in C57BL/6J and DBA/2J Mouse Striatum Using RNA-Seq and Microarrays

library("getDEE2")
library("DESeq2")
library("limma")
library("edgeR")
library("ABSSeq")
library("plyr")
library("msigdbr")
library("mitch")

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

# filter counts
tx1<-x$GeneCounts[which(rowSums(tx)/ncol(tx)>=(10)),]

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(tx1), colData = s, design = ~ s)
res <- DESeq(dds)
res_deseq2<- DESeq2::results(res)
res_deseq2$pvalue[is.na(res_deseq2$pvalue)] <- 1
res_deseq2$log2FoldChange[is.na(res_deseq2$log2FoldChange)] <- 0
rnk_deseq2<-as.data.frame( sign(res_deseq2$log2FoldChange) * -log10(res_deseq2$pvalue))
rownames(rnk_deseq2)<-rownames(res_deseq2)
colnames(rnk_deseq2)<-"DESeq2"
rnk_deseq2$geneID<-rownames(rnk_deseq2)

# edgeR GLMRT
design<-model.matrix(~s$s)
rownames(design)<-s$s
y<-DGEList(counts=tx1)
y<-calcNormFactors(y)
y <- estimateDisp(y, design,robust=TRUE)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
res_edger_glmrt<-as.data.frame(topTags(lrt,n=Inf))
rnk_edger_glmrt<-as.data.frame( sign(res_edger_glmrt$logFC) * -log10(res_edger_glmrt$PValue))
rownames(rnk_edger_glmrt)<-rownames(res_edger_glmrt)
colnames(rnk_edger_glmrt)<-"edgeR_GLMRT"
rnk_edger_glmrt$geneID<-rownames(rnk_edger_glmrt)

# edgeR QL
design<-model.matrix(~s$s)
rownames(design)<-s$s
y<-DGEList(counts=tx1)
y<-calcNormFactors(y)
y <- estimateDisp(y, design,robust=TRUE)
fit<-glmQLFit(y, design)
lrt<-glmQLFTest(fit)
res_edger_ql<-as.data.frame(topTags(lrt,n=Inf))
rnk_edger_ql<-as.data.frame( sign(res_edger_ql$logFC) * -log10(res_edger_ql$PValue))
rownames(rnk_edger_ql)<-rownames(res_edger_ql)
colnames(rnk_edger_ql)<-"edgeR_QL"
rnk_edger_ql$geneID<-rownames(rnk_edger_ql)

# Voom-Limma
z<-DGEList(counts=tx1)
z <- calcNormFactors(z)
v <- voom(z,design,plot=F)
fit <- lmFit(v, design)
fit.de <- eBayes(fit, robust=TRUE)
res_voom_limma<-topTable(fit.de,n=Inf)
rnk_voom_limma<-as.data.frame( sign(res_voom_limma$logFC) * -log10(res_voom_limma$P.Value))
rownames(rnk_voom_limma)<-rownames(res_voom_limma)
colnames(rnk_voom_limma)<-"voom-limma"
rnk_voom_limma$geneID<-rownames(rnk_voom_limma)

# ABSseq
obj<-ABSDataSet(tx1, factor(s$s))  #default normalisation is qtotal
obj<-ABSSeq(obj)
dge<- as.data.frame(cbind(obj$Amean,obj$Bmean,obj$foldChange,obj$pvalue,obj$adj.pvalue))
colnames(dge)=c("Amean","Bmean","logFC","PValue","FDR")
res_absseq<-dge[order(dge$PValue),]
rnk_absseq<-as.data.frame( sign(res_absseq$logFC) * -log10(res_absseq$PValue))
rownames(rnk_absseq)<-rownames(res_absseq)
colnames(rnk_absseq)<-"ABSSeq"
rnk_absseq$geneID<-rownames(rnk_absseq)

# genetable
gt<-data.frame(as.character(x$GeneInfo$GeneSymbol))
gt$geneID<-rownames(x$GeneInfo)
colnames(gt)<-c("symbol","geneID")

# join the tables together
xx<-join_all(list("symbol"=gt,"DESeq2"=rnk_deseq2,
  "edgeR_GLMRT"=rnk_edger_glmrt,
  "edgeR_QL"=rnk_edger_ql,"voom-limma"=rnk_voom_limma,
  "ABSSeq"=rnk_absseq),by="geneID")


# aggregate
xx$geneID=NULL
xx<-aggregate(. ~ symbol,xx,sum)
rownames(xx)<-xx$symbol
xx$symbol=NULL


# gene sets from my website but this needs to be updated
genesets<-gmt_import("reactome.v5.2.symbols_mouse.gmt")

# run mitch analysis
mitch_res<-mitch_calc(xx,genesets,resrows=25,bootstraps=100,priority="effect")

