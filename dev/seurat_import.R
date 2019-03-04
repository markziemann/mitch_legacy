COLDUCT<-read.table("COLDUCT.WCvWD.tsv")
MACROPHAGE<-read.table("MACROPHAGE.WCvWD.tsv")
MESANGIAL<-read.table("MESANGIAL.WCvWD.tsv")
NK_TLYMPH<-read.table("NK_TLYMPH.WCvWD.tsv")


COLDUCT$COLDUCT<-sign(COLDUCT$avg_logFC)*-log10(COLDUCT$p_val)
COLDUCT$rn<-rownames(COLDUCT)

MACROPHAGE$MACROPHAGE<-sign(MACROPHAGE$avg_logFC)*-log10(MACROPHAGE$p_val)
MACROPHAGE$rn<-rownames(MACROPHAGE)

MESANGIAL$MESANGIAL<-sign(MESANGIAL$avg_logFC)*-log10(MESANGIAL$p_val)
MESANGIAL$rn<-rownames(MESANGIAL)

NK_TLYMPH$NK_LYMPH<-sign(NK_TLYMPH$avg_logFC)*-log10(NK_TLYMPH$p_val)
NK_TLYMPH$rn<-rownames(NK_TLYMPH)

x<-list("COLDUCT"=COLDUCT,"MACROPHAGE"=MACROPHAGE,"MESANGIAL"=MESANGIAL,"NK_TLYMPH"=NK_TLYMPH)

xxx<-join_all(x,by = "rn", type = 'full',match='first')

xxx<-xxx[,c(1,7:ncol(xxx))]

rownames(xxx)<-sapply(strsplit(xxx$rn,"_"),"[")[1,]

xxx$rn=NULL

library("biomaRt")
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
gt<-getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), values=ids,mart=ensembl)
x4<-merge(xxx,gt,by.x=0,by.y="ensembl_gene_id")
row.names(x4)<-x4$external_gene_name
 x4$Row.names=x4$external_gene_name=NULL

source("../mitch.R")
genesets<-gmt_import("../reactome.v5.2.symbols_mouse.gmt")



# remove genes with >1 NA value
x5<-x4[which(apply(x4, 1, function(x) sum(is.na(x)))<2),]

res<-mitch_calc(x5,genesets,resrows=25,bootstraps=100,priority="significance",minsetsize=5)

mitch_plots(res,outfile="seurat_mitch_plots.pdf")

mitch_report(res,"seurat_mitch_report.html")


