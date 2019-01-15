source("nDrich.R")

#Read in tsv file
x<-as.matrix(read.table("rna_LGHGvHGHGV_jn.rnk",header=T))

#read in GMT file
genesets<-GMT2DF("ReactomePathways.gmt")

#run the analysis
res<-endrich(x,genesets)
write.table(res, "rna_LGHGvHGHGV_ReactomeResults.txt",sep="\t",col.names=T,row.names=F, quote=F)
pdf("mRNA_LGvHGvHGV_ReactomeResults.pdf") ; plot2DSets(ss,reactome,res,resrows=1:200) ; dev.off()
save.image("EnDrichAnalysis.Rdata")



