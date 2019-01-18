source("nDrich.R")

#Read in tsv file
x<-as.matrix(read.table("rna_LGHGvHGHGV_jn.rnk",header=T))

#read in GMT file
genesets<-GMT2DF("ReactomePathways.gmt")

#run the analysis
res<-endrich(x,genesets)
pdf("mRNA_LGvHGvHGV_ReactomeResults.pdf") ; plot2DSets(res,resrows=1:10) ; dev.off()
render_report(res,"outres.html")


