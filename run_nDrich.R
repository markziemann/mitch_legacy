source("nDrich.R")

#Read in tsv file
x<-as.matrix(read.table("rna_LGHGvHGHGV_jn.rnk",header=T))

#read in GMT file
genesets<-GMT2DF("ReactomePathways.gmt")

#run the analysis
res<-endrich(x,genesets)
pdf("mRNA_LGvHGvHGV_ReactomeResults.pdf") ; plot2DSets(res,resrows=1:10) ; dev.off()
render_report(res,"outres.html")

# run an example shiny document
rmarkdown::run("test_shiny.Rmd")


#randomise rownames test
xx<-x
rownames(xx)<-sample(rownames(x))
head(xx)
res<-endrich(xx,genesets)
render_report(res,"randres_names.html")

#randomise all data
xx<-x
xx[,1]<-sample(x[,1])
xx[,2]<-sample(x[,2])
head(xx)
res<-endrich(xx,genesets)
render_report(res,"randres_xy.html")


