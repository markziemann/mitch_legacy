source("nDrich.R")

#Read in tsv file
x<-as.matrix(read.table("rna_LGHGvHGHGV_jn.rnk",header=T))

#read in GMT file
#genesets<-GMT2DF("ReactomePathways.gmt")
genesets<-gmt_import("ReactomePathways.gmt")

#run the analysis
res<-endrich(x,genesets,resrows=25)
plot2DSets(res,outfile="outres.pdf")
render_report(res,"outres.html")

# run an example shiny document
rmarkdown::run("test_shiny.Rmd")
