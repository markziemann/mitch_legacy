source("nDrich.R")

#in case the data have not been ranked, this function can 
#v<-list("ami1_edger"=ami1_edger,"ami5_edger"=ami5_edger)
#x<-ndrich_import(w,DEtype="edger",geneIDcol="Row.names",geneTable=gt)


#Read in tsv file
x<-as.matrix(read.table("rna_LGHGvHGHGV_jn.rnk",header=T))

#read in GMT file
genesets<-gmt_import("ReactomePathways.gmt")

#run the analysis
res<-endrich(x,genesets,resrows=25,bootstraps=100,priority="effect")
plotSets(res,outfile="outres.pdf")
render_report(res,"outres.html")

# run an example shiny document
rmarkdown::run("test_shiny.Rmd")
