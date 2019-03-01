source("nDrich.R")

#Read in tsv file
x<-as.matrix(read.table("mh_exercise.tsv",header=T))

#read in GMT file
genesets<-gmt_import("ReactomePathways.gmt")

#run the analysis
res<-mitch_calc(x,genesets,resrows=25)
mitch_plots(res,outfile="outres.pdf")
mitch_report(res,"outres.html")

# run an example shiny document
rmarkdown::run("test_shiny.Rmd")
