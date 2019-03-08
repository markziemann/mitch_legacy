install.packages("devtools")
library("devtools")

devtools::install_github("markziemann/Mitch")
library("mitch")

download.file("https://raw.githubusercontent.com/markziemann/Mitch/master/rna_LGHGvHGHGV_jn.rnk",destfile="rna_LGHGvHGHGV_jn.rnk")

x<-as.matrix(read.table("rna_LGHGvHGHGV_jn.rnk",header=T))

download.file("https://raw.githubusercontent.com/markziemann/Mitch/master/getReactomeGMT.sh", destfile="getReactomeGMT.sh")

system("chmod +x getReactomeGMT.sh")
system("./getReactomeGMT.sh")

genesets<-gmt_import("ReactomePathways.gmt")

#run the analysis
res<-mitch_calc(x,genesets,resrows=25,bootstraps=100,priority="effect")
mitch_plots(res,outfile="outres.pdf")
mitch_report(res,"outres.html")

