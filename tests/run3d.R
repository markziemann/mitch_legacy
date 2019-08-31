detach("package:mitch", unload=TRUE)
remove.packages("mitch")
install.packages("pkg/mitch/",repos = NULL, type="source")
library(mitch)

rna<-read.table("Galaxy121-edgeR_STAR_f.xls",header=T)
k9a<-read.table("Galaxy124_H3K9K14ac_f.xls",header=T)
k36a<-read.table("Galaxy128_H3K36ac_f.xls",header=T)
x<-list("rna"=rna,"k9a"=k9a,"k36a"=k36a)

y<-mitch_import(x,"edger",geneIDcol="Name")

#download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip", destfile="ReactomePathways.gmt.zip")
#unzip("ReactomePathways.gmt.zip")
genesets<-gmt_import("ReactomePathways.gmt")

res<-mitch_calc(y,genesets,resrows=16,priority="significance")

system.time(mitch_plots(res,outfile="3dcharts.pdf"))

unlink("3dreport.html")

mitch_report(res,"3dreport.html")

detach("package:mitch", unload=TRUE)

