detach("package:mitch", unload=TRUE)
remove.packages("mitch")
install.packages("pkg/mitch/",repos = NULL, type="source")
library(mitch)

rna<-read.table("Galaxy121-edgeR_STAR_f.xls",header=T)

y<-mitch_import(rna,"edger",geneIDcol="Name")

#download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip", destfile="ReactomePathways.gmt.zip")
#unzip("ReactomePathways.gmt.zip")
genesets<-gmt_import("ReactomePathways.gmt")

res<-mitch_calc(y,genesets,resrows=50,priority="significance")

system.time(mitch_plots(res,outfile="1dcharts.pdf"))

unlink("1dreport.html")

mitch_report(res,"1dreport.html")

detach("package:mitch", unload=TRUE)

