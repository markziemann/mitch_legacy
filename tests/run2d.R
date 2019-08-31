detach("package:mitch", unload=TRUE)
remove.packages("mitch")
install.packages("pkg/mitch/",repos = NULL, type="source")
library(mitch)

rna<-read.table("Galaxy121-edgeR_STAR_f.xls",header=T)
k9a<-read.table("Galaxy124_H3K9K14ac_f.xls",header=T)
x<-list("rna"=rna,"k9a"=k9a)

y<-mitch_import(x,"edger",geneIDcol="Name")

#download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip", destfile="ReactomePathways.gmt.zip")
#unzip("ReactomePathways.gmt.zip")
genesets<-gmt_import("ReactomePathways.gmt")

res<-mitch_calc(y,genesets,resrows=16,priority="significance")

system.time(mitch_plots(res,outfile="2dcharts.pdf"))

unlink("2dreport.html")

mitch_report(res,"2dreport.html")

detach("package:mitch", unload=TRUE)

