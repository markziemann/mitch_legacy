library(mitch)

rna<-read.table(system.file("rna.tsv", package = "mitch"),header=T)
k9a<-read.table(system.file("k9a.tsv", package = "mitch"),header=T)
x<-list("rna"=rna,"k9a"=k9a)

y<-mitch_import(x,"edger",geneIDcol="Name")

genesets<-gmt_import(system.file("sample_genesets.gmt",package="mitch"))

res<-mitch_calc(y,genesets,resrows=5,priority="effect",cores=2)

mitch_plots(res,outfile="2dcharts.pdf",cores=2)

unlink("2dreport.html")

mitch_report(res,"2dreport.html")

detach("package:mitch", unload=TRUE)

