library("mitch")
library("mdgsa")

x<-read.table("https://raw.githubusercontent.com/markziemann/Mitch/master/ex/rna_LGHGvHGHGV_jn.rnk",row.names=1,header=T)

download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip",destfile="ReactomePathways.gmt.zip")

unzip("ReactomePathways.gmt.zip")

#reactome<-read.table("ReactomePathways.gmt")
reactome<-gmt_import("ReactomePathways.gmt")

res.md <- mdGsa (x, reactome)

head(res.md[order(res.md$padj.I),])

