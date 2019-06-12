library("mitch")
# get reactome gene sets
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip",destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets<-gmt_import("ReactomePathways.gmt")
head(genesets)

# get deseq toptables for Ami1 and Ami5
ami1<-read.table("ami1.tsv",header=T,row.names=1)
ami5<-read.table("ami5.tsv",header=T,row.names=1)
head(ami1)
head(ami5)

# get gt file
gt<-read.table("gt.tsv",header=T)

# mitch import
x<-list("ami1"=ami1,"ami5"=ami5)
y<-mitch_import(x,"deseq2",geneTable=gt)
head(y)
