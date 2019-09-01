install.packages("..",repos = NULL, type="source")
library(mitch)

rna<-read.table("rna.tsv",header=T)
rna<-rna[,1:6]

k9a<-read.table("k9a.tsv",header=T)
k9a<-k9a[,1:6]

k36a<-read.table("k36a.tsv",header=T)
k36a<-k36a[,1:6]

myList<-list("rna"=rna,"k9a"=k9a,"k36a"=k36a)

myImportedData<-mitch_import(myList,"edger",geneIDcol="Name")
myImportedData<-head(myImportedData,2000)

rna<-rna[rna$Name %in% rownames(myImportedData),]
k9a<-k9a[k9a$Name %in% rownames(myImportedData),]
k36a<-k36a[k36a$Name %in% rownames(myImportedData),]

genesetsExample<-gmt_import("ReactomePathways.gmt")

resExample<-mitch_calc(myImportedData,genesetsExample,resrows=5,priority="significance")

save(genesetsExample, file="../data/genesetsExample.RData",compress="xz")
save(k36a, file="../data/k36a.RData",compress="xz")
save(k9a, file="../data/k9a.RData",compress="xz")
save(rna, file="../data/rna.RData",compress="xz")
save(myList, file="../data/myList.RData",compress="xz")
save(myImportedData, file="../data/myImportedData.RData",compress="xz")
save(resExample, file="../data/resExample.RData",compress="xz")

write.table(rna,file="../inst/extdata/rna.tsv",sep="\t",quote=F)
write.table(k9a,file="../inst/extdata/k9a.tsv",sep="\t",quote=F)
write.table(k36a,file="../inst/extdata/k36a.tsv",sep="\t",quote=F)

system("head -100 ReactomePathways.gmt > ../inst/extdata/sample_genesets.gmt")
