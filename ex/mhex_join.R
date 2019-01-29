rna<-read.table("Galaxy121-edgeR_STAR_f.xls.rnk",header=T,row.names=1)
k9ac<-read.table("Galaxy124_H3K9K14ac_f.xls.rnk",header=T,row.names=1)
k36ac<-read.table("Galaxy128_H3K36ac_f.xls.rnk",header=T,row.names=1)

x<-merge(rna,k9ac,by=0)
y<-merge(x,k36ac,by.x="Row.names",by.y=0)
rownames(y)<-y$Row.names
y$Row.names=NULL
colnames(y)=c("rna","k9ac","k36ac")

write.table("mh_exercise.tsv",sep="\t",quote=F)



