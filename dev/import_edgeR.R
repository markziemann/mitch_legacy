
edger2ndrich<-function(x) {
library(plyr)

ndrichScore<-function(y) {
  z<-as.data.frame(sign(y$logFC)*-log10(y$PValue))
  colnames(z)<-nm(y)
  z$Row.names<-y$Row.names
  z
}

mynames<-nm(x)         
xx<-lapply(x,ndrichScore)
names(xx)<-mynames
xxx<-join_all(xx,by = 'Row.names', type = 'inner',match='first')
rownames(xxx)<-xxx$Row.names
xxx$Row.names=NULL
colnames(xxx)<-names(x)
return(xxx)
}

x<-list("ami1_edger"=ami1_edger,"ami5_edger"=ami5_edger)

yy<-edger2ndrich(x)
