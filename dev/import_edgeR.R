
edger2ndrich<-function(x) {

  ndrichScore<-function(y) {
    y<-y[order(y$Row.names),]
    z<-as.data.frame(sign(y$logFC)*-log10(y$PValue))
    row.names(z)=y$Row.names
    
  }
  xx<-lapply(x,ndrichScore)
  xxx<-as.data.frame(do.call(rbind, xx))
  colnames(xxx)=names(x)
  return(xxx) 
}


edger2ndrich(list("ami1_eder"=ami1_edger,"ami5_edger"=ami5_edger))
