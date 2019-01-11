source("nDrich.R")

#Read in tsv file
zdat<-as.matrix(read.table("rna_LGHGvHGHGV_jn.rnk",header=T))
ss<-apply(zdat,2,rank)

#read in GMT file
reactome<-GMT2DF("ReactomePathways.gmt")

#run the analysis
res<-EnDrichMANOVA(zdat, reactome)
res$p.adjustMANOVA<-p.adjust(res$pMANOVA,"fdr")
res<- res[order(res$pMANOVA),]
res$minAbsS<-apply(res[,4:6],1, function(zz){min(abs(zz))})
write.table(res, "rna_LGHGvHGHGV_ReactomeResults.txt",sep="\t",col.names=T,row.names=F, quote=F)
pdf("mRNA_LGvHGvHGV_ReactomeResults.pdf") ; plot2DSets(ss,reactome,res,resrows=1:200) ; dev.off()
save.image("EnDrichAnalysis.Rdata")



cex<--log10(res$p.adjustMANOVA)/8
plot(res$s.rna_LGvHG,res$s.rna_HGvHGV,pch=19,col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),cex=cex,xlab=colnames(zdat)[1],ylab=colnames(zdat)[2])
abline(h=0,lty=2)
abline(v=0,lty=2)
