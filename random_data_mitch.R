source("nDrich.R")

#Read in tsv file
x<-as.matrix(read.table("rna_LGHGvHGHGV_jn.rnk",header=T))

#read in GMT file
genesets<-gmt_import("ReactomePathways.gmt")

#run the analysis
res<-mitch_calc(x,genesets,resrows=25)
mitch_plots(res,outfile="outres.pdf")
mitch_report(res,"outres.html")

# run an example shiny document
rmarkdown::run("test_shiny.Rmd")


#randomise rownames test
xx<-x
rownames(xx)<-sample(rownames(x))
head(xx)
res<-mitch_calc(xx,genesets)
mitch_report(res,"randres_names.html")

#randomise all data
xx<-x
xx[,1]<-sample(x[,1])
xx[,2]<-sample(x[,2])
head(xx)
res<-mitch_calc(xx,genesets)
mitch_report(res,"randres_xy.html")


#run some permutes
PERMUTES=100
numsig1=NULL
print(paste("Running",PERMUTES,"permutes"))
for (i in 1:PERMUTES) {
set.seed(i+100)
xx<-x
rownames(xx)<-sample(rownames(x))
head(xx)
res<-mitch_calc(xx,genesets)
n<-length(which(res$manova_result$p.adjustMANOVA<0.05))
numsig1=c(numsig1,n)
}
numsig1<-as.data.frame(numsig1)
numsig1

numsig2=NULL
print(paste("Running",PERMUTES,"permutes"))
for (i in 1:PERMUTES) {
xx<-x
set.seed(i+200)
xx[,1]<-sample(x[,1])
set.seed(i+300)
xx[,2]<-sample(x[,2])
res<-mitch_calc(xx,genesets)
n<-length(which(res$manova_result$p.adjustMANOVA<0.05))
numsig2=c(numsig2,n)
}
numsig2<-as.data.frame(numsig2)
numsig2

MAX=max(c(max(numsig1),max(numsig2)))

pdf("randres.pdf")
par(mfrow=c(2,1))
plot(numsig1$numsig1,main="Gene name randomisation",ylab="No. FDR MANOVA<0.05 sets",xlab="Run",pch=19,ylim=c(0,MAX))
plot(numsig2$numsig2,main="Profile randomisation",ylab="No. FDR MANOVA<0.05 sets",xlab="Run",pch=19,ylim=c(0,MAX))
dev.off()


