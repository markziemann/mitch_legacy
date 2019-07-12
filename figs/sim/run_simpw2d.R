library("tidyverse")
library("parallel")
library("topconfects")
library("edgeR")
library("DESeq2")
library("limma")
library("ABSSeq")
library("stringi")

source("simpw_func.R")

# obtain count data
a<-countData()

# generate some gene sets
gsets<-randomGeneSets(a)

# create some random data with two contraste
N_REPS=5 ; SUM_COUNT=40000000 ; VARIANCE=0.2 ; FRAC_DE=0.1 ; FC=1 ; DGE_FUNC="deseq2"

x<-simrna2d(a,REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,gsets)

# run the DESeq2 DE analysis
x<-deseq2(x)

# run mitch
x<-run_mitch(x,DGE_FUNC,gsets, N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC)

# look at the results
str(x)

###############################################
# run simulations over a range of parameters
###############################################
SIMS=10
unlink("simpw_res_running.tsv")
res=NULL
for ( FRAC_DE in c(0.2)) {
  PDFNAME=paste(FRAC_DE,"_pw.pdf",sep="")
  pdf(file=PDFNAME,width=11.7,height=6.9)
  for (FC in c(1)) {
    par(mfrow=c(3,3))
    for (N_REPS in c(3,5,10)) {
      for (DGE_FUNC in c("deseq2")) {
        for ( SUM_COUNT in c(10000000,40000000,100000000)) {
          for  ( VARIANCE in c(0,0.2,0.3,0.4,0.5)) {
            x<-agg_dge(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,gsets)
            x<-as.data.frame(do.call(rbind, x))
            write.table(x,file="simpw_res_running.tsv",quote=F,sep='\t',append=T)
            res=c(res,x)
          }
        }
      }
    }
  }
}
