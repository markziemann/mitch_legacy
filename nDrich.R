
# v - 2 value vector or 2 col matrix
#s 


diffRankDistance2<-function(v, test, nh) {


	#see stackoverflow.com/questions/22231773/
	h1dist<-sqrt(colSums( (test-t(v) )^2   ))
	h0dist<-sqrt(colSums( (nh-t(v))^2   ))
	rr<-rank( h1dist - h0dist )
	return(rr)

}

proj2DdirMag<-function(v,s=c(1,1) ) {

	qr<- (as.vector(v %*% s) / (s %*% s))
	rrr <-cbind(qr*s[1],qr*s[2])
	return((rrr[,1]^2 + rrr[,2]^2)^0.5  * sign(rrr[,2]))
	

}

# Laod a GMT into dataframe
GMT2DF<-function(gmtfile) {

	rr<-read.table(pipe( paste(" cat ",gmtfile,
		" | awk '{print \"#\"$0}'| tr '\\t' '\\n'|  awk ' /^#/{s=$0; getline;getline} !/^#/ {print s\"\\t\"$1}'|sed 's/^#//'")),
		 header=F)
	return(rr)

}



GSTSets<-function(testset, setsTable) {
	
	#Remove items from set unkown by msigdb
	ntestset = testset[ testset %in% setsTable[,2]  ]
	universesize = length(unique(setsTable[,2]))
	nDrawn<-length(ntestset)
	#Hit vector
	setsTable$hit = setsTable[,2] %in% ntestset
		
	ini<-aggregate(msigdb$hit, list(msigdb[,1]),sum)
	colnames(ini)<-c("Set", "nFound")	
	ini$inCat<- aggregate(msigdb$hit, list(msigdb[,1]),length)[,2]
	ini$NotinCat<- universesize - ini$inCat

	ini$p.value<-apply(ini, 1, function(x){phyper(as.numeric(x[2]),as.numeric(x[3]),as.numeric(x[4]), nDrawn, lower.tail=FALSE )})
	ini$adj.p.value <-p.adjust(ini$p.value,"fdr")
	ini<-ini[order(ini$p.value),]
	return(ini)

}



# x - matrix of 2 columns to be tested
# table of genesets
EnDrichProject<-function(x, geneset, topfig=1) {
	library(limma)
	library(parallel)
	library(plyr)

	xCor<-proj2DdirMag(x)
	xAntiCor<-proj2DdirMag(x,c(1,-1))

	medxCor<-median(xCor)
	medxAntiCor<-median(xAntiCor)

	sets<-unique(geneset[,1])

	res<-mclapply(sets,function(set) { 
		lx<- row.names(x) %in% geneset[geneset[,1]==set,2]
		#cat(sum(lx))
		pCor      <- wilcoxGST(lx, xCor)
		pAntiCor  <- wilcoxGST(lx, xAntiCor)
		
		if( sum(lx[ xCor < medxCor ]) > sum(lx[ xCor > medxCor ])  ) {pCor<- pCor*-1}
		if( sum(lx[ xAntiCor < medxAntiCor ]) > sum(lx[ xAntiCor > medxAntiCor ])  ) {pAntiCor<- pAntiCor*-1}

		return(data.frame(set,pCor,pAntiCor))
	},
	mc.cores=8 )

	return(ldply(res, data.frame))
	
}



# x - matrix of 2 columns to be tested
# table of genesets
EnDrichDist<-function(x, geneset) {
	library(limma)
	library(parallel)
	library(plyr)

#y=x
	p_1_1   = c( max(x[,1]), max(x[,2])        )
	p_n1_n1 = c( min(x[,1]), min(x[,2])         )

#x axis
	p_1_0   = c( max(x[,1]), median(x[,2])      )
	p_n1_0   = c( min(x[,1]), median(x[,2])      )

#y= -x
	p_n1_1  = c( min(x[,1]), max(x[,2])      )
	p_1_n1  = c( max(x[,1]), min(x[,2])      )

	xCorUp      <- diffRankDistance2(x,p_1_1,p_1_0)
	xCorDn      <- diffRankDistance2(x,p_n1_n1,p_n1_0)
	xAntiCorUp  <- diffRankDistance2(x,p_1_n1,p_1_0)
	xAntiCorDn  <- diffRankDistance2(x,p_n1_1,p_n1_0)

	sets<-unique(as.character(geneset[,1]))
	cat(sets)
	res<-mclapply(sets,function(set) { 
#		cat("FSADSS")
		lx<- row.names(x) %in% geneset[geneset[,1]==set,2]
#cat("FddddddSADSS")
#		cat(sum(lx))
#		pCorUp     <- wilcoxGST(lx, xCorUp,alternative="down")
#		pCorDn     <- wilcoxGST(lx, xCorDn,alternative="down")
#		pAntiCorUp <- wilcoxGST(lx, xAntiCorUp,alternative="down")
#		pAntiCorDn <- wilcoxGST(lx, xAntiCorDn,alternative="down")
		pCorUp     <- wilcoxGST(lx, xCorUp)
		pCorDn     <- wilcoxGST(lx, xCorDn)
		pAntiCorUp <- wilcoxGST(lx, xAntiCorUp)
		pAntiCorDn <- wilcoxGST(lx, xAntiCorDn)
#		pCorUp     <- wilcox.test(xCorUp[lx], xCorUp[!lx], alternative="down")
#		pCorDn     <- wilcox.test(xCorDn[lx], xCorDn[!lx],alternative="down")
#		pAntiCorUp <- wilcox.test(xAntiCorUp[lx], xAntiCorUp[!lx],  alternative="down")
#		pAntiCorDn <- wilcox.test(xAntiCorDn[lx], xAntiCorDn[!lx],  alternative="down")

#		cat("FFF")
		return(data.frame(set,setSize=sum(lx),pCorUp, pCorDn, pAntiCorUp, pAntiCorDn))
	},
	mc.cores=8 )

	return(ldply(res, data.frame))
	
}





# x - matrix of 2 columns to be tested
# table of genesets
EnDrichDist3<-function(x, geneset, topfig=1) {
	library(limma)
	library(parallel)
	library(plyr)

#y=x
	p_1_1   = c( max(x[,1]), max(x[,2]) ,max(x[,3])       )
	#p_n1_n1 = c( min(x[,1]), min(x[,2])         )

#x axis
	p_1_0   = c( max(x[,1]), median(x[,2])      )
	p_n1_0   = c( min(x[,1]), median(x[,2])      )

#y= -x
	p_n1_1  = c( min(x[,1]), max(x[,2])      )
	p_1_n1  = c( max(x[,1]), min(x[,2])      )

	xCorUp      <- diffRankDistance2(x,p_1_1,p_1_0)
	xCorDn      <- diffRankDistance2(x,p_n1_n1,p_n1_0)
	xAntiCorUp  <- diffRankDistance2(x,p_1_n1,p_1_0)
	xAntiCorDn  <- diffRankDistance2(x,p_n1_1,p_n1_0)

	sets<-unique(geneset[,1])

	res<-mclapply(sets,function(set) { 
		lx<- row.names(x) %in% geneset[geneset[,1]==set,2]
		#cat(sum(lx))
		pCorUp     <- wilcoxGST(lx, xCorUp)
		pCorDn     <- wilcoxGST(lx, xCorDn)
		pAntiCorUp <- wilcoxGST(lx, xAntiCorUp)
		pAntiCorDn <- wilcoxGST(lx, xAntiCorDn)
		return(data.frame(set,setSize=sum(lx),pCorUp, pCorDn, pAntiCorUp, pAntiCorDn))
	},
	mc.cores=8 )

	return(ldply(res, data.frame))
	
}



#TODO does not work with neg numbers !!
EnDrichMANOVA<-function(x,geneset, minsetsize=10) {
	library(parallel)
	library(plyr)

	#Just in case you did not rank
	x<-apply(x,2,rank)

	sets<-unique(geneset[,1])

	res<-mclapply(sets,function(set) { 

		inset<- row.names(x) %in% geneset[geneset[,1]==set,2]

		fit<- manova(x ~ inset)
		sumMANOVA <- summary.manova(fit)
		sumAOV    <- summary.aov(fit)
		
		pMANOVA <- sumMANOVA$stats[1,"Pr(>F)"]

		raov<-sapply(sumAOV, function(zz) {zz[1,"Pr(>F)"]})
		names(raov)<-gsub("^ Response ","p",names(raov))
		
		#S coordinates
		scord<-apply(zdat,2,function(zz){2*(mean(zz[inset])-mean(zz[!inset]))/length(inset)})
		#sDist<-dist(rbind(rep(0,length(scord)),scord))[1]
		names(scord)<-paste0("s",names(scord))


		return(data.frame(set,setSize=sum(inset),pMANOVA,t(scord),t(raov) ))
	},
	mc.cores=8 )
	fres<-ldply(res, data.frame)
	fres<-fres[fres$setSize >=minsetsize,]
	return(fres)

}





plot2DSets <- function(dat, setdb, restable,  resrows=1:20) {

	 jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
	#Contour of all the data
	ss<-as.data.frame(dat)
	k<-MASS:::kde2d(ss[,1],ss[,2])

	filled.contour(k, color = jet.colors,

# TODO fix 	
#			plot.title={
#				title(
#					xlab=paste( colnames(ss)[1],"rank"),
#					mtext(paste( colnames(ss)[2],"rank"),2,cex=2,line=3,las=1),
#					main=paste("All Data"))
#				)
#			}
	)

	for(i in resrows) {
		ll<-restable[i,]
		size<-ll$setSize
		ss<-as.data.frame(dat[rownames(dat) %in% setdb[setdb[,1]== ll$set,2],])
		k<-MASS:::kde2d(ss[,1],ss[,2])

		filled.contour(k, color = jet.colors,
	
			plot.title={
				title(
					xlab=paste( colnames(ss)[1],"rank"),
					mtext(paste( colnames(ss)[2],"rank"),2,cex=2,line=3,las=1),
					main=paste(ll$set,"\n(",size,")",ll$variable,format(ll$value,digits=3))
				)
				   
			}
		)

	}
}


RankRankBinPlot<-function(x, binsize=500) {
	library(ggplot2)
	 bin=floor(x[,1]/binsize)
	 yy<- aggregate(x[,2],list(bin), function(zz){quantile(zz,c(0.25,0.5,0.75))})
	yy[,1]<-yy[,1]*binsize	
	zz<-cbind(xpos=yy$Group.1,yy$x)
	 colnames(zz)<-gsub("%","",colnames(zz))
         zz<-data.frame(zz)
	 ggplot(zz, aes(xpos)) + 
	geom_line(aes(y=X50) , size=1.4) + 
	theme_bw() + geom_ribbon(aes(ymin=X25, ymax=X75), alpha=0.2) +
	xlab(colnames(x)[1]) + ylab(colnames(x)[2])

}

