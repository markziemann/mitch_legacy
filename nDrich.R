# v - 2 value vector or 2 col matrix
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

# Load a GMT into dataframe
GMT2DF<-function(gmtfile) {
	rr<-read.table(pipe( paste(" cat ",gmtfile,
	" | awk '{print \"#\"$0}'| tr '\\t' '\\n'|  awk ' /^#/{s=$0; getline;getline} !/^#/ {print s\"\\t\"$1}'|sed 's/^#//'")),
	 header=F)
	return(rr)
}

gmt_import<-function(gmtfile){
    pathwayLines <- strsplit(readLines(gmtfile), "\t")
    pathways <- lapply(pathwayLines, utils::tail, -2)
    names(pathways) <- sapply(pathwayLines, head, 1)
    pathways
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
EnDrichProject<-function(x, genesets, topfig=1) {
	library(limma)
	library(parallel)
	library(plyr)

	xCor<-proj2DdirMag(x)
	xAntiCor<-proj2DdirMag(x,c(1,-1))
	medxCor<-median(xCor)
	medxAntiCor<-median(xAntiCor)
	sets<-unique(genesets[,1])
	res<-mclapply(sets,function(set) { 
		lx<- row.names(x) %in% genesets[genesets[,1]==set,2]
		#cat(sum(lx))
		pCor      <- wilcoxGST(lx, xCor)
		pAntiCor  <- wilcoxGST(lx, xAntiCor)
		if( sum(lx[ xCor < medxCor ]) > sum(lx[ xCor > medxCor ])  ) {pCor<- pCor*-1}
		if( sum(lx[ xAntiCor < medxAntiCor ]) > sum(lx[ xAntiCor > medxAntiCor ])  ) {pAntiCor<- pAntiCor*-1}
		return(data.frame(set,pCor,pAntiCor))
	},
	mc.cores=detectCores()-1 )
	return(ldply(res, data.frame))
}


# x - matrix of 2 columns to be tested
# table of genesets
EnDrichDist<-function(x, genesets) {
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

	sets<-unique(as.character(genesets[,1]))
	cat(sets)
	res<-mclapply(sets,function(set) { 
		lx<- row.names(x) %in% genesets[genesets[,1]==set,2]
		pCorUp     <- wilcoxGST(lx, xCorUp)
		pCorDn     <- wilcoxGST(lx, xCorDn)
		pAntiCorUp <- wilcoxGST(lx, xAntiCorUp)
		pAntiCorDn <- wilcoxGST(lx, xAntiCorDn)
		return(data.frame(set,setSize=sum(lx),pCorUp, pCorDn, pAntiCorUp, pAntiCorDn))
	},
	mc.cores=detectCores()-1 )

	return(ldply(res, data.frame))
}


# x - matrix of 2 columns to be tested
# table of genesets
EnDrichDist3<-function(x, genesets, topfig=1) {
	library(limma)
	library(parallel)
	library(plyr)

	p_1_1   = c( max(x[,1]), max(x[,2]) ,max(x[,3])       )

	p_1_0   = c( max(x[,1]), median(x[,2])      )
	p_n1_0   = c( min(x[,1]), median(x[,2])      )

	p_n1_1  = c( min(x[,1]), max(x[,2])      )
	p_1_n1  = c( max(x[,1]), min(x[,2])      )

	xCorUp      <- diffRankDistance2(x,p_1_1,p_1_0)
	xCorDn      <- diffRankDistance2(x,p_n1_n1,p_n1_0)
	xAntiCorUp  <- diffRankDistance2(x,p_1_n1,p_1_0)
	xAntiCorDn  <- diffRankDistance2(x,p_n1_1,p_n1_0)

	sets<-unique(genesets[,1])

	res<-mclapply(sets,function(set) { 
		lx<- row.names(x) %in% genesets[genesets[,1]==set,2]
		#cat(sum(lx))
		pCorUp     <- wilcoxGST(lx, xCorUp)
		pCorDn     <- wilcoxGST(lx, xCorDn)
		pAntiCorUp <- wilcoxGST(lx, xAntiCorUp)
		pAntiCorDn <- wilcoxGST(lx, xAntiCorDn)
		return(data.frame(set,setSize=sum(lx),pCorUp, pCorDn, pAntiCorUp, pAntiCorDn))
	},
	mc.cores=detectCores()-1 )

	return(ldply(res, data.frame))	
}



#TODO does not work with neg numbers !!
EnDrichMANOVA<-function(x,genesets, minsetsize=10, cores=detectCores()-1) {
	library(parallel)
	library(plyr)

	sets<-names(genesets)

	res<-mclapply(sets,function(set){
		inset<-rownames(x) %in% as.character(unlist(genesets[set]))
#		inset<- row.names(x) %in% genesets[genesets[,1]==set,2]
	        fit<- manova(x ~ inset)
        	sumMANOVA <- summary.manova(fit)
	        sumAOV    <- summary.aov(fit)
        	pMANOVA <- sumMANOVA$stats[1,"Pr(>F)"]
	        raov<-sapply(sumAOV, function(zz) {zz[1,"Pr(>F)"]})
        	names(raov)<-gsub("^ Response ","p",names(raov))
	        #S coordinates
	        scord<-apply(x,2,function(zz){2*(mean(zz[inset])-mean(zz[!inset]))/length(inset)})
	        names(scord)<-paste0("s-",names(scord))
        	return(data.frame(set,setSize=sum(inset),pMANOVA,t(scord),t(raov) ))
	},mc.cores=cores )
	fres<-ldply(res, data.frame)
	fres<-fres[fres$setSize >=minsetsize,]
        fres$p.adjustMANOVA<-p.adjust(fres$pMANOVA,"fdr")
        fres$minAbsS<-apply(fres[,4:6],1, function(zz){min(abs(zz))})
	fres<- fres[order(fres$pMANOVA),]
	return(fres)
}





manova_analysis_metrics_calc<-function(x, genesets, manova_result, minsetsize=10 ) {
	library(dplyr)
	num_genesets=length(genesets)
	included_genesets<-nrow(manova_result)
        geneset_counts<-as.data.frame(as.vector(unlist(lapply(genesets,function(set){ length(which(as.vector(unlist(set)) %in% rownames(x))) } ))))
	rownames(geneset_counts)<-names(genesets)
	colnames(geneset_counts)="count"
	genesets_excluded=names(genesets)[which(geneset_counts$count<minsetsize)]
        genesets_included=names(genesets)[which(geneset_counts$count>=minsetsize)]
        num_genesets_excluded=length(genesets_excluded)
        num_genesets_included=length(genesets_included)
	num_genes_in_genesets=length(unique(as.vector(unlist(genesets))))
	num_genes_in_profile=length(unique(rownames(x)))
	duplicated_genes_present=length(rownames(x))>num_genes_in_profile
	num_profile_genes_in_sets=length(which(rownames(x) %in% as.vector(unlist(genesets))))
	num_profile_genes_not_in_sets=num_genes_in_profile - num_profile_genes_in_sets
	profile_pearson_correl=cor(x,method="p")[2,1]
	profile_spearman_correl=cor(x,method="s")[2,1]
	num_sets_significant=nrow( manova_result[which(manova_result$p.adjustMANOVA<0.05),] )

	#genes in each quadrant
        g1=length( which(x[,1]>0 & x[,2]>0) )
	g2=length( which(x[,1]>0 & x[,2]<0) )
        g3=length( which(x[,1]<0 & x[,2]<0) )
        g2=length( which(x[,1]<0 & x[,2]>0) )

	#genesets in each quadrant
	ns1=nrow( subset(manova_result,p.adjustMANOVA<0.05 & manova_result[,4]>0 & manova_result[,5]>0) )
        ns2=nrow( subset(manova_result,p.adjustMANOVA<0.05 & manova_result[,4]>0 & manova_result[,5]<0) )
        ns3=nrow( subset(manova_result,p.adjustMANOVA<0.05 & manova_result[,4]<0 & manova_result[,5]<0) )
        ns4=nrow( subset(manova_result,p.adjustMANOVA<0.05 & manova_result[,4]<0 & manova_result[,5]>0) )
	num_sets_significant_by_quadrant=paste(ns1,ns2,ns3,ns4,sep=",")

	dat <- list("num_genesets" = num_genesets, 
		"num_genes_in_profile" = num_genes_in_profile,
		"duplicated_genes_present" = duplicated_genes_present,
		"num_profile_genes_in_sets" = num_profile_genes_in_sets,
                "num_profile_genes_not_in_sets" = num_profile_genes_not_in_sets,
		"num_genesets_excluded" = num_genesets_excluded,
		"num_genesets_included" = num_genesets_included,
                "num_genes_in_genesets" = num_genes_in_genesets,
		"genesets_excluded" = genesets_excluded,
		"genesets_included" = genesets_included,
		"profile_pearson_correl" = profile_pearson_correl,
		"profile_spearman_correl" = profile_spearman_correl,
		"num_sets_significant" = num_sets_significant,
		"num_sets_significant_by_quadrant" = num_sets_significant_by_quadrant,
		"geneset_counts" = geneset_counts)
	dat
}
#hist(geneset_counts$count,200,xlim=c(0,500))


endrichrank<-function(x) {
  input_profile<-x
  ranked_profile<-apply(x,2,rank)

  x_num_neg=length( which( input_profile[,1]<0 ) )
  x_num_zero=length( which( input_profile[,1]==0 ) )
  x_num_adj=x_num_neg+(x_num_zero/2)

  y_num_neg=length( which( input_profile[,2]<0 ) )
  y_num_zero=length( which( input_profile[,2]==0 ) )
  y_num_adj=y_num_neg+(y_num_zero/2)

  adj_x<-ranked_profile[,1]-x_num_adj
  adj_y<-ranked_profile[,2]-y_num_adj
  adj<-cbind(adj_x,adj_y)

  colnames(adj)=colnames(x)
  adj
}


detailed_sets<-function(res,  resrows=50) {
  #collect ranked genelist of each genest
  ss<-res$ranked_profile
  mykeys <- as.character(res$manova_result[1:resrows,1])
  dat <- vector(mode="list", length=resrows)
  names(dat) <- mykeys

  for(i in 1:resrows) {
    ll<-res$manova_result[i,]
    size<-ll$setSize
    setindex=as.numeric(ll[1])
    sss<-ss[which(rownames(ss) %in% genesets[[as.numeric(ll[1])]]),]

    dat[[i]]<-sss
  }
  dat
}


endrich<-function(x,genesets, minsetsize=10, cores=detectCores()-1 , resrows=50) {
	input_profile<-x

        input_genesets<-genesets

	ranked_profile<-endrichrank(input_profile)

	manova_result<-EnDrichMANOVA(ranked_profile, genesets, minsetsize=minsetsize, cores=cores)

	manova_analysis_metrics<-manova_analysis_metrics_calc(x,genesets,manova_result)

	dat <- list("input_profile" = input_profile,
		"input_genesets" = input_genesets,
		"ranked_profile" = ranked_profile,
		"manova_result" = manova_result,
		"manova_analysis_metrics" = manova_analysis_metrics)

        dat$detailed_sets<-detailed_sets(dat,resrows)

	dat
}


plot2DSets <- function(res,outfile="Rplots.pdf") {
  palette <- colorRampPalette(c("white", "yellow","orange" ,"red","darkred","black"))

  resrows=length(res$detailed_sets)
  #Contour of all the data
  ss<-res$ranked_profile

  xmin=min(ss[,1])
  xmax=max(ss[,1])
  ymin=min(ss[,2])
  ymax=max(ss[,2])

  k<-MASS:::kde2d(ss[,1],ss[,2])
  X_AXIS=paste("Rank in contrast",colnames(ss)[1])
  Y_AXIS=paste("Rank in contrast",colnames(ss)[2])

  pdf(outfile)
  filled.contour(k, xlim=c(xmin,xmax),ylim=c(ymin,ymax),
    color=palette , 
    plot.title={ abline(v=0,h=0,lty=2,lwd=2,col="blue")
       title( main="Rank-rank plot of all genes",xlab=X_AXIS,ylab=Y_AXIS )
    }
  )

  for(i in 1:resrows) {
    ll<-res$manova_result[i,]
    size<-ll$setSize
    sss<-res$detailed_sets[[i]] 
#    sss<-as.data.frame(ss[rownames(ss) %in% res$input_genesets[res$input_genesets[,1]== ll$set,2],])
    k<-MASS:::kde2d(sss[,1],sss[,2])
    filled.contour( k, color = palette, xlim=c(xmin,xmax),ylim=c(ymin,ymax),
        plot.title={ abline(v=0,h=0,lty=2,lwd=2,col="blue")
        title( main=paste(ll$set,"\n(",size,")",ll$variable,format(ll$value,digits=3)),
          xlab=X_AXIS,ylab=Y_AXIS
        )
      }
    )
  }
  dev.off()
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

render_report<-function(res,out) {
  library(knitr)
  library(markdown)
  library(rmarkdown)

  DATANAME<-gsub(".html$",".RData",out)
  save.image(DATANAME)
  MYMESSAGE=paste("Dataset saved as \"",DATANAME,"\".")
  message(MYMESSAGE)

  knitrenv <- new.env()
  assign("DATANAME", DATANAME, knitrenv)
  assign("res",res,knitrenv)
  HTMLNAME=paste(out,".html",sep="")
#  knit2html("nDrich.Rmd", envir=knitrenv , output=out)
  rmarkdown::render("nDrich.Rmd",output_file=out)
}
