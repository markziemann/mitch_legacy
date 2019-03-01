library("plyr")
library("parallel")
library("pbmcapply")
library("Rmisc")


ndrich_import<-function(x , DEtype, geneIDcol=NULL, geneTable=NULL ) {

if ( !is.list(x) ){
  stop("Error: Input (x) must be a LIST of dataframes.")
}

if ( is.null(names(x)) ){
  stop("Error: Input (x) must be a NAMED list of dataframes.")
}

if (length(which(FALSE==unname(unlist(lapply(x,is.data.frame))))) >0 ) {
  stop("Error: Input (x) must be a named list of DATAFRAMES only.")
}

if ( !is.null(geneTable) && !is.data.frame(geneTable) ) {
  stop("Error: the geneTable needs to be a dataframe.")
}

if ( !is.null(geneTable) &&  (ncol(geneTable)<2 || ncol(geneTable)>2 ) ) {
  stop("Error: the geneTable needs to be a dataframe of 2 columns.")
}

#Accession number to gene ID mapping
mapGeneIds<-function(y,z) {
  if ( !is.null(attributes(y)$geneTable) ) {
    gt<-attributes(y)$geneTable
    col1<-length(which (z$geneidentifiers %in% gt[,1]))
    col2<-length(which (z$geneidentifiers %in% gt[,2]))

    if ( col1 + col2 < (nrow(y)/2) ) {
      stop("Error it looks as if the Gene IDs in the profile don't match the geneTable")
    }

    if ( col1 > col2 ) {
      colnames(gt)=c("geneidentifiers","GeneSymbol")
      z<-merge(gt,z,by="geneidentifiers")
    } else {
      colnames(gt)=c("GeneSymbol","geneidentifiers")
      z<-merge(gt,z,by="geneidentifiers")
    }
    z<-aggregate(. ~ GeneSymbol,z,sum)
    z$geneidentifiers=NULL
    colnames(z)=c("geneidentifiers","y")
  }
  z
}

#the geneIDcol should be an attribute added to each list item
for (i in 1:length(x) ) {
  if ( !is.null(geneIDcol) ) {
    LEN=length( which(names(x1[[i]]) %in% geneIDcol) )
    if (LEN<1) { stop("Error: the specified geneIDcol doesn't seem to exist") }
    if (LEN>1) { stop("Error: there are multiple matches for the  specified geneIDcol") }
    attributes(x[[i]])$geneIDcol<-which( names(x1[[i]]) %in% geneIDcol )
  } else {
    attributes(x[[i]])$geneIDcol<-NULL
  }
  if ( !is.null(geneTable) ) {
    if ( !is.data.frame(geneTable) ) { stop("Error: geneTable is not a data frame.") }
    if ( ncol(geneTable)!=2 ) { stop("Error: geneTable must be a 2 column dataframe.") }
    attributes(x[[i]])$geneTable<-geneTable
  }
}

edger_score<-function(y) {

  NCOL=ncol(y)
  if (NCOL<2){ stop("Error: there are <2 columns in the input, 'PValue' and 'logFC' are required ") }

  PCOL=length(which(names(y)=="PValue"))
  if (PCOL>1){ stop("Error, there is more than 1 column named 'PValue' in the input") }
  if (PCOL<1){ stop("Error, there is no column named 'PValue' in the input") }

  FCCOL=length(which(names(y)=="logFC"))
  if (FCCOL>1){ stop("Error, there is more than 1 column named 'logFC' in the input") }
  if (FCCOL<1){ stop("Error, there is no column named 'logFC' in the input") }

  z<-as.data.frame(sign(y$logFC)*-log10(y$PValue))
  colnames(z)<-"y"
  if ( !is.null(attributes(y)$geneIDcol) ) {
    z$geneidentifiers<-y[,attributes(y)$geneIDcol]
  } else {
    z$geneidentifiers<-rownames(y)
  }
  z<-mapGeneIds(y,z)
  z
}

deseq2_score<-function(y) {

  NCOL=ncol(y)
  if (NCOL<2){ stop("Error: there are <2 columns in the input, 'pvalue' and 'log2FoldChange' are required ") }

  PCOL=length(which(names(y)=="pvalue"))
  if (PCOL>1){ stop("Error, there is more than 1 column named 'pvalue' in the input") }
  if (PCOL<1){ stop("Error, there is no column named 'pvalue' in the input") }

  FCCOL=length(which(names(y)=="log2FoldChange"))
  if (FCCOL>1){ stop("Error, there is more than 1 column named 'log2FoldChange' in the input") }
  if (FCCOL<1){ stop("Error, there is no column named 'log2FoldChange' in the input") }

  z<-as.data.frame(sign(y$log2FoldChange)*-log10(y$pvalue))
  colnames(z)<-"y"
  COL<-attributes(y)$geneIDcol
  if ( !is.null(attributes(y)$geneIDcol) ) {
    z$geneidentifiers<-y[,attributes(y)$geneIDcol]
  } else {
    z$geneidentifiers<-rownames(y)
  }
  z<-mapGeneIds(y,z)
  z
}

limma_score<-function(y) {

  NCOL=ncol(y)
  if (NCOL<2){ stop("Error: there are <2 columns in the input, 'P.Value' and 'logFC' are required ") }

  PCOL=length(which(names(y)=="P.Value"))
  if (PCOL>1){ stop("Error, there is more than 1 column named 'P.Value' in the input") }
  if (PCOL<1){ stop("Error, there is no column named 'P.Value' in the input") }

  FCCOL=length(which(names(y)=="logFC"))
  if (FCCOL>1){ stop("Error, there is more than 1 column named 'logFC' in the input") }
  if (FCCOL<1){ stop("Error, there is no column named 'logFC' in the input") }

  z<-as.data.frame(sign(y$logFC)*-log10(y$P.Value))
  colnames(z)<-"y"
  COL<-attributes(y)$geneIDcol
  if ( !is.null(attributes(y)$geneIDcol) ) {
    z$geneidentifiers<-y[,attributes(y)$geneIDcol]
  } else {
    z$geneidentifiers<-rownames(y)
  }
  z<-mapGeneIds(y,z)
  z
}

absseq_score<-function(y) {

  NCOL=ncol(y)
  if (NCOL<2){ stop("Error: there are <2 columns in the input, 'pvalue' and 'foldChange' are required ") }

  PCOL=length(which(names(y)=="pvalue"))
  if (PCOL>1){ stop("Error, there is more than 1 column named 'pvalue' in the input") }
  if (PCOL<1){ stop("Error, there is no column named 'pvalue' in the input") }

  FCCOL=length(which(names(y)=="foldChange"))
  if (FCCOL>1){ stop("Error, there is more than 1 column named 'foldChange' in the input") }
  if (FCCOL<1){ stop("Error, there is no column named 'foldChange' in the input") }

  z<-as.data.frame(sign(y$foldChange)*-log10(y$pvalue))
  colnames(z)<-"y"
  COL<-attributes(y)$geneIDcol
  if ( !is.null(attributes(y)$geneIDcol) ) {
    z$geneidentifiers<-y[,attributes(y)$geneIDcol]
  } else {
    z$geneidentifiers<-rownames(y)
  }
  z<-mapGeneIds(y,z)
  z
}

sleuth_score<-function(y) {

  NCOL=ncol(y)
  if (NCOL<2){ stop("Error: there are <2 columns in the input, 'pval' and 'b' are required ") }

  PCOL=length(which(names(y)=="pval"))
  if (PCOL>1){ stop("Error, there is more than 1 column named 'pval' in the input") }
  if (PCOL<1){ stop("Error, there is no column named 'pval' in the input") }

  FCCOL=length(which(names(y)=="b"))
  if (FCCOL>1){ stop("Error, there is more than 1 column named 'b' in the input") }
  if (FCCOL<1){ stop("Error, there is no column named 'b' in the input") }

  z<-as.data.frame(sign(y$b)*-log10(y$pval))
  colnames(z)<-"y"
  COL<-attributes(y)$geneIDcol
  if ( !is.null(attributes(y)$geneIDcol) ) {
    z$geneidentifiers<-y[,attributes(y)$geneIDcol]
  } else {
    z$geneidentifiers<-rownames(y)
  }
  z<-mapGeneIds(y,z)
  z
}

topconfect_score<-function(y) {

  FCCOL=length(which(names(y)=="confect"))
  if (FCCOL>1){ stop("Error, there is more than 1 column named 'confect' in the input") }
  if (FCCOL<1){ stop("Error, there is no column named 'confect' in the input") }

  z<-as.data.frame(y$confect)
  z[is.na(z)] <- 0
  colnames(z)<-"y"
  COL<-attributes(y)$geneIDcol
  if ( !is.null(attributes(y)$geneIDcol) ) {
    z$geneidentifiers<-y[,attributes(y)$geneIDcol]
  } else {
    z$geneidentifiers<-rownames(y)
  }
  z<-mapGeneIds(y,z)
  z
}

if ( DEtype == "edger" ) {
  xx<-lapply(x,edger_score)
} else if ( DEtype == "deseq2" ) {
  xx<-lapply(x,deseq2_score)
} else if ( DEtype == "limma" ) {
  xx<-lapply(x,limma_score)
} else if ( DEtype == "absseq" ) {
  xx<-lapply(x,absseq_score)
} else if ( DEtype == "sleuth" ) {
  xx<-lapply(x,sleuth_score)
} else if ( DEtype == "topconfects" ) {
  xx<-lapply(x,topconfect_score)
}

xxx<-join_all(xx,by = 'geneidentifiers', type = 'inner',match='first')
rownames(xxx)<-xxx$geneidentifiers
xxx$geneidentifiers=NULL
colnames(xxx)<-names(x)

MEAN_N_GENES_IN=mean(unlist(lapply(x,nrow)))
N_GENES_OUT=nrow(xxx)
PROP=signif(N_GENES_OUT/MEAN_N_GENES_IN,3)
message(paste("Note: Mean no. genes in input =",MEAN_N_GENES_IN))
message(paste("Note: no. genes in output =",N_GENES_OUT))
if (PROP<0.05) {
  warning("Warning: less than half of the input genes are also in the output")
} else {
  message(paste("Note: estimated proportion of input genes in output =",PROP))
}
return(xxx)
}


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

gmt_import<-function(gmtfile){
    genesetLines <- strsplit(readLines(gmtfile), "\t")
    genesets <- lapply(genesetLines, utils::tail, -2)
    names(genesets) <- sapply(genesetLines, head, 1)
    attributes(genesets)$originfile<-gmtfile
    genesets
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
	library("limma")

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
	library("limma")

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
	library("limma")

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
EnDrichMANOVA<-function(x,genesets, minsetsize=10, cores=detectCores()-1, priority=NULL, bootstraps=0) {

sets<-names(genesets)

if (  is.null(priority) ) {
  priority="confidence"
}

if (priority!="significance" && priority!="effect" && priority!="confidence") {
  stop("Error: Parameter 'priority' must be either 'confidence'(the default),'significance' or 'effect'.")
}

hypotenuse <- function(x){ sqrt(sum(unlist(lapply(x,function(x) {x^2} )))) }

#calculate the hypotenuse for downstream use
HYPOT=hypotenuse(apply(x,2,length))

res<-pbmclapply(sets,function(set){

  resample<-function(x,set){
    sss<-x[which (rownames(x) %in% as.character(unlist(genesets[set]))),]
    mysample<-sss[sample(nrow(sss),nrow(sss),replace=T),]
    colMeans(mysample)
  }

  bootstrap<-function(x,n,set){
    xx<-as.data.frame(t(replicate(n,resample(x,set))),stringsAsFactors=F)
    NOTINSET<-colMeans(x[!inset,])
    NROW=nrow(x)
    xxx<- ( 2* (xx - NOTINSET ) ) / NROW
    b<-apply(xxx,1,hypotenuse)
    return(b)
  }

  inset<-rownames(x) %in% as.character(unlist(genesets[set]))

  NROW=nrow(x)

  if ( length(which(inset)) >= minsetsize ) {
    fit<- manova(x ~ inset)
    sumMANOVA <- summary.manova(fit)
    sumAOV    <- summary.aov(fit)
    pMANOVA <- sumMANOVA$stats[1,"Pr(>F)"]
    raov<-sapply(sumAOV, function(zz) {zz[1,"Pr(>F)"]})
    names(raov)<-gsub("^ Response ","p.",names(raov))
    #S coordinates
    #scord<-apply(x,2,function(zz){ 2 * ( mean(zz[inset]) - mean(zz[!inset]) ) / length(inset) } )

    NOTINSET<-colMeans(x[!inset,])
    scord<- ( 2*( colMeans(x[inset,]) - NOTINSET ) ) / NROW 
    names(scord)<-paste0("s-",names(scord))
    #calculate the hypotenuse length of s scores
    sdist<-hypotenuse(scord)
    names(sdist)="s.dist"

    #confidence interval calc by resampling
    if (bootstraps > 0) {

      STRAPSDONE=0
      BAND_SIZE=1
      CHUNK=50
      b=NULL
      # use an approach to stop bootstrapping un
      while ( BAND_SIZE>0.1 && STRAPSDONE<bootstraps) {
        bb<-bootstrap(x,CHUNK,set)
        b<-c(b,bb)
        STRAPSDONE=STRAPSDONE+CHUNK
        BAND_SIZE<-unname((CI(b)[1]-CI(b)[3] ) /CI(b)[2])
      }
      #b<-bootstrap(x,bootstraps,set)
      confESp<-sdist-quantile(abs(sdist-b),0.05)
      names(confESp)="confESp"
    } else {
      confESp<-NA
      names(confESp)="confESp"
    }

    return(data.frame(set,setSize=sum(inset), pMANOVA, t(scord), t(raov), t(sdist), t(confESp) , stringsAsFactors=F ))
  }

},mc.cores=cores )
fres<-ldply(res, data.frame)
#fres<-fres[fres$setSize >=minsetsize,]
fres$p.adjustMANOVA<-p.adjust(fres$pMANOVA,"fdr")
#fres$minAbsS<-apply(fres[,4:6],1, function(zz){min(abs(zz))})

#prioritisation
if (priority=="significance") {
  fres<-fres[order(fres$pMANOVA),]
  message("Note: When prioritising by significance (ie: small p-values), large effect sizes might be missed.")
}
if (priority=="effect") {
  fres<-fres[order(-fres$s.dist),]
  message("Note: Enrichments with large effect sizes may not be statistically significant.")
}
if (priority=="confidence") {
  fres<-fres[order(-fres$confES),]
  message("Note: default gene set prioritisation is by confidence interval (confES). Alternatives are 'significance' and 'effect'.")
}
attributes(fres)$priority<-priority
return(fres)
}


manova_analysis_metrics_calc<-function(x, genesets, manova_result, minsetsize=10 ) {
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
# may implement some type of jitter in future
#  jitter_df<-function(x){
#    rand<-matrix(0.001*rnorm(ncol(x)*nrow(x), mean = 0, sd = 1),ncol=ncol(x))
#    x+rand
#  }
#  x<-jitter_df(x)

for ( i in 1:ncol(x)) {
  LEN=length(x[,i])
  UNIQLEN=length(unique(x[,i]))
  if ( UNIQLEN/LEN<0.99 ) { stop("Error: >99% of genes have the same score. More granular measurements needed for rank based enrichment analysis.") }
  if ( UNIQLEN/LEN<0.4 ) { warning("Warning: >60% of genes have the same score. This isn't optimal for rank based enrichment analysis.") }
}

  rank_adj<-function(x){
    xx<-rank(x)
    num_neg=length( which( x<0 ) )
    num_zero=length( which( x==0 ) )
    num_adj=num_neg+(num_zero/2)
    adj<-xx-num_adj
    adj
  }
  adj<-apply(x,2,rank_adj)
  adj
}


detailed_sets<-function(res,  resrows=50) {
  #collect ranked genelist of each genest
  ss<-res$ranked_profile
  mykeys <- as.character(res$manova_result[1:resrows,1])
  dat <- vector(mode="list", length=resrows)
  names(dat) <- mykeys

  for(i in 1:resrows) {
    sss<-ss[which(rownames(ss) %in% genesets[[which(names(genesets) %in% as.character(res$manova_result[i,1]))]]),]
    dat[[i]]<-sss
  }
  dat
}


endrich<-function(x,genesets, minsetsize=10, cores=detectCores()-1 , resrows=50, priority=NULL, bootstraps=0) {
	input_profile<-x

        input_genesets<-genesets

	ranked_profile<-endrichrank(input_profile)

	manova_result<-EnDrichMANOVA(ranked_profile, genesets, minsetsize=minsetsize, cores=cores, priority=priority, bootstraps=bootstraps)

	manova_analysis_metrics<-manova_analysis_metrics_calc(x,genesets,manova_result)

	dat <- list("input_profile" = input_profile,
		"input_genesets" = input_genesets,
		"ranked_profile" = ranked_profile,
		"manova_result" = manova_result,
		"manova_analysis_metrics" = manova_analysis_metrics)

        dat$detailed_sets<-detailed_sets(dat,resrows)

	dat
}


plotSets <- function(res,outfile="Rplots.pdf") {
  library("GGally")
  library("vioplot")
  library("gridExtra")
  palette <- colorRampPalette(c("white", "yellow","orange" ,"red","darkred","black"))

  resrows=length(res$detailed_sets)

  ss<-res$ranked_profile

  xmin=min(ss[,1])
  xmax=max(ss[,1])
  ymin=min(ss[,2])
  ymax=max(ss[,2])

  if ( ncol(ss)<3 ) {

    pdf(outfile)
    k<-MASS:::kde2d(ss[,1],ss[,2])
    X_AXIS=paste("Rank in contrast",colnames(ss)[1])
    Y_AXIS=paste("Rank in contrast",colnames(ss)[2])

    plot(res$input_profile , pch=19, col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),
      main="Scatterplot of all genes" )
    abline(v=0,h=0,lty=2,lwd=2,col="blue")

    filled.contour(k, xlim=c(xmin,xmax),ylim=c(ymin,ymax),
      color=palette , 
      plot.title={ abline(v=0,h=0,lty=2,lwd=2,col="blue")
        title( main="Rank-rank plot of all genes",xlab=X_AXIS,ylab=Y_AXIS )
      }
    )

    #barcghart of gene locations by quadrant
    uu=length(which(res$input_profile[,1]>0 & res$input_profile[,2]>0))
    ud=length(which(res$input_profile[,1]>0 & res$input_profile[,2]<0))
    dd=length(which(res$input_profile[,1]<0 & res$input_profile[,2]<0))
    du=length(which(res$input_profile[,1]<0 & res$input_profile[,2]>0))
    a<-as.data.frame(c(uu,ud,dd,du))
    rownames(a)=c("top-right","bottom-right","bottom-left","top-left")
    colnames(a)="a"
    barplot(a$a,names.arg=rownames(a),main="number of genes in each quadrant")

    #histograms of gene set counts
    geneset_counts<-res$manova_analysis_metrics$geneset_counts
    boxplot(geneset_counts$count,horizontal=T,frame=F,main="Gene set size",xlab="number of member genes included in profile")
    hist(geneset_counts$count,100,xlab="geneset size",main="Histogram of geneset size")
    hist(geneset_counts$count,100,xlim=c(0,500),xlab="geneset size",main="Trimmed histogram of geneset size")

    #barcghart of gene set locations by quadrant
    a<-res$manova_analysis_metrics[14]
    a<-as.data.frame(as.numeric(unlist(strsplit(as.character(a),','))),stringsAsFactors=F)
    rownames(a)=c("top-right","bottom-right","bottom-left","top-left")
    colnames(a)="a"
    barplot(a$a,names.arg=rownames(a),main="number of genesets FDR<0.05")

    sig<-subset(res$manova_result , p.adjustMANOVA<0.05)
    plot(res$manova_result[,4:5] , pch=19, col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),
      main="Scatterplot of all gene sets; FDR<0.05 in red" )
    abline(v=0,h=0,lty=2,lwd=2,col="blue")
    points(sig[,4:5] , pch=19, col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5))

    top<-head(res$manova_result ,resrows)
    plot(res$manova_result[,4:5] , pch=19, col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),
      main=paste("Scatterplot of all gene sets; top",resrows,"in red") )
    abline(v=0,h=0,lty=2,lwd=2,col="blue")
    points(top[,4:5] , pch=19, col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5))

    # plot effect size versus significance 
    plot(res$manova_result$s.dist,-log(res$manova_result$p.adjustMANOVA), 
      xlab="s.dist (effect size)",ylab="-log(p.adjustMANOVA) (significance)",
      pch=19, col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2), 
      main="effect size versus statistical significance")


   for(i in 1:resrows) {
      ll<-res$manova_result[i,]
      size<-ll$setSize
      sss<-res$detailed_sets[[i]] 

      k<-MASS:::kde2d(sss[,1],sss[,2])
      filled.contour( k, color = palette, xlim=c(xmin,xmax),ylim=c(ymin,ymax),
          plot.title={ abline(v=0,h=0,lty=2,lwd=2,col="blue")
          title( main=ll$set , xlab=X_AXIS,ylab=Y_AXIS  )
        }
      )

      plot(sss, pch=19, col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2),
        main=ll$set ,
        xlim=c(xmin,xmax),ylim=c(ymin,ymax),
        xlab=X_AXIS,ylab=Y_AXIS
      )
      abline(v=0,h=0,lty=2,lwd=2,col="blue")

      do.call(vioplot,c(unname(as.data.frame(sss)),col='gray',drawRect=T,names=list(names(as.data.frame(sss)))))
      grid()
      abline(h=0,lty=2,lwd=2,col="blue")
      title(main = ll[,1] , ylab = "Position in rank")

    }
    dev.off()
  } else {

  pdf(outfile)
  #pairs points plot for genes
  ggpairs_points_plot <- function(data ,mapping, ...){
    p <- ggplot(data = data, mapping = mapping) +
      geom_point(alpha=0.05) +
      geom_vline(xintercept=0,linetype="dashed") +
      geom_hline(yintercept=0,linetype="dashed")
  }

  p<-ggpairs(as.data.frame(x), title="Scatterplot of all genes" , lower  = list(continuous = ggpairs_points_plot ))
  print( p +  theme_bw() ) 

  #pairs contour plot function
  ggpairs_func <- function(data, mapping, ...){
    p <- ggplot(data = data, mapping = mapping) +
      stat_density2d(aes(fill=..density..), geom="tile", contour = FALSE) +
      geom_vline(xintercept=0,linetype="dashed") +
      geom_hline(yintercept=0,linetype="dashed") +
      scale_fill_gradientn(colours=palette(25))
    p
  }

  #pairs contour plot
  p<-ggpairs(as.data.frame(ss), title="Contour plot of all genes after ranking" , lower=list(continuous=ggpairs_func),
    diag=list(continuous=wrap("barDiag", binwidth=nrow(ss)/100)))
  print( p + theme_bw() )

  #subset contour plot
  ggpairs_contour_limit_range <- function(data ,mapping, ...){
    p <- ggplot(data = data, mapping = mapping) +
      stat_density2d(aes(fill=..density..), geom="tile", contour = FALSE) +
      geom_vline(xintercept=0,linetype="dashed") +
      geom_hline(yintercept=0,linetype="dashed") +
      scale_fill_gradientn(colours=palette(25)) +
      scale_x_continuous( limits = range(min(ss[,gsub("~","",as.character(mapping[1]))]),max(ss[,gsub("~","",as.character(mapping[1]))])) ) +
      scale_y_continuous( limits = range(min(ss[,gsub("~","",as.character(mapping[2]))]),max(ss[,gsub("~","",as.character(mapping[2]))])) )
    p
  }

  #a table of gene location by sector
  d=ncol(ss)
  sig<-sign(ss)
  sector_count<-aggregate(1:nrow(sig) ~ ., sig, FUN = length)
  colnames(sector_count)[ncol(sector_count)]<-"Number of genes in each sector"
  grid.table(sector_count)

  #histograms of gene set counts
  par(mfrow=c(3,1))
  geneset_counts<-res$manova_analysis_metrics$geneset_counts
  boxplot(geneset_counts$count,horizontal=T,frame=F,main="Gene set size",xlab="number of member genes included in profile")
  hist(geneset_counts$count,100,xlab="geneset size",main="Histogram of geneset size")
  hist(geneset_counts$count,100,xlim=c(0,500),xlab="geneset size",main="Trimmed histogram of geneset size")

  #a table of geneset location by sector
  sig<-sign(res$manova_result[which(res$manova_result$p.adjustMANOVA<0.05),4:(4+d-1)])
  sector_count<-aggregate(1:nrow(sig) ~ ., sig, FUN = length)
  colnames(sector_count)[ncol(sector_count)]<-"Number of gene sets in each sector"
  grid.table(sector_count)








  DIMS=ncol(ss)
  #pairs points plot for gene sets
  p<-ggpairs(res$manova_result[,4:(3+DIMS)] , title="Scatterplot of all genessets; FDR<0.05 in red" , lower  = list(continuous = ggpairs_points_plot ))
  print( p +  theme_bw() )

  #subset points plot function
  ggpairs_points_limit_range <- function(data ,mapping, ...){
    p <- ggplot(data = data, mapping = mapping) +
      geom_point(alpha=0.1) +
      geom_vline(xintercept=0,linetype="dashed") +
      geom_hline(yintercept=0,linetype="dashed") +
      scale_x_continuous( limits = range(min(ss[,gsub("~","",as.character(mapping[1]))]),max(ss[,gsub("~","",as.character(mapping[1]))])) ) +
      scale_y_continuous( limits = range(min(ss[,gsub("~","",as.character(mapping[2]))]),max(ss[,gsub("~","",as.character(mapping[2]))])) )
    p
  }

  # plot effect size versus significance 
  plot(res$manova_result$s.dist,-log(res$manova_result$p.adjustMANOVA), 
    xlab="s.dist (effect size)",ylab="-log(p.adjustMANOVA) (significance)",
    pch=19, col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2), 
    main="effect size versus statistical significance")


  for(i in 1:resrows) {
    ll<-res$manova_result[i,]
    size<-ll$setSize
    sss<-res$detailed_sets[[i]]

    p<-ggpairs(as.data.frame(sss), title=ll[,1], lower=list(continuous=ggpairs_contour_limit_range),
      diag=list(continuous=wrap("barDiag", binwidth=nrow(ss)/10)) )
    print( p + theme_bw() )

    p<-ggpairs(as.data.frame(sss), title=ll[,1], lower= list(continuous = ggpairs_points_limit_range ),
          diag=list(continuous=wrap("barDiag", binwidth=nrow(ss)/10)))
    print( p + theme_bw() )

    do.call(vioplot,c(unname(as.data.frame(sss)),col='gray',drawRect=T,names=list(names(as.data.frame(sss)))))
    grid()
    abline(h=0,lty=2,lwd=2)
    title(main = ll[,1] , ylab = "Position in rank")

  }
  dev.off()

  }

}


RankRankBinPlot<-function(x, binsize=500) {
	library("ggplot2")
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
  library("knitr")
  library("markdown")
  library("rmarkdown")

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
