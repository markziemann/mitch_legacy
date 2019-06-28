#' mitch_import
#'
#' This function imports differential omics data from common differential tools like edgeR, limma, DESeq2, Sleuth, 
#' ABSSeq and Seurat. It calculates a summarised differential expression metric by multiplying the sign of the log
#' fold change by the -log10 of the p-value. If this behaviour is not desired, mitch_import can be bypassed in 
#' favour of another scoring metric.
#' @param x a list of differential expression tables
#' @param DEtype the program that generated the differential expression table
#' @param geneIDcol the column containing gene names. If gene names are specified as row names, then geneIDcol=NULL.
#' @param geneTable a 2 column table mapping gene identifiers in the profile to gene identifiers in the gene sets. 
#' @return a multi-column table compatible with mitch_calc analysis.
#' @keywords import mitch
#' @export
#' @examples
#' # first step is to create a list of differential 
#' w<-list("edger1"=edger1,"edger2"=edger2)
#' # import as edgeR table with gene accessions in column named "GeneAccession" and "gt" mapping gene names.
#' x<-mitch_import(w,DEtype="edger",geneIDcol="GeneAccession",geneTable=gt)
mitch_import<-function(x , DEtype, geneIDcol=NULL, geneTable=NULL ) {

library("plyr")

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
    z<-aggregate(. ~ GeneSymbol,z,function(x) sum(as.numeric(as.character(x))))
#    z<-aggregate(. ~ GeneSymbol,z,sum)
    z$geneidentifiers=NULL
    colnames(z)=c("geneidentifiers","y")
  }
  z
}

#the geneIDcol should be an attribute added to each list item
for (i in 1:length(x) ) {
  if ( !is.null(geneIDcol) ) {
    LEN=length( which(names(x[[i]]) %in% geneIDcol) )
    if (LEN<1) { stop("Error: the specified geneIDcol doesn't seem to exist") }
    if (LEN>1) { stop("Error: there are multiple matches for the  specified geneIDcol") }
    attributes(x[[i]])$geneIDcol<-which( names(x[[i]]) %in% geneIDcol )
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

  s<-sign(y$logFC)*-log10(y$PValue)

  if ( !is.null(attributes(y)$geneIDcol) ) {
    g<-y[,attributes(y)$geneIDcol]
  } else {
    g<-rownames(y)
  }
  z<-data.frame(g,s,stringsAsFactors=F)
  colnames(z)<-c("geneidentifiers","y")
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

  s<-sign(y$log2FoldChange)*-log10(y$pvalue)

  if ( !is.null(attributes(y)$geneIDcol) ) {
    g<-y[,attributes(y)$geneIDcol]
  } else {
    g<-rownames(y)
  }
  z<-data.frame(g,s,stringsAsFactors=F)
  colnames(z)<-c("geneidentifiers","y")
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

  s<-sign(y$logFC)*-log10(y$P.Value)

  if ( !is.null(attributes(y)$geneIDcol) ) {
    g<-y[,attributes(y)$geneIDcol]
  } else {
    g<-rownames(y)
  }
  z<-data.frame(g,s,stringsAsFactors=F)
  colnames(z)<-c("geneidentifiers","y")
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

  s<-sign(y$foldChange)*-log10(y$pvalue)

  if ( !is.null(attributes(y)$geneIDcol) ) {
    g<-y[,attributes(y)$geneIDcol]
  } else {
    g<-rownames(y)
  }
  z<-data.frame(g,s,stringsAsFactors=F)
  colnames(z)<-c("geneidentifiers","y")
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

  s<-sign(y$b)*-log10(y$pval)

  if ( !is.null(attributes(y)$geneIDcol) ) {
    g<-y[,attributes(y)$geneIDcol]
  } else {
    g<-rownames(y)
  }
  z<-data.frame(g,s,stringsAsFactors=F)
  colnames(z)<-c("geneidentifiers","y")
  z<-mapGeneIds(y,z)
  z
}

topconfect_score<-function(y) {

  FCCOL=length(which(names(y)=="confect"))
  if (FCCOL>1){ stop("Error, there is more than 1 column named 'confect' in the input") }
  if (FCCOL<1){ stop("Error, there is no column named 'confect' in the input") }

  s<-y$confect
  s[is.na(s)] <- 0

  if ( !is.null(attributes(y)$geneIDcol) ) {
    g<-y[,attributes(y)$geneIDcol]
  } else {
    g<-rownames(y)
  }
  z<-data.frame(g,s,stringsAsFactors=F)
  colnames(z)<-c("geneidentifiers","y")
  z<-mapGeneIds(y,z)
  z
}

seurat_score<-function(y) {

  NCOL=ncol(y)
  if (NCOL<2){ stop("Error: there are <2 columns in the input, 'PValue' and 'logFC' are required ") }

  PCOL=length(which(names(y)=="p_val"))
  if (PCOL>1){ stop("Error, there is more than 1 column named 'PValue' in the input") }
  if (PCOL<1){ stop("Error, there is no column named 'PValue' in the input") }

  FCCOL=length(which(names(y)=="avg_logFC"))
  if (FCCOL>1){ stop("Error, there is more than 1 column named 'logFC' in the input") }
  if (FCCOL<1){ stop("Error, there is no column named 'logFC' in the input") }

  s<-sign(y$avg_logFC)*-log10(y$p_val)

  if ( !is.null(attributes(y)$geneIDcol) ) {
    g<-y[,attributes(y)$geneIDcol]
  } else {
    g<-rownames(y)
  }
  z<-data.frame(g,s,stringsAsFactors=F)
  colnames(z)<-c("geneidentifiers","y")
  z<-mapGeneIds(y,z)

  z$y[is.infinite(z$y) & z$y < 0] <- min(z$y[!is.infinite(z$y)])-.01
  z$y[is.infinite(z$y) & z$y > 0] <- max(z$y[!is.infinite(z$y)])+.01
  z<-z[,c(2,1)]
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
} else if ( DEtype == "seurat" ) {
  xx<-lapply(x,seurat_score)
} else if ( DEtype == "topconfects" ) {
  xx<-lapply(x,topconfect_score)
}


# give the colums a unique name otherwise join_all will fail
for (i in 1:length(xx) )  {
  colnames(xx[[i]])<-c("geneidentifiers" , paste("y",i,sep=""))
}

if ( DEtype == "seurat" ) {
  xxx<-join_all(xx,by = "geneidentifiers",type="full")
} else {
  xxx<-join_all(xx,by = 'geneidentifiers', type = 'inner')
}
rownames(xxx)<-xxx$geneidentifiers
xxx$geneidentifiers=NULL
colnames(xxx)<-names(x)

STARTSWITHNUM=length(grep("^[0-9]",colnames(xxx)))
if (STARTSWITHNUM>0) {
  stop("Error: it looks like one or more column names starts with a number. This is incompatible with downstream analysis. Please modify")
}



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

#' gmt_import
#'
#' This function imports GMT files into a list of character vectors for mitch analysis. GMT files are a commonly used
#' format for lists of genes used in pathway enrichment analysis. GMT files can be obtained from Reactome, MSigDB, etc.
#' @param gmtfile a gmt file.
#' @return a list of gene sets.
#' @keywords import genesets
#' @export
#' @examples
#' # Import some gene sets
#' genesets<-gmt_import("MyGeneSets.gmt")

gmt_import<-function(gmtfile){
    genesetLines <- strsplit(readLines(gmtfile), "\t")
    genesets <- lapply(genesetLines, utils::tail, -2)
    names(genesets) <- sapply(genesetLines, head, 1)
    attributes(genesets)$originfile<-gmtfile
    genesets
}

#' MANOVA
#'
#' This function performs the multivariate analysis of variance for each set of genes. It wraps around the existing
#' manova function in R. This function is not meant to be used directly.
#' @param x a multicolumn numerical table with each column containing differential expression scores for a contrast.
#' Rownames must match genesets.
#' @param genesets lists of genes imported by the gmt_imprt function or similar.
#' @param minsetsize the minimum number of genes required in a set for it to be included in the statistical analysis.
#' @param cores the number of parallel threads for computation. Defaults to the number of cores present minus 1.
#' @param priority the prioritisation metric to selecting top gene sets. Valid options are "significance", 
#' "effect" and "confidence"
#' @param bootstraps the number of bootstraps to recalculate in estimate the confidence intervals.
#' @return a mitch results object 
#' @keywords mitch MANOVA
#' @export
#' @examples
#' #This function is not designed to be used directly
MANOVA<-function(x,genesets, minsetsize=10, cores=detectCores()-1, priority=NULL, bootstraps=0) {

STARTSWITHNUM=length(grep("^[0-9]",colnames(x)))
if (STARTSWITHNUM>0) {
  stop("Error: it looks like one or more column names starts with a number. This is incompatible with downstream analysis. Please modify")
}

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
      while ( STRAPSDONE<bootstraps) {
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

if (nrow(fres)<1) {

  message("Warning: No results found. Check that the gene names in the profile match the gene sets and consider loosening the minsetsize parameter.")

  } else {
  fres$p.adjustMANOVA<-p.adjust(fres$pMANOVA,"fdr")

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
}

#' mitch_metrics_calc
#'
#' This function collects some metrics from the manova analysis. This function is not meant to be used directly.
#' @param x a multicolumn numerical table with each column containing differential expression scores for a contrast.
#' Rownames must match genesets.
#' @param genesets lists of genes imported by the gmt_imprt function or similar.
#' @param manova_result a valid result of the MANOVA function
#' @param minsetsize the minimum number of genes required in a set for it to be included in the statistical analysis.
#' @return a list of metrics
#' @keywords mitch metrics
#' @export
#' @examples
#' #This function is not designed to be used directly
mitch_metrics_calc<-function(x, genesets, manova_result, minsetsize=10 ) {

if (!is.null(manova_result)){

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
}


#' mitch_rank
#'
#' This function performs zero-centred ranking of differential contrast omics data. This function is not meant to 
#' be used directly.
#' @param x a multicolumn numerical table with each column containing differential expression scores for a contrast.
#' @return a ranked table differential expression data.
#' @keywords mitch rank ranking
#' @export
#' @examples
#' #This function is not designed to be used directly
mitch_rank<-function(x) {

for ( i in 1:ncol(x)) {
  LEN=length(x[,i])
  UNIQLEN=length(unique(x[,i]))
#  if ( UNIQLEN/LEN<0.1 ) { stop("Error: >90% of genes have the same score. More granular measurements needed for rank based enrichment analysis.") }
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

#' detailed_sets
#'
#' This function fetches differential expression scores for members of sets of genes that were top ranked. This 
#' function is not meant to be used directly.
#' @param res a mitch results object
#' @param resrows an integer representing the number of top genesets for which a detailed report is to be generated.
#' Default is 50.
#' @return mitch res object with dataframes of the high priority gene sets attached.
#' @keywords mitch detailed sets
#' @export
#' @examples
#' #This function is not designed to be used directly
detailed_sets<-function(res,  resrows=50) {
  #collect ranked genelist of each genest
  genesets<-res$input_genesets
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

#' mitch_calc
#'
#' This function performs multivariate gene set enrichment analysis. 

#' @param x a multicolumn numerical table with each column containing differential expression scores for a contrast.
#' Rownames must match genesets.
#' @param genesets lists of genes imported by the gmt_imprt function or similar.
#' @param minsetsize the minimum number of genes required in a set for it to be included in the statistical analysis.
#' Default is 10.
#' @param cores the number of parallel threads for computation. Defaults to the number of cores present minus 1.
#' @param resrows an integer value representing the number of top genesets for which a detailed report is to be 
#' generated. Default is 50.
#' @param priority the prioritisation metric used to selecting top gene sets. Valid options are "significance", 
#' "effect" and "confidence". If using "confidence", then bootstraps must be set >0.
#' @param bootstraps the number of bootstraps to recalculate in estimate the confidence intervals.
#' @return mitch res object with the following parts:
#'  $input_profile: the supplied input differential profile
#'  $input_genesets: the supplied input gene sets
#'  $ranked_profile: the differential profile after ranking
#'  $manova_result: the table of MANOVA enrichment results for each gene set
#'  $manova_analysis_metrics:  several metrics that are important to the interpretation of the results
#'  $detailed_sets: a list of dataframes containing ranks of members of prioritised gene sets.
#' @keywords mitch calc calculate manova 
#' @export
#' @examples
#' # An example of using mitch to calculate multivariate enrichments and prioritise based on confidence intervals 
#' # generated from 100 bootstraps.
#' res<-mitch_calc(x,genesets,resrows=25,bootstraps=100,priority="confidence")
mitch_calc<-function(x,genesets, minsetsize=10, cores=detectCores()-1 , resrows=50, priority=NULL, bootstraps=0) {
library("plyr")
library("parallel")
library("pbmcapply")
library("Rmisc")

        colnames(x)<-sub("-","_",colnames(x))
	input_profile<-x

        input_genesets<-genesets

	ranked_profile<-mitch_rank(input_profile)

	manova_result<-MANOVA(ranked_profile, genesets, minsetsize=minsetsize, cores=cores, priority=priority, bootstraps=bootstraps)

	if (!is.null(manova_result)) {
		mitch_metrics <-mitch_metrics_calc(x,genesets,manova_result)

		dat <- list("input_profile" = input_profile,
			"input_genesets" = input_genesets,
			"ranked_profile" = ranked_profile,
			"manova_result" = manova_result,
			"manova_analysis_metrics" = mitch_metrics)

	        if ( nrow(manova_result) < resrows ) { resrows<-nrow(manova_result) }

	        dat$detailed_sets<-detailed_sets(dat,resrows)

                attr(dat, 'profile_dimensions') <- colnames(dat$input_profile)
	
		dat

	}

}

#' mitch_plots
#'
#' This function generates several plots of multivariate gene set enrichment in high resolution PDF format.
#' The number of detailed sets to generate is dictated by the resrows set in the mitch_calc command.
#' @param res a mitch results object.
#' @param outfile the destination file for the plots in PDF format. should contain "pdf" suffix. Defaults to 
#' "Rplots.pdf"
#' @return generates a PDF file containing enrichment plots.
#' @keywords mitch plot plots pdf 
#' @export
#' @examples
#' # render enrichment plots in high res pdf
#' mitch_plots(res,outfile="outres.pdf")
mitch_plots <- function(res,outfile="Rplots.pdf") {
  library("gplots")
  library("plyr")
  library("reshape2")
  library("GGally")
  library("ggplot2")
#  library("vioplot")
  library("grid")
  library("gridExtra")
  palette <- colorRampPalette(c("white", "yellow","orange" ,"red","darkred","black"))

  mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 0.5)),
    colhead = list(fg_params=list(cex = 0.7)),
    rowhead = list(fg_params=list(cex = 0.7)))

  resrows=length(res$detailed_sets)

  ss<-res$ranked_profile

  xmin=min(ss[,1])
  xmax=max(ss[,1])
  ymin=min(ss[,2])
  ymax=max(ss[,2])

  # d is the number of dimensions - this is very important and refered to frequently
  d=ncol(ss)

  if ( d<3 ) {

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

    #barchart of gene locations by quadrant
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

    #barchart of gene set locations by quadrant
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

    # A heatmap of s values for the resrows sets
    hmapx<-head( res$manova_result[,4:(4+d-1)] ,resrows)
    rownames(hmapx)<-head(res$manova_result$set,resrows)
    colnames(hmapx)<-gsub("^s.","",colnames(hmapx))
    my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
    heatmap.2(as.matrix(hmapx),scale="none",margin=c(10, 25),cexRow=0.8,trace="none",cexCol=0.8,col=my_palette)

    # plot effect size versus significance 
    plot(res$manova_result$s.dist,-log(res$manova_result$p.adjustMANOVA), 
      xlab="s.dist (effect size)",ylab="-log(p.adjustMANOVA) (significance)",
      pch=19, col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2), 
      main="effect size versus statistical significance")

  ss_long<-melt(ss)

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

  sss_long<-melt(sss)

  p<-ggplot(ss_long,aes(Var2,value)) + 
 geom_violin(data=ss_long,fill = "grey", colour = "grey") +
    geom_boxplot(data=ss_long,width=0.9,fill="grey",outlier.shape = NA,coef = 0) +
    geom_violin(data=sss_long,fill = "black", colour = "black") +
    geom_boxplot(data=sss_long,width=0.1,outlier.shape = NA) +
    labs(y = "Position in rank",title = ll[,1] )

  print(
    p + 
    theme_bw() + 
    theme( axis.text=element_text(size=14),
    axis.title=element_text(size=15),
    plot.title = element_text(size = 20)))

#      do.call(vioplot,c(unname(as.data.frame(sss)),col='gray',drawRect=T,names=list(names(as.data.frame(sss)))))
#      grid()
#      abline(h=0,lty=2,lwd=2,col="blue")
#      title(main = ll[,1] , ylab = "Position in rank")

    }
    dev.off()
  } else {

  pdf(outfile)

  # if working with >5 dimensions, then substitute the dimension (colnames) names with a number
  if ( d>5 ) {
    mydims<-data.frame( attributes(res)$profile_dimensions )
    colnames(mydims)<-"dimensions"
    grid.newpage()
    grid.table(mydims,theme=mytheme)
    colnames(res$input_profile)<- paste("d",1:ncol(res$input_profile),sep="")
    colnames(res$ranked_profile)<- paste("d",1:ncol(res$ranked_profile),sep="")
    ss<-res$ranked_profile
  }

  #pairs points plot for genes
  ggpairs_points_plot <- function(data ,mapping, ...){
    p <- ggplot(data = data, mapping = mapping) +
      geom_point(alpha=0.05) +
      geom_vline(xintercept=0,linetype="dashed") +
      geom_hline(yintercept=0,linetype="dashed")
  }

  p<-ggpairs(as.data.frame(res$input_profile), title="Scatterplot of all genes" , lower  = list(continuous = ggpairs_points_plot ))
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
  sig<-sign(ss)
  sector_count<-aggregate(1:nrow(sig) ~ ., sig, FUN = length)
  colnames(sector_count)[ncol(sector_count)]<-"Number of genes in each sector"
  grid.newpage()
  grid.table(sector_count,theme=mytheme)

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
  grid.newpage()
  grid.table(sector_count,theme=mytheme)

  #pairs points plot for gene sets
  manova_result_clipped<-res$manova_result[,4:(3+d)] 
  colnames(manova_result_clipped)<-colnames(res$input_profile)
  p<-ggpairs(manova_result_clipped , title="Scatterplot of all genessets; FDR<0.05 in red" , lower  = list(continuous = ggpairs_points_plot ))
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

  # A heatmap of s values for the resrows sets
  hmapx<-head( res$manova_result[,4:(4+d-1)] ,resrows)
  rownames(hmapx)<-head(res$manova_result$set,resrows)
  colnames(hmapx)<-gsub("^s.","",colnames(hmapx))
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
  heatmap.2(as.matrix(hmapx),scale="none",margin=c(10, 25),cexRow=0.8,trace="none",cexCol=0.8,col=my_palette)

  # plot effect size versus significance 
  par(mfrow=c(1,1))
  plot(res$manova_result$s.dist,-log(res$manova_result$p.adjustMANOVA), 
    xlab="s.dist (effect size)",ylab="-log(p.adjustMANOVA) (significance)",
    pch=19, col=rgb(red = 0, green = 0, blue = 0, alpha = 0.2), 
    main="effect size versus statistical significance")

  ss_long<-melt(ss)

  for(i in 1:resrows) {
    ll<-res$manova_result[i,]
    size<-ll$setSize
    sss<-res$detailed_sets[[i]]

    if ( d>5 ) {
      colnames(sss)<- paste("d",1:ncol(res$input_profile),sep="")
    }

    p<-ggpairs(as.data.frame(sss), title=ll[,1], lower=list(continuous=ggpairs_contour_limit_range),
      diag=list(continuous=wrap("barDiag", binwidth=nrow(ss)/10 )) )
    print( p + theme_bw() )

    p<-ggpairs(as.data.frame(sss), title=ll[,1], lower= list(continuous = ggpairs_points_limit_range ),
          diag=list(continuous=wrap("barDiag", binwidth=nrow(ss)/10)))
    print( p + theme_bw() )

  sss_long<-melt(sss)

  p<-ggplot(ss_long,aes(Var2,value)) +
    geom_violin(data=ss_long,fill = "grey", colour = "grey") +
    geom_boxplot(data=ss_long,width=0.9,fill="grey",outlier.shape = NA,coef = 0) +
    geom_violin(data=sss_long,fill = "black", colour = "black") +
    geom_boxplot(data=sss_long,width=0.1,outlier.shape = NA) +
    labs(y = "Position in rank",title = ll[,1] )

  print(
    p +
    theme_bw() +
    theme( axis.text=element_text(size=14),
    axis.title=element_text(size=15),
    plot.title = element_text(size = 20)))

#    do.call(vioplot,c(unname(as.data.frame(sss)),col='gray',drawRect=T,names=list(names(as.data.frame(sss)))))
#    grid()
#    abline(h=0,lty=2,lwd=2)
#    title(main = ll[,1] , ylab = "Position in rank")

  }
  dev.off()

  }

}


#' mitch_report
#'
#' This function generates an R markdown based html report containing tables and several plots of mitch results 
#' The plots are in png format, so are not as high in resolution as compared to the PDF generated by mitch_plots 
#' function. The number of detailed sets to generate is dictated by the resrows set in the mitch_calc command.
#' @param res a mitch results object.
#' @param outfile the destination file for the html report. should contain "html" suffix. Defaults to 
#' "Rplots.pdf"
#' @return generates a HTML file containing enrichment plots.
#' @keywords mitch report html markdown knitr
#' @export
#' @examples
#' # render mitch results in the form of a HTML report
#' mitch_report(res,"outres.html")
mitch_report<-function(res,out) {
  library("plyr")
  library("knitr")
  library("markdown")
  library("rmarkdown")

  HTMLNAME<-paste(out,".html",sep="")
  HTMLNAME<-gsub(".html.html",".html",HTMLNAME)
  HTMLNAME<-paste(getwd(),HTMLNAME,sep="/")
  rmd_tmpdir<-tempdir()
  rmd_tmpfile<-paste(rmd_tmpdir,"/mitch.Rmd",sep="")
  html_tmp<-paste(paste(rmd_tmpdir,"/mitch_report.html",sep=""))
  download.file("https://raw.githubusercontent.com/markziemann/Mitch/master/mitch.Rmd",destfile=rmd_tmpfile)

  DATANAME<-gsub(".html$",".RData",out)
  DATANAME<-paste(rmd_tmpdir,"/",DATANAME,sep="")
  save.image(DATANAME)
  MYMESSAGE=paste("Dataset saved as \"",DATANAME,"\".")
  message(MYMESSAGE)

  knitrenv <- new.env()
  assign("DATANAME", DATANAME, knitrenv)
  assign("res",res,knitrenv)

  download.file("https://raw.githubusercontent.com/markziemann/Mitch/master/mitch.Rmd",destfile=rmd_tmpfile)
  rmarkdown::render(rmd_tmpfile,output_file=html_tmp)
  file.copy(html_tmp,HTMLNAME)
}
