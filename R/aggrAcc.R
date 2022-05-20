#' Aggregate peak counts to gene features
#'
#' Use a peak-count matrix (rows=peaks, columns=samples) and as csAnno object from ChIPseeker
#' to aggregate counts from all or specific features for each gene
#'
#' @param counts count matrix where rows=peaks and columns=individual samples
#' @param anno csAnno object (from ChIPseeker annotatePeak()) containing the annotations for only and all of the peaks in 'counts'
#' @param annoCol column name in anno@anno corresponding to the rownames of counts. Defaults to "V4".
#' @param feature Type of feature to aggregate. 'all" will be all features. Other options include "promoter","intron","intergenic", etc. Default is "all".
#' @param method method to use for aggregating counts. One of "mean","median","geoMean", "max", or "sum". Defaults to "mean"
#' @param raw Boolean indicating if exported object should be whole integers (i.e. raw counts). Default=FALSE.
#' @return Data.frame with rows=gene symbols, columns= samples, and aggregated gene counts for each gene in each sample.
#' @export

aggrAcc<-function(counts, anno, annoCol="V4", feature="all", method=c("mean","median","geoMean", "max", "sum"), raw=F, geneCol="SYMBOL"){
  colGeoMean<-function (a){
    x <- NULL
    for (i in 1:ncol(a)) {
      x <- c(x, prod(a[,i]^(1/length(a[,i]))))
    }
    return(x)
  }
  colMax<-function(a){
    x<-NULL
    for(i in 1:ncol(a)){
      #b<-abs(a[,i]) #Convert to absolute maximums
      #x<-c(x, a[which(b==max(b)),i]) #Take the corresponding value to the absolute max (that way, if it's negative, it will still be negative)
      x<-c(x, max(a[,i]))
    }
    return(x)
  }
  #Calculate an aggregate accessibility for each gene using normalized ATAC-seq counts
  anno<-as.data.frame(anno@anno)
  if(feature!="all"){
    anno<-anno[grep(feature, anno$annotation, ignore.case=T),]
  }
  genes<-unique(anno[,which(colnames(anno) %in% geneCol)])
  message(length(genes)," genes identified from dataset for ", feature," features.")
  x<-list()
  pb<-txtProgressBar(min=0, max=length(genes), style=3)
  for(i in 1:length(genes)){
    rows<-anno[,which(colnames(anno) %in% annoCol)][which(anno[,which(colnames(anno) %in% geneCol)] %in% genes[i])]
    sub<-counts[which(rownames(counts) %in% rows),]
    if(method[1]=="mean"){
      x[[i]]<-colMeans(sub)
    }
    if(method[1]=="median"){
      x[[i]]<-colMedians(sub)
    }
    if(method[1]=="geoMean"){
      x[[i]]<-colGeoMean(sub)
    }
    if(method[1]=="max"){
      x[[i]]<-colMax(sub)
    }
    if(method[1]=="sum"){
      if(!is.null(nrow(sub))){
        x[[i]]<-colSums(sub)
      } else {
        x[[i]]<-sub
      }
    }
    names(x)[i]<-genes[i]
    names(x[[i]])<-colnames(counts)
    setTxtProgressBar(pb, i)
  }

  res<-do.call("rbind", x)
  if(isTRUE(raw)){
    return(round(res))
  } else {
    return(res)
  }
}
