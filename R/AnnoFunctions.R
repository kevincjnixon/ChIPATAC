#' Subset and extract genes from ChIPseeker annotations
#'
#' @param anno ChIPseeker annotation (csAnno) object
#' @param feature string indicating the genomic features to subset e.g. "promoter", "intergenic". Default is "all" and will not subset features.
#' @param unique Boolean indicating if only unique genes should be returned. Default is TRUE.
#' @param geneCol string indicating the column name in anno@anno containing the gene identifiers to be returned. Default is "SYMBOL" for gene symbols.
#' @return Character vector of gene IDs.
#' @export


getAnnoGenes<-function(anno, feature="all", unique=T, geneCol="SYMBOL"){
  anno<-as.data.frame(anno@anno)
  if(feature!="all"){
    anno<-anno[grep(feature, anno$annotation, ignore.case=T),]
  }
  genes<-anno[,which(colnames(anno) %in% geneCol)]
  if(isTRUE(unique)){
    genes<-unique(genes)
  }
  return(genes)
}

#' Return genomic regions associated with a specific gene or set of genes
#'
#' @param anno ChIPseeker annotation (csAnno) object
#' @param symbol String or vector of strings indicating the gene symbol(s) to search
#' @return data frame containing all information about all peaks associated with gene(s) of interest
#' @export

retRegion<-function(anno, symbol){
  anno<-as.data.frame(anno@anno)
  return(anno[which(anno$SYMBOL %in% symbol),])
}
