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


custAnnoSum<-function(x, order=c("TSS","Gene_body","TTS","Intergenic"), TTSdist=3000){
  require(dplyr, quietly=T)
  #Gene body = Exons and Introns
  #Intergenic = intergenic
  #TTS = 3'-UTR and 3kb downstream
  #TSS = 5'-UTR and 3kb upstream (run annotatePeak with tssRegion=c(-3000,0))
  x<-as.data.frame(x@anno)
  #to_remove<-grep("intron", x$annotation, ignore.case=T)
  #to_remove<-c(to_remove,grep("intergenic", x$annotaion, ignore.case=T))
  #message("removing ",length(to_remove)," peaks not in annotation regions")
  #x<-x[-to_remove,]

  x$stat<-sapply(strsplit(x$annotation, split="(", T),'[[',1)
  x$stat[grep("exon", x$annotation, ignore.case=T)]<-"Gene_body"
  x$stat[grep("intron",x$annotation, ignore.case=T)]<-"Gene_body"
  x$stat[grep("promoter", x$annotation, ignore.case=T)]<-"TSS"
  x$stat[grep("intergenic", x$annotation, ignore.case=T)]<-"Intergenic"
  x$stat[grep("5'",x$annotation, ignore.case=T)]<-"TSS"
  x$stat[grep("3'",x$annotation, ignore.case=T)]<-"TTS"
  scdist<-x$geneEnd-x$end #Get distance from end of gene to see if peak is within 3kp of stop codon
  scdist2<-x$geneEnd-x$start
  #Negative scdist means peak starts/ends after geneEnd - these are TTS
  x$stat[which(scdist>= -TTSdist & scdist <= 0)]<-"TTS"
  x$stat[which(scdist2>= -TTSdist & scdist2 <=0 )]<-"TTS"
  y<-as.data.frame(table(x$stat))
  colnames(y)<-c("Feature", "Total")
  y$Frequency<-y$Total/sum(y$Total)*100
  y<-y[match(order, y$Feature),]
  y$Feature<-as.factor(y$Feature)
  y$Feature<-forcats::fct_relevel(y$Feature, order)
  #y<-y[complete.cases(y),]
  return(y)
}

# custAnnoStat<-list(Control=custAnnoSum(C_anno),shRNF8=custAnnoSum(R_anno), Random=custAnnoSum(rand_anno))
#
# #Paper colours
# cols<-c(rgb(228,108,10,255, maxColorValue=255),rgb(139,163,89,255, maxColorValue=255),
#         rgb(196,189,151,255, maxColorValue=255),rgb(235,63,8,255, maxColorValue=255),
#         rgb(96,74,123,255, maxColorValue=255))
#
# test<-plyr::ldply(custAnnoStat, data.frame)
# test<-test[,c("Feature","Frequency","Total",".id")]
#
# custPlotAnnoBar(C_stat, title="Updated Annotations", categoryColumn = ".id", col=rev(cols))


custPlotAnnoBar<-function (custAnnoList, xlab = "", ylab = "Percentage(%)", title = "Feature Distribution",
                           col="Dark2")
{
  if(is.null(col)){
    col<-c(rgb(228,108,10,255, maxColorValue=255),rgb(139,163,89,255, maxColorValue=255),
                          rgb(196,189,151,255, maxColorValue=255),rgb(235,63,8,255, maxColorValue=255),
                          rgb(96,74,123,255, maxColorValue=255))
  }
  anno.df<-plyr::ldply(custAnnoList, data.frame)
  #print(anno.df)
  anno.df<-anno.df[,c("Feature","Frequency","Total",".id")]
  categoryColumn=".id"
  anno.df$Feature <- factor(anno.df$Feature, levels = rev(levels(anno.df$Feature)))
  p <- ggplot2::ggplot(anno.df, ggplot2::aes_string(x = categoryColumn, fill = "Feature",
                                                    y = "Frequency"))
  p <- p + ggplot2::geom_bar(stat = "identity") + ggplot2::coord_flip() + ggplot2::theme_bw()
  p <- p + ggplot2::ylab(ylab) + ggplot2::xlab(xlab) + ggplot2::ggtitle(title)
  if (categoryColumn == 1) {
    p <- p + ggplot2::scale_x_continuous(breaks = NULL)
    if(is.null(col)){
      p <- p + ggplot2::scale_fill_manual(values = rev(getCols(nrow(anno.df))),
                                          guide = ggplot2::guide_legend(reverse = TRUE))
    } else {
      p<-p+ggplot2::scale_fil_manual(values=BinfTools::colPal(col)[c(1:nrow(anno.df))],
                                     guide= ggplot2::guide_legend(reverse=T))
    }
  }
  else {
    if(is.null(col)){
      p <- p + ggplot2::scale_fill_manual(values = rev(getCols(length(unique(anno.df$Feature)))),
                                          guide = ggplot2::guide_legend(reverse = TRUE))
    } else {
      p<-p+ggplot2::scale_fill_manual(values=BinfTools::colPal(col)[c(1:length(unique(anno.df$Feature)))],
                                      guide=ggplot2::guide_legend(reverse=F))
    }
  }
  return(p)
}
