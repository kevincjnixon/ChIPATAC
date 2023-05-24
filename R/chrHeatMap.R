#'Average signal across the genome in bins for a list of GRanges objects
#'
#'Takes a named list of GRanges objects containing a score column to be averaged in bins across the genome
#'
#'@param GRList Named list of GRanges objects - each GRanges object should be from the same organism/genome build and should contain seqlengths
#'@param scoreCol Name of column in GRanges object containing the metric to be averaged across each bin
#'@param binSize Size of each genomic bin. Default is 100kb to help save memory for larger datasets
#'@return a list containing a data frame ('mat') of the binned results and a GenomicRanges ('GR') object containing the binned results.
#'@export

chrHeatMat<-function(GRList, scoreCol, binSize=100000){
  #Divide each chromosome into bins of binSize

  bins<-GenomicRanges::tileGenome(seqlengths(GRList[[1]]),
                                  tilewidth = binSize,
                                  cut.last.tile.in.chrom = T)
  #Now, we calculate the scores for each GR object
  score<-lapply(GRList, GenomicRanges::coverage, weight=scoreCol)
  #Now get the average score for each bin
  binned_data<-lapply(score, function(x){return(GenomicRanges::binnedAverage(bins, x, "average_score"))})
  #Now collapse all of the results into one GRanges object and data frame
  binned.df<-lapply(binned_data, as.data.frame)
  res.df<-binned.df[[1]]
  colnames(res.df)[ncol(res.df)]<-names(binned.df)[1]
  for(i in 2:length(binned.df)){
    res.df<-cbind(res.df, binned.df[[i]]$average_score)
    colnames(res.df)[ncol(res.df)]<-names(binned.df)[i]
  }
  res.gr<-GenomicRanges::makeGRangesFromDataFrame(res.df, keep.extra.columns=T)
  return(list(mat=res.df, GR=res.gr))
}

#'Plot a heatmap of binned signal across the whole genome
#'
#'@param chrHeatMat output of chrHeatMat() after binning. Not that smaller binSize may cause memory errors.
#'@param clusterSamples Boolean indicating if samples (columns) should be hierarchically clustered (TRUE) or remain in order (FALSE). Can be a data frame for grouping samples with row.names=names(GRList) and columns indicating annotations for groups.
#'@param showColNames Boolean indicating if sample names (column names) should be shown. Default is TRUE.
#'@param title character containing title for heatmap
#'@param scaleValue character containing the title for the heatmap scale- Default is 'score'
#'@return Prints a ComplexHeatmap with binned chromosomes as rows, and samples as columns
#'@export

chrHeatMap<-function(chrHeatMat, clusterSamples = TRUE, showColNames = TRUE, title = "", scaleValue = "score"){
  mat<-chrHeatMat$mat
  annodf<-data.frame(row.names=rownames(mat),
                     chr=mat$seqnames)
  anno_cols<-list(chr=rep(c("lightgrey","black"), length(unique(mat$seqnames)))[1:length(unique(mat$seqnames))])
  names(anno_cols$chr)<-as.character(unique(mat$seqnames))

  chr_anno<-ComplexHeatmap::rowAnnotation(name1=ComplexHeatmap::anno_empty(border=F, width=ComplexHeatmap::max_text_width(unique(annodf$chr))),
                                          name2=ComplexHeatmap::anno_empty(border=F, width=ComplexHeatmap::max_text_width(unique(annodf$chr))),
                                          chr=annodf$chr, show_annotation_name=F,col=anno_cols, show_legend=F)

  columnLabs<-colnames(mat[,c(6:ncol(mat))])
  if(isFALSE(showColNames)){
    columnLabs<-rep("", length(6:ncol(mat)))
  }
  t<-NA
  if(isTRUE(clusterSamples)){
    t<-ComplexHeatmap::Heatmap(as.matrix(mat[,c(6:ncol(mat))]), cluster_rows=F, row_labels = rep("",nrow(mat)),
                               left_annotation = chr_anno, name=scaleValue, row_title="Chromosome", row_title_rot=90,
                               split=annodf$chr, row_gap=grid::unit(0,"mm"), column_labels=columnLabs, column_title=title)
    print(t)
  }

  if(isFALSE(clusterSamples)){
    t<-ComplexHeatmap::Heatmap(as.matrix(mat[,c(6:ncol(mat))]), cluster_rows=F, row_labels = rep("",nrow(mat)),
                               left_annotation = chr_anno, name=scaleValue, row_title="Chromosome", row_title_rot=90,
                               split=annodf$chr, row_gap=grid::unit(0,"mm"), cluster_columns=F, column_labels=columnLabs,
                               column_title=title)
    print(t)
  }
  if(is.data.frame(clusterSamples)){
    colAnno<-ComplexHeatmap::HeatmapAnnotation(df=clusterSamples, show_annotation_name=F, show_legend=F)
    t<-ComplexHeatmap::Heatmap(as.matrix(mat[,c(6:ncol(mat))]), cluster_rows=F, row_labels=rep("",nrow(mat)),
                               left_annotation = chr_anno, name=scaleValue, row_title="Chromosome", row_title_rot=90,
                               split=annodf$chr, row_gap=grid::unit(0,"mm"), column_title=title, column_labels=columnLabs, top_annotation = colAnno, cluster_columns=F,
                               column_split=clusterSamples[,1])
    print(t)
    for(i in 1:length(unique(clusterSamples[,1]))){
      ComplexHeatmap::decorate_annotation(colnames(clusterSamples)[1], slice=i, {
        grid::grid.text(unique(clusterSamples[,1])[i], just="center",gp=grid::gpar(col="white"))
      })
    }
  }
  for(i in 1:length(unique(annodf$chr))){
    if(i%%2!=0){
      ComplexHeatmap::decorate_annotation("name1", slice=i,{
        grid::grid.text(unique(annodf$chr)[i], just="right")})
    } else {
      ComplexHeatmap::decorate_annotation("name2", slice=i,{
        grid::grid.text(unique(annodf$chr)[i], just="right")})
    }
  }

}
