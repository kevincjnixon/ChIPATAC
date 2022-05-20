
# readBed<-function(x){
#   bed<-read.delim(x, header=F)
#   colnames(bed)[c(1:3)]<-c("seqnames","start","end")
#   bed<-GenomicRanges::makeGRangesFromDataFrame(bed, keep.extra.columns = T)
#   return(bed)
# }

#' Function to export GenomicRanges object for use with Cistrome BETA
#'
#' @param x GenomicRanges object to export
#' @param peakID String or String vector to use to identify peak names (will be concatenated to the row number of each peak)
#' @param filename String of path to export results (Should end in .bed)
#' @param scorCol String indicating column name in GenomicRanges object that should be used for scoring peaks in BETA analysis. Default="V7"
#' @export

writeBETA<-function(x, peakID, filename, scoreCol="V7"){
  x<-as.data.frame(x)
  y<-data.frame(chrom=x$seqnames, start=x$start, end=x$end,
                name=paste(peakID,1:nrow(x),sep="_"), score=x[,which(colnames(x) %in% scoreCol)])
  write.table(y, file=filename, quote=F, col.names=F, row.names=F, sep="\t")
}

#'Function to export GenomicRanges to GTF format
#'
#'@param x GenomicRanges object to export
#'@param source String to identify the source of peaks in GTF file. Default="ATAC-peakset"
#'@param peakID String or string vector to identify unique peaks. Will be set as 'gene_id' in gtf file and concatenated with row number of each peak.
#'@param filename String of path to export gtf file (should end in .gtf)
#'@export

writeGTF<-function(x, source="ATAC-peakset", peakID, filename){
  x<-as.data.frame(x)
  y<-data.frame(chrom=x$seqnames, source=rep(source, nrow(x)), feature=rep("exon", nrow(x)),
                start=x$start, end=x$end, blank1=rep(".", nrow(x)), blank2=rep(".", nrow(x)),
                blank3=rep(".", nrow(x)), gene_id=paste0("gene_id \"",peakID,".",1:nrow(x),"\";"))
  if(length(peakID)>0){
    y$gene_id<-paste0("gene_id \"",peakID,"\";")
  }
  write.table(y, file=filename, quote=F, col.names=F, row.names=F, sep="\t")
}

#'Function to read tab-delimited bed file to GenomicRanges object
#'
#'@param filename String of path to bed file to read in
#'@param header Boolean indicating if there is a header in the bed file. Default=FALSE.
#'@return GenomicRanges object
#'@export

readBed<-function(filename, header=F){
  x<-read.delim(filename, header)
  colnames(x)[1:3]<-c("seqnames","start","end")
  x<-GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = T)
  return(x)
}
