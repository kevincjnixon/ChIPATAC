###Feature Enrichment
#Function to take a given bed file, annotate it
#And then bootstrap N number of random peaks using bedtools shuffle,
#Annotate those, and use values to look for significance and enrichment
#Of each genomic feature

#Calculate empirical p-value
#p=(number of times the number of random peaks in a given feature is > number of raw peaks in feature)/N
#I guess this assumes

#' Calculate empirical p-value for genomic feature enrichment
#'
#' @param bedFile String with path to bedfile to analyze. Must be a bed file
#' @param species string indicating genus/species of genome. One of "hsapiens" or "mmusculus"
#' @param N numeric indicating number of iterations
#' @param tss numeric vector of length 2 indicating the promoter region around the tss site (default= c(-3000,3000))
#' @param sys string indicating opertating system. - Not working at the moment, currently only works for windows running windows subsystem for linux with bedtools installed
#' @return Data frame with p-values and adjusted p-values (BH corrected) for each genomic feature
#' @export

FeatureEnrichment<-function(bedFile, species="hsapiens", N=1000, tss=c(-3000,3000), sys=c("win","")){
  require(TxDb.Hsapiens.UCSC.hg38.knownGene, quietly = T)
  require(TxDb.Mmusculus.UCSC.mm10.knownGene, quietly= T)
  require(ChIPseeker, quietly= T)
  require(org.Hs.eg.db, quietly= T)
  require(org.Mm.eg.db, quietly= T)
  require(data.table, quietly= T)
  TxDb<-NULL
  org<-NULL
  genome<-NULL
  if(species=="hsapiens"){
    TxDb<-TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    org<-"org.Hs.eg.db"
    genome<-"hg38.chrNameLength.txt"
  }
  if(species=="mmusculus"){
    TxDb<-TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    org<-"org.Mm.eg.db"
    genome<-"mm10.chrNameLength.txt"
  }
  #Annotate the original bedFile
  raw_Anno<-ChIPseeker::annotatePeak(bedFile, tssRegion = tss, TxDb = TxDb, annoDb=org)
  #Get the number of peaks annotated to each feature:
  raw_nums<-data.frame(NPeak=(raw_Anno@annoStat$Frequency/100)*nrow(as.data.frame(raw_Anno@anno)),
                       row.names=raw_Anno@annoStat$Feature)
  #Now to set up a data frame for calculating the empricial p-value
  p_frame<-matrix(data=0, nrow=N, ncol=nrow(raw_Anno@annoStat))
  colnames(p_frame)<-raw_Anno@annoStat$Feature
  features<-c("Promoter","5UTR","3UTR","Exon","Intron","Downstream","Intergenic")
  test_feat<-c("Promoter", "5'", "3'", "Exon", "Intron", "Downstream", "Intergenic")
  tmp<-as.character(raw_Anno@annoStat$Feature)
  to_remove<-c()
  for(i in 1:length(features)){
    if(length(grep(test_feat[i],tmp))==0){
      to_remove<-c(to_remove,i)
    }
  }
  if(length(to_remove)>0){
    features<-features[-to_remove]
    features<-c(features,test_feat[to_remove])
    message("Investigating only these features:")
    print(features)
  }
  #Now to loop and boostrap random peak files and annotating them
  pb<-txtProgressBar(min=0, max=N, style=3)
  for(i in 1:N){
    x<-GenomicRanges::makeGRangesFromDataFrame(data.table::fread(
      paste0("bash.exe -c 'bedtools shuffle -i ",bedFile," -g ",genome," -chrom'")),
      seqnames.field = "V1", start.field = "V2",end.field = "V3")
    #Annotate it
    #print(head(x))
    anno<-ChIPseeker::annotatePeak(x, tssRegion=tss, TxDb=TxDb, annoDb=org,
                                   genomicAnnotationPriority=features)
    #print(anno@annoStat)
    #Get the peak numbers
    nums<-data.frame(NPeak=(anno@annoStat$Frequency/100)*nrow(as.data.frame(anno@anno)),
                     row.names=anno@annoStat$Feature)
    #print(head(nums))
    #Now compare the numbers and add to p_frame
    for(k in 1:nrow(nums)){
      if(nums[k,1]>raw_nums[k,1]){
        p_frame[i,k]<-1
      }
    }
    setTxtProgressBar(pb, i)
  }
  #Now calculate the empirical p-value for each feature:
  p<-colSums(p_frame)/N
  padj<-p.adjust(p, method="BH")
  return(data.frame(p=p, padj=padj, row.names=colnames(p_frame)))
}
