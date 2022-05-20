getGenomeCov<-function(x, TxDb, annoDb, promoter, tssRegion, gs=NULL){
  entrez<-NULL
  if(class(annoDb)[1]!="data.frame"){
    entrez<-suppressMessages(AnnotationDbi::select(annoDb, keys=x, keytype="SYMBOL", columns=c("ENTREZID"))$ENTREZID)
  } else {
    entrez<-annoDb$ENTREZID[which(annoDb$SYMBOL %in% x)]
  }
  if(length(entrez)<1){
    message("Gene symbols could not be mapped to Entrez IDs")
    return(NA)
  }
  g<-NULL
  if(class(TxDb)=="TxDb"){
    if(isFALSE(promoter)){
      g<-suppressMessages(GenomicFeatures::genes(TxDb))
    }
    if(isTRUE(promoter)){
      g<-suppressWarnings(GenomicFeatures::promoters(TxDb, columns="gene_id", upstream = abs(tssRegion[1]),
                                                     downstream=abs(tssRegion[2])))
    }
  }else{
    g<-TxDb
  }
  #message("g")
  g<-tryCatch({g[(S4Vectors::elementMetadata(g)[,"gene_id"] %in% entrez)]},
              error=function(e){
                message("Cannot subset S4. Trying alternative method...")
                tmp<-as.data.frame(S4Vectors::elementMetadata(g))
                g<-g[as.character(tmp[,"gene_id"]) %in% entrez]
                return(g)
              })
  #message("filteredg")
  #remove alternative chromosomes from 'g'
  if(isTRUE(promoter)){
    to_remove<-grep("_", as.data.frame(g, row.names=c(1:length(g)))$seqnames)
    if(length(to_remove)>0){
      g<-g[-to_remove]
    }
  }
  g2<-g
  if(isFALSE(promoter)){
    g2<-GenomicRanges::resize(g, width(g) + 5000, fix="end")
    g2<-GenomicRanges::resize(g2, width(g2) + 1000, fix="start")
  }
  #message("g2")
  if(is.null(gs)){
    clen<-as.data.frame(GenomeInfoDb::seqinfo(g2))
    to_remove<-grep("_", rownames(clen))
    if(length(to_remove)>0){
      clen<-clen[-to_remove,]
    }
    gs<-sum(clen$seqlengths)
  } else {
    message(gs)
  }
  #message("gs")
  g2<-GenomicRanges::reduce(g2)
  if(is.logical(gs)){
    gs<-sum(as.data.frame(IRanges::ranges(g2)) $width)
    return(gs)
  } else {
    tot.length<-sum(as.data.frame(IRanges::ranges(g2))$width)
    return(tot.length/gs)
  }
}

getFreqHits<-function(genes, annodf){
  return(nrow(annodf[which(annodf$SYMBOL %in% genes),]))
}

#' Function to perform GREAT analysis using custom gene sets
#'
#' This function performs the binomial test on a given set of genomic regions and a custom gene set file
#'
#' @param hits GenomicRanges object of regions to query, or string indicating the path to a bed file
#' @param gmt Named list of gene sets of string indicating path/address to custom gene set (in gmt format)
#' @param TxDb TxDb object corresponding to the genome build specific to the hits
#' @param annoDb String indiciating the annotation DB object, i.e. "org.hs.eg.db"
#' @param promoter Boolean indicating if analysis should be limited to hits annotated to only promoter regions (TRUE) or not (FALSE). Default is FALSE.
#' @param tssRegion Numeric vector of length 2 indicating the size of the promoter region around the transcription start site (default is c(-1000,1000))
#' @param parallel Boolean indicating if genomic coverage process should be run in parallel (requires packages 'parallel' and 'doparallel'). Only works for Windows at the moment. Default is TRUE.
#' @param returnCov Boolean indicating if background coverage of each term in gmt should be returned. Default is FALSE. If you plan on running the same gmt file for multiple sets of regions, it is recommended you set returnCov=T for the first iteration and use the returned object in the genCov argument to save time.
#' @param genCov Genomic coverage for each term in gmt file. Object returned if 'returnCov=T'. Default=NULL. Recommended to use if using same gmt file for multiple sets of regions.
#' @param gsName String providing identifier for gene set name in results. Default="customGMT"
#' @param FDR Boolean indicating if reported p-values should be FDR (BH) corrected. Default is TRUE.
#' @param enr String indicating type of enrichment. Default="pos" for positive enrichment. Other options: "neg" for negative enrichment or "ts" for two-sided.
#' @param significant Boolean indicating if only significant (p<0.05) results should be returned. Default=TRUE. FALSE returns all results.
#' @return Data frame of results in the same format as BinfTools::GOGEM() results for compatibility with BinfTools.
#' @export

customGREAT<-function(hits, gmt, TxDb, annoDb, promoter=F, tssRegion=c(-1000,1000), parallel=T,
                      returnCov=F, genCov=NULL, gsName="customGMT", FDR=T, enr="pos", significant=T){
  anno<-ChIPseeker::annotatePeak(hits, tssRegion, TxDb, annoDb=annoDb)
  anno<-as.data.frame(anno@anno)
  if(isTRUE(promoter)){
    anno<-anno[grep("promoter", anno$annotation, ignore.case=T),]
  }
  if(is.character(gmt)){
    message("Reading in gmt...")
    gmt<-BinfTools::read.gmt(gmt)
  }
  require(annoDb, character.only = TRUE)
  annoDb <- eval(parse(text = annoDb))
  keys<-keys(annoDb, keytype="SYMBOL")
  annoDb<-select(annoDb, keys=keys, columns=c("ENTREZID"), keytype="SYMBOL")

  if(isFALSE(promoter)){
    TxDb<-suppressMessages(GenomicFeatures::genes(TxDb))
  }
  if(isTRUE(promoter)){
    TxDb<-suppressWarnings(GenomicFeatures::promoters(TxDb, columns="gene_id", upstream = abs(tssRegion[1]),
                                                      downstream=abs(tssRegion[2])))
  }
  if(is.null(genCov)){
    message("Calculating Background Genome Size...")
    gs<-getGenomeCov(x=unique(unlist(gmt)), TxDb=TxDb, annoDb=annoDb, promoter=promoter, tssRegion=tssRegion, gs=T)
    message("Background Genome Size calculated as: ", gs,".")
    message("Genome Coverage for gene set not provided. Caculating...")
    genCov<-c()
    if(isFALSE(parallel)){
      pb<-txtProgressBar(min=0, max=length(gmt), style=3)
      for(i in 1:length(gmt)){
        genCov<-c(genCov,getGenomeCov(gmt[[i]], TxDb, annoDb, promoter, tssRegion, gs))
        names(genCov)[i]<-names(gmt)[i]
        setTxtProgressBar(pb, i)
      }
    }
    if(isTRUE(parallel)){
      message("Calculating genome coverage in parallel using half available cores...")
      no_cores<-parallel::detectCores(logical=T)
      c1<-parallel::makeCluster(no_cores/2)
      doParallel::registerDoParallel(c1)
      parallel::clusterExport(c1, list('getGenomeCov','TxDb','annoDb','promoter','tssRegion', 'gs'), envir=environment())
      genCov<-c(parallel::parLapply(c1, gmt, fun=getGenomeCov, annoDb=annoDb, TxDb=TxDb, promoter=promoter, tssRegion=tssRegion, gs=gs))
      parallel::stopCluster(c1)
    }
    genCov<-unlist(genCov[!is.na(unlist(genCov))])
    if(isTRUE(returnCov)){
      return(genCov)
    }
  }
  genCov<-genCov[which(genCov>0)]
  genCov<-genCov[which(genCov<=1)]
  gmt<-gmt[which(names(gmt) %in% names(genCov))]
  message("Obtaining hit frequency...")
  hits<-unlist(lapply(gmt, getFreqHits, annodf=anno))

  res<-data.frame(query=rep("query_1", length(gmt)),
                  significant=rep("FALSE", length(gmt)),
                  p_value=rep(1, length(gmt)),
                  term_size=lengths(gmt),
                  term_cov=genCov,
                  query_size=rep(nrow(anno), length(gmt)),
                  intersection_size=hits,
                  precision=rep(0, length(gmt)),
                  recall=rep(0, length(gmt)),
                  #term_id=names(gmt),
                  source=rep(gsName, length(gmt)),
                  term_name=names(gmt),
                  #effective_domain_size=length(bg),
                  #intersection=unlist(lapply(OL2, toString)),
                  enrichment=rep(0, length(gmt)))

  vectorBinom<-function(intersection_size, query_size, term_cov, adj=T, alt="greater"){
    pvals<-c()
    #pb<-txtProgressBar(0,1,style=3)
    for(i in 1:length(intersection_size)){
      pvals<-c(pvals, binom.test(intersection_size[i], query_size[i], term_cov[i], alternative=alt)$p.value)
      #setTxtProgressBar(pb, i)
    }
    if(isTRUE(adj)){
      pvals<-p.adjust(pvals, "BH")
    }
    return(pvals)
  }

  if(enr=="pos"){
    message("Returning only positively enriched pathways...")
    enr<-"greater"
  } else {
    if(enr=="neg"){
      message("Returning only negatively enriched pathways...")
      enr<-"less"
    } else {
      message("Running two-sided binomial analysis...")
      enr<-"two.sided"
    }
  }
  #print(dim(res))
  message("Running binomial test...")
  res$p_value<-vectorBinom(res$intersection_size, res$query_size, res$term_cov, adj=FDR, alt=enr)
  res$significant<-ifelse(res$p_value<0.05, TRUE, FALSE)
  res$precision<-res$intersection_size/res$query_size
  res$recall<-res$intersection_size/res$term_cov
  res$enrichment<-(res$intersection_size/res$query_size)/(res$term_cov)
  res<-res[order(res$enrichment, decreasing=T),]

  if(isTRUE(significant)){
    if(nrow(subset(res, significant=="TRUE"))<1){
      message("No significant results. Returning NA.")
      return(NA)
    } else {
      return(subset(res, significant=="TRUE"))
    }
  } else {
    return(res)
  }
}
