options(connectionObserver = NULL)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ATACseqQC)
library(MotifDb)
library(motifStack)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)


#' Read a bed file to genomicRanges Object with option for seqlengths
#'
#' @param filename Character of .bed file name and location (assumes no header)
#' @param HOCO Boolean indicaitng if .bed file is HOCOMOCO TFBS bed file. Default is FALSE
#' @param species One of "human" or "mouse" to use already specified BSgenome objects. Leave NULL to use custom BSgenome.
#' @param BSg BSgenome object containign seqlengths for genome of interest. Use if you want to specify specific BSgenome (especially if not mouse or human). Or leave NULL along with 'species' to not include seqlengths.
#' @return GRanges object of bedfile regions
#' @export

readBed2<-function(filename, HOCO=F, species="human", BSg=NULL){
  bed<-read.delim(filename, header=F)
  if(isTRUE(HOCO)){
    colnames(bed)<-c("seqnames","start","end","motif","score","strand")
  } else {
    colnames(bed)[c(1:3)]<-c("seqnames","start","end")
  }
  bed<-GenomicRanges::makeGRangesFromDataFrame(bed, keep.extra.columns = T)
  seqs<-names(seqlengths(bed))
  if(species=="human"){
    seqlengths(bed)<-seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[seqs]
  }
  if(species=="mouse"){
    seqlengths(bed)<-seqlengths(BSgenome.Mmusculus.UCSC.mm10::Mmusculus)[seqs]
  }
  if(is.null(species)){
    if(!is.null(BSg)){
      seqlengths(bed)<-seqlengths(BSg)[seps]
    }
  }
  return(bed)

}

roundup<-function(motif=CTCF){
  for(i in 1:ncol(motif)){
    if(round(colSums(motif)[i], digits=4) != 1){
      diff<- 1-round(colSums(motif)[i], digits=4)
      motif[which(motif[,i]==max(motif[,i])),i]<-motif[which(motif[,i]==max(motif[,i])),i]+diff
    }
  }
  return(motif)
}

factorFootprints2 <- function(bamfiles, index=bamfiles, pfm, genome,
                              min.score="95%", bindingSites,
                              seqlev=paste0("chr", c(1:22, "X", "Y")),
                              upstream=100, downstream=100,
                              maxSiteNum=1e6, anchor="cut site",
                              group="strand", deleteNoCoverageBindingSites=T,
                              ...){
  #stopifnot(length(bamfiles)==4)
  if(!missing(genome)) stopifnot(is(genome, "BSgenome"))
  stopifnot(all(round(colSums(pfm), digits=4)==1))
  stopifnot(upstream>10 && downstream>10)
  stopifnot(is.numeric(maxSiteNum))
  seqlev<-seqlevels(bindingSites)
  stopifnot(all(seqlev %in% names(getBamTargets(bamfiles[1], index[1]))))
  maxSiteNum <- ceiling(maxSiteNum[1])
  stopifnot(maxSiteNum>1)
  anchor <- match.arg(anchor, choices = c("cut site", "fragment center"))
  if(length(group) == 1){
    if(group!="strand"){
      stop("If length of group is 1, it must be strand.")
    }
    groupFlag <- TRUE
  }else{
    stopifnot(length(group)==length(bamfiles))
    group <- as.factor(group)
    if(length(levels(group))!=2){
      stop("The length of levels of group must be 2.")
    }
    groupFlag <- FALSE
  }
  if(anchor=="fragment center"){
    null <- mapply(function(.bam, .index){
      suppressMessages(pe <- testPairedEndBam(.bam, .index))
      if(!pe){
        stop("If anchor is fragment center, the bamfiles must be paired end reads.")
      }
    }, bamfiles, index)
  }
  if(missing(bindingSites)){
    pwm <- motifStack::pfm2pwm(pfm)
    maxS <- maxScore(pwm)
    if(!is.numeric(min.score)){
      if(!is.character(min.score)){
        stop("'min.score' must be a single number or string")
      }else{
        nc <- nchar(min.score)
        if (substr(min.score, nc, nc) == "%"){
          min.score <- substr(min.score, 1L, nc - 1L)
        }else{
          stop("'min.score' can be given as a character string containing a percentage",
               "(e.g. '85%') of the highest possible score")
        }
        min.score <- maxS * as.double(min.score)/100
      }
    }else{
      min.score <- min.score[1]
    }
    predefined.score <- maxS * as.double(0.85)
    suppressWarnings({
      mt <- tryCatch(matchPWM(pwm, genome, min.score = min(predefined.score, min.score),
                              with.score=TRUE,
                              exclude=paste0("^", names(genome)[!names(genome) %in% seqlev], "$")),
                     error=function(e){
                       message(e)
                       stop("No predicted binding sites available by giving seqlev. ")
                     })
    })
    if (min.score <= predefined.score){
      mt$userdefined <- TRUE
    } else {
      mt$userdefined <- FALSE
      mt$userdefined[mt$score >= min.score] <- TRUE
    }
    if(length(mt)>maxSiteNum){## subsample
      mt$oid <- seq_along(mt)
      mt <- mt[order(mt$score, decreasing = TRUE)]
      mt <- mt[seq.int(maxSiteNum)]
      mt <- mt[order(mt$oid)]
      mt$oid <- NULL
    }
    subsampleNum <- min(10000, maxSiteNum)
    if(length(mt)>subsampleNum && sum(mt$userdefined)<subsampleNum){## subsample
      set.seed(seed = 1)
      mt.keep <- seq_along(mt)[!mt$userdefined]
      n <- subsampleNum-sum(mt$userdefined)
      if(length(mt.keep)>n){
        mt.keep <- mt.keep[order(mt[mt.keep]$score, decreasing = TRUE)]
        mt.keep <- mt.keep[seq.int(n)]
        mt.keep <- sort(mt.keep)
        mt.keep <- seq_along(mt) %in% mt.keep
        mt <- mt[mt$userdefined | mt.keep]
      }
    }
  }else{
    stopifnot(is(bindingSites, "GRanges"))
    stopifnot(all(!is.na(seqlengths(bindingSites))))
    stopifnot(length(bindingSites)>1)
    stopifnot(length(bindingSites$score)==length(bindingSites))
    mt <- bindingSites
    mt$userdefined <- TRUE
  }
  if(sum(mt$userdefined)<2){
    stop("less than 2 binding sites by input min.score")
  }
  seqlevelsStyle(mt) <- checkBamSeqStyle(bamfiles[1], index[1])[1]
  wid <- ncol(pfm)
  #mt <- mt[seqnames(mt) %in% seqlev]
  seqlevels(mt) <- seqlev
  seqinfo(mt) <- Seqinfo(seqlev, seqlengths = seqlengths(mt))
  ## read in bam file with input seqlev specified by users
  which <- as(seqinfo(mt), "GRanges")
  #param <- ScanBamParam(which=which)
  if(anchor=="cut site"){
    bamIn <- mapply(function(.b, .i) {
      seqlevelsStyle(which) <- checkBamSeqStyle(.b, .i)[1]
      param <- Rsamtools::ScanBamParam(which=which)
      GenomicAlignments::readGAlignments(.b, .i, param = param)
    }, bamfiles, index, SIMPLIFY = FALSE)
  }else{
    bamIn <- mapply(function(.b, .i) {
      seqlevelsStyle(which) <- checkBamSeqStyle(.b, .i)[1]
      param <- Rsamtools::ScanBamParam(which=which)
      GenomicAlignments::readGAlignmentPairs(.b, .i, param = param)
    }, bamfiles, index, SIMPLIFY = FALSE)
  }
  bamIn <- lapply(bamIn, as, Class = "GRanges")
  bamIn <- split(bamIn, group)
  bamIn <- lapply(bamIn, function(.ele){
    if(!is(.ele, "GRangesList")) .ele <- GRangesList(.ele)
    .ele <- unlist(.ele)
    #seqlevelsStyle(.ele) <- seqlevelsStyle(genome)[1]
    if(anchor=="cut site"){
      ## keep 5'end as cutting sites
      promoters(.ele, upstream=0, downstream=1)
    }else{
      ## keep fragment center
      ChIPpeakAnno::reCenterPeaks(.ele, width=1)
    }
  })
  libSize <- lengths(bamIn)
  if(groupFlag){
    libFactor <- mapply(bamIn, libSize, FUN=function(.ele, .libSize){
      coverageSize <- sum(as.numeric(width(reduce(.ele, ignore.strand=TRUE))))
      .libSize / coverageSize
    })
  }else{
    libFactor <- libSize/mean(libSize)
  }


  if(groupFlag){
    ## split into positive strand and negative strand
    bamIn <- split(bamIn[[1]], strand(bamIn[[1]]))
    ## get coverage
    cvglist <- sapply(bamIn, coverage)
    cvglist <- cvglist[c("+", "-")]
    cvglist <- lapply(cvglist, function(.ele)
      .ele[names(.ele) %in% seqlev])
    ## coverage of mt, must be filtered, otherwise too much
    cvgSum <- cvglist[["+"]] + cvglist[["-"]] ## used for filter mt
  }else{
    ## it was splitted by groups
    ## get coverage
    cvglist <- sapply(bamIn, coverage)
    cvglist <- lapply(cvglist, function(.ele)
      .ele[names(.ele) %in% seqlev])
    ## coverage of mt, must be filtered, otherwise too much
    cvgSum <- cvglist[[1]] + cvglist[[2]] ## used for filter mt
  }

  mt.s <- split(mt, seqnames(mt))
  seqlev <- intersect(names(cvgSum), names(mt.s))
  cvgSum <- cvgSum[seqlev]
  mt.s <- mt.s[seqlev]
  ## too much if use upstream and downstream, just use 3*wid maybe better.
  mt.s.ext <- promoters(mt.s, upstream=wid, downstream=wid+wid)
  stopifnot(all(lengths(mt.s.ext)==lengths(mt.s)))
  mt.v <- Views(cvgSum, mt.s.ext)
  if(deleteNoCoverageBindingSites){
    mt.s <- mt.s[viewSums(mt.v)>0]
  }
  mt <- unlist(mt.s)
  mt.ids <- promoters(ChIPpeakAnno::reCenterPeaks(mt, width=1),
                      upstream=upstream+floor(wid/2),
                      downstream=downstream+ceiling(wid/2)+1)
  mt.ids <- paste0(as.character(seqnames(mt.ids)), ":", start(mt.ids), "-", end(mt.ids))
  sigs <- ChIPpeakAnno::featureAlignedSignal(cvglists=cvglist,
                                             feature.gr=ChIPpeakAnno::reCenterPeaks(mt, width=1),
                                             upstream=upstream+floor(wid/2),
                                             downstream=downstream+ceiling(wid/2),
                                             n.tile=upstream+downstream+wid)
  mt <- mt[match(rownames(sigs[[1]]), mt.ids)]
  cor <- lapply(sigs, function(sig){
    sig.colMeans <- colMeans(sig)
    ## calculate correlation of footprinting and binding score
    windows <- slidingWindows(IRanges(1, ncol(sig)), width = wid, step = 1)[[1]]
    # remove the windows with overlaps of motif binding region
    windows <- windows[end(windows)<=upstream | start(windows)>=upstream+wid]
    sig.windowMeans <- viewMeans(Views(sig.colMeans, windows))
    windows.sel <- windows[which.max(sig.windowMeans)][1]
    highest.sig.windows <-
      rowMeans(sig[, start(windows.sel):end(windows.sel)])
    predictedBindingSiteScore <- mt$score
    if(length(predictedBindingSiteScore) == length(highest.sig.windows)){
      suppressWarnings({
        cor <- cor.test(x = predictedBindingSiteScore,
                        y = highest.sig.windows,
                        method = "spearman")
      })
    }else{
      cor <- NA
    }
    cor
  })
  sigs <- lapply(sigs, function(.ele) .ele[mt$userdefined, ])
  mt <- mt[mt$userdefined]
  mt$userdefined <- NULL
  ## segmentation the signals
  if(groupFlag){
    ## x2 because stranded.
    Profile <- lapply(sigs, function(.ele) colMeans(.ele, na.rm = TRUE)*2/libFactor[1])
  }else{
    Profile <- mapply(sigs, libFactor, FUN = function(.ele, .libFac){
      colMeans(.ele, na.rm = TRUE)/.libFac
    }, SIMPLIFY = FALSE)
  }
  ## upstream + wid + downstream
  Profile.split <- lapply(Profile, function(.ele){
    list(upstream=.ele[seq.int(upstream)],
         binding=.ele[upstream+seq.int(wid)],
         downstream=.ele[upstream+wid+seq.int(downstream)])
  })
  optimalSegmentation <- function(.ele){
    .l <- length(.ele)
    short_abun <- cumsum(.ele)/seq.int(.l)
    long_abun <- cumsum(rev(.ele))/seq.int(.l)
    long_abun <- rev(long_abun)
    short_abun <- short_abun[-length(short_abun)]
    long_abun <- long_abun[-1]
    ##long_abun should always greater than short_abun
    long_abun <- long_abun - short_abun
    long_abun[long_abun<0] <- 0
    cov_diff <- numeric(length(short_abun))
    for(i in seq_along(.ele)){
      cov_diff_tmp <- .ele
      cov_diff_tmp <- cov_diff_tmp-short_abun[i]
      cov_diff_tmp[-seq.int(i)] <- cov_diff_tmp[-seq.int(i)] - long_abun[i]
      cov_diff[i] <- mean(cov_diff_tmp^2)
    }
    .ids <- which(cov_diff==min(cov_diff, na.rm = TRUE))
    data.frame(pos=.ids, short_abun=short_abun[.ids], long_abun=long_abun[.ids])
  }
  Profile.seg <- lapply(Profile.split, function(.ele){
    ups <- optimalSegmentation(.ele$upstream)
    downs <- optimalSegmentation(rev(.ele$downstream))
    ## find the nearest pair
    .min <- c(max(rbind(ups, downs)), 0, 0)
    for(i in seq.int(nrow(ups))){
      for(j in seq.int(nrow(downs))){
        tmp <- sum(abs(ups[i, -1] - downs[j, -1]))
        if(tmp < .min[1]){
          .min <- c(tmp, i, j)
        }
      }
    }
    c(colMeans(rbind(ups[.min[2], ], downs[.min[3], ])), binding=mean(.ele$binding, na.rm=TRUE))
  })
  Profile.seg <- colMeans(do.call(rbind, Profile.seg))
  Profile.seg[3] <- Profile.seg[2]+Profile.seg[3]
  names(Profile.seg)[2:3] <- c("distal_abun", "proximal_abun")
  tryCatch({ ## try to avoid the error when ploting.
    args <- list(...)
    if(groupFlag){
      args$Profile <- c(Profile[["+"]], Profile[["-"]])
      args$legLabels <- c("For. strand", "Rev. strand")
    }else{
      args$Profile <- c(Profile[[1]], Profile[[2]])
      args$legLabels <- names(Profile)
    }
    args$ylab <- ifelse(anchor=="cut site", "Cut-site probability", "reads density (arbitrary unit)")
    args$Mlen <- wid
    args$motif <- pwm2pfm(pfm)
    args$segmentation <- Profile.seg
    do.call(plotFootprints, args = args)
  }, error=function(e){
    message(e)
  })
  return(list(signal=sigs,
              spearman.correlation=cor,
              bindingSites=mt,
              Mlen=wid,
              estLibSize=libSize,
              Profile.segmentation=Profile.seg))
}

pwm2pfm <- function(pfm, name="motif"){
  if(!all(round(colSums(pfm), digits=4)==1)){
    return(NULL)
  }
  new("pfm", mat=as.matrix(pfm), name=name)
}

getBamTargets <- function(bamfile, index){
  header <- Rsamtools::scanBamHeader(bamfile, index=index)
  header[[1]]$targets
}
checkBamSeqStyle <- function(bamfile, index){
  which <- getBamTargets(bamfile, index)
  seqlevelsStyle(names(which))
}


listFFp2<-function(bs, bam, pfm, delete=T){
  bs<-readBed(bs, HOCO=T)
  factorFootprints2(bamfiles=bam, index=bam, bindingSites=bs, pfm=pfm, genome=Hsapiens, min.score="95%",
                    deleteNoCoverageBindingSites = delete, seqlev=names(seqlengths(bs)))
}

bsAnno<-function(bs, species="human"){
  names(bs)<-c(1:length(bs))
  Tx<-TxDb.Hsapiens.UCSC.hg38.knownGene
  adb<-"org.Hs.eg.db"
  if(species!="human"){
    Tx<-TxDb.Mmusculus.UCSC.mm10.knownGene
    adb<-"org.Mm.eg.db"
  }
  anno<-annotatePeak(bs, tssRegion=c(-1000,1000), TxDb=Tx,
                     annoDb=adb)
  return(as.data.frame(anno@anno))
}

plotFootprints2<-function (Profile, Mlen = 0, xlab = "Dist. to motif (bp)", ylab = "Cut-site probability",
                           legTitle, newpage = TRUE, motif, cols="Dark2",
                           anno=NULL, feature="all", genes=NULL, legend=T, ymax=NULL, refSamp=NULL)
{
  stopifnot(is(motif, "pfm"))
  if (newpage)
    grid.newpage()
  maxval<-c()
  if(!is.null(refSamp)){
    rows<-which(rowSums(Profile[[refSamp]][[1]])>0)
    message("Keeping ",length(rows)," from ", names(Profile)[refSamp],"...")
    for(i in 1:length(Profile)){
      Profile[[i]][[1]]<-Profile[[i]][[1]][rows,]
      Profile[[i]][[2]]<-Profile[[i]][[2]][rows,]
      if(!is.null(anno)){
        anno[[i]]<-anno[[i]][rows,]
      }
    }
  }
  for(i in 1:length(Profile)){
    rows<-c(1:nrow(Profile[[i]][[1]]))
    if(!is.null(anno)){
      if(!is.null(genes) && feature=="all"){
        rows<-which(anno[[i]]$SYMBOL %in% genes)
        message("Subsetting to ",length(rows)," annotated to ",length(genes)," genes...")
      }
      if(feature!="all"){
        if(!is.null(genes)){
          rows<-which(anno[[i]]$SYMBOL %in% genes)
          rows2<-grep(feature, anno[[i]]$annotation)
          rows<-rows[which(rows %in% rows2)]
          message("Subsetting to ", length(rows), "annotated to ", length(genes)," genes and ", feature," features...")
        } else {
          rows<-grep(feature, anno[[i]]$annotation)
          message("Subsetting to ",length(rows)," ",feature,"features...")
        }
      }
    }
    for(k in 1:length(Profile[[i]])){
      Profile[[i]][[k]]<-Profile[[i]][[k]][rows,]
    }
    message("Profile rows: ", nrow(Profile[[i]][[1]]))
    Profile[[i]]<-colMeans(do.call(cbind, Profile[[i]]))
    maxval<-c(maxval, max(Profile[[i]]))
  }
  if(!is.null(ymax)){
    maxval<-ymax
  }

  S <- length(Profile[[1]])
  W <- ((S/2) - Mlen)/2
  vp <- plotViewport(name = "plotRegion")
  pushViewport(vp)
  vp1 <- viewport(y = 0.4, height = 0.8, xscale = c(0, S/2 +
                                                      1), yscale = c(0, max(maxval) * 1.12), name = "footprints")
  pushViewport(vp1)

  for(i in 1:length(Profile)){
    avg<-colMeans(rbind(Profile[[i]][1:(S/2)], Profile[[i]][(S/2+1):S]))
    grid.lines(x=1:(S/2), y=avg, default.units="native", gp=gpar(lwd=2, col=BinfTools::colPal(cols)[i]))
  }
  grid.xaxis(at = c(seq(1, W, length.out = 3), W + seq(1,
                                                       Mlen), W + Mlen + seq(1, W, length.out = 3)), label = c(-(W +
                                                                                                                   1 - seq(1, W + 1, length.out = 3)), rep("", Mlen), seq(0,
                                                                                                                                                                          W, len = 3)))
  grid.yaxis()
  grid.lines(x = c(W, W, 0), y = c(0, max(maxval), max(maxval) *
                                     1.12), default.units = "native", gp = gpar(lty = 2))
  grid.lines(x = c(W + Mlen + 1, W + Mlen + 1, S/2), y = c(0,
                                                           max(maxval), max(maxval) * 1.12), default.units = "native",
             gp = gpar(lty = 2))
  upViewport()
  vp2 <- viewport(y = 0.9, height = 0.2, xscale = c(0, S/2 +
                                                      1), name = "motif")
  pushViewport(vp2)
  motifStack::plotMotifLogoA(motif)
  upViewport()
  upViewport()
  grid.text(xlab, y = unit(1, "lines"))
  grid.text(ylab, x = unit(1, "line"), rot = 90)
  if(isTRUE(legend)){
    if (missing(legTitle)) {
      legvp <- viewport(x = unit(1, "npc") - convertX(unit(1,
                                                           "lines"), unitTo = "npc"), y = unit(1, "npc") -
                          convertY(unit(1, "lines"), unitTo = "npc"), width = convertX(unit(14,
                                                                                            "lines"), unitTo = "npc"), height = convertY(unit(3,
                                                                                                                                              "lines"), unitTo = "npc"), just = c("right", "top"),
                        name = "legendWraper")
      pushViewport(legvp)
      # grid.legend(labels = c("For. strand", "Rev. strand", "Avg"),
      #             gp = gpar(lwd = 2, lty = 1, col = c("darkblue",
      #                                                 "darkred",
      #                                                 "green")))
      grid.legend(labels=names(Profile), gp=gpar(lwd=2, lty=1, col=BinfTools::colPal(cols)))
      upViewport()
    }
    else {
      grid.text(legTitle, y = unit(1, "npc") - convertY(unit(1,
                                                             "lines"), unitTo = "npc"), gp = gpar(cex = 1.2,
                                                                                                  fontface = "bold"))
    }
  }
  return(invisible())
}
