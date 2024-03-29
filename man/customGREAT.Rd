% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GREAT.R
\name{customGREAT}
\alias{customGREAT}
\title{Function to perform GREAT analysis using custom gene sets}
\usage{
customGREAT(
  hits,
  gmt,
  TxDb,
  annoDb,
  promoter = F,
  tssRegion = c(-1000, 1000),
  parallel = T,
  returnCov = F,
  genCov = NULL,
  gsName = "customGMT",
  FDR = T,
  enr = "pos",
  significant = T
)
}
\arguments{
\item{hits}{GenomicRanges object of regions to query, or string indicating the path to a bed file}

\item{gmt}{Named list of gene sets of string indicating path/address to custom gene set (in gmt format)}

\item{TxDb}{TxDb object corresponding to the genome build specific to the hits}

\item{annoDb}{String indiciating the annotation DB object, i.e. "org.hs.eg.db"}

\item{promoter}{Boolean indicating if analysis should be limited to hits annotated to only promoter regions (TRUE) or not (FALSE). Default is FALSE.}

\item{tssRegion}{Numeric vector of length 2 indicating the size of the promoter region around the transcription start site (default is c(-1000,1000))}

\item{parallel}{Boolean indicating if genomic coverage process should be run in parallel (requires packages 'parallel' and 'doparallel'). Only works for Windows at the moment. Default is TRUE.}

\item{returnCov}{Boolean indicating if background coverage of each term in gmt should be returned. Default is FALSE. If you plan on running the same gmt file for multiple sets of regions, it is recommended you set returnCov=T for the first iteration and use the returned object in the genCov argument to save time.}

\item{genCov}{Genomic coverage for each term in gmt file. Object returned if 'returnCov=T'. Default=NULL. Recommended to use if using same gmt file for multiple sets of regions.}

\item{gsName}{String providing identifier for gene set name in results. Default="customGMT"}

\item{FDR}{Boolean indicating if reported p-values should be FDR (BH) corrected. Default is TRUE.}

\item{enr}{String indicating type of enrichment. Default="pos" for positive enrichment. Other options: "neg" for negative enrichment or "ts" for two-sided.}

\item{significant}{Boolean indicating if only significant (p<0.05) results should be returned. Default=TRUE. FALSE returns all results.}
}
\value{
Data frame of results in the same format as BinfTools::GOGEM() results for compatibility with BinfTools.
}
\description{
This function performs the binomial test on a given set of genomic regions and a custom gene set file
}
