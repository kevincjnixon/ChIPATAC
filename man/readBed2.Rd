% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FactorFootprints.R
\name{readBed2}
\alias{readBed2}
\title{Read a bed file to genomicRanges Object with option for seqlengths}
\usage{
readBed2(filename, HOCO = F, species = "human", BSg = NULL)
}
\arguments{
\item{filename}{Character of .bed file name and location (assumes no header)}

\item{HOCO}{Boolean indicaitng if .bed file is HOCOMOCO TFBS bed file. Default is FALSE}

\item{species}{One of "human" or "mouse" to use already specified BSgenome objects. Leave NULL to use custom BSgenome.}

\item{BSg}{BSgenome object containign seqlengths for genome of interest. Use if you want to specify specific BSgenome (especially if not mouse or human). Or leave NULL along with 'species' to not include seqlengths.}
}
\value{
GRanges object of bedfile regions
}
\description{
Read a bed file to genomicRanges Object with option for seqlengths
}
