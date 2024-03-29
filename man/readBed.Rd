% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BedFunctions.R
\name{readBed}
\alias{readBed}
\title{Function to read tab-delimited bed file to GenomicRanges object}
\usage{
readBed(filename, header = F)
}
\arguments{
\item{filename}{String of path to bed file to read in}

\item{header}{Boolean indicating if there is a header in the bed file. Default=FALSE.}
}
\value{
GenomicRanges object
}
\description{
Function to read tab-delimited bed file to GenomicRanges object
}
