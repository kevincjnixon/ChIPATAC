% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BedFunctions.R
\name{writeGTF}
\alias{writeGTF}
\title{Function to export GenomicRanges to GTF format}
\usage{
writeGTF(x, source = "ATAC-peakset", peakID, filename)
}
\arguments{
\item{x}{GenomicRanges object to export}

\item{source}{String to identify the source of peaks in GTF file. Default="ATAC-peakset"}

\item{peakID}{String or string vector to identify unique peaks. Will be set as 'gene_id' in gtf file and concatenated with row number of each peak.}

\item{filename}{String of path to export gtf file (should end in .gtf)}
}
\description{
Function to export GenomicRanges to GTF format
}
