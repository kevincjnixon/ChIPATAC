% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FeatureEnrichment.R
\name{FeatureEnrichment}
\alias{FeatureEnrichment}
\title{Calculate empirical p-value for genomic feature enrichment}
\usage{
FeatureEnrichment(
  bedFile,
  species = "hsapiens",
  N = 1000,
  tss = c(-3000, 3000),
  sys = c("win", "")
)
}
\arguments{
\item{bedFile}{String with path to bedfile to analyze. Must be a bed file}

\item{species}{string indicating genus/species of genome. One of "hsapiens" or "mmusculus"}

\item{N}{numeric indicating number of iterations}

\item{tss}{numeric vector of length 2 indicating the promoter region around the tss site (default= c(-3000,3000))}

\item{sys}{string indicating opertating system. - Not working at the moment, currently only works for windows running windows subsystem for linux with bedtools installed}
}
\value{
Data frame with p-values and adjusted p-values (BH corrected) for each genomic feature
}
\description{
Calculate empirical p-value for genomic feature enrichment
}
