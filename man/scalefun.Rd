% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotMetagene.R
\name{scalefun}
\alias{scalefun}
\title{Scale the upstream/downstream regions of a coverage matrix to a certain number of columns}
\usage{
scalefun(df, cols = c(0, 100, 300, 700), b = 100, a = 100)
}
\arguments{
\item{df}{data.frame object of a single coverage matrix}

\item{cols}{numeric vector of length 4 indicating landmark columns (start=0, end of upstream region, end of main region, ncol(df))}

\item{b}{numeric indicating the number of columns the upstream region (cols[1]-cols[2]) should be scaled down to. Default=100.}

\item{a}{numeric indicating the number of columns the downstream region (cols[3]-cols[4]) should be scaled down to. Default=100.}
}
\value{
Data.frame of a coverage matrix with upstream and downstream regions scaled down (i.e. reduced columns) for clearer plotting
}
\description{
Scale the upstream/downstream regions of a coverage matrix to a certain number of columns
}
