% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions_BSS_SGL.R
\name{block.finder}
\alias{block.finder}
\title{cluster the points by neighborhood size a_n}
\usage{
block.finder(pts, an)
}
\arguments{
\item{pts}{vector of candidate change points}

\item{an}{radius of the cluster}
}
\value{
A list of change points clusters
}
\description{
helper function for determining the clusters
}
\keyword{internal}
