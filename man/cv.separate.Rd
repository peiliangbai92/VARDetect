% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions_LpS.R
\name{cv.separate}
\alias{cv.separate}
\title{cross-validation index function, separate train and test sets}
\usage{
cv.separate(period, nfold = 5)
}
\arguments{
\item{period}{an input dataset with size of n by p}

\item{nfold}{an positive integer, indicates the number of folds in cv}
}
\value{
a list object, includes
\describe{
    \item{train}{training dataset}
    \item{test}{test dataset}
}
}
\description{
a function to setup training and testing datasets, select data point in test every nfold time points
}
\keyword{internal}
