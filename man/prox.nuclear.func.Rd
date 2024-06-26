% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions_LpS.R
\name{prox.nuclear.func}
\alias{prox.nuclear.func}
\title{Proximal function with nuclear norm penalty updating}
\usage{
prox.nuclear.func(w1, y, A, b, L, lambda, AtA, Atb)
}
\arguments{
\item{w1}{previously updated model parameter}

\item{y}{updated model parameter}

\item{A}{design matrix}

\item{b}{correspond vector, or matrix}

\item{L}{learning rate}

\item{lambda}{tuning parameter for low-rank component}

\item{AtA}{Gram matrix of design matrix A}

\item{Atb}{inner product of design matrix A and correspond vector b}
}
\value{
Value of proximal function with nuclear norm penalty
}
\description{
Proximal function with nuclear norm
}
\keyword{internal}
