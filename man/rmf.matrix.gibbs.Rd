% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmf.matrix.gibbs.R
\name{rmf.matrix.gibbs}
\alias{rmf.matrix.gibbs}
\title{Gibbs Sampling for the Matrix-variate von Mises-Fisher Distribution}
\usage{
rmf.matrix.gibbs(M, X, rscol = NULL)
}
\arguments{
\item{M}{a matrix.}

\item{X}{the current value of the random orthonormal matrix.}

\item{rscol}{the number of columns to update simultaneously.}
}
\value{
a new value of the matrix \code{X} obtained by Gibbs sampling.
}
\description{
Simulate a random orthonormal matrix from the matrix von Mises-Fisher
distribution using Gibbs sampling.
}
\note{
This provides one Gibbs scan. The function should be used iteratively.
}
\examples{

Z<-matrix(rnorm(10*5),10,5) 

U<-rmf.matrix(Z)
U<-rmf.matrix.gibbs(Z,U)


## The function is currently defined as
function (M, X, rscol = NULL) 
{
    if (is.null(rscol)) {
        rscol <- max(2, min(round(log(dim(M)[1])), dim(M)[2]))
    }
    sM <- svd(M)
    H <- sM$u \%*\% diag(sM$d)
    Y <- X \%*\% sM$v
    m <- dim(H)[1]
    R <- dim(H)[2]
    for (iter in 1:round(R/rscol)) {
        r <- sample(seq(1, R, length = R), rscol)
        N <- NullC(Y[, -r])
        y <- rmf.matrix(t(N) \%*\% H[, r])
        Y[, r] <- N \%*\% y
    }
    Y \%*\% t(sM$v)
  }

}
\references{
Hoff(2009)
}
\author{
Peter Hoff
}
