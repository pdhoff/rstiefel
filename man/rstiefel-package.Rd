\name{rstiefel-package}
\alias{rstiefel-package}
\alias{rstiefel}
\docType{package}
\title{
Random Orthonormal Matrix Generation on the Stiefel Manifold
}
\description{
This package simulates random orthonormal matrices from linear and quadratic exponential family distributions on the Stiefel manifold. The most general type of distribution covered is the matrix-variate  Bingham-von Mises-Fisher distribution. Most of the simulation methods are presented in Hoff(2009) "Simulation of the Matrix Bingham-von Mises-Fisher Distribution, With Applications to Multivariate and Relational Data."
}
\details{
\tabular{ll}{
Package: \tab rstiefel\cr
Type: \tab Package\cr
Version: \tab 0.10\cr
Date: \tab 2014-11-24\cr
License: \tab GPL-3\cr
}

}
\author{
Peter Hoff 

Maintainer: Peter Hoff <pdhoff@uw.edu>
}
\references{
Hoff(2009)
}
\keyword{ package }
\examples{
Z<-matrix(rnorm(10*5),10,5) ; A<-t(Z)\%*\%Z 
B<-diag(sort(rexp(5),decreasing=TRUE))
U<-rbing.Op(A,B)
U<-rbing.matrix.gibbs(A,B,U)

U<-rmf.matrix(Z)
U<-rmf.matrix.gibbs(Z,U)
}
