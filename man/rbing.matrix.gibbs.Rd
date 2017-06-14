\name{rbing.matrix.gibbs}
\alias{rbing.matrix.gibbs}
\title{Gibbs Sampling for the Matrix-variate Bingham Distribution
}
\description{
Simulate a random orthonormal matrix from the Bingham distribution using Gibbs sampling. 
}
\usage{
rbing.matrix.gibbs(A, B, X)
}
\arguments{
  \item{A}{
a symmetric matrix. 
}
  \item{B}{
a diagonal matrix with decreasing entries.
}
  \item{X}{
the current value of the random orthonormal matrix. 
}
}
\value{
a new value of the matrix \code{X} obtained by Gibbs sampling. 
}
\note{This provides one Gibbs scan. The function should be used iteratively.}
\references{
Hoff(2009)
}
\author{
Peter Hoff
}

\examples{
Z<-matrix(rnorm(10*5),10,5) ; A<-t(Z)\%*\%Z
B<-diag(sort(rexp(5),decreasing=TRUE))
U<-rbing.Op(A,B)
U<-rbing.matrix.gibbs(A,B,U)

## The function is currently defined as
function (A, B, X) 
{
    m <- dim(X)[1]
    R <- dim(X)[2]
    if (m > R) {
        for (r in sample(seq(1, R, length = R))) {
            N <- NullC(X[, -r])
            An <- B[r, r] * t(N) \%*\% (A) \%*\% N
            X[, r] <- N \%*\% rbing.vector.gibbs(An, t(N) \%*\% X[, 
                r])
        }
    }
    if (m == R) {
        for (s in seq(1, R, length = R)) {
            r <- sort(sample(seq(1, R, length = R), 2))
            N <- NullC(X[, -r])
            An <- t(N) \%*\% A \%*\% N
            X[, r] <- N \%*\% rbing.Op(An, B[r, r])
        }
    }
    X
  }
}
