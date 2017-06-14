\name{rbmf.matrix.gibbs}
\alias{rbmf.matrix.gibbs}
\title{
Gibbs Sampling for the Matrix-variate Bingham-von Mises-Fisher Distribution.
}
\description{
Simulate a random orthonormal matrix from the Bingham distribution using Gibbs sampling.
}
\usage{
rbmf.matrix.gibbs(A, B, C, X)
}
\arguments{
  \item{A}{
a symmetric matrix.
}
  \item{B}{
a diagonal matrix with decreasing entries.
}
  \item{C}{
a matrix with the same dimension as X. 
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
## The function is currently defined as
function (A, B, C, X) 
{
    m <- dim(X)[1]
    R <- dim(X)[2]
    if (m > R) {
        for (r in sample(seq(1, R, length = R))) {
            N <- NullC(X[, -r])
            An <- B[r, r] * t(N) \%*\% (A) \%*\% N
            cn <- t(N) \%*\% C[, r]
            X[, r] <- N \%*\% rbmf.vector.gibbs(An, cn, t(N) \%*\% 
                X[, r])
        }
    }
    if (m == R) {
        for (s in seq(1, R, length = R)) {
            r <- sort(sample(seq(1, R, length = R), 2))
            N <- NullC(X[, -r])
            An <- t(N) \%*\% A \%*\% N
            Cn <- t(N) \%*\% C[, r]
            X[, r] <- N \%*\% rbmf.O2(An, B[r, r], Cn)
        }
    }
    X
  }
}
