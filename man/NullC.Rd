\name{NullC}
\alias{NullC}
\title{
Null Space of a Matrix
}
\description{
Given a matrix \code{M}, find a matrix \code{N} giving a basis for the null
space. This is a modified version of Null from the package MASS.
}
\usage{
NullC(M)
}
\arguments{
  \item{M}{input matrix. 
}
}
\value{
an orthonormal matrix such that \code{t(N)\%*\%M} is a matrix of zeros. }
\note{The MASS function \code{Null(matrix(0,4,2))}  returns a 4*2 matrix, whereas \code{NullC(matrix(0,4,2))} returns \code{diag(4)}. 
}
\author{
Peter Hoff
}

\examples{
NullC(matrix(0,4,2))

## The function is currently defined as
function (M) 
{
    tmp <- qr(M)
    set <- if (tmp$rank == 0L) 
        1L:nrow(M)
    else -(1L:tmp$rank)
    qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
  }
}
