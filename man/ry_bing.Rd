\name{ry_bing}
\alias{ry_bing}
\title{
Helper Function for Sampling a Bingham-distributed Vector
}
\description{
C interface to perform a Gibbs update of 
\code{y} with invariant distribution  proportional to 
\code{exp( sum(l*y^2) }                 
with respect to the uniform measure on the sphere. 
}
\usage{
ry_bing(y, l)
}
\arguments{
  \item{y}{a normal vector. 
}
  \item{l}{a vector.
}
}
\value{
a normal vector. 
}
\references{
Hoff(2009)
}
\author{
Peter Hoff
}
\examples{
## The function is currently defined as
function (y, l) 
{
    .C("ry_bing", y = as.double(y), l = as.double(l), n = as.integer(length(y)))$y
  }
}
