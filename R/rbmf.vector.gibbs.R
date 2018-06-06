#' Gibbs Sampling for the Vector-variate Bingham-von Mises-Fisher Distribution
#' 
#' Simulate a random normal vector from the Bingham-von Mises-Fisher
#' distribution using Gibbs sampling.
#' 
#' 
#' @param A a symmetric matrix.
#' @param c a vector with the same length as \code{x}.
#' @param x the current value of the random normal vector.
#' @return a new value of the vector \code{x} obtained by Gibbs sampling.
#' @note This provides one Gibbs scan. The function should be used iteratively.
#' @author Peter Hoff
#' @references Hoff(2009)
#' @examples
#' 
#' ## The function is currently defined as
#' function (A, c, x) 
#' {
#'     evdA <- eigen(A)
#'     E <- evdA$vec
#'     l <- evdA$val
#'     y <- t(E) %*% x
#'     d <- t(E) %*% c
#'     x <- E %*% ry_bmf(y, l, d)
#'     x/sqrt(sum(x^2))
#'   }
#' 
#' @export rbmf.vector.gibbs
rbmf.vector.gibbs <-
function(A,c,x)
{
  #simulate from the vector bmf distribution as described in Hoff(2009) 
  #this is one Gibbs step, and must be used iteratively
  evdA<-eigen(A)
  E<-evdA$vec
  l<-evdA$val

  y<-t(E)%*%x
  d<-t(E)%*%c
  x<-E%*%ry_bmf(y,l,d)
  x/sqrt(sum(x^2))
}
