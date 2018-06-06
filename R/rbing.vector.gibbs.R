#' Gibbs Sampling for the Vector-variate Bingham Distribution
#' 
#' Simulate a random normal vector from the Bingham distribution using Gibbs
#' sampling.
#' 
#' 
#' @param A a symmetric matrix.
#' @param x the current value of the random normal vector.
#' @return a new value of the vector \code{x} obtained by Gibbs sampling.
#' @note This provides one Gibbs scan. The function should be used iteratively.
#' @author Peter Hoff
#' @references Hoff(2009)
#' @examples
#' 
#' ## The function is currently defined as
#' rbing.vector.gibbs <-
#' function(A,x)
#' {
#'   #simulate from the vector bmf distribution as described in Hoff(2009) 
#'   #this is one Gibbs step, and must be used iteratively
#'   evdA<-eigen(A,symmetric=TRUE)
#'   E<-evdA$vec
#'   l<-evdA$val
#' 
#'   y<-t(E)%*%x
#'   x<-E%*%ry_bing(y,l)
#'   x/sqrt(sum(x^2))
#'   #One improvement might be a rejection sampler 
#'   #based on a mixture of vector mf distributions. 
#'   #The difficulty is finding the max of the ratio.
#' }
#' 
#' 
#' @export rbing.vector.gibbs
rbing.vector.gibbs <-
function(A,x)
{
  #simulate from the vector bmf distribution as described in Hoff(2009) 
  #this is one Gibbs step, and must be used iteratively
  evdA<-eigen(A,symmetric=TRUE)
  E<-evdA$vec
  l<-evdA$val

  y<-t(E)%*%x
  x<-E%*%ry_bing(y,l)
  x/sqrt(sum(x^2))
  #One improvement might be a rejection sampler 
  #based on a mixture of vector mf distributions. 
  #The difficulty is finding the max of the ratio.
}
