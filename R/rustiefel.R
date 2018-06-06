#' Siumlate a Uniformly Distributed Random Orthonormal Matrix
#' 
#' Siumlate a random orthonormal matrix from the uniform distribution on the
#' Stiefel manifold.
#' 
#' 
#' @param m the length of each column vector.
#' @param R the number of column vectors.
#' @return an \code{m*R} orthonormal matrix.
#' @author Peter Hoff
#' @references Hoff(2007)
#' @examples
#' 
#' ## The function is currently defined as
#' function (m, R) 
#' {
#'     X <- matrix(rnorm(m * R), m, R)
#'     tmp <- eigen(t(X) %*% X)
#'     X %*% (tmp$vec %*% sqrt(diag(1/tmp$val, nrow = R)) %*% t(tmp$vec))
#'   }
#' 
#' @export rustiefel
rustiefel <-
function(m,R)
{
  #simulate uniformly from V_{R,m}  
  #see Chikuse 
  #note R is given second, the result is an m*R matrix
  X<-matrix(rnorm(m*R),m,R)
  tmp<-eigen(t(X)%*%X)
  X%*%( tmp$vec%*%sqrt(diag(1/tmp$val,nrow=R))%*%t(tmp$vec) )
}
