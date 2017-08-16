#' @title Simulate a Uniformly Distributed Random Orthonormal Matrix
#'
#' @description
#' Simulate uniformly from \deqn{V(R, m)} (See Chikuse, 1994).
#' @param m The length of each column vector. 
#' @param R The number of column vectors. 
#'
#' @return An \code{m x R}  uniformly distributed random semi-orthogonal matrix
#' @export
#'
#' @examples
#' rustiefel(10, 2)
rustiefel <- function(m, R) {
  
  # note R is given second, the result is an m * R matrix
  X <- matrix(rnorm(m*R), m, R)
  tmp <- eigen(t(X) %*% X)
  X %*% (tmp$vec %*% sqrt(diag(1/tmp$val, nrow=R)) %*%t(tmp$vec))
}
