#' Given a matrix M, find a matrix N giving a basis for the null space. This is a modified version of Null from the package MASS.
#'
#' @param M An semi-orthogonal matrix
#'
#' @return A semi-orthogonal matrix N such that t(N) %*% M is a matrix of zeros.
#' @export
#' @note 
#' The MASS function Null(matrix(0,4,2)) returns a 4*2 matrix, whereas NullC(matrix(0,4,2)) returns diag(4).
#' @examples
#' NullC(matrix(0,4,2))
NullC <-
function(M)
{
  #modified from package "MASS" 
  #MASS version : Null(matrix(0,4,2))  returns a 4*2 matrix
  #this version : NullC(matrix(0,4,2)) returns diag(4)

  tmp <- qr(M)
  set <- if (tmp$rank == 0L)
      1L:nrow(M)
  else -(1L:tmp$rank)
  qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
}
