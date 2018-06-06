#' Gibbs Sampling for the Matrix-variate von Mises-Fisher Distribution
#' 
#' Simulate a random orthonormal matrix from the matrix von Mises-Fisher
#' distribution using Gibbs sampling.
#' 
#' 
#' @param M a matrix.
#' @param X the current value of the random orthonormal matrix.
#' @param rscol the number of columns to update simultaneously.
#' @return a new value of the matrix \code{X} obtained by Gibbs sampling.
#' @note This provides one Gibbs scan. The function should be used iteratively.
#' @author Peter Hoff
#' @references Hoff(2009)
#' @examples
#' 
#' Z<-matrix(rnorm(10*5),10,5) 
#' 
#' U<-rmf.matrix(Z)
#' U<-rmf.matrix.gibbs(Z,U)
#' 
#' 
#' ## The function is currently defined as
#' function (M, X, rscol = NULL) 
#' {
#'     if (is.null(rscol)) {
#'         rscol <- max(2, min(round(log(dim(M)[1])), dim(M)[2]))
#'     }
#'     sM <- svd(M)
#'     H <- sM$u %*% diag(sM$d)
#'     Y <- X %*% sM$v
#'     m <- dim(H)[1]
#'     R <- dim(H)[2]
#'     for (iter in 1:round(R/rscol)) {
#'         r <- sample(seq(1, R, length = R), rscol)
#'         N <- NullC(Y[, -r])
#'         y <- rmf.matrix(t(N) %*% H[, r])
#'         Y[, r] <- N %*% y
#'     }
#'     Y %*% t(sM$v)
#'   }
#' 
#' @export rmf.matrix.gibbs
rmf.matrix.gibbs <-
function(M,X,rscol=NULL)
{
  #simulate from the matrix mf distribution as described in Hoff(2009) 
  #this is one Gibbs step, and must be used iteratively

  #The number of columns to replace should be small enough so that 
  #rmf.matrix will be quick, but not so small that you are not taking 
  #advantage of rmf.matrix. In particular, for square matrices you must 
  #be replacing two or more columns at a time, or else the chain is 
  #not irreducible. 

  if(is.null(rscol)){rscol<-max(2,min( round(log(dim(M)[1])),dim(M)[2])) }

  sM<-svd(M)
  H<-sM$u%*%diag(sM$d)
  Y<-X%*%sM$v

  m<-dim(H)[1] ; R<-dim(H)[2]
  for(iter in 1:round(R/rscol))
  {
    r<-sample(seq(1,R,length=R),rscol)
    N<-NullC(Y[,-r])
    y<-rmf.matrix(t(N)%*%H[,r])
    Y[,r]<- N%*%y
  }
  Y%*%t(sM$v)
}
