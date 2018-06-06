#' Gibbs Sampling for the Matrix-variate Bingham-von Mises-Fisher Distribution.
#' 
#' Simulate a random orthonormal matrix from the Bingham distribution using
#' Gibbs sampling.
#' 
#' 
#' @param A a symmetric matrix.
#' @param B a diagonal matrix with decreasing entries.
#' @param C a matrix with the same dimension as X.
#' @param X the current value of the random orthonormal matrix.
#' @return a new value of the matrix \code{X} obtained by Gibbs sampling.
#' @note This provides one Gibbs scan. The function should be used iteratively.
#' @author Peter Hoff
#' @references Hoff(2009)
#' @examples
#' 
#' ## The function is currently defined as
#' function (A, B, C, X) 
#' {
#'     m <- dim(X)[1]
#'     R <- dim(X)[2]
#'     if (m > R) {
#'         for (r in sample(seq(1, R, length = R))) {
#'             N <- NullC(X[, -r])
#'             An <- B[r, r] * t(N) %*% (A) %*% N
#'             cn <- t(N) %*% C[, r]
#'             X[, r] <- N %*% rbmf.vector.gibbs(An, cn, t(N) %*% 
#'                 X[, r])
#'         }
#'     }
#'     if (m == R) {
#'         for (s in seq(1, R, length = R)) {
#'             r <- sort(sample(seq(1, R, length = R), 2))
#'             N <- NullC(X[, -r])
#'             An <- t(N) %*% A %*% N
#'             Cn <- t(N) %*% C[, r]
#'             X[, r] <- N %*% rbmf.O2(An, B[r, r], Cn)
#'         }
#'     }
#'     X
#'   }
#' 
#' @export rbmf.matrix.gibbs
rbmf.matrix.gibbs <-
function(A,B,C,X)
{
  #simulate from the matrix bmf distribution as described in Hoff(2009) 
  #this is one Gibbs step, and must be used iteratively

  #warning - do not start X at the eigenvectors of A:
  #The relative weights on A then can become infinity. 
  #Instead, you can start close to A, e.g. 
  # X<-rmf.matrix(UA[,1:R]*m)

  m<-dim(X)[1] ;  R<-dim(X)[2]
  if(m>R)
  {
    for(r in sample( seq(1,R,length=R)))
    {
      N<-NullC(X[,-r])
      An<-B[r,r]*t(N)%*%(A)%*%N ; cn<- t(N)%*%C[,r]
      X[,r]<-N%*%rbmf.vector.gibbs(An,cn,t(N)%*%X[,r])
    }
  }

  if(m==R)
  {
    for(s in seq(1,R,length=R))
    {
      r<-sort(sample(seq(1,R,length=R),2))
      N<-NullC( X[,-r]  )
      An<- t(N)%*%A%*%N
      Cn<- t(N)%*%C[,r]
      X[,r]<-N%*%rbmf.O2(An,B[r,r],Cn)
    }
  }

  X
}
