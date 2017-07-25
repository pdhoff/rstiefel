#' Gibbs Sampling for the Vector-variate Bingham-von Mises-Fisher Distribution
#'
#' Simulate a random normal vector from the Bingham-von Mises-Fisher distribution using Gibbs sampling as described in Hoff(2009) as described in Hoff(2009) as described in Hoff(2009) .
#' @param A A symmetric matrix
#' @param c A vector with the same length as \code{x}
#' @param x The current value of the random normal vector.
#' @note This is one Gibbs step, and must be used iteratively.
#'
#' @return A new value of the vector \code{x} obtained by Gibbs sampling.
#' @export
#'
#' @examples
rbmf.vector.gibbs <- function(A, c, x) {
    

    evdA<-eigen(A)
    E<-evdA$vec
    l<-evdA$val
    
    y<-t(E)%*%x
    d<-t(E)%*%c
    x<-E%*%ry_bmf(y,l,d)
    x/sqrt(sum(x^2))
  }


#' Gibbs Sampling for the Matrix-variate Bingham-von Mises-Fisher Distribution.
#' Simulate from the matrix bmf distribution as described in Hoff(2009).  
#'
#' @param A A symmetric matrix.
#' @param B A diagonal matrix with decreasing entries.
#' @param C A matrix with the same dimension as X. 
#' @param X The current value of the random orthonormal matrix.
#'
#' @return A new value of the matrix \code{X} obtained by Gibbs sampling.
#' @export
#'
#' @note This provides one Gibbs scan. The function should be used iteratively.    
#' WARNING: - do not start X at the eigenvectors of A: The relative weights on A then can become infinity. 
#'
#' @examples
rbmf.matrix.gibbs <- function(A, B, C, X) {

    m <- dim(X)[1];  R <- dim(X)[2]
    if(m > R)
    {
      for(r in sample( seq(1,R,length=R)))
      {
        N <- NullC(X[,-r])
        An <- B[r,r]*t(N) %*% (A) %*% N; cn <- t(N) %*% C[,r]
        X[,r] <- N %*% rbmf.vector.gibbs(An, cn, t(N) %*% X[,r])
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

#' Simulate a \code{2 x 2} Orthogonal Random Matrix
#'
#' @param A A symmetric matrix.
#' @param B A diagonal matrix with decreasing entries.
#' @param C A \code{2 x 2} matrix.
#' @param env 
#'
#' @return A random \code{2 x 2} orthogonal matrix simulated from the Bingham-von Mises-Fisher distribution.
#' @export
#'
#' @examples
rbmf.O2 <- function(A, B, C, env=FALSE) {

    sC<-svd(C)
    d1<-sum(sC$d)
    eA<-eigen(A) 
    ab<-sum(eA$val*diag(B)) 
    
    ### if Bingham part dominates, use Bingham envelope
    if(d1<=ab | env=="bingham")   
    {
      lrmx<-sum(sC$d) ; lr<- -Inf
      while(lr<log(runif(1)))
      {      
        X<-rbing.O2(A,B,a=(eA$val[1]-eA$val[2])*(B[1,1]-B[2,2]),E=eA$vec)
        lr<-sum(diag(t(X)%*%C)) - lrmx
      }
    }   
    ###
    
    ### if MF part dominates, use MF envelope
    if(d1>ab | env=="mf")  
    {
      lrmx<-sum(eA$val*sort(diag(B),decreasing=TRUE)) ; lr<- -Inf  
      while(lr< log(runif(1)))
      {
        X<-rmf.matrix(C)
        lr<-sum(diag(B%*%t(X)%*%A%*%X)) - lrmx
      }
    }
    ###
    
    X
}

#' Helper Function for Sampling a Bingham-von Mises-Fisher-distributed Vector
#'
#' @param y A normal vector.
#' @param l A vector
#' @param d A vector 
#'
#' @return A normal vector
#' @export
#'
#' @examples
ry_bmf <- function(y, l, d) {
  
    .C("ry_bmf",y=as.double(y), l=as.double(l), d=as.double(d),
      n=as.integer(length(y)))$y
}
