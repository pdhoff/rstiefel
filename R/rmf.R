#' Simulate a Random Orthonormal Matrix
#'
#' Simulate a random orthonormal matrix from the von Mises-Fisher distribution. 
#' @param M A matrix parameterizing the von Mises-Fisher, e.g. p(X) ~ etr(M^TX)
#'
#' @return An orthonormal matrix of the same dimension as \code{M}. 
#' @export
#'
#' @examples
#' Z <- matrix(rnorm(10*5), 10, 5) 
#'
#' U <- rmf.matrix(Z)
#' U <- rmf.matrix.gibbs(Z, U)
rmf.matrix <- function(M) {
  
  if(dim(M)[2]==1) { X<-rmf.vector(M) } 
  if(dim(M)[2]>1) 
  {
    #simulate from the matrix mf distribution using the rejection 
    #sampler as described in Hoff(2009)
    svdM<-svd(M)
    H<-svdM$u%*%diag(svdM$d)
    m<-dim(H)[1] ; R<-dim(H)[2]
    
    cmet<-FALSE
    rej<-0
    while(!cmet)
    {
      U<-matrix(0,m,R)
      U[,1]<-rmf.vector(H[,1])
      lr<-0
      
      for(j in seq(2,R,length=R-1))
      {
        N<-NullC(U[,seq(1,j-1,length=j-1)])
        x<-rmf.vector(t(N)%*%H[,j])
        U[,j]<- N%*%x
        
        if(svdM$d[j]>0) 
        {
          xn<- sqrt(sum( (t(N)%*%H[,j])^2))
          xd<- sqrt(sum( H[,j]^2 ))
          lbr<-  log(besselI(xn, .5*(m-j-1),expon.scaled=TRUE))-
            log(besselI(xd, .5*(m-j-1),expon.scaled=TRUE))
          if(is.na(lbr)){lbr<- .5*(log(xd) - log(xn)) }
          lr<- lr+ lbr + (xn-xd) + .5*(m-j-1)*( log(xd)-log(xn) )
        }
      }
      
      cmet<- (log(runif(1)) <  lr ) ; rej<-rej+(1-1*cmet)
    }
    X<-U%*%t(svd(M)$v)
  }
  X
}

#' Gibbs Sampling for the Matrix-variate Bingham-von Mises-Fisher Distribution.
#'
#' Simulate a random orthonormal matrix from the matrix von Mises-Fisher distribution using Gibbs sampling, Hoff(2009). 
#'
#' @param M A matrix
#' @param X The current value of the random orthonormal matrix.
#' @param rscol The number of columns to update simultaneously. 
#'
#' @return A new value of the matrix \code{X} obtained by Gibbs sampling.
#' @export
#'
#' @note This provides one Gibbs scan. The function should be used iteratively. The number of columns to replace should be small enough so that 
#' rmf.matrix will be quick, but not so small that you are not taking 
#' advantage of rmf.matrix. In particular, for square matrices you must 
#' be replacing two or more columns at a time, or else the chain is 
#' not irreducible. 
#' @examples
#' Z <- matrix(rnorm(10*5), 10, 5) 
#'
#' U <- rmf.matrix(Z)
#' U <- rmf.matrix.gibbs(Z, U)
rmf.matrix.gibbs <- function(M, X, rscol=NULL) {
    
  
    
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


#' Simulate a Random Normal Vector
#'
#' @param kmu A vector. 
#'
#' @return A vector. 
#' @export
#'
#' @references Wood(1994), Hoff(2009)
#' @examples
#' N <- 100
#' kmu <- rustiefel(N, 1)
#' rmf.vector(kmu)
rmf.vector <- function(kmu) {
  
    #simulate from the vector mf distribution as described in Wood(1994)
    kap<-sqrt(sum(kmu^2)) ; mu<-kmu/kap ; m<-length(mu)
    if(kap==0){ u<-rnorm(length(kmu)) ; u<-matrix(u/sqrt(sum(u^2)),m,1) }
    if(kap>0)
    {
      if(m==1){ u<- (-1)^rbinom( 1,1,1/(1+exp(2*kap*mu))) }
      if(m>1)
      {
        W<-rW(kap,m)
        V<-rnorm(m-1) ;  V<-V/sqrt(sum(V^2))
        x<-c((1-W^2)^.5*t(V),W)
        u<-cbind( NullC(mu),mu)%*%x
      }
    }
    u
  }

#' Simulate \code{W} as Described in Wood(1994)
#' 
#' Auxilliary variable simulation for rejection sampling of \code{rmf.vector}, 
#' as described in Wood(1994).
#' @param kap A positive scalar.
#' @param m A positive integer.
#' 
#' @return A number between zero and one. 
#' @export
#'
#' @examples
#' rW(pi, 4)
rW <- function(kap, m) {
  
    #simulate W as described in Wood(1994)
    .C("rW",kap=as.double(kap),m=as.integer(m),w=double(1))$w
}
