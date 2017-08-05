#' Gibbs Sampling for the Matrix-variate Bingham Distribution
#'
#'Simulate a random orthonormal matrix from the Bingham distribution using Gibbs sampling.
#' @useDynLib rstiefel
#' @param A A symmetric matrix 
#' @param B A diagonal matrix with decreasing entries
#' @param X The current value of the random orthonormal Matrix
#'
#' @return A new value the matrix X obtained by Gibbs sampling
#' @export
#'
#' @examples
#' N <- 5; P <- 2
#' X <- rustiefel(N, P)
#' U <- rustiefel(N, N)
#' A <- U %*% t(U)
#' B <- diag(5:1)
#' for(i in 1:100) {
#'    x <- rbing.matrix.gibbs(diag(N:1), diag(N:1), X)
#' }
rbing.matrix.gibbs <- function(A, B, X) {
    
    #simulate from the matrix bmf distribution as described in Hoff(2009) 
    #this is one Gibbs step, and must be used iteratively
    
    ### assumes B is a diagonal matrix with *decreasing* entries 
    
    m<-dim(X)[1] ;  R<-dim(X)[2]
    if(m>R)
    {
      for(r in sample( seq(1,R,length=R)))
      {
        N<-NullC(X[,-r])
        An<-B[r,r]*t(N)%*%(A)%*%N 
        X[,r]<-N%*%rbing.vector.gibbs(An,t(N)%*%X[,r])
      }
    }
    
    #If m=R then the fc of one vector given all the others is 
    #just +-1 times the vector in the null space. In this case, 
    #the matrix needs to be updated at least two columns at a 
    #time. 
    if(m==R)
    {
      for(s in seq(1,R,length=R))
      {
        r<-sort( sample(seq(1,R,length=R),2) )
        N<-NullC( X[,-r]  )
        An<- t(N)%*%A%*%N
        #X[,r]<-N%*%rbing.O2(An,B[r,r]) 
        X[,r]<-N%*%rbing.Op(An,B[r,r]) 
      }
    }
    X
  }

#' Gibbs Sampling for the Vector-variate Bingham Distribution
#' 
#' Simulate a random normal vector from the Bingham distribution using Gibbs sampling.
#' @param A A symmetric matrix.
#' @param x The current value of the normal random vector.
#'
#' @return A new value of the vector \code{x} obtained by Gibbs sampling.
#' 
#' @note This provides one Gibbs scan. The function should be used iteratively.
#' @export
#'
#' @examples
#' N <- 5
#' x <- rustiefel(N, 1)
#' U <- rustiefel(N, N)
#' A <- U %*% t(U)
#' for(i in 1:100) {
#'    X <- rbing.vector.gibbs(A, x)
#' }
rbing.vector.gibbs <- function(A, x) {
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


#' Efficiently simulate a 2x2 Orthogonal Random Matrix
#' 
#' Simulate a 2 x 2 random orthogonal matrix from the Bingham distribution using a rejection sampler.
#' @param A A symmetric matrix.
#' @param B A diagonal matrix with decreasing entries.
#' @param a Sum of the eigenvalues of A, multiplied by the difference in B-values
#' @param E Eigenvectors of A. 
#'
#' @return A random 2x2 orthogonal matrix simulated from the Bingham distribution. 
#' @export
#' 
#'
#' @examples
#' N <- 2;
#' U <- rustiefel(N, N)
#' A <- U %*% t(U)
#' B <- diag(N:1)
#' X <- rbing.O2(A, B)
rbing.O2 <- function(A, B, a=NULL, E=NULL) {
    
    ### assumes B is a diagonal matrix with *decreasing* entries 
    if(is.null(a))
    {
      trA<-A[1,1]+A[2,2]
      lA<-2*sqrt(trA^2/4-A[1,1]*A[2,2]+A[1,2]^2 )
      a<-lA*(B[1,1]-B[2,2])
      E<-diag(2) ; if(A[1,2]!=0)
      {
        E<-cbind(c(.5*(trA+lA)-A[2,2],A[1,2]), c(.5*(trA-lA)-A[2,2],A[1,2]) )
        E[,1]<-E[,1]/sqrt(sum(E[,1]^2)) ; E[,2]<-E[,2]/sqrt(sum(E[,2]^2))
      }
    }
    
    b<-min(1/a^2,.5) ; beta<- .5-b
    lrmx<- a
    if(beta>0) { lrmx<-lrmx + beta*(log(beta/a)-1) }
    lr<- -Inf
    while(lr< log(runif(1)))
    {
      w<-rbeta(1,.5,b)
      lr<- a*w + beta*log(1-w)   - lrmx
    }
    u<-c(sqrt(w),sqrt(1-w) )*(-1)^rbinom(2,1,.5)
    x1<-E%*%u ; x2<-(x1[2:1]*c(-1,1)*(-1)^rbinom(1,1,.5))
    cbind(x1,x2)
  }

#' Simulate a \code{p x p} Orthogonal Random Matrix from the Bingham distribution
#' 
#' Simulate a \code{p x p} random orthogonal matrix from the Bingham distribution on O(p) having density proportional to \deqn{etr(BU^TAU)} using the rejection sampler described in Hoff(2009)
#' This is a rejection sampler and only works for small matrices, otherwise the sampler will reject too frequently    
#'
#' @param A A symmetric matrix.
#' @param B A diagonal matrix with decreasing entries
#'
#' @return A random \code{p x p} orthogonal matrix U simulated from the Bingham distribution.
#' @export
#'
#' @examples
#' N <- 5;
#' U <- rustiefel(N, N)
#' A <- U %*% t(U)
#' B <- diag(N:1)
#' X <- rbing.Op(A, B)
rbing.Op <- function(A, B) {
    
    ### assumes B is a diagonal matrix with *decreasing* entries 
    
    b <-diag(B); bmx <- max(b); bmn <- min(b)  
    if(bmx>bmn)
    { 
      A<-A*(bmx-bmn) ; b<-(b-bmn)/(bmx -bmn)
      vlA<-eigen(A)$val  
      diag(A)<-diag(A)-vlA[1]
      vlA<-eigen(A)$val  
      
      nu<- max(dim(A)[1]+1,round(-vlA[length(vlA)]))
      del<- nu/2
      M<- solve( diag(del,nrow=dim(A)[1] ) - A )/2
      
      rej<-TRUE
      cholM<-chol(M)
      nrej<-0
      while(rej)
      {
        Z<-matrix(rnorm(nu*dim(M)[1]),nrow=nu,ncol=dim(M)[1])
        Y<-Z%*%cholM ; tmp<-eigen(t(Y)%*%Y)
        U<-tmp$vec%*%diag((-1)^rbinom(dim(A)[1],1,.5)) ; L<-diag(tmp$val)
        D<-diag(b)-L
        lrr<- sum(diag(( D%*%t(U)%*%A%*%U)) ) - sum( -sort(diag(-D))*vlA)
        rej<- ( log(runif(1))> lrr )
        nrej<-nrej+1
      }
    }
    if(bmx==bmn) { U<-rustiefel(dim(A)[1],dim(A)[1]) } 
    U
}

#' Helper Function for Sampling a Bingham-distributed Vector
#' 
#' C interface to perform a Gibbs update of 
#' \code{y} with invariant distribution  proportional to 
#' \deqn{\exp(sum(l*y^2))} with respect to the uniform measure on the sphere. 
#' @param y A normal vector
#' @param l A vector
#'
#' @return A normal vector
ry_bing <- function(y, l) {
    .C("ry_bing",y=as.double(y),l=as.double(l),n=as.integer(length(y)))$y
}
