#' @title Optimize a function on the Stiefel manifold
#' 
#' @param F A function V(P, S) -> \code{R^1}
#' @param dF A function to compute the gradient of F
#' @param Vinit The starting point on the stiefel manifold for the optimization
#' @param method: "bb" or curvilinear 
#' @return A stationary point of F on the Stiefel manifold.
#' @references (Wen and Yin, 2013)
#' @examples
#' Find the first eigenspace spanned by the first P eigenvectors for a large matrix M
#' library(rstiefel)
#' 
#' N <- 1000
#' P <- 3
#' Lam <- diag(c(10, 5, 3, rep(1, N-P)))
#' U <- rustiefel(N, N)
#' M <- U %*% Lam %*% t(U)
#' 
#' F <- function(V) { - sum(diag(t(V) %*% M %*% V)) }
#' dF <- function(V) { - 2*M %*% V }
#' V = optStiefel(F, dF, Vinit=rustiefel(N, P),
#'                method="curvilinear",
#'                searchParams=list(rho1=0.1, rho2=0.9, tau=1))
#'                
#' print(sprintf("Sum of first %d eigenvalues is %f", P, -F(V)))
#' 
#' @export
optStiefel <- function(F, dF, Vinit, method="bb",
                       searchParams=NULL,
                       reltol=sqrt(.Machine$double.eps),
                       maxIters=100, verbose=FALSE, 
                       maxLineSearchIters=100) {

    P <- nrow(Vinit)
    S <- ncol(Vinit)
    
    if (method == "bb") {

        if(is.null(searchParams)) {

            ## Default parameters 
            rho <- 0.1
            eta <- 0.2

        } else {
            passedParams <- c("rho", "eta") %in% names(searchParams)
            if(!all(passedParams)) {
                stop(sprintf("For linesearch with BB updates need to specify rho and eta.  Missing %s.",
                             paste(names(passedParams)[!passedParams], collapse=", ")))
            }

            rho = searchParams$rho
            eta = searchParams$eta
        }
        
        Ccur <- Qcur <- 1
        Vprev <- rustiefel(P, S)

    } else if (method == "curvilinear") {

        if(is.null(searchParams)) {

            ## Default parameters
            rho1 <- 0.1
            rho2 <- 0.9
            tau <- 1

        } else {
            passedParams <- c("rho1", "rho2", "tau") %in% names(searchParams)
            if(!all(passedParams)) {
                stop(sprintf("Curvilinear linesearch requires specification of rho1, rho2 and initial stepsize tau.  Missing %s.",
                             paste(names(passedParams)[!passedParams], collapse = ", ")))
            }

            rho1 = searchParams$rho1
            rho2 = searchParams$rho2
            tau = searchParams$tau
        }

        
    }

    V <- Vinit
    Fcur <- F(V)
    Fprev <- Inf
    Gcur <- Gprev <- dF(V)
    iter <- 1
    
    while((Fprev - Fcur) > reltol * (abs(Fcur) + reltol) & iter < maxIters) {
        
        Fprev <- Fcur
        if ( method == "bb") {

            newV <- lineSearchBB(F, V, Vprev, Gcur, Gprev, rho, Ccur, maxIters=maxLineSearchIters)

            Vprev <- V
            V <- newV
            Fcur <- F(V)
            Gprev <- Gcur
            Gcur <- dF(V)

            Qprev <- Qcur
            Qcur <- eta*Qcur + 1
            Ccur <- (eta*Qprev*Ccur + Fcur) / Qcur
            
        } else if ( method == "curvilinear" ) { 

            V <- lineSearch(F, V, Gcur, rho1, rho2, tau, maxIters=maxLineSearchIters)

            Fcur <- F(V)
            Gprev <- Gcur
            Gcur <- dF(V)
        }

        if(verbose) {
            print(sprintf("Iteration %i: %f", iter, Fcur))
        }

        iter <- iter + 1
    }

    V
}


#' A curvilinear search on the Stiefel manifold (Wen and Yin 2013, Algo 1)
#' 
#' @param X an n x p semi-orthogonal matrix (starting point)
#' @param F A function V(n, p) -> \code{R^1}
#' @param G_x an n x p matrix with \code{(G_x)_ij = dF(X)/dX_ij}
#' @return A semi-orthogonal matrix, Ytau, which satisfies Armijo-Wolfe conditions
#'  
#' @references (Wen and Yin, 2013)
#' @export
lineSearch <- function(F, X, G_x, rho1, rho2, tauStart, maxIters=100) {

    n <- nrow(X)
    p <- ncol(X)
    
    reached <- FALSE
    tau <- tauStart
    
    A <- G_x %*% t(X) - X %*% t(G_x)
    U <- cbind(G_x, X)
    V <- cbind(X, -1*G_x)

    ## If tau is too large condition number is too large
    ## and matrix can't be inverted, so reduce the value of initTau until its invertible
    while(kappa(diag(2*p) + tau/2*t(V) %*% U) > 1e12 ) {
        tau <- tau/2
    }
    H <- solve(diag(2*p) + tau/2*t(V) %*% U)
    
    Ytau <- X - tau * U %*% (H %*% t(V) %*% X)
    FprimeY0 <- tr( t(G_x) %*% -A %*% X )
    B <- diag(n) - tau/2*U %*% H %*% t(V)
    FprimeYtau <- tr( t(G_x) %*% -B %*% A %*% (X + Ytau) / 2 )

    ## Check Armijo-Wolfe conditions
    iter <- 0

    while(F(Ytau) > (F(X) + rho1*tau*FprimeY0) | FprimeYtau < rho2*FprimeY0) {

        if(iter > maxIters) {
            print("Reached maximum iterations in line search.")
            break
        }

        tau <- tau/2

        HV <- solve(diag(2*p) + tau/2 * t(V) %*% U) %*% t(V)

        Ytau <- X - tau * U %*% (HV %*% X)
        FprimeY0 <- tr( t(G_x) %*% -A %*% X )
        B <- diag(n) - tau/2*U %*% HV
        FprimeYtau <- tr( t(G_x) %*% -B %*% A %*% (X + Ytau)/2 )
        iter <- iter + 1

    }

    Ytau 
}

#' A curvilinear search on the Stiefel manifold with BB steps (Wen and Yin 2013, Algo 2)
#' This is based on the line search algorithm described in (Zhang and Hager, 2004)
#' 
#' @param X an n x p semi-orthogonal matrix
#' @param F A function V(n, p) -> R
#' @param G_x an n x p matrix with (G_x)_ij = dF(X)/dX_ij
#' @return A semi-orthogonal matrix Ytau which satisfies convergence criter (Eqn 29 in Wen & Yin '13)
#' @references (Wen and Yin, 2013) and (Zhang and Hager, 2004)
#' @export
lineSearchBB <- function(F, X, Xprev, G_x, G_xprev, rho, C, maxIters=100) {

    n <- nrow(X)
    p <- ncol(X)

    A <- G_x %*% t(X) - X %*% t(G_x)
    U <- cbind(G_x, X)
    V <- cbind(X, -1*G_x)

    Sk <- X - Xprev
    Mk <- (G_x - X %*% t(G_x) %*% X) - (G_xprev - Xprev %*% t(G_xprev) %*% Xprev)

    tau <- tr(t(Sk) %*% Sk) / abs(tr(t(Sk) %*%  Mk))
    
    ## If tau is too large condition number is too large
    ## and matrix can't be inverted, so reduce
    while(kappa(diag(2*p) + tau/2*t(V) %*% U) > 1e12 ) {
        tau <- tau/2
    }


    HV <- solve(diag(2*p) + tau/2 * t(V) %*% U) %*% t(V)
    Ytau <- X - tau * U %*% (HV %*% X)
    FprimeY0 <- sum(diag(t(G_x) %*% -A %*% X))

    iter <- 1

    while(F(Ytau) > C + rho*tau*FprimeY0) {

        tau <- tau/2

        if(iter > maxIters) {
            print("Reached maximum iterations in line search.")
            break
        }

        HV <- solve(diag(2*p) + tau/2 * t(V) %*% U) %*% t(V)

        Ytau <- X - tau * U %*% (HV %*% X)
        FprimeY0 <- sum(diag(t(G_x) %*% -A %*% X))
        B <- diag(n) - tau/2*U %*% HV
        FprimeYtau <- tr( t(G_x) %*% -B %*% A %*% (X + Ytau)/2 )

        iter <- iter + 1

    }

    Ytau

}

#' Compute the trace of a matrix
#' @export
tr <- function(X) { sum(diag(X)) }

