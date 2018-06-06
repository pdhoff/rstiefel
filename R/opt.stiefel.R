#' @title Optimize a function on the Stiefel manifold
#'
#' @description
#' Find a local optimum of a function defined on the stiefel manifold using algorithms described in Wen and Yin (2013).
#' @param F A function V(P, S) -> \code{R^1}
#' @param dF A function to compute the gradient of F.  Returns a \code{P x S} matrix with \code{dF(X)_ij  = d(F(X))/dX_ij}
#' @param Vinit The starting point on the stiefel manifold for the optimization
#' @param method Line search type: "bb" or curvilinear
#' @param searchParams List of parameters for the line search algorithm.  If the line search algorithm is the standard curvilinear search than the search parameters are rho1 and rho2.  If the line search algorithm is "bb" then the parameters are rho and eta.
#' @param tol Convergence tolerance.  Optimization stops when Fprime < abs(tol), an approximate stationary point.  
#' @param maxIters Maximum iterations for each gradient step
#' @param maxLineSearchIters Maximum iterations for for each line search (one step in the gradient descent algorithm)
#' @param verbose Boolean indicating whether to print function value and iteration number at each step.
#' 
#' @return A stationary point of F on the Stiefel manifold.
#' @author Alexander Franks
#' @references (Wen and Yin, 2013)
#' @examples
#' ## Find the first eigenspace spanned by the first P eigenvectors for a 
#' ## matrix M. The size of the matrix has been kept small and the tolerance 
#' ## has been raised to keep the runtime 
#' ## of this example below the CRAN submission threshold. 
#' 
#' N <- 500
#' P <- 3
#' Lam <- diag(c(10, 5, 3, rep(1, N-P)))
#' U <- rustiefel(N, N)
#' M <- U %*% Lam %*% t(U)
#' 
#' F <- function(V) { - sum(diag(t(V) %*% M %*% V)) }
#' dF <- function(V) { - 2*M %*% V }
#' V = optStiefel(F, dF, Vinit=rustiefel(N, P),
#'                method="curvilinear",
#'                searchParams=list(rho1=0.1, rho2=0.9, tau=1),tol=1e-4)
#'                
#' print(sprintf("Sum of first %d eigenvalues is %f", P, -F(V)))
#' 
#' @export
optStiefel <- function(F, dF, Vinit, method="bb",
                       searchParams=NULL,
                       tol=1e-5,
                       maxIters=100, verbose=FALSE, 
                       maxLineSearchIters=20) {

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
        
        Qcur <- 1
        Ccur <- F(Vinit)
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
    Fprime <- Inf

    ## While ||gradF(V)|| > eps
    while(abs(Fprime) > tol & iter < maxIters) {
        
        Fprev <- Fcur
        if ( method == "bb") {

            res <- lineSearchBB(F, V, Vprev, Gcur, Gprev, rho, Ccur, maxIters=maxLineSearchIters)

            Vprev <- V
            V <- res$Y
            Fcur <- F(V)
            Gprev <- Gcur
            Gcur <- dF(V)

            Qprev <- Qcur
            Qcur <- eta*Qcur + 1
            Ccur <- (eta*Qprev*Ccur + Fcur) / Qcur
            Gcur <- dF(V)
            
        } else if ( method == "curvilinear" ) { 

            res <- lineSearch(F, dF, V, rho1, rho2, tau, maxIters=maxLineSearchIters)
            V <- res$Y
            tau <- res$tau

            Fcur <- F(V)
            Gprev <- Gcur
            Gcur <- dF(V)
        }

        if(verbose) {
            print(sprintf("Iteration %i: F = %f, dF = %f", iter, Fcur, Fprime))
        }

        ## compute ||gradF(V)||
        A <- Gcur %*% t(V) - V %*% t(Gcur)
        Fprime <- tr( t(Gcur) %*% -A %*% V )
        
        iter <- iter + 1
    }

    V
}


#' A curvilinear search on the Stiefel manifold (Wen and Yin 2013, Algo 1)
#' 
#' @param X an n x p semi-orthogonal matrix (starting point)
#' @param dF A function V(n, p) -> \code{R^1}
#' @param F A function V(n, p) -> \code{R^1}
#' @param rho1 Parameter for Armijo condition.  Between 0 and 1 and usually small, e.g < 0.1
#' @param rho2 Parameter for Wolfe condition Between 0 and 1 usually large, > 0.9
#' @param tauStart Initial step size
#' @param maxIters Maximum number of iterations
#' 
#' @return A list containing: Y, the semi-orthogonal matrix satisfying the Armijo-Wolfe conditions and tau: the stepsize satisfying these conditions
#'
#' @examples
#' N <- 10
#' P <- 2
#' M <- diag(10:1)
#' F <- function(V) { - sum(diag(t(V) %*% M %*% V)) }
#' dF <- function(V) { - 2*M %*% V }
#' X <- rustiefel(N, P)
#' res <- lineSearch(F, dF, X, rho1=0.1, rho2=0.9, tauStart=1)
#' 
#' @author Alexander Franks
#' @references (Wen and Yin, 2013)
#' @export
lineSearch <- function(F, dF, X, rho1, rho2, tauStart, maxIters=20) {

    G_x <- dF(X)
    
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

    HV <- solve(diag(2*p) + tau/2 * t(V) %*% U) %*% t(V)
    
    Ytau <- X - tau * U %*% (HV %*% X)
    FprimeY0 <- tr( t(G_x) %*% -A %*% X )

    B <- diag(n) - tau/2*U %*% HV
    FprimeYtau <- tr( t(dF(Ytau)) %*% -B %*% A %*% (X + Ytau) / 2 )

    ## Check Armijo-Wolfe conditions
    iter <- 0

    lower <- 0
    upper <- Inf
    
    Armijo <- F(Ytau) > (F(X) + rho1*tau*FprimeY0)
    Wolfe <- FprimeYtau < rho2*FprimeY0

    while(Armijo | Wolfe) {

        if(iter > maxIters) {
            print("Reached maximum iterations in line search.")
            break
        }

        if(Armijo) {
            upper <- tau
            tau <- (lower + upper) / 2
        } else if(Wolfe) {
            lower <- tau
            if(upper == Inf)
                tau <- 2 * lower
            else
                tau <- (lower + upper) / 2
        }

        HV <- solve(diag(2*p) + tau/2 * t(V) %*% U) %*% t(V)
        Ytau <- X - tau * U %*% (HV %*% X)
        ## B <- diag(n) - tau/2*U %*% HV
        B <- solve(diag(n) + tau/2 * A)
        FprimeYtau <- tr( t(dF(Ytau)) %*% -B %*% A %*% (X + Ytau) / 2 )
        
        Armijo <- F(Ytau) > (F(X) + rho1*tau*FprimeY0)
        Wolfe <- FprimeYtau < rho2*FprimeY0
        iter <- iter + 1
        ## print(F(X))

    }

    list(Y=Ytau, tau=tau)
}

#' A curvilinear search on the Stiefel manifold with BB steps (Wen and Yin 2013, Algo 2)
#' This is based on the line search algorithm described in (Zhang and Hager, 2004)
#' @param F A function V(n, p) -> R
#' @param X an n x p semi-orthogonal matrix (the current )
#' @param Xprev an n x p semi-orthogonal matrix (the previous)
#' @param G_x an n x p matrix with (G_x)_ij = dF(X)/dX_ij
#' @param G_xprev an n x p matrix with (G_xprev)_ij = dF(X_prev)/dX_prev_ij
#' @param rho Convergence parameter, usually small (e.g. 0.1)
#' @param C C_t+1 = (etaQ_t + F(X_t+1))/Q_t+1 See section 3.2 in Wen and Yin, 2013
#' @param maxIters Maximum number of iterations
#' 
#' @return A list containing Y: a semi-orthogonal matrix Ytau which satisfies convergence criteria (Eqn 29 in Wen & Yin '13), and tau: the stepsize satisfying these criteria
#'
#' @examples
#' N <- 10
#' P <- 2
#' M <- diag(10:1)
#' F <- function(V) { - sum(diag(t(V) %*% M %*% V)) }
#' dF <- function(V) { - 2*M %*% V }
#' Xprev <- rustiefel(N, P)
#' G_xprev <- dF(Xprev)
#' X <- rustiefel(N, P)
#' G_x <- dF(X)
#' Xprev <- dF(X)
#' res <- lineSearchBB(F, X, Xprev, G_x, G_xprev, rho=0.1, C=F(X))
#' 
#' @references (Wen and Yin, 2013) and (Zhang and Hager, 2004)
#' @author Alexander Franks
#' @export
lineSearchBB <- function(F, X, Xprev, G_x, G_xprev, rho, C, maxIters=20) {

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

        iter <- iter + 1

    }

    list(Y=Ytau, tau=tau)

}

#' @title Compute the trace of a matrix
#' @description compute the trace of a square matrix
#' @param X Square matrix
#' @export
tr <- function(X) { sum(diag(X)) }

