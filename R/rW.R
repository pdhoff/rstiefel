#' Simulate \code{W} as Described in Wood(1994)
#' 
#' Auxilliary variable simulation for rejection sampling of \code{rmf.vector},
#' as described in Wood(1994).
#' 
#' 
#' @param kap a positive scalar.
#' @param m a positive integer.
#' @return a number between zero and one.
#' @author Peter Hoff
#' @examples
#' 
#' rW(pi,4)
#' 
#' ## The function is currently defined as
#' function (kap, m) 
#' {
#'     .C("rW", kap = as.double(kap), m = as.integer(m), w = double(1))$w
#'   }
#' 
#' @export rW
rW <-
function(kap,m)
{
  #simulate W as described in Wood(1994)
  .C("rWc",kap=as.double(kap),m=as.integer(m),w=double(1))$w
}
