#' Helper Function for Sampling a Bingham-von Mises-Fisher-distributed Vector
#' 
#' C interface to perform a Gibbs update of \code{y} with invariant
#' distribution proportional to \code{exp( sum(l*y^2+y*d) } with respect to the
#' uniform measure on the sphere.
#' 
#' 
#' @param y a normal vector.
#' @param l a vector.
#' @param d a vector.
#' @return a normal vector
#' @author Peter Hoff
#' @references Hoff(2009)
#' @examples
#' 
#' ## The function is currently defined as
#' function (y, l, d) 
#' {
#'     .C("ry_bmf", PACKAGE="rstiefel", y = as.double(y), l = as.double(l), d = as.double(d), 
#'         n = as.integer(length(y)))$y
#'   }
#' 
#' @export ry_bmf
ry_bmf <-
function(y,l,d)
{
  .C("ry_bmfc",PACKAGE="rstiefel",y=as.double(y),l=as.double(l),d=as.double(d),
               n=as.integer(length(y)))$y
}
