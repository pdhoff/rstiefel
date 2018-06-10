#' Helper Function for Sampling a Bingham-distributed Vector
#' 
#' C interface to perform a Gibbs update of \code{y} with invariant
#' distribution proportional to \code{exp( sum(l*y^2) } with respect to the
#' uniform measure on the sphere.
#' 
#' 
#' @param y a normal vector.
#' @param l a vector.
#' @return a normal vector.
#' @author Peter Hoff
#' @references Hoff(2009)
#' @examples
#' 
#' ## The function is currently defined as
#' function (y, l) 
#' {
#'     .C("ry_bing", y = as.double(y), l = as.double(l), n = as.integer(length(y)))$y
#'   }
#' 
#' @export ry_bing
ry_bing <-
function(y,l)
{
  #C function to perform a Gibbs update of 
  #y with invariant distribution 
  #p(y|l) \propto \exp{ sum(l*y^2) } 
  #with respect to the uniform measure on the sphere
  .C("ry_bingc",y=as.double(y),l=as.double(l),n=as.integer(length(y)))$y
}
