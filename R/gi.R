

#' @title Generation interval as a Gamma distribution
#' 
#' @param m Numerical. Mean.
#' @param v Numerical. Variance
#' @param tmax Numerical. Maximum value (support)
#' 
#' @return A numerical vector representing the probability density (sums to 1)
#' 
gi_gamma <- function(m, v, tmax) {
  y = dgamma(0:tmax, shape=m^2/v, scale=v/m)
  res = y / sum(y)
  return(res)
}

#' @title Generation interval distribution definition.
#'
#' @param prms List of distribution parameters.
#' The names of elements depends on \code{type}. 
#' If \code{type = "manual"}, then only one element named \code{values} is 
#' expected. The values must sum to 1. 
#'  If \code{type = "gamma"}, then elements named \code{m, v} 
#'  and \code{tmax} are expected to respectively define the mean, variance and maximum
#'  values of the Gamma distribution. 
#' @param type String. Distribution family. Implemented: 
#' \code{manual} to define the distribution manually, 
#' \code{gamma} for a Gamma shape distribution
#'
#' @return Numerical vector representing the probability densities of
#' the intrinsic generation interval distribution.
#' 
#' @export
#'
#' @examples
#' 
#' # Defined as Gamma
#' gg = get_gi(type = 'gamma', prms = list(m=5, v=2, tmax = 12)) 
#' plot(gg, typ = 'b')
#' 
#' # Defined manually
#' y = c(0.01, 0.1, 0.5, 0.2, 0.15, 0.04)
#' gm = get_gi(type = 'manual', prms = list(values = y/sum(y)))
#' plot(gm, typ = 'b')
#' 
#' 
#' 
get_gi <- function(prms = list(m=3, v=1, tmax=10), type = 'gamma') {
  
  if(type == 'manual')  g = prms[['values']]
  if(type == 'gamma') 
    g = gi_gamma(m=prms$m, v=prms$v, tmax=prms$tmax)
  
  if( abs(sum(g)-1) > 1e-3 ){
    stop('The generation interval distribution does not sum to 1. Aborting!')
  }
  return(g)
}

