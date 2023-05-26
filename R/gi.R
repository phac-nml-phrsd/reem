

#' Generation interval distribution
#' 
gi_gamma <- function(m, v, tmax) {
  y = dgamma(0:tmax, shape=m^2/v, scale=v/m)
  res = y / sum(y)
  return(res)
}

#' Title
#'
#' @param prms List of distribution parameters 
#' @param type String. Distribution family (e.g., gamma)
#'
#' @return Numerical vector representing 
#' the intrinsic generation interval distribution.
#' @export
#'
#' @examples
get_gi <- function(prms = list(m=3, v=1, tmax=10), type = 'gamma') {
  
  if(type == 'manual')  g = prms[['values']]
  if(type == 'gamma') 
    g = gi_gamma(m=prms$m, v=prms$v, tmax=prms$tmax)
  
  if( abs(sum(g)-1) > 1e-3 ){
    stop('The generation interval distribution does not sum to 1. Aborting!')
  }
  return(g)
}

