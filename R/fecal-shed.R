
#' Title
#'
#' @param prms List of distribution parameters
#'
#' @return Numerical vector representing the fecal shedding distribution
#' 
#' @export
#'
get_fecalshed <- function(prms = list(m=5, v=3, tmax=14)) {
  y = stats::dgamma(
    x     = 1:prms$tmax, 
    shape = prms$m^2 / prms$v, 
    scale = prms$v/prms$m)
  fec = y/sum(y)
  return(fec)
}
