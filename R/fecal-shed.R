
#' Title
#'
#' @param prms 
#'
#' @return
#' @export
#'
#' @examples
get_fecalshed <- function(prms = list(m=5, v=3, tmax=14)) {
  y = dgamma(x=1:prms$tmax, 
             shape = prms$m^2 / prms$v, 
             scale = prms$v/prms$m)
  fec = y/sum(y)
  return(fec)
}
