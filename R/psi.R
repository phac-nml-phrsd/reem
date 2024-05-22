#' @title Pathogen decay in wastewater
#'
#' @param prm List of distribution parameters 
#'
#' @return Numerical vector representing  the distribution
#' of viral decay in wastewater.
#' @export
#'
#' @examples
#' 
get_psi <- function(prm) {
  return(c(0.85, 0.10, 0.05))
}
