###
###   VARIOUS HELPER FUNCTIONS
###



#' Timestamp as a short string.
#'
#' @return String representing the time now. 
#' @export
#' @importFrom magrittr %>%
#' @examples
timestamp_short <- function() {
  stamp = lubridate::now() %>%
    stringr::str_replace_all('\\s','a') %>%
    stringr::str_replace_all('\\:','\\_')
  return(stamp)
}
