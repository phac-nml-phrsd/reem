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
    stringr::str_replace_all('\\:','\\_') %>% 
    stringr::str_remove('.\\d+$')
  return(stamp)
}



#' Translate the initial PROPORTION 
#' of infected individuals into a
#' NuMBER of individuals initially infected.
#'
#' @param prm List. Model parameters.
#'
#' @return
#'
#' @examples
set_I_init <- function(prm) {
  i0val = round(prm$N * 10^prm$i0prop)
  prm$I.init <- rep(i0val, prm$lag)
  return(prm)
}


#' @title Mean incidence from the renewal equation.
#'
#' @param t Numerical. Time.
#' @param R0 Numerical. Basic reproduction number.
#' @param B Numerical vector. Multiplicative factor changing the contact rate.
#' (if no change, set \code{B[t]=1} for all t)
#' @param S Numerical for the number of susceptible.
#' @param N Numerical. Total population size
#' @param alpha Numerical. Parameter representing mixing heterogeneity.
#' \code{alpha = 0} means homogeneous mixing.
#' @param g Numerical vector representing the intrinsic generation interval distribution.
#' @param I Numerical vector. Incidence at each time step.
#'
#' @return
#'
#' @examples
#' 
mean_inc <- function(t, R0, B, S, N, alpha, g, I) {
  
  #stopifnot(t>1)
  
  # work with existing data
  n = sum(!is.na(I))
  #stopifnot(n>0)
  
  tmp1 = R0 * B[t] * (S[t-1]/N)^(exp(alpha))  
  
  revI = I[n:1]  # faster than using `rev()`
  n2   = min(n,length(g))
  tmp2 = g[1:n2] * revI[1:n2] 
  
  m = tmp1 * sum(tmp2)
  return(m)
}

calc_I <- function(lambda, deterministic){
  res = lambda
  if(!deterministic) res = rpois(n=1, lambda = lambda)
  return(res)
}

calc_Y <- function(lambdaY,n, deterministic){
  res = lambdaY
  if(!deterministic) res = rpois(n=n, lambda = lambdaY)
  return(res)
}



#' Check if the parameter should be an integer
#'
#' @param prm.name String. Name of the parameter
#'
#' @return Logical
#' 
is_prm_integer <- function(prm.name) {
  res = FALSE
  namelist = c('start.delta', 'foo')
  if(prm.name %in% namelist) 
    res = TRUE
  return(res)
}



#' Aggregate values across time 
#'
#' @param df Dataframe oftemporal data. 
#' Must have a variable named 
#' either \code{date} or \`code{time}. 
#' @param dt.aggr Vector of integer (time) or dates (date)
#' were the aggregation is performed.
#' @param var.name String. Name of the variable to be aggregated.
#'
#' @return A dataframe where the values of the 
#' variable \code{var.name} has been aggregated (summed).
#' 
#' @export
#'
#' @examples
#' 
aggregate_time <- function(df, 
                           dt.aggr, 
                           var.name,
                           aggreg.name = NULL) {
  
  if(0){ # DEBUG
    dt.aggr = obs.cl$date
    dt.aggr = obs.cl$t
    var.name = 'Y'
  }
  
  # Determine if we work with times or dates
  dt = 't'
  if(class(dt.aggr) == 'Date') dt = 'date'
  
  # Generic variable names
  df$zzz <- df[[var.name]]
  df$dt  <- df[[dt]]
  
  # Calculate aggregation
  
  # define the indices of each group
  idx.tmp = df$dt %in% dt.aggr 
  if(sum(idx.tmp)==0) {
    stop('The aggregation times/dates are not found in the simulation times.\n',
         '(clinical observations dates and/or start date are likely mispecified)')
  }
  idx     = cumsum(idx.tmp) + 1   # `+1` because we don't want index `0` in vectors 
  # need to shift by one 
  # to aggregate until `dt.agrr`
  idx = c(1, idx[-length(idx)])
 
  # aggregation
  tmp = data.frame(dt = df$dt, zzz = df$zzz) %>% 
    dplyr::arrange(dt) %>% 
    dplyr::mutate(group = idx) %>%
    dplyr::group_by(group) %>% 
    dplyr::summarise(
      aggregation = sum(zzz), 
      dt = max(dt))
  
  # Rename appropriately
  res = dplyr::select(tmp, dt, aggregation) 
  names(res)[1] <- dt
  
  if(!is.null(aggreg.name)) names(res)[2] <- aggreg.name
  
  return(res) 
}

#' Helper function
#' 
aggcl <- function(df, dt.aggr, vars ) {
  tmp = list()
  for(v in vars){
    a = aggregate_time(df = df, 
                       dt.aggr = dt.aggr, 
                       var.name = v, 
                       aggreg.name = v) 
    d = a$date
    tmp[[v]] =  select(a, -date)
  }
  res = tmp %>% 
    bind_cols( ) %>% 
    mutate(date = d) 
  return(res)
}
