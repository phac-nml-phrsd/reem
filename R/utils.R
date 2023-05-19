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



#' Title
#'
#' @param prm 
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
  
  revI = rev(I[1:n])
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



#' Error function for a fit algorithm
#'
#' @param cl.i Dataframe of simulated clinical observations
#' @param ww.i Dataframe of simulated wastewater observations
#' @param obs.cl Dataframe of clinical observations
#' @param obs.ww Dataframe of wastewater observations
#' @param use.cl Logical. Use clinical data in fit?
#' @param use.ww Logical. Use wastewater data in fit?
#'
#' @return Numerical value representing the distance 
#' of the simulation from the observed data.
#' 
err_fct <- function(cl.i, ww.i, 
                    obs.cl, obs.ww,
                    use.cl, use.ww,
                    err.type = 'rel',
                    do.plot = FALSE) {
  
  # --- Checks  
  
  if(!all(cl.i$date %in% obs.cl$date)){
    print('observed cl:')
    print(obs.cl)
    print('simulated cl (ABC internal):')
    print(cl.i)
    stop('ERROR : date mismatch for clinical data')
  }
  if(0){ # TODO: reactivate once ww included
    if(!all(ww.i$date %in% obs.ww$date)){
      print('observed ww:')
      print(obs.ww)
      print('simulated ww (ABC internal):')
      print(ww.i)
      stop('ERROR : date mismatch for ww')
    }
  }
  
  # --- Adjust to fitting dates
  obs.cla = obs.cl
  obs.wwa = obs.ww
  if(nrow(cl.i) < nrow(obs.cl)){
    obs.cla = obs.cl[obs.cl$date %in% cl.i$date,]
  }
  if(nrow(ww.i) < nrow(obs.ww)){
    obs.wwa = obs.ww[obs.ww$date %in% ww.i$date,]
  }
  
  
  if(do.plot){
    g = obs.cl %>% 
      ggplot(aes(x=date))+
      geom_point(aes(y=Y))+
      geom_point(data = obs.cla, aes(y=Y), color='red1', size=2)+
      geom_line(data = obs.cla, aes(y=Y), color='red1',linetype='dashed')+
      geom_line(data=cl.i, aes(y=Ym)) 
    plot(g)
  }
  
  err.cl = 0
  err.ww = 0
  
  if(err.type == 'L2'){
    if(use.cl) err.cl = sqrt(sum( (obs.cla$obs - cl.i$Ym )^2 ))
    if(use.ww) err.ww = sqrt(sum( (obs.wwa$obs - ww.i$Wm )^2 ))
  }
  if(err.type == 'rel'){
    if(use.cl) err.cl = sqrt(sum( (cl.i$Ym/obs.cla$Y - 1 )^2 )) #TODO: handle div by 0
    if(use.ww) err.ww = sqrt(sum( (ww.i$Wm/obs.wwa$W - 1 )^2 ))
  }
  
  # WARNING -- FIXME??
  # The errors calculated above have a number of elements
  # equal to nrow(cl.i), which varies as delta.start changes.
  # Hence, the error may not be consistent across all priors.
  
  err =  err.cl + err.ww
  return(err)
}
