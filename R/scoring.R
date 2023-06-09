
# Retrieve the forecasted value at a 
# single date for a single simulation
extract_fcst_value <- function(i, fcst, d, var) {
  tmp = fcst$simfwd[[i]] %>% 
    select(date, !!var) %>% 
    filter(date == d) %>% 
    select(!!var) %>% 
    as.numeric()
  return(tmp)
}



#' Calculate the density of forecasted values 
#' at a single date 
#'
#' @param fcst.vals 
#' @param density.n 
#' @param density.adjust 
#'
#' @return
#'
calc_density_one <- function(fcst.vals,
                             density.n, 
                             density.adjust) {
  # Density for all these forecasted values
  a = density(x      = fcst.vals, 
              n      = density.n, 
              adjust = density.adjust)
  
  # Filter only positive values
  # (the density function expand the 
  # range of values and may allow 
  # negative values, which we know do not 
  # make sense from incidence & concentrations)
  idx = (a$x > 0)
  b   = data.frame(x = a$x[idx], y = a$y[idx])
  b$y = b$y / sum(b$y)   # must normalize after filtering out.
  
  return(b)
}


#' Calculate the forecasting score at 
#' a given date of new observations
#'
#' @param i 
#' @param var 
#' @param obs.new 
#' @param fcst 
#' @param density.n 
#' @param density.adjust 
#'
#' @return
calc_score_one <- function(i, var, 
                           obs.new, 
                           fcst,
                           density.n, 
                           density.adjust) {
  
  # All forecasted values at a given 
  # date across all simulations
  fcst.vals = sapply(X    = 1:length(fcst$simfwd), 
                     FUN  = extract_fcst_value, 
                     fcst = fcst, 
                     d    = obs.new$date[i], 
                     var  = var)
  
  # Density for all these forecasted values
  b = calc_density_one(fcst.vals = fcst.vals, 
                       density.n = density.n, 
                       density.adjust = density.adjust) 
  
  # Retrieve the x-value step of the density
  dx = diff(b$x)[1]
  
  # Identify the density value 
  # at the new observation 
  xx    = which( abs(b$x - obs.new$obs[i]) < dx/2 )[1]
  score = -log(b$y[xx])
  
  # Notes: 
  #  - smaller score is better.
  #  - if observed value is outside range of forecasts, 
  #    then returned score is `NA`.
  
  return(score)  
}



#' Forecasting scores given new observations
#'
#' @param var 
#' @param obs.new 
#' @param fcst 
#' @param density.n 
#' @param density.adjust 
#'
#' @return
#'
reem_calc_scores <- function(var, 
                             obs.new, 
                             fcst,
                             density.n = 100, 
                             density.adjust = 0.4) {
  
  scores = sapply(X       = 1:nrow(obs.new), 
                  FUN     = calc_score_one,
                  var     = var, 
                  obs.new = obs.new, 
                  fcst    = fcst,
                  density.n      = density.n,
                  density.adjust = density.adjust)
  return(scores)    
}



#' Densities of the forecasted values
#' given new observations.
#'
#' @param var 
#' @param obs.new 
#' @param fcst 
#' @param density.n 
#' @param density.adjust 
#'
#' @return
#'
reem_forecast_densities <- function(var, 
                                    obs.new, 
                                    fcst,
                                    density.n , 
                                    density.adjust) {
  res = list()
  
  # Loop through all dates
  
  for(i in 1:nrow(obs.new)){
    
    d =  obs.new$date[i]
    
    # All forecasted values at a single 
    # date across all simulations
    fcst.vals = sapply(X    = 1:length(fcst$simfwd), 
                       FUN  = extract_fcst_value, 
                       fcst = fcst, 
                       d    = d, 
                       var  = var)
    
    tmp  = calc_density_one(fcst.vals = fcst.vals,
                              density.n = density.n, 
                              density.adjust = density.adjust)
    
    res[[i]] = tmp  %>% 
      mutate(date = d, 
             obs  = obs.new$obs[i])
  }
  return(res)
}

