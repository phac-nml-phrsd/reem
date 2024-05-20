
# Retrieve the forecasted value at a 
# single date for a single simulation
extract_fcst_value <- function(i, fcst, d, var) {
  
  is.aggr = (var == 'Y.aggr')
  
  if(is.aggr)  a = fcst$simfwd.aggr$Y.aggr[[i]]
  if(!is.aggr) a = fcst$simfwd[[i]]
  
  tmp = a %>%
    select(date, !!var) %>% 
    ungroup() %>%
    filter(date == d) 
  
  if(nrow(tmp)==0){
    stop('Date ',d,' not found in the simulated forecasts.\n',
         'Cannot extract the forecasted value. ',
         'Maybe your forcast horizon is earlier than the new observation dates?\n ',
         'These are the available dates:\n',
         paste(fcst$simfwd.aggr$Y.aggr[[i]]$date, sep= ' ; '),
         '\nAborting!\n')
  }
  
  tmp = tmp %>% 
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
  
  # Add two artificial points such 
  # that it is pretty if it is plotted
  n = nrow(b)
  dx = diff(b$x)[1]
  bb = data.frame( 
    x = c(max(0, b$x[1] - dx), b$x, b$x[n] + dx ),
    y = c(0, b$y, 0))
  return(bb)
}


#' @title Helper function that extract forecasted values
#' with the option to aggregate them (e.g. clinical cases) 
#'
#' @param i 
#' @param fcst 
#' @param aggr.window 
#' @param obs.new 
#' @param var 
#'
#' @return Numerical vector.
#' 
#' @keywords internal
#' 
extract_helper <- function(i, fcst, aggr.window, obs.new, var) {
  
  # Number of simulations informing the forecast
  nfs = length(fcst$simfwd)
  
  if(is.null(aggr.window)){
    # All forecasted values at a given 
    # date across all simulations
    res = sapply(X    = 1:length(fcst$simfwd), 
                 FUN  = extract_fcst_value, 
                 fcst = fcst, 
                 d    = obs.new$date[i], 
                 var  = var)
  }
  
  if(!is.null(aggr.window)){
    tmp = matrix(nrow = aggr.window, ncol = nfs)
    for(k in 1:aggr.window){
      tmp[k,] = sapply(X    = 1:nfs, 
                       FUN  = extract_fcst_value, 
                       fcst = fcst, 
                       # TODO FIXME
                       # 
                       # The 1st line (commented) works with variables other than `Y.aggr`
                       # The 2nd line (uncommented) works with `Y.aggr` 
                       #d    = obs.new$date[i] - k+1, # `+1` bc we want `-k+1=0` 
                       d    = obs.new$date[i],  
                       var  = var)
    }
    res = apply(tmp, MARGIN = 2, FUN = sum)
  }
  return(res) 
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
#' @param aggr.window 
#'
#' @return
calc_score_one <- function(i, var, 
                           obs.new, 
                           fcst,
                           density.n, 
                           density.adjust,
                           aggr.window,
                           pb) {
  # print(i)
  
  setTxtProgressBar(pb, value = i)
  
  tmp.fv = extract_helper(i = i, 
                          fcst = fcst, 
                          aggr.window = aggr.window,
                          obs.new = obs.new, 
                          var = var)
  fcst.vals = tmp.fv[!is.na(tmp.fv)]
  
  if(length(fcst.vals) == 0) return(NA)
    
  # Density for all these forecasted values
  b = calc_density_one(fcst.vals = fcst.vals, 
                       density.n = density.n, 
                       density.adjust = density.adjust) 
  
  # Retrieve the x-value step of the density
  dx = diff(b$x)[2]
  
  # Identify the density value 
  # at the new observation 
  xx    = which( abs(b$x - obs.new$obs[i]) < dx )[1]
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
                             density.adjust = 0.4,
                             aggr.window = NULL) {
  n = nrow(obs.new)
  pb = txtProgressBar(min = 0, max = n, style = 3, width = 30)
  scores = sapply(X       = 1:n, 
                  FUN     = calc_score_one,
                  var     = var, 
                  obs.new = obs.new, 
                  fcst    = fcst,
                  density.n      = density.n,
                  density.adjust = density.adjust,
                  aggr.window    = aggr.window,
                  pb = pb)
  return(scores)    
}



#' Assess if the forecasted credible interval contains
#' the actual (future) observed value.
#'
#' @param var String. Name of the variable.
#' @param obs.new Dataframe of the new/actual observations of this variable.
#' @param ci.width Numerical. Width of the centered credible interval.
#' @param aggr.window Numerical. Width of the aggregation window for the variable. 
#' \code{NULL} if no aggregation.
#' @param fcst Forecast object of the \code{reem} object.
#'
#' @return A logical vector indicating if the actual observed value
#' is inside the credible interval. Each element corresponds to an 
#' observation date. (Hence the vector has the same length as the 
#' number of rows of \code{obs.new}.)
#'
reem_inside_CI <- function(var, obs.new, ci.width, aggr.window, fcst) {
  
  if(0){
    ci.width = 0.80
  }
  
  n = nrow(obs.new)
  inside = logical(n)
  for(i in 1:n){
    message('inside CI: ',i,'/',n)
    tmp.fv = extract_helper(i = i, 
                            fcst = fcst, 
                            aggr.window = aggr.window,
                            obs.new = obs.new, 
                            var = var)
    fcst.vals = tmp.fv[!is.na(tmp.fv)]
    
    q = quantile(fcst.vals, probs = 0.5 + c(-1,1) * ci.width/2)
    inside[i] = (q[1] <= obs.new$obs[i]) & (obs.new$obs[i] <= q[2]) 
  }
  return(inside)
}


#' Densities of the forecasted values
#' given new observations.
#'
#' @param var 
#' @param obs.new 
#' @param fcst 
#' @param density.n 
#' @param density.adjust 
#' @param aggr.window 
#'
#' @return
#'
reem_forecast_densities <- function(var, 
                                    obs.new, 
                                    fcst,
                                    density.n , 
                                    density.adjust,
                                    aggr.window) {
  res = list()
  
  # Loop through all dates
  
  for(i in 1:nrow(obs.new)){
   # print(i) # DEBUG
    fcst.vals = extract_helper(
      i = i, 
      fcst = fcst, 
      aggr.window = aggr.window, 
      obs.new = obs.new, 
      var = var)
    
    tmp  = calc_density_one(fcst.vals = fcst.vals,
                            density.n = density.n, 
                            density.adjust = density.adjust)
    
    res[[i]] = tmp  %>% 
      mutate(date = obs.new$date[i], 
             obs  = obs.new$obs[i])
  }
  return(res)
}



#' Title
#'
#' @param obj 
#' @param obs.new 
#' @param var 
#' @param aggr.window 
#'
#' @return
#' @export
#'
#' @examples
plot_forecast_scores <- function( obj,  
                                  obs.new,
                                  var, 
                                  aggr.window = NULL) {
  
  # Forecasting score of the new observations
  scores = obj$calc_scores(var         = var, 
                           obs.new     = obs.new, 
                           aggr.window = aggr.window) 
  
  # Attach scores to associated observation
  obs.new$score <- scores
  
  # Retrieve the forecast densities
  fd  = obj$forecast_densities(var         = var, 
                               obs.new     = obs.new, 
                               aggr.window = aggr.window)
  dfd = bind_rows(fd)
  
  # normalize for pretty figure
  dfdn = dfd %>% 
    group_by(date) %>% 
    mutate(yn = y / max(y))
  
  # === Plot the new observations
  # === and the associated scores:
  
  # Make sure the height of the 
  # densities do not overlap
  dt = min(as.numeric(diff(obs.new$date)))
  md = dt * 0.6
  
  g = obj$plot_forecast()
  
  g.tmp = g$cl
  if(var == 'Wr') g.tmp = g$ww
  
  g.new = g.tmp  +
    geom_text(data = obs.new, 
              color = 'indianred2', size=3,
              mapping = aes(x=date, y=-3,
                            label = round(score,2)))+
    # Densities at new observation dates
    geom_polygon(data = dfdn, 
                 aes(x = date + md*yn, y = x,
                     group = date),
                 fill  = 'mediumpurple3',
                 color = 'mediumpurple1',
                 alpha     = 0.2,
                 linewidth = 0.2) + 
    geom_point(data = obs.new, color='indianred2',
               shape = 7,size = 3,stroke = 0.8, 
               mapping = aes(x=date, y=obs)) + 
    labs(title = paste('Forecast scoring:', var))
  g.new
  
}


