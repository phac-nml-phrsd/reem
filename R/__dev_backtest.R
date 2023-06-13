if(0){
  
  suppressMessages({
    library(tidyr)
    library(dplyr)
    library(ggplot2) ; theme_set(theme_bw())
    library(lubridate)
    library(stringr)
    library(patchwork)
  })
  
  # Retrieve a fitted object
  obj = readRDS('obj-dev.rds')
  
  # Set forecast parameters
  prm.fcst = list(
    asof         = ymd('2022-03-01') ,
    horizon.fcst = ymd('2022-06-01'),
    use.fit.post = TRUE,
    n.resample   = 15,
    ci           = 0.99
  )
  
  # Perform the forecast
  fcst = obj$forecast(prm = prm.fcst, verbose = 0)
  
  # Create new (future) observations
  obs.new = data.frame(
    date = ymd('2022-03-14') + c(0,14,30),#c(0:2)*14,
    # obs  = c(50, 15,4)
    obs = c(300,100,10)
  )
  var = 'Y' #Wr'
  
  aggr.window = ifelse(var=='Y', 7, NULL)
   
  g.fs = plot_forecast_scores(obj = obj, 
                              obs.new = obs.new, 
                              var = var, 
                              aggr.window = aggr.window)
  plot(g.fs)
}

