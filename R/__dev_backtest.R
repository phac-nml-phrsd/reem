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
  
  
  # Fri Jun  9 12:08:35 2023 ------------------------------
  # STOPPED HERE
  # implement this as a class method...
 
  # Create new (future) observations
  obs.new = data.frame(
    date = c(ymd('2022-03-14'), ymd('2022-04-01')),
    obs  = c(50, 15)
  )
  var = 'Wr'
  scores = obj$calc_scores(var, obs.new) 
  scores
  obs.new$score <- scores
  
  # Retrieve the forecast densities
  fd = obj$forecast_densities(var, obs.new)
  dfd = bind_rows(fd)
  
  # Plot the new observations
  # and the associated scores:
  md = 200
  g = obj$plot_forecast()
  g.new = g$ww  +
    geom_text(data = obs.new, 
              color = 'red2', size=3,
              mapping = aes(x=date, y=-3,
                            label = round(score,2)))+
    # Densities at new observation dates
    geom_polygon(data = dfd, 
               aes(x = date + md*y, y = x,
                   group = date),
               fill = 'mediumpurple3',
               color = 'mediumpurple1',
               alpha = 0.2,
               linewidth = 0.2) + 
    geom_point(data = obs.new, color='red2',
               shape = 7,size = 3,stroke = 0.8, 
               mapping = aes(x=date, y=obs)) + 
    labs(title = paste('Backtesting', var))
  plot(g.new)
  
  
  
  
}

