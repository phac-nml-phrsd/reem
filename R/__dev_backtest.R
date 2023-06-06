if(0){
  
  obj = readRDS('obj-dev.rds')
  asof = ymd('2022-03-01') 
  
  
  prm.fcst = list(
    asof         = asof,
    horizon.fcst = ymd('2022-06-01'),
    use.fit.post = TRUE,
    n.resample   = 15,
    ci           = 0.95
  )
  
  fcst = obj$forecast(prm = prm.fcst, verbose = 1)
  
  obj$obs.ww
  
  obs.new = data.frame(
    date = c(ymd('2022-03-14'), ymd('2022-04-01')),
    obs  = c(25, 8)
  )
  
  g = obj$plot_forecast()
  
  g.new = g$ww + geom_point(data = obs.new, color='red2',
                            shape = 15,size = 2, 
                            mapping = aes(x=date, y=obs))
  plot(g.new)
  
  dates.new = obs.new$date
  
  # Tue Jun  6 16:45:30 2023 ------------------------------
  ## STOPPED HERE
  ## implement this as a class method
  
  i = 1
  d = dates.new[1]
  
  extract_fcst_value <- function(i, fcst, d, var) {
    tmp = fcst$simfwd[[i]] %>% 
      select(date, !!var) %>% 
      filter(date == d) %>% 
      select(!!var) %>% 
      as.numeric()
    return(tmp)
  }
  
  foo = sapply(1:length(fcst$simfwd), extract_fcst_value, 
               fcst = fcst, d = dates.new[1], var = 'Wr')

  
  a = density(foo, n = 100, adjust = 0.4)
  plot(a);grid()
  idx = (a$x > 0)
  b = data.frame(x = a$x[idx], y = a$y[idx])
  b$y = b$y / sum(b$y)
  
  dx = diff(b$x)[1]
  xx = which( abs(b$x - 10) < dx/2)[1]
  score = -log(b$y[xx])
  score  # smaller = better
}
