###
### DEVELOPMENT SCRIPT : forecast
###

if(0){
  suppressMessages({
    library(tidyr)
    library(dplyr)
    library(ggplot2) ; theme_set(theme_bw())
    library(lubridate)
    library(stringr)
    library(patchwork)
  })
  devtools::load_all()
  
  
  obj = readRDS('debug-fit.rds')
  asof = ymd('2022-03-21') 
  
  prm.fcst = list(
    asof         = asof,
    horizon.fcst = ymd('2022-06-01'),
    use.fit.post = TRUE,
    n.resample   = 20,
    vars.to.fcst = c('Y', 'Wr', 'H'),
    ci           = seq(0.1,0.9, by = 0.1)
  )
  
  
  fcst = obj$forecast(prm = prm.fcst, verbose = 1)
  
  g.fcst = obj$plot_forecast(date_breaks = '1 month')
  g      = patchwork::wrap_plots(g.fcst, nrow=1)
  g
 
  pk = obj$forecast_peak(var = 'H.aggr')
  pk
  mean(pk$peak.value > 10000) 
  
  g.peak.ha = obj$plot_peak(var = 'H.aggr', logscale = 0)
  g.peak.cl = obj$plot_peak(var = 'Y.aggr', logscale = 0)
  g.peak.cl | g.peak.ha
  
  pdf(paste0('plot-fcst-', reem::timestamp_short(),'.pdf'))
  plot(g)
  dev.off()
  
  var = 'Y'  # H.aggr   Y
  date.lower = ymd('2022-04-01')
  date.upper = ymd('2022-09-12')
  val.lower = 2500
  val.upper = 25000
  
  a = obj$proba_box(var        = var, 
                    date.lower = date.lower, 
                    date.upper = date.upper,
                    val.lower  = val.lower, 
                    val.upper  = val.upper)
  print(a)
  
  
  
  
  fcst.obj = obj$fcst.obj 
}
