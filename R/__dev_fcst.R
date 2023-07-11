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
  theme_set(theme_bw())
  obj = readRDS('obj-dev.rds')
  asof = ymd('2022-03-01') 
  
  prm.fcst = list(
    asof         = asof,
    horizon.fcst = ymd('2022-06-01'),
    use.fit.post = TRUE,
    n.resample   = 20,
    ci           = seq(0.1,0.9, by = 0.1)
  )
  
  fcst = obj$forecast(prm = prm.fcst, verbose = 1)
  
  g.fcst = obj$plot_forecast(date_breaks = '1 month')
  g      = patchwork::wrap_plots(g.fcst, ncol=1)
  pdf(paste0('plot-fcst-', reem::timestamp_short(),'.pdf'))
  plot(g)
  dev.off()
  
  var = 'Y.aggr'
  date.lower = ymd('2022-03-10')
  date.upper = ymd('2099-01-01')
  val.lower = 700
  val.upper = 9900
  
  a = obj$proba_box(var        = var, 
                    date.lower = date.lower, 
                    date.upper = date.upper,
                    val.lower  = val.lower, 
                    val.upper  = val.upper)
  print(a)
  
  fcst.obj = obj$fcst.obj 
}
