
if(0){
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(reem)
  # devtools::load_all() 
  
  date.start = ymd('2022-01-01')
  asof       = ymd('2022-03-01') 
  
  prms = list(
    horizon = 300,  # horizon of the simulation
    last.obs = 299,  # last observation time (must be < horizon)
    B       = rep(1,300), # Behavior change
    freq.obs.ww = 3, # average frequency of ww observation
    t.obs.cl = seq(7,280, by = 7),
    t.obs.ww = seq(3,200, by=3),
    i0prop  = 1e-3,
    date.start = date.start,
    start.delta = 0, 
    R0      = 1.5, # Basic reproduction number
    N       = 1e4, # population size
    alpha   = 0.2, # transmission heterogeneity (alpha=0: homogeneous)
    I.init  = c(1,1,3,5), # initial incidence (overwritten in fit ABC)
    lag     = 7,   # Aggregation lag for clinical reports
    rho     = 0.1, # mean reporting ratio
    g       = get_gi(), # Generation interval distribution
    fec     = get_fecalshed(), # fecal shedding kinetics
    kappa   = 0.18, # decay in ww
    psi     = get_psi(), # plug flow simulation,
    shed.mult = 0.2 # deposited fecal shedding multiplier  
  )
  
  obj = new('reem', 
            name = 'foo', 
            prms = prms, 
            is.fitted = FALSE)
  
  obj$print_prms()
  
  simepi  = obj$simulate_epi(deterministic = FALSE)
  
  simepi$sim %>% ggplot(aes(x=date, y = Wd)) + geom_line()
  
} 
