
if(0){
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(reem)
  # devtools::load_all() 
  
  date.start = ymd('2022-01-01')
  asof       = ymd('2022-03-01') 
  
  prms0 = list(
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
    N       = 9999, # population size
    alpha   = 0.2, # transmission heterogeneity (alpha=0: homogeneous)
    I.init  = c(1,1,3,5), # initial incidence (overwritten in fit ABC)
    lag     = 7,   # Aggregation lag for clinical reports
    rho     = 0.1, # mean reporting ratio
    g       = get_gi(), # Generation interval distribution
    fec     = get_fecalshed(), # fecal shedding kinetics
    kappa   = 0.18, # decay in ww
    psi     = get_psi()   # plug flow simulation
  )
  
  obj0 = new('reem', 
             name = 'foo', 
             prms = prms0, 
             is.fitted = FALSE)
  
  obj0$print_prms()
  
  simepi  = obj0$simulate_epi(deterministic = FALSE)
  
  obs.cl = filter(simepi$obs.cl, date <= asof)
  obs.ww = filter(simepi$obs.ww, date <= asof)
  
  prms$R0 <- 2.5
  
  # Attached simulated data to new `reem` object:
  
  prms = prms0
  prms$t.obs.cl <- NULL
  prms$t.obs.ww <- NULL
  obj  = new('reem', 
             name = 'foo2', 
             prms = prms, 
             obs.cl = obs.cl,
             obs.ww = obs.ww,
             is.fitted = FALSE)
  
  # ---- Fit ----
  
  prm.abc = list(
    n.abc = 2e3,
    n.sim = 0,     #`0` for deterministic, else`8` should be enough
    p.abc = 0.01, #1e-2,
    n.cores = min(12, parallel::detectCores() - 1),
    use.cl = 1, 
    use.ww = 1,
    err.type = 'L2'
  )
  
  prms.to.fit = list(
    R0    = c(1.2, 3),
    alpha = c(0, 2),
    i0prop = c(-5,-2),
    start.delta = c(-7,7)  
  )
  
  system.time({
    foo = obj$fit_abc(prm.abc, prms.to.fit)  
  })
  
  gg = obj$plot_fit()
  
  fname = paste0('plot-',timestamp_short(),'.pdf')
  
  pdf(file = fname, width = 25, height = 10)
  plot(gg$all)
  dev.off()
  
  
  d = obj$fit.obj$all.distances
  fit.obj = obj$fit.obj
  
  saveRDS(object = obj, file = 'obj-dev.rds')
  
  

  # - - - - - - - - - -  - - - - - - - - 
  # --- other stuff  ----
  if(0){
    g = simepi$sim %>% 
      ggplot(aes(x=date))+
      geom_point(aes(y=I))+
      geom_step(data = obs.cl, aes(y=obs), 
                color='chartreuse3',
                linewidth = 1)+
      geom_step(data = obs.ww, aes(y=obs), 
                color='chocolate2',
                linewidth = 1)+
      geom_line(data = a$sim, aes(y=I), 
                color = 'steelblue2', alpha = 0.2)
    plot(g + scale_y_log10())
    # plot(g)
    
    
    # --- by reference
    
    a = new('reem', 
            name = 'foo', 
            prms = prms, 
            is.fitted = FALSE)
    
    a$is.fitted
    
    b = a        # copy by reference
    q = a$copy() # copy by value
    
    a$is.fitted <- TRUE  
    
    b$is.fitted
    q$is.fitted
    
    
    a$show()
  }
}


