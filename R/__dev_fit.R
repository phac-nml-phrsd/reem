
if(0){
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(reem)
  devtools::load_all()
  
  set.seed(1234)
  date.start = ymd('2022-01-01')
  asof       = ymd('2022-03-10') #ymd('2022-02-21') 
  hz = 220
 
  B.date = date.start + c(0:(hz+1))
  date.break = asof + 30 
  B = data.frame(date = B.date, mult = rep(1,length(B.date))) |>
    mutate(mult = if_else(date >= date.break, 2, mult))
   
  prms0 = list(
    horizon = hz,  # horizon of the simulation
    last.obs = hz-1,  # last observation time (must be < horizon)
    B       = B, # Behavior change
    freq.obs.ww = 3, # average frequency of ww observation
    date.obs.cl = date.start + seq(7,hz-2, by = 7),
    date.obs.ha = date.start + seq(12,hz-2, by = 10),
    date.obs.ww = date.start + seq(3,hz-2, by=14),
    i0prop  = 1e-3,
    date.start = date.start,
    start.delta = 0, 
    R0      = 1.9, # Basic reproduction number
    N       = 1e6, # population size
    alpha   = 0.2, # transmission heterogeneity (alpha=0: homogeneous)
    I.init  = c(1,1,3,5)*20, # initial incidence (overwritten in fit ABC)
    lag     = 7,   # Aggregation lag for clinical reports
    rho     = 0.1, # mean reporting ratio
    g       = get_gi(), # Generation interval distribution
    fec     = get_fecalshed(), # fecal shedding kinetics
    h.prop  = 0.05, # total proportion hospitalized for one cohort
    h.lags  = c(rep(0,3), 1, 2, 2, 1, 0), # Lag infection-hospitalization
    kappa   = 0.18, # decay in ww
    psi     = get_psi(),   # plug flow simulation,
    shed.mult = 1e-3
  )
  
  obj0 = new('reem', 
             name = 'foo', 
             prms = prms0, 
             is.fitted = FALSE)
  
  obj0$print_prms()
  
  simepi  = obj0$simulate_epi(deterministic = 0)
  g.epi = plot_epi(simepi) 
  g.epi |> patchwork::wrap_plots(ncol=1)
  g.epi$populations + geom_vline(xintercept = date.break, linetype = 'dashed'  )
  
  obs.cl = filter(simepi$obs.cl, date <= asof)
  obs.ha = filter(simepi$obs.ha, date <= asof)
  obs.ww = filter(simepi$obs.ww, date <= asof)
  
  
  # Attached simulated data to new `reem` object:
  prms = prms0
  prms$t.obs.cl <- NULL
  prms$t.obs.ha <- NULL
  prms$t.obs.ww <- NULL
  obj  = new('reem', 
             name = 'foo2', 
             prms = prms, 
             obs.cl = obs.cl,
             obs.ha = obs.ha,
             obs.ww = obs.ww,
             is.fitted = FALSE)
  
  g.obs = plot_obs(obj)
  g.obs
  
  # ---- Fit ----
  
  prm.abc = list(
    n.abc = 1e3,
    n.sim = 0,     #`0` for deterministic, else`8` should be enough
    p.abc = 0.01, #1e-2,
    n.cores = 1, #min(12, parallel::detectCores() - 1),
    use.cl = 1, 
    use.ha = 1, 
    use.ww = 1,
    err.type = 'normlarge'   # normlarge, L2
  )
  
  prms.to.fit = list(
    R0          = list('gamma', 1.5, 0.251),
    alpha       = list('normp', 2, 1),
    start.delta = list('unif_int', -7, 7)  ,
    i0prop      = list('unif', -5.8, -2),
    h.prop      = list('unif', 0.001, 0.2),
    B20220201   = list('normp', 1.0, 0.3),
    B20220215   = list('normp', 1.0, 0.6)
  )
  
  foo = obj$fit_abc(prm.abc, prms.to.fit)  
  
  saveRDS(obj, file = 'debug-fit.rds')
  
  gg = obj$plot_fit()
  gg$all
  
  fname = paste0('plot-',timestamp_short(),'.pdf')
  
  plot(gg$traj.cl)
  plot(gg$traj.ha)
  plot(gg$traj.ww)
  plot(gg$post.prms)
  plot(gg$post.prms.2d)
  
  
  pdf(file = fname, width = 25, height = 10)
  plot(gg$all)
  dev.off()
  
  d = obj$fit.obj$all.distances
  fit.obj = obj$fit.obj
  
  saveRDS(object = obj, file = 'obj-dev.rds')
  
  message('\nDONE.\n')
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


