if(FALSE){
  
  # prms = list(R0=1.4, alpha = 0.2, bar = 3.14)
  # obj = new('reem', name = 'testmodel', prms = prms, is.fitted = FALSE)
  # 
  # obj$name
  # obj$print_prms()
  
  prms = list(
    horizon = 300,  # horizon of the simulation
    last.obs = 299,  # last observation time (must be < horizon)
    B       = rep(1,1e3), # Behavior change
    freq.obs.ww = 3, # average frequency of ww observation
    t.obs.ww = seq(3,75, by=3),
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
  
  obj2 = new('reem', name = 'foo', prms = prms, is.fitted = FALSE)
  obj2
  sim  = obj2$simulate(deterministic = FALSE)
  
  
  library(dplyr)
  library(ggplot2)
  sim.obs = sim %>% 
    mutate(tmp = t %% prms$lag) %>% 
    filter(tmp==0) %>% 
    select(t, Y)
  
  
  g = sim %>% 
    ggplot(aes(x=t))+
    geom_step(aes(y=I))+
    geom_point(data = sim.obs, aes(y=Y), color='blue')+
    geom_line(aes(y=Wd), color='chocolate2', linewidth=1)
  plot(g)
}
