
if(0){
  
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  # devtools::load_all() 
  
  prms = list(
    horizon = 300,  # horizon of the simulation
    last.obs = 299,  # last observation time (must be < horizon)
    B       = rep(1,1e3), # Behavior change
    freq.obs.ww = 3, # average frequency of ww observation
    t.obs.cl = seq(7,280, by = 7),
    t.obs.ww = seq(3,200, by=3),
    i0prop  = 1e-3,
    date.start = ymd('2022-01-01'),
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
            prms = prms, 
            is.fitted = FALSE)
  
  obj0$print_prms()
  
  simepi  = obj0$simulate_epi(deterministic = FALSE)
  
  obs.cl = simepi$obs.cl
  obs.ww = simepi$obs.ww
  
  prms$R0 <- 2.5
  
  # Attached simulated data to new `reem` object:
  obj  = new('reem', 
             name = 'foo2', 
             prms = prms, 
             obs.cl = obs.cl,
             obs.ww = obs.ww,
             is.fitted = FALSE)
  
  a = obj$traj_dist_obs(use.cl = 1, 
                        use.ww = 1, 
                        deterministic = FALSE,
                        err.type = 'L2',
                        verbose = TRUE)
  
  print(a$distance)

  
  # Tue May 23 11:31:34 2023 ------------------------------
  
  prm.abc = list(
    n.abc = 20,
    n.sim = 8,     # `8` should be enough
    p.abc = 0.20,
    n.cores = 4,  # parallel::detectCores() - 1,
    use.cl = 1, 
    use.ww = 1,
    err.type = 'L2'
  )
  
  prms.to.fit = list(
    R0    = c(0.90, 1.99),
    alpha = c(0, 9),
    i0prop = c(-6,-2),
    start.delta = c(-7,7)  
  )
  
  foo = obj$fit_abc(prm.abc, prms.to.fit)  
  
  obj$is.fitted 
  
  gp = foo$posteriors %>% 
    select(-abc.err) %>% 
    tidyr::pivot_longer(cols = -abc.index ) %>%
    ggplot(aes(x=value)) + 
    geom_histogram() + 
    facet_wrap(~name, scales = 'free')
    
  gp  
  
  
  
  
  g = simepi$sim %>% 
    ggplot(aes(x=t))+
    geom_point(aes(y=I))+
    geom_step(data = obs.cl, aes(y=obs), color='chartreuse3',
              linewidth = 1)+
    geom_step(data = obs.ww, aes(y=obs), color='chocolate2',
              linewidth = 1)+
    geom_line(data = a$sim, aes(y=I, group = n.sim), 
              color = 'steelblue2', alpha = 0.2)
  plot(g + scale_y_log10())
  # plot(g)
}
