
if(0){
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(reem)
  devtools::load_all()
  
  date.start = ymd('2022-01-01')
  horizon = 120
  
  B.date = date.start + c(0:(horizon+1))
  B = data.frame(date = B.date, mult = rep(1,length(B.date))) |>
    mutate(mult = if_else(date >= ymd('2022-03-01'), 0.3, mult))
  
  prms = list(
    horizon  = horizon,  # horizon of the simulation
    last.obs = horizon,  # last observation time (must be < horizon)
    B        = B,  # Behavior change
    freq.obs.ww = 3, # average frequency of ww observation
    # Define the schedule of `observations`
    date.obs.cl = date.start + seq(7,280, by = 7),
    date.obs.ha = date.start + seq(12,280, by = 20),
    date.obs.ww = date.start + seq(3,200, by=15),
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
    h.prop  = 0.05, # total proportion hospitalized for one cohort
    h.lags  = c(rep(0,3), 1, 2, 2, 1, 0), # Lag infection-hospitalization
    kappa   = 0.18, # decay in ww
    psi     = get_psi(), # plug flow simulation,
    shed.mult = 0.2 # deposited fecal shedding multiplier  
  )
  
  obj = new('reem', 
            name = 'foo', 
            prms = prms, 
            is.fitted = FALSE)
  
  obj$print_prms()
  head(obj$prms$B)
  
  simepi = obj$simulate_epi(deterministic = TRUE)
  
  obj$print_prms()
  head(obj$prms$B)
  g = plot_epi(simepi) 
  patchwork::wrap_plots(g, ncol = 1)
  
  
  # impact of `alpha`
  obj$prms$alpha <- 0
  s0  = obj$simulate_epi(deterministic = TRUE)
  obj$prms$alpha <- 1
  s2  = obj$simulate_epi(deterministic = TRUE)
  
  df = rbind(
    mutate(s0$sim, alpha = 'alpha = 0'),
    mutate(s2$sim, alpha = 'alpha = 1'))
    
  
  g.alpha = df %>% 
    ggplot(aes(x=date, y = I, color = alpha)) + 
    geom_line(linewidth = 2) + scale_y_log10() + 
    coord_cartesian(ylim=c(1,1e3))+
    theme(panel.grid.minor = element_blank()) +
    labs(title = 'impact of alpha', x='', y='incidence')
  g.alpha
  
  # impact of `R0`
  obj$prms$alpha <- 0
  obj$prms$R0    <- 1.5
  
  s0  = obj$simulate_epi(deterministic = TRUE)
  obj$prms$R0 <- 2
  s2  = obj$simulate_epi(deterministic = TRUE)
  
  df = rbind(
    mutate(s0$sim, R0 = 'R0 = 1.5'),
    mutate(s2$sim, R0 = 'R0 = 2'))
  
  
  g.R0 = df %>% 
    ggplot(aes(x=date, y = I, color = R0)) + 
    geom_line(linewidth = 2) + scale_y_log10()+
    coord_cartesian(ylim=c(1,1e3))+
    theme(panel.grid.minor = element_blank()) +
    scale_color_brewer(palette = 'Dark2')+
    labs(title = 'impact of Ro', x='', y='incidence')
  g.R0
  
} 
