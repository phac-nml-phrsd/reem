
#' Simulate an epidemic with a REEM.
#'
#' @param deterministic 
#' @param prms 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
reem_simulate <- function(prms, deterministic) {
  # Unpack parameters
  R0      = prms$R0
  B       = prms$B
  N       = prms$N
  alpha   = prms$alpha
  I.init  = prms$I.init
  horizon = prms$horizon
  rho     = prms$rho
  lag     = prms$lag
  g       = prms$g
  fec     = prms$fec
  kappa   = prms$kappa
  psi     = prms$psi
  t.obs.ww = prms$t.obs.ww
  shed.mult = prms$shed.mult
  
  ni = length(I.init)
  
  m = numeric(horizon)
  I = rep(NA, horizon)
  S = rep(NA, horizon)
  A = rep(NA, horizon)
  Y = rep(NA, horizon)
  Wd = rep(NA, horizon)
  
  # Initial period when incidence is known:
  m[1:ni] = I.init
  I[1:ni] = I.init
  S[1:ni] = N - cumsum(I.init)
  
  for(t in 1:ni){
    tlag = max(1, t-lag)
    A[t] = sum(I[t:tlag])
  }
  
  if(!deterministic) rho = rnorm(n=horizon, mean = rho, sd = 0.05)
  rho[rho <= 0] <- 0.001
  
  for(t in (ni+1):horizon){ # t = ni+1
    
    m[t] = mean_inc(t, R0, B, S, N, alpha, g, I)
    I[t] = calc_I(lambda = m[t], deterministic)
    S[t] = max(0, S[t-1] - I[t])
    
    tlag = max(1, t-lag)
    A[t] = sum(I[t:tlag])
  }
  
  # Observed aggregated incidence (Y):
  lambdaY = rho*A
  lambdaY[lambdaY==0] <- 1e-3
  lambdaY[is.na(lambdaY)] <- 1e-3
  Y = calc_Y(lambdaY, n=length(A), deterministic)
  
  # --- Wastewater ---
  
  # -- deposited
  
  nf = length(fec)
  w.m = numeric(horizon)
  for(t in 1:horizon){
    idx = 1:min(t-1,nf)
    w.m[t] = shed.mult * sum(fec[idx]*I[t-idx])
  }
  Wd =  w.m
  if(!deterministic) Wd = rnorm(n = horizon, mean = w.m, sd = w.m * 0.1)
  
  # -- present at sampling site
  
  n.psi= length(psi)
  wp.m = numeric(horizon)
  
  for(t in 1:horizon){
    idx = 1:min(t-1,n.psi)
    wp.m[t] = sum(psi[idx] * Wd[t-idx] * exp(-kappa * idx))
  }
  Wp = wp.m
  if(!deterministic) Wp = rnorm(n=horizon, 
                                mean = wp.m, 
                                sd = wp.m * 0.1)
  
  # -- observed ("reported") at sampling site
  
  t.obs.ww = t.obs.ww[t.obs.ww < horizon]
  n.ww = length(t.obs.ww)
  wr.m = numeric(n.ww)
  
  for(i in seq_along(t.obs.ww)) 
    wr.m[i] = Wp[t.obs.ww[i]]
  
  Wr = wr.m
  if(!deterministic) Wr = rnorm(n=n.ww, mean = wr.m, sd = wr.m * 0.2)
  
  df = data.frame(
    t = 1:horizon, 
    m = m, 
    I = I,
    S = S,
    A = A,
    Y = Y,
    Wd = Wd, 
    Wp = Wp)
  
  # This is equivalent as but quicker than a `left_join()`
  # and we want this code to be as fast as possible!
  df$Wr <- NA
  idx = df$t %in% t.obs.ww
  df$Wr[idx] <- Wr
   
  return(df)  
}



#' Simulate a full epidemic with a REEM
#' including clinical and wastewater "observations".
#'
#' @param variables 
#'
#' @return
#' @export
#'
#' @examples
reem_simulate_epi <- function(obj, 
                              deterministic) {
  
  prms = obj$prms
  
  # == Build times ==
   
  # When working with real data, 
  # we usually are in "date mode"
  has.date.start = !is.null(prms$date.start)
  
  if(!has.date.start){
    prms$date.start <- lubridate::ymd('2020-01-01')
    
    warning('The start date (prms$date.start) of ',
            'the REEM model has not been specified.',
            ' Setting it to 2020-01-01. ')
  }
  
  # -- Observation times 
  
  # If observation times/dates not specified:
  # - if `obs.xx` exists, set obs date/times to the ones of this object.
  # - if `obs.xx` does not exist, assume observation at every time step.
  
  # Check if dates or times are missing as input 
  miss.obs.dt.cl = is.null(prms$t.obs.cl) & is.null(prms$date.obs.cl)
  miss.obs.dt.ww = is.null(prms$t.obs.ww) & is.null(prms$date.obs.ww)
  
  if(miss.obs.dt.ww){
    # If not specified, then assumed 
    # observed at all time steps:
    date.obs.ww  = prms$date.start + 1:prms$horizon
    t.obs.ww     = 1:prms$horizon
    prms = c(prms, 
             t.obs.ww    = list(t.obs.ww),
             date.obs.ww = list(date.obs.ww))
  }
  
  if(miss.obs.dt.cl){
    # same times as when clinical is observed:
    date.obs.cl = obj$obs.cl$date
    
    # if clinical observation times are not defined
    times.cl.def = 't' %in% names(obj$obs.cl)
    if(! times.cl.def){
      t.obs.cl = as.integer(obj$obs.cl$date - prms$date.start)
    }
    if(times.cl.def) t.obs.cl = obj$obs.cl$t
  }
  
  # For simulations (TODO: change that...)
  if(!is.null(prms$t.obs.cl)) {
    t.obs.cl    = prms$t.obs.cl
    date.obs.cl = prms$date.start + t.obs.cl
  }
  
   # Simulate epidemic to generate data
  sim = reem_simulate(prms, deterministic)
  
  sim = dplyr::mutate(sim, date = prms$date.start + t)
  
  # Create two dataframes of observations only:
  #  - clinical data
  #  - wastewater data
  # (they may not be observed on the same schedule)
  
  
  # Aggregate clinical reports
  sim.obs.cl = aggregate_time(df       = sim, 
                              dt.aggr  = date.obs.cl, 
                              var.name = 'Y') %>% 
    dplyr::transmute(
      date, 
      t = as.integer(date - prms$date.start),
      obs  = aggregation) 
  
  # Extract wastewater observations
  sim.obs.ww = sim %>% 
    dplyr::select(t, Wr) %>% 
    dplyr::filter(t <= prms$last.obs, 
                  !is.na(Wr)) %>%   # filtering out NA keeps observed times only.
    dplyr::rename(obs = Wr) %>%
    dplyr::mutate(date = prms$date.start + t) %>%
    select(t, date, obs)
  
  return(list(
    sim    = sim, 
    obs.cl = sim.obs.cl,
    obs.ww = sim.obs.ww
  ))
}


