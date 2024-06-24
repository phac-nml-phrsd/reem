

#' Check (and potentially set) the simulation start date
#'
#' @param obj 
#'
#' @return
#' 
check_date_start <- function(obj) {
  
  prms = obj$prms
  
  # When working with real data, 
  # we usually are in "date mode"
  has.date.start = !is.null(prms$date.start)
  
  if(!has.date.start){
    prms$date.start <- lubridate::ymd('2020-01-01')
    
    warning('The start date (`prms$date.start`) of ',
            'the REEM model has not been specified.',
            ' Setting it to 2020-01-01. ')
  }
  return(prms)
}

#' Helper function
#' @keywords internal
#' 
create_obs_if_missing <- function(obj, type, verbose) {
  
  # DEBUG:  type = 'cl'
  
  thetype = ifelse(type == 'cl', 'obs.cl', 'obs.ha')
  
  theobs = obj[[thetype]]
  n  = nrow(theobs)
  
  # If no observations were attached, 
  # create an empty one.
  if(n == 0){
    
    if(verbose) message('no existing ',type,' observations, creating empty one...')
    
    h = prms$horizon
    obs_type = data.frame(
      t    = 0:(h-1), 
      date = lubridate::ymd('2020-01-01') + 0:(h-1), 
      inc  = NA)
    prms[[thetype]] = obs_type
    
    warning('No ', 
            ifelse(type == 'cl', 'clinical reports', 'hospital admission'),
            ' observations (`prms$',
            ifelse(type == 'cl', 'obs.cl', 'obs.ha'),
            '`) attached to the REEM object. ',
            'Setting an empty one. ')
  }
  
  # If observations attached, use this one
  if(n > 0){
    if(verbose) message('Using existing ',type,' observations for schedule...')
    obs_type = theobs
  }
  
  # if observation _times_ are not defined
  times.def = 't' %in% names(obs_type)
  if(! times.def){
    t.obs.type = as.integer(obs_type$date - prms$date.start)
  }
  if(times.def) t.obs.type = obs_type$t
  
  if(type == 'cl')  prms = c(prms, 
                             t.obs.cl    = list(t.obs.type),
                             date.obs.cl = list(obs_type$date))
  
  if(type == 'ha')  prms = c(prms, 
                             t.obs.ha    = list(t.obs.type),
                             date.obs.ha = list(obs_type$date))
  
  if(verbose){
    message('type = ', type, '  ==> prms$t.obs.cl = ', prms$t.obs.cl)
    message('type = ', type, '  ==> prms$t.obs.ha = ', prms$t.obs.ha)
  }
  
  return(prms)
}


#' Check observations (data) times
#'
#' @param obj 
#'
#' @return
#'
check_obs_schedule <- function(obj, verbose = 0) {
  
  prms = obj$prms
  
  # If observation times/dates not specified:
  # - if `[t/date].obs.xx` exists, set obs date/times to the ones already existing of this object.
  # - if `[t/date].obs.xx` does not exist, assume observation at every time step.
  
  # Check if dates or times are missing as input 
  miss.obs.dt.cl = is.null(prms$t.obs.cl) & is.null(prms$date.obs.cl)
  miss.obs.dt.ha = is.null(prms$t.obs.ha) & is.null(prms$date.obs.ha)
  miss.obs.dt.ww = is.null(prms$t.obs.ww) & is.null(prms$date.obs.ww)
  
  if(verbose){ # DEBUG 
    message('miss.obs.dt.cl = ', miss.obs.dt.cl)
    message('miss.obs.dt.ha = ', miss.obs.dt.ha)
    message('miss.obs.dt.ww = ', miss.obs.dt.ww)
  }
  
  if(is.null(prms$date.start)) stop('`date.start` must be specified. ABORTING!')
  if(is.null(prms$horizon)) stop('`horizon` must be specified. ABORTING!')
  
  # -- Wastewater data
  
  if(miss.obs.dt.ww){
    
    has.obs.ww = !is.null(obj$obs.ww)
    
    if(verbose) message('has.obs.ww = ', has.obs.ww)
    
    # if observations specified, use their times and dates
    if(has.obs.ww){
      prms = c(prms, 
               t.obs.ww    = list(obj$obs.ww$t),
               date.obs.ww = list(obj$obs.ww$date))
      
      if(verbose) message('prms$t.obs.ww = ', prms$t.obs.ww)
    }
    
    # If not specified, then assumed 
    # observed at all time steps:
    if(! has.obs.ww){
      date.obs.ww = prms$date.start + 1:prms$horizon
      t.obs.ww    = 1:prms$horizon
      
      prms = c(prms, 
               t.obs.ww    = list(t.obs.ww),
               date.obs.ww = list(date.obs.ww))
      
      warning('Observation times/dates for wastewater data not specified: ',
              'assuming ww observation schedule at every time step.')
    }
  }
  
  # -- Clinical data
  
  if(miss.obs.dt.cl) prms = create_obs_if_missing(obj, type = 'cl', verbose) 
  if(miss.obs.dt.ha) prms = create_obs_if_missing(obj, type = 'ha', verbose) 
  
  
  return(prms)
}


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
  R0      = prms[['R0']]
  B       = prms[['B']]
  N       = prms[['N']]
  alpha   = prms[['alpha']]
  I.init  = prms[['I.init']]
  horizon = prms[['horizon']]
  rho     = prms[['rho']]
  lag     = prms[['lag']]
  g       = prms[['g']]
  fec     = prms[['fec']]
  h.prop  = prms[['h.prop']]
  h.lags  = prms[['h.lags']]
  kappa   = prms[['kappa']]
  psi     = prms[['psi']]
  t.obs.ww = prms[['t.obs.ww']]
  shed.mult = prms[['shed.mult']]
  
  ni = length(I.init)
  
  m  = numeric(horizon)
  I  = rep(NA, horizon)   # daily incidence
  S  = rep(NA, horizon)   # daily number of susceptible
  A  = rep(NA, horizon)   # rolling sum of daily incidence aggregated over `lag`
  Y  = rep(NA, horizon)   # observed aggregated incidence (:stochastic fraction of `A`)
  H  = rep(0, horizon)    # daily hospital admissions
  Wd = rep(NA, horizon)   # wastewater concentration deposited
  
  # -- Initial period when incidence is known:
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
  
  # -- Observed aggregated incidence (Y):
  lambdaY = rho*A
  lambdaY[lambdaY==0] <- 1e-3
  lambdaY[is.na(lambdaY)] <- 1e-3
  Y = calc_Y(lambdaY, n=length(A), deterministic)
  
  # -- Hospital admissions (H):
  
  
  # calculate `h` from lags and total proportion
  h.undefined = is.null(h.lags) | is.null(h.prop)
  if(h.undefined) {
    h = rep(0,2)
  }
  if(!h.undefined){
    h = h.lags / sum(h.lags) * h.prop
  }
  
  nh = length(h)
  for(t in 2:horizon){
    H[t] = 0
    upperidx = min(nh, t-1)
    for(k in 1:upperidx){
      H[t] = H[t] + h[k] * I[t-k]
    }
    H[t] = round(H[t])
  }
  
  # --- Wastewater ---
  
  # -- deposited
  
  nf = length(fec)
  w.m = numeric(horizon)
  for(t in 1:horizon){
    idx = 1:min(t-1,nf)
    w.m[t] = shed.mult * sum(fec[idx]*I[t-idx])
  }
  Wd =  w.m
  if(!deterministic) Wd = rnorm(n    = horizon, 
                                mean = w.m, 
                                sd   = w.m * 0.1)
  
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
  if(!deterministic) Wr = rnorm(n    = n.ww, 
                                mean = wr.m, 
                                sd   = wr.m * 0.2)
  
  # Ending
  
  df = data.frame(
    t = 1:horizon, 
    m = m, 
    I = I,
    S = S,
    A = A,
    Y = Y,
    H = H,
    Wd = Wd, 
    Wp = Wp)
  
  # This is equivalent as, but quicker than, a `left_join()`
  # because we want this code to be as fast as possible!
  df$Wr <- NA
  idx = df$t %in% t.obs.ww
  df$Wr[idx] <- Wr
  
  return(df)  
}

#' Helper function to aggregate simulated observations
#' @keywords internal
#' 
helper_aggreg <- function(sim, type, dateobs, prms) {
  
  if(type == 'cl') v = 'Y'
  if(type == 'ha') v = 'H'
  
  a = aggregate_time(
    df       = sim, 
    dt.aggr  = dateobs, 
    var.name = v) %>% 
    dplyr::transmute(
      date, 
      t    = as.integer(date - prms$date.start),
      obs  = as.integer(aggregation))
  
  return(a)
}

#' Simulate a full epidemic with a REEM
#' including clinical and wastewater "observations".
#'
#' @param obj A object of class `reem`
#' @param deterministic Logical. Is the simulation deterministic? (or stochastic) 
#'
#' @return A list
#'
reem_simulate_epi <- function(obj, 
                              deterministic) {
  
  # -- Checks 
  obj$prms = check_date_start(obj)
  obj$prms = check_obs_schedule(obj, verbose = 1)
  
  prms = obj$prms
  
  # Simulate epidemic to generate data
  sim = reem_simulate(prms, deterministic)
  
  sim = dplyr::mutate(sim, date = prms$date.start + t)
  
  # Create dataframes of observations only:
  #  - clinical data
  #  - hospital admissions data
  #  - wastewater data
  # (they may not be observed on the same schedule)
  
  date.obs.cl = prms$date.start + prms$t.obs.cl
  date.obs.ha = prms$date.start + prms$t.obs.ha
  
  # Aggregate clinical reports
  sim.obs.cl = helper_aggreg(sim = sim, 
                             type = 'cl', 
                             dateobs = date.obs.cl, 
                             prms = prms)
  
  # Aggregate hospital admissions
  sim.obs.ha = helper_aggreg(sim = sim, 
                             type = 'ha', 
                             dateobs = date.obs.ha, 
                             prms = prms)
  
  # Extract wastewater observations
  sim.obs.ww = sim %>% 
    dplyr::select(t, Wr) %>% 
    dplyr::filter(t <= prms$last.obs, 
                  !is.na(Wr)) %>%   # filtering out NA keeps observed times only.
    dplyr::rename(obs = Wr) %>%
    dplyr::mutate(date = prms$date.start + t) %>%
    dplyr::select(date, t, obs)
  
  return(list(
    sim    = sim, 
    obs.cl = sim.obs.cl,
    obs.ha = sim.obs.ha,
    obs.ww = sim.obs.ww
  ))
}


