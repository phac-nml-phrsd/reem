

#' Check (and potentially set) the simulation start date
#'
#' @param obj 
#'
#' @return
#' 
check_date_start <- function(obj) {
  has.date.start = !is.null(obj$prms$date.start)
  if(!has.date.start){
    stop('The start date of ',
         'the REEM model has not been specified.',
         ' It MUST be specified in `obj$prms$date.start`\nABORTING.')
  }
}

check_length_timeseries <- function(x, horizon) {
  n = rlang::as_label(eval(parse(text=dplyr::enquo(x)))[2])
  if(length(x) != horizon) 
    message('The size of vector ',n, ' (',
         length(x),
         ') is different from `horizon` (=',
         horizon,')')  
}

#' Set the schedule of times and dates for the observations (if any).
#' Handles `date.start` changes.
#'
#' @param type string. variable type
#' @param obj reem object
#'
#' @return a reem object
#' @keywords internal
#'
set_obs_schedule <- function(type, obj) {
  # DEBUG:  type = 'cl'
  
  vtype    = paste0('obs.',type)
  obs.type = obj[[vtype]]
  n        = nrow(obs.type)
  arg.d    = paste0('date.obs.', type)
  
  # If there are observations attached to 
  # the `reem` object, simply use their schedules
  if(n > 0){
    d = obs.type$date
    t = as.integer(obs.type$date - obj$prms$date.start)
    idx = which(t >= 0)
    if(length(idx)==0){
      stop('`prms$date.start` is larger than the last `',
           type,'` observation date.\nABORTING!')
    }
    obj$prms[[arg.d]] <- d[idx]
  }
  
  # If there are no observations attached
  # and no the dates defined in `prms`, 
  # create an ad-hoc one
  if(n == 0){
    has.d.type = !is.null(obj$prms[[arg.d]])
    if(!has.d.type) {
      t.adhoc = 0:obj$prms$horizon
      obj$prms[[arg.d]] <- obj$prms$date.start + t.adhoc
    }
  }
  return(obj)
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
  nB = length(B)
  
  if(horizon <= ni){
    warning("reem: horizon (=",horizon,") for simulation is shorter than `I.init` (=",ni,"). Cropping vector `I.init`!")
    I.init = I.init[1:(horizon-1)]  # "-1" for at least one simulation step
    ni = horizon - 1
  }
  
  if(nB < horizon){
    B2 = B
    B2 = c(B, rep(B[nB], horizon - nB))
    B  = B2
    warning('Length of parameter B is too short (',nB,').',
         ' It must be at least as long as horizon (',horizon,').\n',
         'B has been extended to match the required length',
         ' duplicating its last element to fill in.\n',
         'You may want to define B with a longer length when building the `reem` object.')
  }
  
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
  
  Wr = rnorm(n    = length(Wp), 
             mean = Wp, 
             sd   = Wp * 0.2)
  
  # Ending
  
  check_length_timeseries(m, horizon)
  check_length_timeseries(I, horizon)
  check_length_timeseries(S, horizon)
  check_length_timeseries(A, horizon)
  check_length_timeseries(Y, horizon)
  check_length_timeseries(H, horizon)
  check_length_timeseries(Wd, horizon)
  check_length_timeseries(Wp, horizon)
  check_length_timeseries(Wr, horizon)
  
  
  df = data.frame(
    t = 1:horizon, 
    m = m, 
    I = I,
    S = S,
    A = A,
    Y = Y,
    H = H,
    Wd = Wd, 
    Wp = Wp, 
    Wr = Wr)
  
  return(df)  
}

#' Helper function to aggregate simulated observations
#' @keywords internal
#' 
helper_aggreg <- function(sim, type, dateobs, prms, var.name = 'obs') {
  
  if(type == 'cl') v = 'I'
  if(type == 'ha') v = 'H'
  
  # Remove simulation data beyond the last observation date
  df = sim[sim$date <= max(dateobs),]
  
  a = aggregate_time(
    df       = df, 
    dt.aggr  = dateobs, 
    var.name = v) %>% 
    dplyr::transmute(
      date, 
      t    = as.integer(date - prms$date.start),
      obs  = as.integer(aggregation))
  
  if(var.name != 'obs'){
    idx = which(names(a) == 'obs')
    names(a)[idx] <- var.name
  } 
  
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
  check_date_start(obj)
  
  # Set observation schedules
  for(vt in c('cl', 'ha', 'ww')) {
    obj = set_obs_schedule(vt,obj)
  }
 
  # DEBUG
  # obj$print_prms()
   
  # Simulate epidemic to generate data
  sim = reem_simulate(prms = obj$prms, deterministic)
  
  # append the dates
  sim = dplyr::mutate(sim, date = obj$prms$date.start + t)
  
  if(0){ # DEBUG 
    print('\nin sim')
    print(head(sim, n=12))
  }
  
  # Create dataframes of "simulated observations" 
  # at the specified observation dates:
  #  - clinical data
  #  - hospital admissions data
  #  - wastewater data
  # (they may not be observed on the same schedule)
  #
  # Note: clinical reports and hospital admissions
  #       are assumed to be AGGREGATED observations
  #       unlike wastewater concentration.
  
  # Clinical reports are already aggregated 
  # if we consider the variable `Y`, so to 
  # calculate the "simulated observations"
  # we just need to retrieve `Y` at the 
  # observation dates defined in `prms`:
  dateobs.cl = obj[['prms']][['date.obs.cl']]
  sim.obs.cl = sim |> 
    dplyr::filter(date %in% dateobs.cl) |> 
    dplyr::transmute(date, t, obs = round(Y))
  
  # Aggregate hospital admissions because the 
  # variable `H` in `sim` is _not_aggregated 
  # (it is a simple proportion of the daily incidence)
  sim.obs.ha = helper_aggreg(sim     = sim, 
                             type    = 'ha', 
                             dateobs = obj[['prms']][['date.obs.ha']], 
                             prms    = obj$prms)
  
  # Extract wastewater observations
  # (wastewater concentration is NOT aggregated in time)
  sim.obs.ww = sim |> 
    dplyr::filter(date %in% obj[['prms']][['date.obs.ww']]) |>
    dplyr::select(date, t, obs = Wr)
  
  return(list(
    sim    = sim, 
    obs.cl = sim.obs.cl,
    obs.ha = sim.obs.ha,
    obs.ww = sim.obs.ww
  ))
}


