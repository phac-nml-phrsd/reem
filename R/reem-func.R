



reem_create <-  function(){
  return(
    new(Class     = 'reem', 
        name      = 'unnamed_reem',
        prms      = list(),
        is.fitted = FALSE,
        obs.cl    = data.frame(),
        obs.ww    = data.frame()))
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
    w.m[t] = sum(fec[idx]*I[t-idx])
  }
  Wd = w.m
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
  
  tmp = data.frame(t = t.obs.ww, Wr = Wr)
  
  df = data.frame(
    t = 1:horizon, 
    m = m, 
    I = I,
    S = S,
    A = A,
    Y = Y,
    Wd = Wd, 
    Wp = Wp)
  
  df = dplyr::left_join(df, tmp, by='t')
  
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
reem_simulate_epi <- function(prms, 
                              deterministic) {
  
  if(is.null(prms$t.obs.ww)){
    # generate times when ww is observed:
    t.obs.ww  = cumsum(1 + rpois(n=prms$horizon, 
                                 lambda= prms$freq.obs.ww-1))
    t.obs.ww  = t.obs.ww[t.obs.ww < prms$horizon]
    prms = c(prms, t.obs.ww = list(t.obs.ww))
  }
  
  # Simulate epidemic to generate data
  sim = reem_simulate(prms, deterministic)
  
  # When working with real data, 
  # we usually are in "date mode"
  has.date.start = !is.null(prms$date.start)
  
  # Build times
  if(!has.date.start){
    prms$date.start <- lubridate::ymd('2020-01-01')
    # TODO: throw a warning..
  }
  
  if(is.null(prms$t.obs.cl)){
    # generate times when clinical is observed:
    tmax = max(sim$t)
    t.obs.cl = prms$lag * c(1:999)
    t.obs.cl = t.obs.cl[t.obs.cl < tmax]
  }
  if(!is.null(prms$t.obs.cl)) t.obs.cl = prms$t.obs.cl
  
  sim = dplyr::mutate(sim, date = prms$date.start + t)
  
  # Create two dataframes of observations only:
  #  - clinical data
  #  - wastewater data
  # (they may not be observed on the same schedule)
  
  
  # Aggregate clinical reports
  sim.obs.cl = aggregate_time(df       = sim, 
                              dt.aggr  = t.obs.cl, 
                              var.name = 'Y') %>% 
    dplyr::transmute(
      t, 
      date = prms$date.start + t,
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





#' Calculate the distance of the 
#' epidemic trajectory from 
#' clinical and wastewater observations
#' for a given set of REEM parameters.
#'
#' @param obs.cl 
#' @param obs.ww 
#' @param use.cl 
#' @param use.ww 
#' @param err.type 
#' @param deterministic 
#' @param n.sim 
#' @param verbose 
#' @param prms 
#'
#' @return
#' @export
#'
#' @examples
#' 
reem_traj_dist_obs <- function(
    obj,
    use.cl, 
    use.ww, 
    err.type,
    deterministic,
    n.sim = 10, 
    verbose = FALSE
) {
  
  prms   = obj$prms
  obs.cl = obj$obs.cl
  obs.ww = obj$obs.ww
  
  if(nrow(obs.cl)==0 & nrow(obs.ww)==0) {
    stop('The REEM object does not have any observation attached.
         Hence, cannot calculate a distance from observation. ABORTING!')
  }
  
  # Adjust the epidemic start time
  # and check if we did not create
  # negative times inadvertently:
  prms$date.start <- prms$date.start + prms$start.delta
  
  if(verbose){
    cat('-- Distance debug\n')
    cat('\nprms$N:', prms$N)
    cat('\nprms$i0prop:', prms$i0prop)
    cat('\nprms$I.init:', paste(prms$I.init, collapse = '; '))
    cat('\nprms$date.start:', as.character(prms$date.start))
  }
  
  if(deterministic){
    s     = obj$simulate_epi(deterministic = TRUE)
    a.cl  = s$obs.cl
    a.ww  = s$obs.ww
    a.sim = s$sim
  }
  
  if(!deterministic){
    # Simulate `n.sim` times with a given set of prior parameters.
    # The ABC distance from the observations will be computed 
    # using the _mean_ value across the `n.sim` simualtions:
    tmp.cl = tmp.ww = tmp.sim = list()
    
    for(k in 1:n.sim){
      s            = obj$simulate_epi(deterministic = FALSE)
      tmp.cl[[k]]  = s$obs.cl %>% mutate(iter = k)
      tmp.ww[[k]]  = s$obs.ww %>% mutate(iter = k)
      tmp.sim[[k]] = s$sim %>% mutate(n.sim = k)
    }
    a.cl  = dplyr::bind_rows(tmp.cl)
    a.ww  = dplyr::bind_rows(tmp.ww)
    a.sim = dplyr::bind_rows(tmp.sim)
  }
  
  # only the "observed" variables are averaged (e.g., not `a.sim`)
  # for comparison with actual observations (see further down)
  cl.i = a.cl %>% 
    dplyr::group_by(date) %>% 
    dplyr::summarize(Ym = mean(obs)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(date <= max(obs.cl$date))
  
  ww.i = a.ww %>% 
    dplyr::group_by(date) %>% 
    dplyr::summarize(Wm = mean(obs)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(date <= max(obs.ww$date))
  
  # Calculate the ABC distance
  res = err_fct(cl.i, ww.i, 
                obs.cl, obs.ww, 
                use.cl, use.ww,
                err.type) 
  return(list(
    distance = res, 
    sim = a.sim
  ))
  
  
}


# Util function for parallel computation
calc_dist_parallel <- function(i,
                               priors,
                               obj,
                               use.cl, 
                               use.ww, 
                               err.type,
                               deterministic,
                               n.sim = 10, 
                               verbose = FALSE ) {
  cat('ABC iteration #',i,'\n') 
  pp = priors[i,]
  
  if('i0prop' %in% names(priors)){
    obj$prms$i0prop <- pp$i0prop
    obj$prms = set_I_init(obj$prms)
  }
  
  obj$prms[names(pp)] <- pp
  
  x = reem_traj_dist_obs(
    obj           = obj,
    use.cl        = use.cl, 
    use.ww        = use.ww, 
    err.type      = err.type,
    deterministic = deterministic,
    n.sim         = n.sim, 
    verbose       = verbose
  )
  return(x)
}

#' Fit model to observed data using the ABC algorithm.
#'
#' @param obj 
#' @param prm.abc 
#' @param prms.to.fit 
#'
#' @return
#' @export
#'
#' @examples
#' 
reem_fit_abc <- function(obj,
                         prm.abc,
                         prms.to.fit,
                         verbose = FALSE) {
  
  # Unpack ABC parameters:
  n.abc   = prm.abc$n.abc
  n.sim   = prm.abc$n.sim
  p.abc   = prm.abc$p.abc
  n.cores = prm.abc$n.cores
  
  use.ww = prm.abc$use.ww
  use.cl = prm.abc$use.cl
  
  d.max = max(obj$obs.cl$date, obj$obs.ww$date)
  hz    = as.integer(d.max - prms$date.start)
  
  obj$prms$horizon <- hz 
  
  message('\n----- ABC FIT -----\n\n',
          'Target data sources :\n',
          '  clinical   = ', ifelse(use.cl,'yes', 'NO'),'\n',
          '  wastewater = ', ifelse(use.ww,'yes', 'NO'),'\n\n',
          'Number of priors     : ', prm.abc$n.abc, '\n',
          'Number of posteriors : ', 
          round(prm.abc$p.abc * prm.abc$n.abc,0),
          ' (accept ratio = ',prm.abc$p.abc,')\n',
          'Number of cores      : ', prm.abc$n.cores,
          ' (', round(prm.abc$n.abc / prm.abc$n.cores,0),
          ' iters per core)\n',
          '\nData horizon : ', hz, '\n',
          '\n---------------------\n\n')
  
  err.type = prm.abc$err.type
  
  deterministic = ifelse(n.sim > 0, FALSE, TRUE)
  
  # Priors 
  priors = generate_priors(prms.to.fit, n.priors = n.abc) 
  
  # calculate distance from observations
  # for each parameter prior value
  snowfall::sfInit(parallel = n.cores > 1, cpus = n.cores)
  snowfall::sfExportAll()
  snowfall::sfLibrary(dplyr)
  
  z = snowfall::sfLapply(
    x         = 1:n.abc, 
    fun       = calc_dist_parallel,
    priors    = priors, 
    obj       = obj,
    use.cl    = use.cl, 
    use.ww    = use.ww,
    n.sim     = n.sim,
    verbose   = verbose,
    err.type  = err.type,
    deterministic = deterministic)
  snowfall::sfStop()
  
  abc.err = sapply(z, '[[', 'distance')
  abc.sim = lapply(z, '[[', 'sim') %>% 
    # Add the simulation number (`index`) 
    # to each iteration:
    purrr::imap(~ dplyr::mutate(.x, index = .y))
  
  # -- Posteriors
  
  df.abc = cbind(priors, abc.err) %>% 
    mutate(abc.index = 1:nrow(priors)) %>%
    arrange(abc.err)
  
  n.post  = round(n.abc*p.abc)
  df.post = df.abc[1:n.post,]
  
  res = list(
    all.distances    = df.abc,
    all.simulations  = abc.sim,
    post.prms        = df.post,
    post.simulations = abc.sim[df.post$abc.index]
  ) 
  
  return(res)
}



summarize_fcst <- function(simfwd, prm) {
  
  if(0){ 
    ci = 0.95
  }
  
  ci = prm$ci
  
  
  res = simfwd %>% 
    bind_rows() %>% 
    group_by(date) %>%
    summarise(
      across(.cols = c(Y,Wr), 
             .fns = mean,
             .names = '{.col}_mean',
             na.rm = TRUE),
      across(.cols = c(Y,Wr), 
             .fns = quantile,
             probs = 0.5 - ci/2, 
             .names = '{.col}_lo',
             na.rm = TRUE),
      across(.cols = c(Y,Wr), 
             .fns = quantile,
             probs = 0.5 + ci/2, 
             .names = '{.col}_hi',
             na.rm = TRUE)
    )
  
  return(res)
 
  # Mon May 29 11:26:27 2023 ------------------------------
  # do not delete. try to find an elegant way 
  # to return multiple CIs (for nice plots)
  if(0){
    quantile_df <- function(x, probs = c(0.25, 0.5, 0.75)) {
      tibble(
        val = quantile(x, probs, na.rm = TRUE),
        quant = probs
      )
    }
    
    tmp = simfwd %>% 
      bind_rows() %>% 
      reframe(
        across(c(Y, Wr), quantile_df, .unpack = TRUE),
        .by = date
      )
    
    tmp %>% 
      ggplot(aes(x=date)) + 
      geom_line(aes(y=Y_mean))+ 
      geom_ribbon(aes(ymin = Y_lo, ymax = Y_hi), alpha=0.2)
    
    
    tmp %>%
      drop_na(starts_with('Wr')) %>%
      ggplot(aes(x=date)) + 
      geom_point(aes(y=Wr_mean))+ 
      geom_ribbon(aes(ymin = Wr_lo, ymax = Wr_hi), alpha=0.2)   
  }
  
}


reem_forecast <- function(obj, prm.fcst) {
  
  if(0){ # DEBUG
    
    prm.fcst = list(
      asof = ymd('2022-03-01'),
      horizon.fcst = ymd('2022-07-01'),
      use.fit.post = TRUE,
      n.resample = 20,
      ci = 0.95
    )
    
  }
  
  a = obj$fit.obj
  
  # Extend the horizon to match the one requested for the forecast
  obj$prms$horizon <- as.integer(prm.fcst$horizon.fcst - obj$prms$date.start)
  
  # In this case, the forecast is simply reusing
  # the simulations from the posterior fits
  # (i.e., no resampling from posterior distributions)
  if(prm.fcst$use.fit.post){
    
    pp = a$post.prms %>% dplyr::select(!tidyr::starts_with('abc'))
   
    # Helper function 
    update_and_simulate <- function(i, pp, obj) {
      cat('Simulating forward with posterior sample #',i,'\n')
      obj$prms[names(pp)] <- pp[i,]
      s = obj$simulate_epi(deterministic = TRUE) #TODO: let user choose
      s$sim$index <- i
      return(s$sim)
    }
    
    npp = nrow(pp)
    ns  = prm.fcst$n.resample
    
    if(ns > npp) {
      warning('
      Number of samples required for forecast (',ns,') is larger
      than the total number of posteriors (',npp,').
      Because option to resample from posterior is selected, 
      only ',npp,' samples will be used for the forecast.')
      
      ii = 1:npp
    }
    
    if(ns <= npp) ii = sample(1:npp, ns, replace = FALSE)
    
    simfwd = lapply(X   = ii,
                    FUN = update_and_simulate, 
                    pp  = pp, 
                    obj = obj)
  }
  
  
  # Here, we resample from the posterior distributions.
  # Assume multidimensional normal distribution
  # to account for correlations between variables
  
  if(! prm.fcst$use.fit.post){
    
    #pp = a$post.prms
    
    # TODO: finish if you think it makes sense 
    # to d it this way...
  }
  
  summary.fcst = summarize_fcst(simfwd, prm.fcst)
  
  return( list(
    simfwd = simfwd, 
    summary.fcst = summary.fcst
  ))
  
}










