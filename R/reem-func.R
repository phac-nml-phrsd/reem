



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



summarize_fcst <- function(simfwd, prm.fcst) {
  
  if(0){ 
    ci = 0.95
  }
  message('\nSummarizing forecasts...',appendLF = FALSE)
  ci = prm.fcst$ci
  
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
  message(' done.')
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


reem_forecast <- function(obj, prm.fcst, verbose ) {
  
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
    update_and_simulate <- function(i, pp, obj, verbose, tpb) {
      
      if(verbose) cat('Simulating forward with posterior sample #',i,'\n')
      setTxtProgressBar(tpb, value = i)
      
      # update fitted parameters 
      # with their posterior values
      obj$prms[names(pp)] <- pp[i,]
      
      # Forward simulations are calculated 
      # for every day in the future (no unobserved date):
      tww = obj$prms$t.obs.ww
      obj$prms$t.obs.ww <- min(tww):max(tww)
      
      # Simulate forward
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
    
    if(ns <= npp) ii = sort(sample(1:npp, ns, replace = FALSE))
    
    message('Sampling ',ns, ' posterior parameter sets out of ', npp,' available.')
    tpb = txtProgressBar(min = 0, max = max(ii), style = 3, width = 50)
    
    simfwd = lapply(X   = ii,
                    FUN = update_and_simulate, 
                    pp  = pp, 
                    obj = obj,
                    verbose = verbose,
                    tpb = tpb)
  }
  
  
  if(! prm.fcst$use.fit.post){
    stop('`use.fit.post = FALSE` is not implemented!')
    # TODO: finish if you think it makes sense 
    # to d it this way...
  }
  
  summary.fcst = summarize_fcst(simfwd, prm.fcst)
  
  return( list(
    simfwd = simfwd, 
    summary.fcst = summary.fcst
  ))
  
}



reem_plot_fit <- function(obj) {
  
  # Prepare dataframes for plotting
 
  fit.obj = obj$fit.obj 
  ps = fit.obj$post.simulations 
  
  ps.cl = lapply(ps, aggregate_time, 
                 dt.aggr = obs.cl$date, 
                 # `Y` is the aggregated clinical reports
                 var.name = 'Y') %>% 
    dplyr::bind_rows() %>% 
    dplyr::group_by(date) %>%
    dplyr::summarise(Y.m = mean(aggregation),
                     Y.lo = min(aggregation),
                     Y.hi = max(aggregation))
  
  ps.ww = ps %>% 
    dplyr::bind_rows() %>% 
    tidyr::drop_na(Wr) %>%
    dplyr::group_by(date) %>%
    # `Wr` is the reported wastewater concentration
    dplyr::summarise(Wr.m = mean(Wr),
                     Wr.lo = min(Wr),
                     Wr.hi = max(Wr))
  
  
  ggplot2::theme_set(ggplot2::theme_bw())
  
  # ---- Trajectories
  
  g.cl = ps.cl %>% 
    ggplot2::ggplot(ggplot2::aes(x=date)) + 
    ggplot2::geom_line(ggplot2::aes(y = Y.m), color = 'chartreuse3')+
    ggplot2::geom_ribbon(ggplot2::aes(ymin=Y.lo, ymax=Y.hi), 
                         alpha=0.2, fill='chartreuse')+
    ggplot2::geom_point(data = obs.cl, ggplot2::aes(y=obs)) +
    labs(title = 'Fit to clinical data', y = 'cases')
  # g.cl
  
  g.ww = ps.ww %>%
    drop_na(starts_with('Wr')) %>%
    ggplot2::ggplot(ggplot2::aes(x=date)) + 
    ggplot2::geom_line(ggplot2::aes(y = Wr.m), 
                       color = 'chocolate3')+
    ggplot2::geom_ribbon(ggplot2::aes(
      ymin = Wr.lo, 
      ymax = Wr.hi), 
      alpha=0.2, fill='chocolate')+
    ggplot2::geom_point(data = obs.ww, ggplot2::aes(y=obs)) +
    labs(title = 'Fit to wastewater data', y = 'concentration')
  # g.ww
  
  # ---- Posterior parameters
  
  # -- 1D density
  
  dp = rbind(mutate(fit.obj$all.distances, type = 'prior'),
             mutate(fit.obj$post.prms, type = 'posterior')) %>%
    select(-starts_with('abc')) %>% 
    pivot_longer(cols = -type)
  
  col.pp =  c(posterior = 'indianred3', 
              prior = 'gray70')
  
  g.post = dp %>% 
    ggplot2::ggplot(ggplot2::aes(x     = value, 
                                 color = type, 
                                 fill  = type)) + 
    ggplot2::geom_density(alpha = 0.3)+
    ggplot2::scale_color_manual(values = col.pp)+
    ggplot2::scale_fill_manual(values = col.pp)+
    ggplot2::facet_wrap(~name, scales = 'free')+
    ggplot2::theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank())+
    ggplot2::labs(title = 'Posterior parameters density')
  
  # -- 2D density
  
  nam = names(fit.obj$post.prms) 
  nam = nam[!grepl('^abc', nam)]
  n = length(nam)
  k = 1 ; gp = list()
  
  for(i in 1:n){
    for(j in 1:n){
      if(i<j){
        tmp =  fit.obj$post.prms[,c(i,j)] 
        names(tmp) = c('x','y')
        gp[[k]] = tmp %>% 
          ggplot2::ggplot(ggplot2::aes(x = x, y = y))+
          ggplot2::geom_density_2d_filled()+
          ggplot2::theme(panel.grid = ggplot2::element_blank())+
          ggplot2::labs(x=nam[i], y=nam[j]) + 
          ggplot2::guides(fill = 'none')
        
        k = k+1
      }
    }
  }
  gpall2d = patchwork::wrap_plots(gp) +
    patchwork::plot_annotation(title = 'Posterior parameters 2D density')
  
  
  # -- Ordered ABC distances
  
  fit.prm = obj$fit.prm
  
  n.post = round(fit.prm$n.abc * fit.prm$p.abc, 0)
  d = fit.obj$all.distances %>%
    dplyr::mutate(i = dplyr::row_number()) %>%
    dplyr::mutate(type = ifelse(i <= n.post,
                                'accepted','rejected'))
  
  g.dist = d %>%
    ggplot2::ggplot(ggplot2::aes(x = 1:nrow(d), 
                                 y = abc.err)) + 
    ggplot2::geom_vline(xintercept = n.post, 
                        linetype = 'dashed')+
    ggplot2::geom_step(ggplot2::aes(color = type), 
                       linewidth = 1 )+
    ggplot2::scale_y_log10() + 
    ggplot2::scale_x_log10() + 
    ggplot2::scale_color_manual(values = c('red2', 'gray80'))+
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) + 
    ggplot2::labs(
      title = 'ABC distances from data',
      x = 'ordered ABC iteration',
      y = 'distance'
    )
  
  
  g.list = list(
    traj.cl = g.cl,
    traj.ww = g.ww,
    post.prms = g.post,
    post.prms.2d = gpall2d,
    dist = g.dist
  )
  g.all = patchwork::wrap_plots(g.list) 
  
  res = c(g.list, list(all = g.all))
  return(res) 
}



#' Plot forecasts
#'
#' @param obj 
#' @param date_breaks 
#' @param date_labels 
#'
#' @return A ggplot object.
#'
reem_plot_forecast <- function(
    obj,
    date_breaks,
    date_labels  ) {
  
  obs.cl = obj$obs.cl
  obs.ww = obj$obs.ww
  fcst.obj = obj$fcst.obj
  fcst.prm = obj$fcst.prm
  
  # - - - Cosmetics  
  
  col.fcst = 'steelblue2'
  col.fit  = 'tan2'
  xaxis = scale_x_date(
    date_breaks = date_breaks, 
    date_labels = date_labels)
  alpha.ribbon = 0.20
  
  # - - - - Retrieve fitted curves
  
  sf = fcst.obj$summary.fcst 
  
  fitsim.ww = obj$fit.obj$post.simulations %>% 
    bind_rows() %>% 
    group_by(date) %>% 
    summarize(m  = mean(Wr), 
              lo = min(Wr), 
              hi = max(Wr)) %>% 
    drop_na(m)
  
  tmp.cl = obj$fit.obj$post.simulations %>% 
    bind_rows() %>% 
    group_by(date) %>% 
    summarize(m = mean(Y), 
              lo = min(Y), 
              hi = max(Y)) %>% 
    drop_na(m)
  
  # Aggregate for clinical reports
  fitsim.cl = tmp.cl %>% 
    aggcl(dt.aggr = obs.cl$date, 
          vars = c('m','lo','hi')) %>%
    filter(date <= max(obs.cl$date))
  
  # set the aggregation dates as 
  # starting from `asof` and a time interval
  # equal to the one of the observations:
  dt = as.numeric(diff(obs.cl$date))[1]
  dt.aggr.fcst = seq.Date(from = prm.fcst$asof , 
                          to = max(sf$date), 
                          by = dt)
  
  # aggregation of clinical reports
  sf.cl = aggcl(df = sf, 
                dt.aggr = dt.aggr.fcst, 
                vars = c('Y_mean','Y_lo','Y_hi')) %>% 
    filter(date >= prm.fcst$asof)
  
  # - - - Plots - - - 
  
  g.cl = ggplot(data = sf.cl,
                aes(x=date))+ 
    geom_line(data = fitsim.cl, aes(y=m), 
              color = col.fit, linetype = 'dashed') + 
    geom_ribbon(data = fitsim.cl, aes(ymin=lo, ymax=hi), 
                fill=col.fit, alpha = alpha.ribbon / 2) + 
    geom_point(data = obs.cl, aes(y=obs))+ 
    geom_line( aes(y = Y_mean), color= col.fcst, 
               linetype = 'dotted') + 
    geom_ribbon(aes(ymin = Y_lo, ymax = Y_hi), 
                alpha = alpha.ribbon, 
                fill= col.fcst,
                color= col.fcst) +
    geom_vline(xintercept = fcst.prm$asof, 
               linetype = 'dashed', 
               color = 'gray50') + 
    annotate(geom = 'text', y=1, x=fcst.prm$asof, 
             label = fcst.prm$asof, size = 2) + 
    xaxis + 
    labs(title = 'Forecast clinical reports', 
         x = '', y = 'cases')
  # g.cl
  
  g.ww = sf %>% 
    drop_na(Wr_mean) %>%
    filter(date >= fcst.prm$asof) %>%
    ggplot(aes(x=date))+ 
    geom_line(data = fitsim.ww, aes(y=m), 
              color = col.fit, linetype = 'dashed') + 
    geom_ribbon(data = fitsim.ww, aes(ymin=lo, ymax=hi), 
                fill=col.fit, alpha = alpha.ribbon / 2) + 
    geom_point(data = obs.ww, aes(y=obs))+ 
    geom_line( aes(y = Wr_mean), color= col.fcst, 
               linetype = 'dotted') + 
    geom_ribbon(aes(ymin = Wr_lo, ymax = Wr_hi), 
                alpha = alpha.ribbon, 
                fill= col.fcst,
                color= col.fcst) +
    geom_vline(xintercept = fcst.prm$asof, 
               linetype = 'dashed', 
               color = 'gray50') + 
    annotate(geom = 'text', y=1, x=fcst.prm$asof, 
             label = fcst.prm$asof, size = 2) + 
    labs(title = 'Forecast wastewater concentration', 
         x = '', y = 'concentration') +
    xaxis 
  # g.ww
  
  return(list(
    cl = g.cl, 
    ww = g.ww
  ))
  
}









