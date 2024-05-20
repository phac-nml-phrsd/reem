
#' @title sqrt sum of squared difference
#'
#' @param x 
#' @param y 
#'
#' @return Numerical. The sum of squared difference
#' @keywords internal
#'  
sss <- function(x,y) {
  return(sqrt( sum( (x-y)^2 ) ))
}


#' Normalization with one of the largest values
#'
#' @param x Vector of positive numerical values.
#'
#' @return Normalized vector.
#'
normlarge <- function(x, largest.rank = 3) {
  if(!all(x>=0)) stop('`normlarge` is implemented for positive vectors only.')
  # Element used for normalization
  idx = min(length(x), largest.rank)
  val = sort(x, decreasing = TRUE)[idx]
  return(val)
}


#' Error function for a fit algorithm
#'
#' @param cl.i Dataframe of simulated clinical observations
#' @param ww.i Dataframe of simulated wastewater observations
#' @param obs.cl Dataframe of clinical observations
#' @param obs.ww Dataframe of wastewater observations
#' @param use.cl Logical. Use clinical data in fit?
#' @param use.ww Logical. Use wastewater data in fit?
#'
#' @return Numerical value representing the distance 
#' of the simulation from the observed data.
#' 
err_fct <- function(cl.i, ww.i, 
                    obs.cl, obs.ww,
                    use.cl, use.ww,
                    err.type = 'rel',
                    do.plot = FALSE) {
  
  # --- Checks  
  
  if(! all(obs.cl$date[-1] %in% cl.i$date)){
    print('observed cl:')
    print(obs.cl$date)
    print('simulated cl (ABC internal):')
    print(cl.i$date)
    stop('Date mismatch: clinical observation dates are missing in simulation.')
  }
  if(!all(obs.ww$date[-1] %in% ww.i$date)){
    print('observed ww:')
    print(obs.ww$date)
    print('simulated ww (ABC internal):')
    print(ww.i$date)
    stop('Date mismatch: wastewater observation dates are missing in simulation.')
  }
  
  # --- Adjust to fitting observation dates
  
  df.cl = dplyr::left_join(cl.i, obs.cl, by='date')
  df.ww = dplyr::left_join(ww.i, obs.ww, by='date') %>% 
    tidyr::drop_na(obs)
  
  err.cl = 0
  err.ww = 0
  
  if(err.type == 'L2'){
    if(use.cl) err.cl = sss(df.cl$obs, df.cl$Ym) 
    if(use.ww) err.ww = sss(df.ww$obs, df.ww$Wm)
  }
  if(err.type == 'rel'){
    if(use.cl) err.cl = sqrt(sum( (df.cl$Ym/df.cl$obs - 1 )^2 )) #TODO: handle div by 0
    if(use.ww) err.ww = sqrt(sum( (df.ww$Wm/df.ww$obs - 1 )^2 ))
  }
  # The vectors are normalized by their respective
  # largest value to put both ww and clinical 
  # on the same scale:
  if(err.type == 'normlarge'){
    largest.rank = 3
    if(use.cl) {
      z.cl = normlarge(df.cl$obs, largest.rank ) 
      err.cl = sss(df.cl$obs/z.cl , df.cl$Ym/z.cl)
    }
    if(use.ww) {
      z.ww = normlarge(df.ww$obs, largest.rank) 
      err.ww = sss(df.ww$obs/z.ww , df.ww$Wm/z.ww)
    }
  }
  
  # WARNING -- FIXME??
  # The errors calculated above have a number of elements
  # equal to nrow(cl.i), which varies as delta.start changes.
  # Hence, the error may not be consistent across all priors.
  
  err =  err.cl + err.ww
  return(err)
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
#' @return List containing:
#' \item{distance: }{Numerical. The distance ofthe simulation from observations}
#' \item{sim: }{Dataframe of the simulations.}
#' 
#' @export
#'
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
  
  # Setting the horizon for each data source
  if(nrow(obs.cl) > 0)  datemax.cl = max(obs.cl$date)
  if(nrow(obs.ww) > 0)  datemax.ww = max(obs.ww$date)
  if(nrow(obs.ww) == 0) datemax.ww = datemax.cl
  if(nrow(obs.cl) == 0) datemax.cl = datemax.ww
  
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
    dplyr::filter(date <= datemax.cl)
  
  ww.i = a.ww %>% 
    dplyr::group_by(date) %>% 
    dplyr::summarize(Wm = mean(obs)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(date <= datemax.ww)
  
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
  if(i%%50 == 0) cat('ABC iteration #',i,'\n') 
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


#' Generate prior samples given 
#' parameter names and ranges.
#'
#' @param prms.to.fit List of parameters name with their 
#' prior distribution and parameters
#' @param n.priors Integer. Number of sample to generate for each parameter.
#'
#' @return Dataframe of prior samples (each column 
#' represents a variable).
#' 
#' @export
#'
generate_priors <- function(prms.to.fit, n.priors) {
  tmp = list()
  
  for(i in seq_along(prms.to.fit)){
    x = prms.to.fit[[i]]
    distrib = x[[1]]
    
    if(distrib == 'unif') 
      tmp[[i]] = runif(n = n.priors, min = x[[2]], max = x[[3]])
    
    if(grepl('^norm', distrib)){
      y = rnorm(n = n.priors, mean = x[[2]], sd = x[[3]])
      if(distrib == 'normp') y[y<0] <- -y[y<0] # force positive samples
      tmp[[i]] = y
    }
    
    if(distrib == 'exp')
      tmp[[i]] = rexp(n = n.priors, rate = 1 / x[[2]])
    
    if(distrib == 'unif_int')
      tmp[[i]] = sample(x = x[[2]]:x[[3]], size = n.priors, replace = TRUE)
    
    if(distrib == 'gamma'){
      # Note: Gamma is assumed parametrized with 
      # mean and variance, not shape and scale.
      m = x[[2]] ; v = x[[3]]
      scale = v / m
      shape = m^2 / v
      tmp[[i]] = rgamma(n = n.priors, shape = shape, scale = scale)
    }
  }
  priors = data.frame(tmp)
  names(priors) = names(prms.to.fit)
  return(priors)
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
  
  d.max = max(obj$obs.cl$date, obj$obs.ww$date) + 1
  hz    = as.integer(d.max - obj$prms$date.start)
  
  if(hz <= 1) {
    stop('Horizon is too short!\n',
         'horizon = ', hz, '\n',
         'max date obs cl & ww : ', d.max,'\n',
         'model start date     : ',obj$prms$date.start,
         '\n ABORTING!\n')
  }
  
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
          '\nData horizon : ', hz, ' (days) \n',
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
  snowfall::sfLibrary(purrr)
  
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
    dplyr::mutate(abc.index = 1:nrow(priors)) %>%
    dplyr::arrange(abc.err)
  
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


reem_plot_fit <- function(obj) {
  
  # Prepare dataframes for plotting
  
  fit.obj = obj$fit.obj 
  ps      = fit.obj$post.simulations 
  obs.cl  = obj$obs.cl
  obs.ww  = obj$obs.ww
  
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
    ggplot2::labs(title = 'Fit to clinical data', y = 'cases')
  # g.cl
  
  g.ww = ps.ww %>%
    tidyr::drop_na(starts_with('Wr')) %>%
    ggplot2::ggplot(ggplot2::aes(x=date)) + 
    ggplot2::geom_line(ggplot2::aes(y = Wr.m), 
                       color = 'chocolate3')+
    ggplot2::geom_ribbon(ggplot2::aes(
      ymin = Wr.lo, 
      ymax = Wr.hi), 
      alpha=0.2, fill='chocolate')+
    ggplot2::geom_point(data = obs.ww, ggplot2::aes(y=obs)) +
    ggplot2::labs(title = 'Fit to wastewater data', y = 'concentration')
  # g.ww
  
  # ---- Posterior parameters
  
  # -- 1D density
  
  dp = rbind(
    dplyr::mutate(fit.obj$all.distances, type = 'prior'),
    dplyr::mutate(fit.obj$post.prms, type = 'posterior')) %>%
    dplyr::select(-starts_with('abc')) %>% 
    tidyr::pivot_longer(cols = -type)
  
  col.pp =  c(posterior = 'indianred3', 
              prior     = 'gray70')
  
  g.post = dp %>% 
    ggplot2::ggplot(ggplot2::aes(x     = value, 
                                 color = type, 
                                 fill  = type)) + 
    ggplot2::geom_density(alpha = 0.3)+
    ggplot2::scale_color_manual(values = col.pp)+
    ggplot2::scale_fill_manual(values = col.pp)+
    ggplot2::facet_wrap(~name, scales = 'free')+
    ggplot2::theme(
      axis.text.y      = ggplot2::element_blank(),
      axis.ticks.y     = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
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
          ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.text  = ggplot2::element_text(size = ggplot2::rel(0.6))
          )+
          ggplot2::labs(x=nam[i], y=nam[j]) + 
          ggplot2::guides(fill = 'none')
        
        k = k+1
      }
    }
  }
  
  
  gpall2d = patchwork::wrap_plots(gp) +
    patchwork::plot_annotation(title = 'Posterior parameters 2D density')
  
  gpall2d
  
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
    ggplot2::scale_color_manual(values = c('indianred2', 'gray90'))+
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) + 
    ggplot2::labs(
      title = 'ABC distances from data',
      x = 'ordered ABC iteration',
      y = 'distance'
    )
  
  g.list = list(
    traj.cl      = g.cl,
    traj.ww      = g.ww,
    post.prms    = g.post,
    post.prms.2d = gpall2d,
    dist         = g.dist
  )
  g.all = patchwork::wrap_plots(g.list) 
  
  res = c(g.list, list(all = g.all))
  return(res) 
}

