
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

#' @title sqrt sum of squared relative difference
#'
#' @param x simulations
#' @param y observations
#'
#' @return Numerical. The sum of squared relative difference
#' @keywords internal
#'  
ssrs <- function(x,y) {
  # deal with zeros:
  idx = which(y!=0)
  yy = y[idx]
  xx = x[idx]
  return(sqrt( sum( (xx/yy-1)^2 ) ))
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


check_dates_cl <- function(obs.cl, cl.i) {
  if(! all(obs.cl$date[-1] %in% cl.i$date)){
    print('observed cl:')
    print(obs.cl$date)
    print('simulated cl (ABC internal):')
    print(cl.i$date)
    stop('Date mismatch: clinical observation dates are missing in simulation.')
  }
}


check_dates_ha <- function(obs.ha, ha.i) {
  if(! all(obs.ha$date[-1] %in% ha.i$date)){
    print('observed hosp. adm.:')
    print(obs.ha$date)
    print('simulated hosp. adm. (ABC internal):')
    print(ha.i$date)
    stop('Date mismatch: hospital admission dates are missing in simulation.')
  }
}

check_dates_ww <- function(obs.ww, ww.i) {
  
  if(!all(obs.ww$date[-1] %in% ww.i$date)){
    print('observed ww:')
    print(obs.ww$date)
    print('simulated ww (ABC internal):')
    print(ww.i$date)
    
    is.in = obs.ww$date[-1] %in% ww.i$date
    is.missing = obs.ww$date[-1][which(!is.in)]
    stop(paste('Date mismatch:',
               'wastewater observation dates are missing in simulation.\n',
               'Here are the dates missing:\n',
               paste(is.missing, sep=' ; ')))
  }
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
err_fct <- function(cl.i, ha.i, ww.i, 
                    obs.cl, obs.ha, obs.ww,
                    use.cl, use.ha, use.ww,
                    err.type = 'rel',
                    do.plot = FALSE) {
  
  # --- Checks  
  
  if(use.cl) check_dates_cl(obs.cl, cl.i)
  if(use.ha) check_dates_ha(obs.ha, ha.i)
  if(use.ww) check_dates_ww(obs.ww, ww.i)
  
  # --- Adjust to fitting observation dates
  
  if(use.cl) df.cl = dplyr::left_join(cl.i, obs.cl, by='date') |> tidyr::drop_na(obs)
  if(use.ha) df.ha = dplyr::left_join(ha.i, obs.ha, by='date') |> tidyr::drop_na(obs)
  if(use.ww) df.ww = dplyr::left_join(ww.i, obs.ww, by='date') |> tidyr::drop_na(obs)
  
  if(0){ # DEBUG 
    message('\nDEBUG err_fct:')
    print((df.cl))
    print((df.ww))
    print((df.ha))
  }
  
  err.cl = 0
  err.ha = 0
  err.ww = 0
  
  if(err.type == 'L2'){
    if(use.cl) err.cl = sss(df.cl$obs, df.cl$Ym) 
    if(use.ha) err.ha = sss(df.ha$obs, df.ha$Hm) 
    if(use.ww) err.ww = sss(df.ww$obs, df.ww$Wm)
  }
  if(err.type == 'rel'){
    #TODO: handle div by 0
    if(use.cl) err.cl = ssrs(df.cl$Ym, df.cl$obs) 
    if(use.ha) err.ha = ssrs(df.ha$Hm, df.ha$obs) 
    if(use.ww) err.ww = ssrs(df.ww$Wm, df.ww$obs) 
  }
  # The vectors are normalized by their respective
  # largest value to have all data sources 
  # on the same scale:
  if(err.type == 'normlarge'){
    largest.rank = 3
    
    if(use.cl) {
      z.cl = normlarge(df.cl$obs, largest.rank ) 
      err.cl = sss(df.cl$obs/z.cl , df.cl$Ym/z.cl)
    }
    
    if(use.ha) {
      z.ha = normlarge(df.ha$obs, largest.rank ) 
      err.ha = sss(df.ha$obs/z.ha , df.ha$Hm/z.ha)
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
  
  err =  err.cl + err.ha + err.ww
  
  res = c(
    total = err, 
    cl = err.cl, 
    ha = err.ha, 
    ww = err.ww
  )
  
  return(res)
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
    use.ha, 
    use.ww, 
    err.type,
    deterministic,
    n.sim = 10, 
    verbose = FALSE
) {
  
  prms = obj$prms
  
  # Adjust the epidemic start time and horizon
  date.start.new  = prms$date.start + prms$start.delta
  horizon.new     = prms$horizon - prms$start.delta  # `-` because we want to keep the same _date_
  
  # inventory of observations
  has.cl = nrow(obj$obs.cl)>0
  has.ha = nrow(obj$obs.ha)>0
  has.ww = nrow(obj$obs.ww)>0
  
  obs.cl = obs.ha = obs.ww = data.frame()
  
  # Crop observations that occurred before the (new) start date
  if(has.cl) obs.cl = obj$obs.cl |> dplyr::filter(date > date.start.new)
  if(has.ha) obs.ha = obj$obs.ha |> dplyr::filter(date > date.start.new)
  if(has.ww) obs.ww = obj$obs.ww |> dplyr::filter(date > date.start.new)
  
  if(!has.cl & !has.ha & !has.ww) {
    stop('The REEM object does not have any observation attached.\n',
         'Hence, cannot calculate a distance from observation.\nABORTING!\n\n')
  }
  
  # Setting the horizon for each data source
  if(has.cl)  datemax.cl = max(obs.cl$date)
  if(has.ha)  datemax.ha = max(obs.ha$date)
  if(has.ww)  datemax.ww = max(obs.ww$date)
  if(!has.cl) datemax.cl = datemax.ww
  if(!has.ha) datemax.ha = datemax.cl
  if(!has.ww) datemax.ww = datemax.cl
  
  if(verbose){
    cat('-- Distance debug\n')
    cat('\nprms$N:', prms$N)
    cat('\nprms$i0prop:', prms$i0prop)
    cat('\nprms$I.init:', paste(prms$I.init, collapse = '; '))
    cat('\nprms$date.start:', as.character(prms$date.start))
  }
  
  obj.x = obj$copy()
  obj.x$obs.cl <- obs.cl
  obj.x$obs.ha <- obs.ha
  obj.x$obs.ww <- obs.ww
  obj.x$prms$date.start <- date.start.new
  obj.x$prms$horizon <- horizon.new
  
  if(deterministic){
    s     = obj.x$simulate_epi(deterministic = TRUE)
    a.cl  = s$obs.cl
    a.ha  = s$obs.ha
    a.ww  = s$obs.ww
    a.sim = s$sim
    if(0){ # DEBUG 
      print('in fit')
      print(head(a.sim))
    }
  }
  
  if(!deterministic){
    # Simulate `n.sim` times with a given set of prior parameters.
    # The ABC distance from the observations will be computed 
    # using the _mean_ value across the `n.sim` simualtions:
    tmp.cl = tmp.ha = tmp.ww = tmp.sim = list()
    
    for(k in 1:n.sim){
      s            = obj.x$simulate_epi(deterministic = FALSE)
      tmp.cl[[k]]  = s$obs.cl %>% dplyr::mutate(iter = k)
      tmp.ha[[k]]  = s$obs.ha %>% dplyr::mutate(iter = k)
      tmp.ww[[k]]  = s$obs.ww %>% dplyr::mutate(iter = k)
      tmp.sim[[k]] = s$sim %>% dplyr::mutate(n.sim = k)
    }
    a.cl  = dplyr::bind_rows(tmp.cl)
    a.ha  = dplyr::bind_rows(tmp.ha)
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
  
  ha.i = a.ha %>% 
    dplyr::group_by(date) %>% 
    dplyr::summarize(Hm = mean(obs)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(date <= datemax.ha)
  
  ww.i = a.ww %>% 
    dplyr::group_by(date) %>% 
    dplyr::summarize(Wm = mean(obs)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(date <= datemax.ww)
  
  # The shift from `date.start` may have 
  # changed the horizon DATE before the 
  # last observations. 
  # Cropping again the tail end
  if(has.cl) obs.cl = dplyr::filter(obs.cl, date <= max(cl.i$date))
  if(has.ha) obs.ha = dplyr::filter(obs.ha, date <= max(ha.i$date))
  if(has.ww) obs.ww = dplyr::filter(obs.ww, date <= max(ww.i$date))
  
  # Plot for debugging
  if(0){
    g.cl = ggplot(cl.i, aes(x=date))+ 
      geom_point(data = obs.cl, aes(y=obs)) + 
      geom_line(aes(y=Ym), color = 'steelblue') + 
      labs(title = 'clinical (Ym)')
    g.cl
  }
  
  
  # Calculate the ABC distance
  res = err_fct(cl.i, ha.i, ww.i, 
                obs.cl, obs.ha, obs.ww, 
                use.cl, use.ha, use.ww,
                err.type) 
  return(list(
    distance    = as.numeric(res['total']),
    distance.cl = as.numeric(res['cl']),
    distance.ha = as.numeric(res['ha']),
    distance.ww = as.numeric(res['ww']),
    sim = a.sim
  ))
}


# Util function for parallel computation
calc_dist_parallel <- function(i,
                               priors,
                               obj,
                               use.cl, 
                               use.ha, 
                               use.ww, 
                               err.type,
                               deterministic,
                               n.sim = 10, 
                               verbose = FALSE ) {
  
  if(i == 1)      cat('ABC iteration # 1 /',nrow(priors),'\n') 
  if(i%%500 == 0) cat('ABC iteration #',i,'/',nrow(priors),'\n') 
  
  if(0){ # DEBUG 
    cat('\n--- DEBUG ABC iter: ',i,'\n')
  }
  
  # Update prior value at this iteration
  pp = priors[i,]
  
  if('i0prop' %in% names(priors)){
    obj$prms$i0prop <- pp$i0prop
    obj$prms = set_I_init(obj$prms)
  }
  
  # overwrite model parameters 
  # with priors values:
  obj$prms[names(pp)] <- pp
  
  x = reem_traj_dist_obs(
    obj           = obj,
    use.cl        = use.cl, 
    use.ha        = use.ha, 
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
  use.ha = prm.abc$use.ha
  use.cl = prm.abc$use.cl
  
  if(is.null(use.ww)) use.ww = 0
  if(is.null(use.ha)) use.ha = 0
  if(is.null(use.cl)) use.cl = 0
  
  d.max = max(obj$obs.cl$date, obj$obs.ha$date, obj$obs.ww$date) + 1
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
          '  hosp. adm. = ', ifelse(use.ha,'yes', 'NO'),'\n',
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
  snowfall::sfLibrary(dplyr)
  snowfall::sfLibrary(purrr)
  snowfall::sfExportAll()
  
  z = snowfall::sfLapply(
    x         = 1:n.abc, 
    fun       = calc_dist_parallel,
    priors    = priors, 
    obj       = obj,
    use.cl    = use.cl, 
    use.ha    = use.ha, 
    use.ww    = use.ww,
    n.sim     = n.sim,
    verbose   = verbose,
    err.type  = err.type,
    deterministic = deterministic)
  snowfall::sfStop()
  
  abc.err = sapply(z, '[[', 'distance')
  abc.err.cl = sapply(z, '[[', 'distance.cl')
  abc.err.ha = sapply(z, '[[', 'distance.ha')
  abc.err.ww = sapply(z, '[[', 'distance.ww')

  abc.sim = lapply(z, '[[', 'sim') %>% 
    # Add the simulation number (`index`) 
    # to each iteration:
    purrr::imap(~ dplyr::mutate(.x, index = .y))
  
  # -- Posteriors
  
  df.abc = cbind(priors, abc.err, abc.err.cl, abc.err.ha, abc.err.ww) %>% 
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

plot_traj <- function(obs, ps, varname, color, title, ylab) {
  colordark = paste0(color,'3')
  varnamem  = paste0(varname,'.m')
  varnamelo = paste0(varname,'.lo')
  varnamehi = paste0(varname,'.hi')
  
  g = ps %>% 
    ggplot2::ggplot(ggplot2::aes(x=date)) + 
    ggplot2::geom_line(ggplot2::aes(y = .data[[varnamem]]),
                       color = colordark)+
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[varnamelo]], 
                                      ymax = .data[[varnamehi]]), 
                         alpha=0.2, fill=color)+
    ggplot2::geom_point(data = obs, ggplot2::aes(y=obs)) +
    ggplot2::labs(title = title, y = ylab)
  return(g)
}


#' Extract posterior variable and aggregate
#'
#' @param obj `reem` object
#' @param type String. Typcal values: `cl`, `ha`. (wastewater is NOT aggregated)
#' @param rename Logical. Rename aggregated variables.
#'
#' @return
#' @keywords internal
#'
extract_fit_aggreg <- function(obj, type, rename = TRUE) {
 
  # type = 'cl'
  
  # Extract posterior simulations 
  ps = obj$fit.obj$post.simulations
  
  vtype = paste0('obs.',type)
  res = lapply(ps, helper_aggreg, 
               type = type, 
               dateobs = obj[[vtype]][['date']], 
               prms= obj$prms) |> 
    dplyr::bind_rows() |> 
    dplyr::group_by(date) |>
    dplyr::summarise(m = mean(obs),
                     lo = min(obs),
                     hi = max(obs),
                     n = n()) |>
    dplyr::filter(date <= max(obj[[vtype]]$date))
  
  if(rename & type == 'cl') res = dplyr::rename(res, 
                                                Y.m=m, Y.lo=lo, Y.hi=hi)
  if(rename & type == 'ha') res = dplyr::rename(res, 
                                                H.m=m, H.lo=lo, H.hi=hi)
  return(res)  
}

#' Helper function to summarize 
#' NON-AGGREGATED posterior trajectories
#'
#' @keywords internal
#' 
summarise_post_traj <- function(obj, obsdata, varname) {
  # obsdata = 'obs.cl'
  # varname = 'Y'
  
  # Retrieve posterior simulations:
  ps = obj[['fit.obj']][['post.simulations']] 
  
  res = ps %>% 
    dplyr::bind_rows() %>% 
    tidyr::drop_na(!!varname) %>%
    # Consider only the dates where the 
    # actual observations were made for 
    # a meaningful comparison:
    dplyr::filter(date %in% obj[[obsdata]][['date']]) |>
    # Summary stats:
    dplyr::group_by(date) %>%
    dplyr::summarise(zz.m = mean(.data[[varname]]),
                     zz.lo = min(.data[[varname]]),
                     zz.hi = max(.data[[varname]])) 
  
  # Rename summary stats using the variable name
  names(res)[names(res)=='zz.m']  <- paste0(varname, '.m')
  names(res)[names(res)=='zz.lo'] <- paste0(varname, '.lo')
  names(res)[names(res)=='zz.hi'] <- paste0(varname, '.hi')

  return(res)
}


#' Plot fit outputs of a (fitted) reem object
#' @param obj fitted reem object
#' 
#' 
reem_plot_fit <- function(obj) {
  
  # Prepare dataframes for plotting
  fit.obj = obj$fit.obj 
  obs.cl  = obj$obs.cl
  obs.ha  = obj$obs.ha
  obs.ww  = obj$obs.ww
  
  # inventory of observations
  has.cl = nrow(obj$obs.cl)>0
  has.ha = nrow(obj$obs.ha)>0
  has.ww = nrow(obj$obs.ww)>0
  
  if(has.cl){
    ps.cl = summarise_post_traj(obj, obsdata = 'obs.cl', varname = 'Y')
  }
  
  if(has.ww){
    ps.ww = summarise_post_traj(obj, obsdata = 'obs.ww', varname = 'Wr')
  } 
  
  if(has.ha) {
    ps.ha = extract_fit_aggreg(obj, type = 'ha')
  }
  
  ggplot2::theme_set(ggplot2::theme_bw())
  
  # ---- Trajectories
  g.cl = g.ha = g.ww = NULL
  
  if(has.cl){
    g.cl = plot_traj(obs     = obs.cl, 
                     ps      = ps.cl, 
                     varname = 'Y',
                     color   = 'chartreuse',
                     title   = 'Fit to clinical data', 
                     ylab    = 'cases')
  }
  
  if(has.ha){
    g.ha = plot_traj(obs     = obs.ha, 
                     ps      = ps.ha, 
                     varname = 'H',
                     color   = 'orchid',
                     title   = 'Fit to hosp.adm.', 
                     ylab    = 'cases')
  }
  
  if(has.ww){
    g.ww = plot_traj(obs     = obs.ww, 
                     ps      = ps.ww, 
                     varname = 'Wr',
                     color   = 'chocolate',
                     title   = 'Fit to wastewater data', 
                     ylab    = 'concentration')
  }
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
  
  
  # -- Ordered ABC distances
  
  fit.prm = obj$fit.prm

  n.post = round(fit.prm$n.abc * fit.prm$p.abc, 0)
  d = fit.obj$all.distances %>%
    dplyr::mutate(i = dplyr::row_number()) %>%
    dplyr::mutate(type = ifelse(i <= n.post,
                                'accepted','rejected'))
  # d.source = dplyr::select(d, dplyr::starts_with('abc.err.'))|>
  #   dplyr::mutate(x = dplyr::row_number()) |> 
  #   tidyr::pivot_longer(cols = dplyr::starts_with('abc') )
  
  alpha.source = 0.7
  col.source =
  
  g.dist = d %>%
    ggplot2::ggplot() + 
    #
    # -- Total ABC errors
    ggplot2::geom_vline(xintercept = n.post, 
                        linetype = 'dashed')+
    ggplot2::geom_step(ggplot2::aes(x = 1:nrow(d),
                                    y = abc.err,
                                    color = type), 
                       linewidth = 1 )+
    #
    # -- Cosmetics
    ggplot2::scale_y_log10() + 
    ggplot2::scale_x_log10() + 
    ggplot2::scale_color_manual(values = c('indianred2', 'gray90'))+
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) + 
    ggplot2::labs(
      title = 'ABC distances from data',
      x = 'ordered ABC iteration',
      y = 'distance'
    )
  
  g.dist2 = d |> 
    dplyr::select(starts_with('abc.err.')) |>
    dplyr::mutate(x = row_number())|>
    tidyr::pivot_longer(cols = dplyr::starts_with('abc.err.')) |> 
    dplyr::filter(x < n.post) |> 
    dplyr::mutate(source = substr(name, 9,10)) |> 
    ggplot2::ggplot() + 
    ggplot2::geom_step(ggplot2::aes(y=value, x=x, color = source), 
              alpha = alpha.source, linewidth = 1) + 
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) + 
    ggplot2::scale_color_manual(values =  c( 'steelblue', 'seagreen3', 'chocolate' )) + 
    ggplot2::labs(title = 'ABC distance by data source',
         subtitle = 'posteriors only',
         x='ordered ABC iteration', y = 'distance')
  # g.dist2
  
  g.list = list(
    traj.cl      = g.cl,
    traj.ha      = g.ha,
    traj.ww      = g.ww,
    post.prms    = g.post,
    post.prms.2d = gpall2d,
    dist         = g.dist,
    dist.source  = g.dist2
  )
  # Remove any NULL element (typically when an observation data set is missing)
  g.list = g.list[!sapply(g.list, is.null)]
  
  g.all = patchwork::wrap_plots(g.list, heights = c(1,2,1)) 
  res = c(g.list, list(all = g.all))
  return(res) 
}

