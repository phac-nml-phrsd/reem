quantile_df <- function(x, probs) {
  tibble::tibble(
    q     = stats::quantile(x, probs, na.rm = TRUE),
    qprob = probs
  )
}


#' @title Calculate quantiles of forecast trajectories
#'
#' @param simfwd Object of forecast trajectories.
#' @param prm.fcst List of forecast parameters. 
#' In particular, must contain an element named `ci`
#' specifying the quantiles probabilities to calculate.
#' @param vars String vector. If not NULL (default), overrides `prm.fcst$vars.to.fcst`.
#' 
#' @return Dataframe of quantiles.
#' @keywords internal
#'
summarize_fcst <- function(simfwd, prm.fcst, vars = NULL) {
  
  if(0){  # --- DEBUG
    simfwd = fcst$simfwd
    ci = c(0.50, 0.80, 0.95)
  }
  
  if(is.null(vars)) vars = prm.fcst$vars.to.fcst
  if(is.null(vars)) vars = c('Y', 'Wr')
  
  message('\nSummarizing forecasts for variables ',
          paste(vars, collapse = ', '), appendLF = FALSE)
  ci    = prm.fcst$ci
  probs = sort(c(0.5 - ci/2, 0.5 + ci/2) )
  
  res = dplyr::bind_rows(simfwd) %>% 
    dplyr::select(date, !!vars) %>% 
    tidyr::pivot_longer(cols = !!vars) %>%
    dplyr::reframe(
      quantile_df(value, probs), 
      mean = mean(value),
      .by = c(name, date))
  
  message(' done.')
  return(res)
}


#' Helper function to aggregate forecasts 
#'
#' @param var.to.aggregate String. Name of the variable to aggregate.
#' @param obj REEM object.
#' @param simfwd Object of forward simulations, 
#' as returned by the function \code{forecast()}.
#'
#' @return A list of dataframes.
#'
aggregate_fcst <- function(var.to.aggregate, obj, simfwd) {
  
  # retrieve the aggregation interval 
  # from the observation data set
  dt = as.numeric(stats::median(diff(obj$obs.cl$date)))
  
  res = list()
  for(i in 1:length(simfwd)){
    # Extract the daily forecasts for the variable
    sf = dplyr::select(simfwd[[i]], date, !!var.to.aggregate)
    
    # define the aggregation schedule
    dateaggr = seq.Date(from = max(obj$obs.cl$date) , 
                        to   = max(sf$date), 
                        by   = dt)
    
    # Calculate the aggregated values
    res[[i]] = sf %>% 
      aggcl(dt.aggr = dateaggr, 
            vars = var.to.aggregate) %>% 
      dplyr::filter(date > obj$fcst.prm$asof) %>%
      dplyr::mutate(index = dplyr::row_number())
    
    # Set proper name
    idx = names(res[[i]])==var.to.aggregate
    names(res[[i]])[idx] <- paste0(var.to.aggregate,'.aggr')
  }
  return(res)
}


summ_aggr_fcst <- function(simfwd, obj, var, prm.fcst) {
  # var = 'ha'
  
  # The aggregation time interval 
  # is defined based on time interval 
  # of the observations (median interval):
  vtype = paste0('obs.',var)
  d     = obj[[vtype]][['date']]
  
  if(!is.null(d)){
    
    ld = length(unique(d))
    
    if(ld > 1)
      dt = median(diff(d)) |> as.integer()
    
    # If there are no observations available,
    # we can still aggregate, but on using an 
    # arbitrary interval (e.g., weekly):   
    if(ld <= 1){
      dt = 7
      message('WARNING: Only one single observation for ', vtype,
              ' => assumming aggregation period to be *weekly*.')
    }
  }
  
  if(is.null(d)){
    dt = 7
    message('WARNING: No observations for ', vtype,
            ' => assumming aggregation period to be *weekly*.')
  }
  
  ci    =   prm.fcst$ci
  probs = sort(c(0.5 - ci/2, 0.5 + ci/2) )
  
  # Create the forward dates schedule (after "asof")
  d.fwd = max(d) + seq(dt, obj$prms$horizon, by = dt)

  # DEBUG
  #message("DEBUG: aggregation schedule: ", paste(d.fwd, collapse = ' ; ') )
  
  # We must consider the values after the last observation date
  # otherwise, will sum from date.start!
  simfwd.crop = purrr::map(simfwd, ~ dplyr::filter(., date > max(d) ))

  # Aggregation
  simfwd.aggr = lapply(simfwd.crop, 
                       helper_aggreg, 
                       type     = var, 
                       dateobs  = d.fwd, 
                       prms     = obj$prms, 
                       var.name = 'value') 
  
  # Summary statistics (quantiles, mean)
  summary.fcst.aggr = simfwd.aggr |>
    dplyr::bind_rows() |>
    dplyr::select(date, value) |> 
    dplyr::reframe(
      quantile_df(value, probs), 
      mean = mean(value),
      .by = c(date))
  
  res = list(
    simfwd.aggr       = simfwd.aggr,
    summary.fcst.aggr = summary.fcst.aggr
  )
  return(res)
}


# Helper function 
update_and_simulate <- function(i, pp, obj, verbose, tpb) {
  
  if(verbose)  cat('Simulating forward with posterior sample #',i,'\n')
  if(!is.null(tpb)) setTxtProgressBar(tpb, value = i)
  
  # update fitted parameters 
  # with their posterior values
  obj$prms[names(pp)] <- pp[i,]
  
  # if Bt was fitted (if not, `obj` is returned unchanged)
  obj = update_Bt(obj, pp[i,])
  
  # Update initial number of infectious individuals `obj$prms$I.init`
  obj$prms = set_I_init(obj$prms)
  
  # Forward simulations are calculated 
  # for every day in the future (no unobserved date):
  obj$prms$t.obs.ww <- 1:obj$fcst.prm$horizon.fcst 
  
  # Simulate forward
  s = obj$simulate_epi(deterministic = TRUE) #TODO: let user choose
  s$sim$index <- i
  return(s$sim)
}


#' @title Forecast a fitted epidemic 
#'
#' @param obj REEM object.
#' @param prm.fcst List. Forecast parameters.
#' @param verbose Logical. If \code{TRUE}, prints details.
#'
#' @return List containing forecast objects, 
#' with the following elements:
#' \item{asof: }{date when the forecast is done.}
#' \item{simfwd: }{object containing the forecasted trajectories} 
#' 
#' @export
#'
reem_forecast <- function(obj, prm.fcst, verbose, progressbar ) {
  
  if(0){   #---  DEBUG
    prm.fcst = list(
      asof = lubridate::ymd('2022-03-01'),
      horizon.fcst = lubridate::ymd('2022-07-01'),
      use.fit.post = TRUE,
      vars.to.fcst = c('Y', 'Wr', 'H'),
      n.resample = 20,
      ci = c(0.50, 0.80, 0.95)
    )
  }
  
  a = obj$fit.obj
  
  # Extend the horizon to match the one requested for the forecast
  obj$prms$horizon <- as.integer(prm.fcst$horizon.fcst - obj$prms$date.start)
  
  # In this case, the forecast is simply reusing
  # the simulations from the posterior fits
  # (i.e., no resampling from posterior distributions)
  if(prm.fcst$use.fit.post){
    
    # Posterior parameters
    pp = a$post.prms %>% dplyr::select(!tidyr::starts_with('abc'))
  
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
    
    message('\nSampling ',ns, ' posterior parameter sets out of ', npp,
            ' available.\n')
    tpb = NULL
    if(progressbar){
      tpb = utils::txtProgressBar(min = 0, max = max(ii), style = 3, width = 25)
    }
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
  
  if(is.null(prm.fcst$vars.to.fcst)) 
    prm.fcst$vars.to.fcst = c('Y', 'Wr')
  
  summary.fcst = summarize_fcst(
    simfwd   = simfwd, 
    prm.fcst = prm.fcst)
  
  # Dealing with aggregated hospital admissions
  simfwd.aggr = list()
  summary.fcst.aggr = list()
  
  has.cl = nrow(obj$obs.cl)>0
  has.ha = nrow(obj$obs.ha)>0
  
  if(has.cl){
    tmp.cl = summ_aggr_fcst(simfwd, obj, var = 'cl', prm.fcst)
    simfwd.aggr[['Y.aggr']]       = tmp.cl$simfwd.aggr
    summary.fcst.aggr[['Y.aggr']] = tmp.cl$summary.fcst.aggr
  }
  if(has.ha){
    tmp.ha = summ_aggr_fcst(simfwd, obj, var = 'ha', prm.fcst)
    simfwd.aggr[['H.aggr']]       = tmp.ha$simfwd.aggr
    summary.fcst.aggr[['H.aggr']] = tmp.ha$summary.fcst.aggr
  }
   
  return( list(
    asof              = prm.fcst$asof,
    simfwd            = simfwd, 
    summary.fcst      = summary.fcst,
    simfwd.aggr       = simfwd.aggr,
    summary.fcst.aggr = summary.fcst.aggr
  ))
}


#' @title Helper function to add ribbons for forecasts.
#'
#' @param g 
#' @param z 
#' @param k 
#' @param col.fcst 
#' @param alpha.ribbon 
#'
#' @return A ggplot object
#' @keywords internal
#'
add_ribbons_quantiles <- function(g, qlist, k,
                                  col.fcst, alpha.ribbon) {
  nq = length(qlist)
  res = g + 
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[[qlist[k]]], 
                   ymax = .data[[qlist[nq - k + 1] ]]), 
      fill = col.fcst, alpha = alpha.ribbon)
  return(res)
}



#' Helper function
#'
#' @param traj.fit 
#' @param traj.fcst 
#' @param obs 
#' @param alpha.ribbon 
#' @param col.fit 
#' @param col.fcst 
#' @param fcst.prm 
#' @param title 
#' @param ylab 
#' @param qlist 
#'
#' @return a ggplot object
#' 
#' @keywords internal
#' 
#' @import ggplot2
#'
plot_fitfcst <- function(traj.fit, traj.fcst, obs, 
                         alpha.ribbon, col.fit, col.fcst, 
                         fcst.prm, 
                         title, ylab, qlist, xaxis) {
  # Last observation
  date.obs.last = max(obs$date)
  obs.last = obs$obs[obs$date == date.obs.last]
  obs.last.plot = lubridate::ymd(date.obs.last)
  
  # Plot
  g = ggplot(data = traj.fcst, aes(x=date))+ 
    # --- Fit
    geom_line(data = traj.fit, aes(y=m), 
              color = col.fit, linetype = 'dashed') + 
    geom_ribbon(data = traj.fit, aes(ymin=lo, ymax=hi), 
                fill = col.fit, color = col.fit, 
                linewidth = 0.1, 
                alpha = alpha.ribbon / 2) + 
    # --- Observations
    geom_point(data = obs, aes(y=obs)) + 
    # --- Forecast
    geom_line( aes(y = mean), color= col.fcst, 
               linetype = 'dotted') + 
    geom_vline(xintercept = fcst.prm$asof, 
               linetype = 'dashed', 
               color = col.fcst) + 
    annotate(geom = 'text', y=1, x = fcst.prm$asof + 1, 
             label = fcst.prm$asof, size = 2, hjust = 0,
             color = col.fcst) + 
    annotate(geom = 'segment', 
             x = date.obs.last, xend = date.obs.last,
             y = 0, yend = obs.last, linetype = 'dotted', color = 'grey75')+
    annotate(geom = 'text', x = date.obs.last, y = 1,
               label = obs.last.plot, hjust = 1, size = 2)+
    xaxis + 
    theme(panel.grid.minor = element_line(color = 'grey97'))+
    labs(title =title, 
         x = '', y = ylab, 
         caption = paste('fit ribbon: min/max\nfcst ribbon: quantiles'))
  
  # Add all the ribbons corresponding to each quantile
  for(k in 1:(length(qlist)/2)) {
    g = add_ribbons_quantiles(g, qlist, k, col.fcst, alpha.ribbon/2)
  }
  # g
  return(g)
}


#' Helper function
#'
#' @param post.sim 
#' @param var 
#'
#' @return dataframe
#' @keywords internal
#'
#' 
sumpost <- function(post.sim, var) {
  res = post.sim %>%  
    dplyr::bind_rows() %>% 
    dplyr::group_by(date) %>% 
    dplyr::summarize(m  = mean(.data[[var]]), 
                     lo = min(.data[[var]]), 
                     hi = max(.data[[var]])) %>% 
    tidyr::drop_na(m)
  return(res)
}

#' Plot forecasts
#'
#' @param obj List. Forecast object as returned by the function \code{forecast()}.
#' @param date_breaks String. Date breaks (as understood by \code{ggplot}).
#' @param date_labels String. Date labels (as understood by \code{ggplot}).
#' @param logscale Logical. Use logarithmic scale for y axis. Default is \code{FALSE}.
#'
#' @return A ggplot object.
#'
reem_plot_forecast <- function(
    obj,
    date_breaks,
    date_labels,
    logscale) {
  
  obs.cl = obj$obs.cl
  obs.ha = obj$obs.ha
  obs.ww = obj$obs.ww
  
  n.cl = nrow(obs.cl)
  n.ha = nrow(obs.ha)
  n.ww = nrow(obs.ww)
  
  fcst.obj = obj$fcst.obj
  fcst.prm = obj$fcst.prm
  
  # Cosmetics  
  
  xaxis = ggplot2::scale_x_date(
    date_breaks = date_breaks, 
    date_labels = date_labels)
  alpha.ribbon = 0.20
  
  col.fcst = "steelblue" ; col.fit = "tan3"
  
  # - - - - Retrieve fitted curves
  
  post.sim = obj$fit.obj$post.simulations
  
  if(n.cl > 0) fitsim.cl = sumpost(post.sim, var = 'Y')
  if(n.ha > 0) fitsim.ha = extract_fit_aggreg(obj, 'ha', rename = F) 
  if(n.ww > 0) fitsim.ww = sumpost(post.sim, var = 'Wr')
  
  # Retrieve the forecast summary
  sf = fcst.obj$summary.fcst 
  
  # set the aggregation dates as 
  # starting from `asof` and a time interval
  # equal to the one of the observations:
  dt = as.numeric(diff(obs.cl$date))[1]
  
  # Reformat to suit ggplot
  sf2 = sf %>% 
    tidyr::pivot_wider(names_from = qprob, values_from = q, 
                names_prefix = 'q_')
  
  z  = names(sf2)
  qlist  = z[grepl('^q_',z)]
  
  # ** WARNING **
  # HERE WE CALCULATE THE SUM OF THE QUANTILE WHICH IS 
  # DIFFERENT FROM THE QUANTILE OF THE SUM (WHAT WE REALLY WANT!)
  # TODO: CHANGE THAT!
  
  if(n.ww > 0) sf.ww = dplyr::filter(sf2, name == 'Wr') %>%
    tidyr::drop_na(mean) %>%
    dplyr::filter(date >= fcst.prm$asof)
  
  if(n.cl > 0) sf.cl = dplyr::filter(sf2, name == 'Y') %>%
    tidyr::drop_na(mean) %>%
    dplyr::filter(date >= fcst.prm$asof)
  
  if(n.ha > 0) sf.ha = fcst.obj$summary.fcst.aggr$H.aggr |> 
    tidyr::pivot_wider(names_from = qprob, 
                       values_from = q, 
                       names_prefix = 'q_')
  
  # - - - Plots - - - 
  
  g.cl = g.ha = g.ww = ggplot2::ggplot()
  
  if(n.cl > 0) g.cl = plot_fitfcst(
    traj.fit = fitsim.cl, 
    traj.fcst = sf.cl, 
    obs = obs.cl, 
    alpha.ribbon = alpha.ribbon,
    col.fit = col.fit, 
    col.fcst = col.fcst, 
    fcst.prm = fcst.prm, 
    title = 'Reported cases', ylab = 'cases',
    qlist = qlist, 
    xaxis = xaxis)
 
  if(n.ha > 0) g.ha = plot_fitfcst(
    traj.fit = fitsim.ha, 
    traj.fcst = sf.ha, 
    obs = obs.ha, 
    alpha.ribbon = alpha.ribbon,
    col.fit = col.fit, 
    col.fcst = col.fcst, 
    fcst.prm = fcst.prm, 
    title = 'Hospital admissions', ylab = 'daily adm',
    qlist = qlist,
    xaxis = xaxis)
  
  if(n.ww > 0) g.ww = plot_fitfcst(
    traj.fit = fitsim.ww, 
    traj.fcst = sf.ww, 
    obs = obs.ww, 
    alpha.ribbon = alpha.ribbon,
    col.fit = col.fit, 
    col.fcst = col.fcst, 
    fcst.prm = fcst.prm, 
    title = 'Wastewater', ylab = 'concentration',
    qlist = qlist,
    xaxis = xaxis)
  
  if(logscale){
    g.cl = g.cl + scale_y_log10()
    g.ha = g.ha + scale_y_log10()
    g.ww = g.ww + scale_y_log10()
  }
  
  return(list(
    cl = g.cl, 
    ha = g.ha, 
    ww = g.ww
  ))
}


#' Plot timing and level of peak forecast 
#'
#' @param var String. Name of the variable
#' @param obj List. Object
#' @param date_labels String. Date labels (as understood by \code{ggplot}).
#' @param logscale Logical. Use logarithmic scale for y axis. Default is \code{FALSE}.
#'
#' @return A ggplot object plotting the time and value of the forecasted peaks.
#'
reem_plot_peak <- function( var,
                            obj         ,
                            date_labels ,
                            logscale    ) {
  
  pk = obj$forecast_peak(var = var)
  
  # Check if all dates or values are the same
  # (typically when the peak has been reached in the past)
  dv = diff(pk$peak.value)
  dd = diff(pk$peak.date)
  if(all(dv==0)) {
    perturb = 0.002 * rnorm(n = nrow(pk))
    pk$peak.value = pk$peak.value * (1 + perturb)
  }
  if(all(dd==0)) {
    perturb =  rnorm(n = nrow(pk)) 
    pk$peak.date = pk$peak.date  + perturb/3
  }
    
  g.pk2d = pk %>% 
    ggplot2::ggplot(ggplot2::aes(x=peak.date , y=peak.value)) +
    ggplot2::geom_point(alpha = 0.5, size = 3) + 
    ggplot2::geom_density_2d_filled(color = 'grey50', alpha = 0.6)+
    ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank(), 
                   panel.border = ggplot2::element_blank(),
                   axis.ticks   = ggplot2::element_blank())  +
    ggplot2::scale_y_continuous(position = "right", 
                                labels = scales::comma_format()) +
    ggplot2::scale_x_date(date_labels = date_labels)+
    ggplot2::labs(title = paste0('Forecast peak for `', var,'`'),
                  x = 'Peak date', y = 'Peak value') + 
    ggplot2::guides(fill = 'none')
  # g.pk2d
  
  g.pk.v = pk |> 
    ggplot2::ggplot(ggplot2::aes(x=peak.value)) + 
    ggplot2::geom_histogram(bins = 12, fill = 'lightgrey', color='darkgrey') + 
    ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank())+
    ggplot2::coord_flip()+
    ggplot2::labs(x='',y='')
  # g.pk.v
  
  g.pk.d = pk |> 
    ggplot2::ggplot(ggplot2::aes(x=peak.date)) + 
    ggplot2::geom_histogram(bins = 12, fill = 'lightgrey', color='darkgrey')+
    ggplot2::scale_x_date(date_labels = date_labels) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = ggplot2::rel(0.6))  )+
    ggplot2::labs(x='',y='')
  # g.pk.d
    
  if(logscale) {
    g.pk2d = g.pk2d + ggplot2::scale_y_log10(position = "right", 
                                             labels = scales::comma_format())
    g.pk.v = g.pk.v + ggplot2::scale_x_log10(labels = scales::comma_format())
  }
  
  g.pk = patchwork::wrap_plots(g.pk2d, g.pk.v, g.pk.d, 
                               heights = c(5,1),
                               widths  = c(13,2),
                               ncol = 2)
  # g.pk
  return(g.pk)
}

#' Extract the forecasted peak values and timing.
#' 
#' @param var String. Name of the model variable. 
#' @param fcst Forecast object `fcst.obj` of the `reem` class object.  
#' @param obs Dataframe. Observations associated with `var`.  
#'
#' @return A dataframe of the forecasted peaks.
#'
reem_forecast_peak <- function(var, fcst, obs) {
  
  is.aggregated = grepl('aggr$', var)
  
  if(is.aggregated){
    df = fcst$simfwd.aggr[[var]]
  }
  if(!is.aggregated){
    df = fcst$simfwd |> 
      purrr::map(~dplyr::mutate(., value = .data[[var]]))
  }
  
  # Evaluate the FUTURE peaks
  # (trajectories AFTER the as of date)
  pkafter = df %>% 
    dplyr::bind_rows(.id = 'post') %>%
    dplyr::group_by(post) %>% 
    dplyr::summarise(
      peak.date  = date[which.max(value)[1]],
      peak.value = max(value,na.rm = TRUE))
  
  # Maximum values (potential peak) from observations
  # If peak reached on multiple dates, return the earliest.
  if(!is.null(obs)){
    pk.obs.value = max(obs$obs)[1]
    pk.obs.date  = obs$date[obs$obs == pk.obs.value][1]
  }
  if(is.null(obs)){
    pk.obs.value = -1
  }
  
  # Compare future peaks with past observations
  # and overwrite future if past peak is larger
  if(!is.null(obs)){
    res = pkafter |> dplyr::mutate(
      peak.value2 = dplyr::if_else(peak.value < pk.obs.value, 
                                   pk.obs.value, peak.value),
      peak.date2  = dplyr::if_else(peak.value < pk.obs.value, 
                                   pk.obs.date, peak.date)) |> 
      dplyr::select(post, peak.date2, peak.value2) |>
      dplyr::rename(peak.date = peak.date2, peak.value = peak.value2)
  }
  if(is.null(obs)){
    res = pkafter
  }
  
  return(res)
}


#' Returns the probability that the forecasted variable trajectory
#' is within a box defined by lower and upper dates and values.
#'
#' @param var String. Name of the variable.
#' @param date.lower Date. Lower date defining the box.
#' @param date.upper Date. Upper date defining the box.
#' @param val.lower Numeric. Lower value defining the box.
#' @param val.upper Numeric. Upper value defining the box.
#' @param fcst Forecast object as returned by the function \code{forecast()}.
#'
#' @return Numeric. The corresponding probability.
#'
reem_proba_box <- function(var, 
                           date.lower, 
                           date.upper,
                           val.lower, 
                           val.upper, 
                           fcst) {
  # --- DEBUG
  # var = 'Y.aggr' ; val.lower = 10 ; val.upper = 70
  # date.lower = ymd('2022-03-01')
  # date.upper = ymd('2022-03-20')
  
  is.aggregated = grepl('\\.aggr$', var)
  
  if(!is.aggregated) {
    fs = fcst$simfwd |> 
      purrr::map(~dplyr::mutate(., value = .data[[var]]))
  }
  if(is.aggregated)  fs = fcst$simfwd.aggr[[var]]
  
  n = length(fs)
  x = logical(n)
  for(i in 1:n){
    val = fs[[i]] %>% 
      dplyr::filter(between(date, date.lower, date.upper)) %>% 
      dplyr::select(date, value) %>% 
      tidyr::drop_na(value) 
    x[i] = any(val.lower <= val$value & val$value <= val.upper)
  }
  return(mean(x))
}



#' @title Retrieve the forecast posterior distribution
#' 
#' @description After the REEM model is fitted to the data,
#' epidemic trajectories are simulated by drawing the model 
#' parameters from their posterior distributions. This generates
#' "posterior trajectories" that are used for forecasting. 
#' This function returns, at given future time "slices", 
#' using the posterior trajectories, the forecasted (posterior) distribution
#' of the variable `var`. 
#'
#' @param var String. Name of the variable forecasted. Typical choices are
#' `Y`, `Y.aggr` and `Wr`
#' @param date.future Vector of dates. Future dates where the variable `var` is forecasted.
#' @param fcst Forecast `reem` object.
#' @param density.n Integer. Number of equally spaced points at which 
#' the density is to be estimated. see \code{?density}.
#' @param density.adjust Numerical. Bandwidth adjustment. see \code{?density}.
#' @param aggr.window Integer or NULL. Time window to aggregate the values 
#' of the variable `var`. No aggregation if NULL. 
#'
#' @return List of dataframes defining the densities as `y=f(x)`.
#' @export
#'
#' @seealso \code{density()}
#'
#' 
reem_get_fcst_density <- function(
    var, 
    date.future, 
    fcst,
    density.n , 
    density.adjust,
    aggr.window) {
  
  res = list()
  
  # Loop through all dates
  df = data.frame(date = date.future)
  
  for(i in seq_along(date.future)){
    # print(i) # DEBUG
    fcst.vals = extract_helper(
      i = i, 
      fcst = fcst, 
      aggr.window = aggr.window, 
      obs.new = df, 
      var = var)
    
    tmp  = calc_density_one(
      fcst.vals = fcst.vals,
      density.n = density.n, 
      density.adjust = density.adjust)
    
    res[[i]] = tmp  %>% 
      dplyr::mutate(date = date.future[i])
  }
  return(res)
}



