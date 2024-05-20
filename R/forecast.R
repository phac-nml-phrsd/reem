

#' @title Calculate quantiles of forecast trajectories
#'
#' @param simfwd Object of forecast trajectories.
#' @param prm.fcst List of forecast parameters. 
#' In particular, must contain an element named `ci`
#' specifying the quantiles probabilities to calculate.
#' @param vars String vector. Names of the variables
#' forecasted to calculate the quantiles.
#'
#' @return Dataframe of quantiles.
#' @keywords internal
#'
summarize_fcst <- function(simfwd, prm.fcst, vars) {
  
  if(0){  # --- DEBUG
    simfwd = fcst$simfwd
    ci = c(0.50, 0.80, 0.95)
    vars = c('Y', 'Wr')
  }
  
  message('\nSummarizing forecasts...',appendLF = FALSE)
  ci    = prm.fcst$ci
  probs = sort(c(0.5 - ci/2, 0.5 + ci/2) )
  
  quantile_df <- function(x, probs) {
    tibble(
      q     = stats::quantile(x, probs, na.rm = TRUE),
      qprob = probs
    )
  }
  
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
  dt = as.numeric(median(diff(obj$obs.cl$date)))
  
  res = list()
  for(i in 1:length(simfwd)){
    # Extract the daily forecasts for the variable
    sf = select(simfwd[[i]], date, !!var.to.aggregate)
    
    # define the aggregation schedule
    dateaggr = seq.Date(from = max(obj$obs.cl$date) , 
                        to   = max(sf$date), 
                        by   = dt)
    
    # Calculate the aggregated values
    res[[i]] = sf %>% 
      aggcl(dt.aggr = dateaggr, 
            vars = var.to.aggregate) %>% 
      filter(date > obj$fcst.prm$asof) %>%
      mutate(index = row_number())
    
    # Set proper name
    idx = names(res[[i]])==var.to.aggregate
    names(res[[i]])[idx] <- paste0(var.to.aggregate,'.aggr')
  }
  return(res)
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
reem_forecast <- function(obj, prm.fcst, verbose ) {
  
  if(0){   #---  DEBUG
    prm.fcst = list(
      asof = ymd('2022-03-01'),
      horizon.fcst = ymd('2022-07-01'),
      use.fit.post = TRUE,
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
    
    pp = a$post.prms %>% dplyr::select(!tidyr::starts_with('abc'))
    
    # Helper function 
    update_and_simulate <- function(i, pp, obj, verbose, tpb) {
      
      if(verbose)  cat('Simulating forward with posterior sample #',i,'\n')
      if(!verbose) setTxtProgressBar(tpb, value = i)
      
      # update fitted parameters 
      # with their posterior values
      obj$prms[names(pp)] <- pp[i,]
      
      # Update initial number of infectious individuals `obj$prms$I.init`
      obj$prms = set_I_init(obj$prms)
      
      # Forward simulations are calculated 
      # for every day in the future (no unobserved date):
      obj$prms$t.obs.ww <- 1:prm.fcst$horizon.fcst 
      
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
    
    message('\nSampling ',ns, ' posterior parameter sets out of ', npp,
            ' available.\n')
    
    tpb = txtProgressBar(min = 0, max = max(ii), style = 3, width = 25)
    
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
  
  summary.fcst = summarize_fcst(
    simfwd   = simfwd, 
    prm.fcst = prm.fcst,
    vars     = c('Y', 'Wr')) # TODO: remove hard code
  
  # Dealing with aggregated incidence
  simfwd.aggr = list()
  summary.fcst.aggr = list()
  
  # TODO: make a function for these 2 function calls!
  simfwd.aggr[['Y.aggr']] = aggregate_fcst(
    var.to.aggregate = 'Y', 
    obj              = obj, 
    simfwd           = simfwd)
  
  summary.fcst.aggr[['Y.aggr']] = summarize_fcst(
    simfwd   = simfwd.aggr, 
    prm.fcst = prm.fcst,
    vars     = c('Y.aggr'))  
  
  return( list(
    asof   = prm.fcst$asof,
    simfwd = simfwd, 
    summary.fcst = summary.fcst,
    simfwd.aggr = simfwd.aggr,
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
add_ribbons_quantiles <- function(g, z, k,
                                  col.fcst, alpha.ribbon) {
  nz = length(z)
  res = g + 
    geom_ribbon(aes(ymin = .data[[z[k]]], 
                    ymax = .data[[z[nz-k+1] ]]), 
                fill = col.fcst, alpha = alpha.ribbon)
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
  obs.ww = obj$obs.ww
  fcst.obj = obj$fcst.obj
  fcst.prm = obj$fcst.prm
  
  # - - - Cosmetics  
  
  col.fcst = 'steelblue2'
  col.fit  = 'tan2'
  xaxis = scale_x_date(
    date_breaks = date_breaks, 
    date_labels = date_labels)
  alpha.ribbon = 0.10
  
  # - - - - Retrieve fitted curves
  
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
  
  # Tue Jul 11 13:49:02 2023 ------------------------------
  # FIXME
  # This should be redundant now that the 
  # forecast object returns aggregated forecasts
  # see function `aggregate_fcst()`
  
  # --- Aggregation of clinical reports ---
  
  # Aggregate clinical reports for the fitted part
  fitsim.cl = tmp.cl %>% 
    aggcl(dt.aggr = obs.cl$date, 
          vars = c('m','lo','hi')) %>%
    filter(date <= max(obs.cl$date))
  
  # Retrieve the forecast summary
  sf = fcst.obj$summary.fcst 
  
  # set the aggregation dates as 
  # starting from `asof` and a time interval
  # equal to the one of the observations:
  dt = as.numeric(diff(obs.cl$date))[1]
  dt.aggr.fcst = seq.Date(from = fcst.prm$asof , 
                          to = max(sf$date), 
                          by = dt)
  
  # Reformat to suit ggplot
  sf2 = sf %>% 
    pivot_wider(names_from = qprob, values_from = q, 
                names_prefix = 'q_')
  
  z  = names(sf2)
  z  = z[grepl('^q_',z)]
  nz = length(z)
  
  # ** WARNING **
  # HERE WE CALCULATE THE SUM OF THE QUANTILE WHICH IS 
  # DIFFERENT FROM THE QUANTILE OF THE SUM (WHAT WE REALLY WANT!)
  # TODO: CHANGE THAT!
  
  sf.cl = sf2 %>% 
    filter(name == 'Y') %>%
    aggcl(dt.aggr = dt.aggr.fcst, 
          vars = c('mean', z)) %>% 
    filter(date > fcst.prm$asof)
  
  # - - - Plots - - - 
  
  g.cl = ggplot(data = sf.cl,
                aes(x=date))+ 
    #
    # --- Fit
    #
    geom_line(data = fitsim.cl, aes(y=m), 
              color = col.fit, linetype = 'dashed') + 
    geom_ribbon(data = fitsim.cl, aes(ymin=lo, ymax=hi), 
                fill=col.fit, alpha = alpha.ribbon / 2) + 
    geom_point(data = obs.cl, aes(y=obs)) + 
    #
    # --- Forecast
    # 
    geom_line( aes(y = mean), color= col.fcst, 
               linetype = 'dotted') + 
    geom_vline(xintercept = fcst.prm$asof, 
               linetype = 'dashed', 
               color = 'gray50') + 
    annotate(geom = 'text', y=1, x=fcst.prm$asof, 
             label = fcst.prm$asof, size = 2) + 
    xaxis + 
    labs(title = 'Forecast infections', 
         x = '', y = 'cases')
  
  # Add all the ribbons corresponding to each quantile
  for(k in 1:(nz/2)) {
    g.cl = add_ribbons_quantiles(g.cl, z, k, 
                                 col.fcst, alpha.ribbon)
  }
  # g.cl
  
  # --- Wastewater 
  
  sf.ww = filter(sf2, name == 'Wr')
  
  g.ww =  sf.ww %>%
    drop_na(mean) %>%
    filter(date >= fcst.prm$asof) %>%
    ggplot(aes(x=date)) + 
    #
    # --- Fit
    # 
    geom_line(data = fitsim.ww, aes(y=m), 
              color = col.fit, linetype = 'dashed') + 
    geom_ribbon(data = fitsim.ww, aes(ymin=lo, ymax=hi), 
                fill=col.fit, alpha = alpha.ribbon / 2) + 
    geom_point(data = obs.ww, aes(y=obs)) + 
    #
    # --- Forecast
    #
    geom_line( aes(y = mean), color= col.fcst, 
               linetype = 'dotted') + 
    geom_vline(xintercept = fcst.prm$asof, 
               linetype = 'dashed', 
               color = 'gray50') + 
    annotate(geom = 'text', y = 0, x = fcst.prm$asof, 
             label = fcst.prm$asof, size = 2) + 
    labs(title = 'Forecast wastewater concentration', 
         x = '', y = 'concentration (gcp/ml)') +
    xaxis 
  
  for(k in 1:(nz/2)) {
    g.ww = add_ribbons_quantiles(g.ww, z, k, 
                                 col.fcst, alpha.ribbon)
  }
  # g.ww
  
  if(logscale){
    g.cl = g.cl + scale_y_log10()
    g.ww = g.ww + scale_y_log10()
  }
  
  return(list(
    cl = g.cl, 
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
#' @return
#'
reem_plot_peak <- function( var,
    obj         ,
    date_labels ,
    logscale    ) {
  
  pk = obj$forecast_peak(var = var)

  g.pk = pk %>% ggplot(aes(x=peak.date , y=peak.value)) +
    geom_density_2d_filled(color = 'grey50', alpha = 0.6)+
    theme(panel.grid.minor.y = element_blank(), 
          panel.border = element_blank(),
          axis.ticks   = element_blank())  +
    scale_x_date(date_labels = date_labels)+
    labs(title = paste0('Forecast peak for `', var,'`'),
         x = 'Peak date', y = 'Peak value') + 
    guides(fill = 'none')
  
  if(logscale) g.pk = g.pk + scale_y_log10()
  return(g.pk)
}

#' Extract the forecasted peak values and timing.
#' 
#' @param var String. Name of the model variable. 
#' @param fcst Forecast object `fcst.obj` of the `reem` class object.  
#'
#' @return A dataframe of the forecasted peaks.
#'
reem_forecast_peak <- function(var, fcst) {
  
  is.aggregated = grepl('aggr$', var)
  
  if(is.aggregated){
    df = fcst$simfwd.aggr[[var]]
  }
  if(!is.aggregated){
    df = fcst$simfwd
  }
  
  res = df %>% 
    bind_rows(.id = 'post') %>%
    group_by(post) %>% 
    summarise(peak.date = date[which.max(.data[[var]])[1]],
              peak.value = max(.data[[var]],na.rm = TRUE))
  
  return(res)
  
  if(FALSE){ # DEBUG
    df %>% 
      bind_rows(.id = 'post') %>% 
      ggplot(aes(x=date, y=.data[[var]])) + 
      geom_line(aes(group = post), alpha = 0.4) + 
      geom_point(data = res, aes(x=peak.date, y=peak.value), size=3)
  }
}


#' Returns the probability that the forecasted variable
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
   
  fs = fcst$simfwd
  if(is.aggregated) fs = fcst$simfwd.aggr[[var]]
  n = length(fs)
  x = logical(n)
  for(i in 1:n){
    val = fs[[i]] %>% 
      filter(between(date, date.lower, date.upper)) %>% 
      select(date, !!var) %>% 
      drop_na(!!var) 
    x[i] = any(val.lower <= val[[var]] & val[[var]] <= val.upper)
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
      mutate(date = date.future[i])
  }
  return(res)
}



