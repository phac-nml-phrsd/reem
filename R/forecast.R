

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
      q     = quantile(x, probs, na.rm = TRUE),
      qprob = probs
    )
  }
  
  res = bind_rows(simfwd) %>% 
    select(date, !!vars) %>% 
    pivot_longer(cols = !!vars) %>%
    reframe(quantile_df(value, probs), 
            mean = mean(value),
            .by = c(name, date))
  
  message(' done.')
  return(res)
}


#' Helper function to aggregate forecasts 
#'
#' @param var.to.aggregate 
#' @param obj 
#' @param simfwd 
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
      filter(date > obj$fcst.prm$asof) 
    
    names(res[[i]])[names(res[[i]])==var.to.aggregate] <- paste0(var.to.aggregate,'.aggr')
    
  }
  return(res)
}

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
  
  summary.fcst = summarize_fcst(simfwd, 
                                prm.fcst,
                                vars = c('Y', 'Wr')) # TODO: remove hard code
  
  # Dealing with aggregated incidence
  
  simfwd.aggr = list()
  summary.fcst.aggr = list()
  
  # TODO: make a function for these 2 function calls!
  simfwd.aggr[['Y.aggr']] = aggregate_fcst(var.to.aggregate = 'Y', 
                                           obj = obj, 
                                           simfwd = simfwd)
  
  summary.fcst.aggr[['Y.aggr']] = summarize_fcst(simfwd.aggr, 
                                                 prm.fcst,
                                                 vars = c('Y.aggr'))  
  
  return( list(
    asof   = prm.fcst$asof,
    simfwd = simfwd, 
    summary.fcst = summary.fcst,
    simfwd.aggr = simfwd.aggr,
    summary.fcst.aggr = summary.fcst.aggr
  ))
  
}


#' Helper function
#'
#' @param g 
#' @param z 
#' @param k 
#' @param col.fcst 
#' @param alpha.ribbon 
#'
#' @return A ggplot object
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
  
  return(list(
    cl = g.cl, 
    ww = g.ww
  ))
  
}


#' Extract the forecasted peak values and timing.
#' 
#' @param var String. Name of the model variable. 
#' @param fcst Forecast object `fcst.obj` of the `reem` class object.  
#'
#' @return A dataframe of the forecasted peaks.
#'
reem_forecast_peak <- function(var, fcst) {
  res = fcst$simfwd %>% 
    bind_rows() %>%
    group_by(index) %>% 
    summarise(peak.date = date[which.max(.data[[var]])[1]],
              peak.value = max(.data[[var]],na.rm = TRUE))
  return(res)
}


#' Returns the probability that the forecasted variable
#' is within a box defined by lower and upper dates and values.
#'
#' @param var 
#' @param date.lower 
#' @param date.upper 
#' @param val.lower 
#' @param val.upper 
#' @param fcst 
#'
#' @return
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
  # aggr.interval = 7
 
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
    x[i] = any(val.lower <= val & val <= val.upper)
  }
  x
  return(mean(x))
  
  # --- DEBUG
  df =  fcst$summary.fcst
  if(is.aggregated) df = fcst$summary.fcst.aggr[[var]]
  g = df %>% 
    filter(name == !!var) %>% 
    ggplot(aes(x=date, y=mean)) + 
    geom_line(color = 'blue') + 
    geom_hline(yintercept = c(val.lower, val.upper), linetype='dashed') +
    geom_vline(xintercept = c(date.lower, date.upper), linetype='dashed')
  plot(g)
    
}



