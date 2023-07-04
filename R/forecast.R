

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
  
  summary.fcst = summarize_fcst(simfwd, prm.fcst)
  
  return( list(
    simfwd = simfwd, 
    summary.fcst = summary.fcst
  ))
  
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
  dt.aggr.fcst = seq.Date(from = fcst.prm$asof , 
                          to = max(sf$date), 
                          by = dt)
  
  # aggregation of clinical reports
  sf.cl = aggcl(df = sf, 
                dt.aggr = dt.aggr.fcst, 
                vars = c('Y_mean','Y_lo','Y_hi')) %>% 
    filter(date >fcst.prm$asof)
  
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


