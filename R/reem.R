
#' Class to represent a Renewal Equation Epidemic Model.
#'
#' @field name Character for the name of the model. 
#' @field prms List of model parameters.
#' @field is.fitted Logical. Is the model fitted to the observation
#' data attached (\code{obs.cl} and/or \code{obs.ww})? 
#' @field obs.cl Dataframe representing the clinical observations. 
#' Must have \code{date} and \code{obs} variables. 
#' @field obs.ww Dataframe representing the wastewater observations. 
#' Must have \code{date} and \code{obs} variables. 
#' @field fit.obj List containing the fitted object. 
#' @field fit.prm List that defines the parameters for the fitting algorithm.
#' 
#' @return
#' 
#' @export
#'
#' 
setRefClass(
  Class = 'reem', 
  fields =  list(
    name      = "character",
    prms      = "list",
    obs.cl    = "data.frame",
    obs.ww    = "data.frame",
    is.fitted = "logical",
    fit.obj   = "list",
    fit.prm   = "list",
    fcst.prm  = "list",
    fcst.obj  = "list"
  ),
  
  # - - - - - - - - - - - - - - - - - -
  # - - - `reem` CLASS METHODS  - - - -
  # - - - - - - - - - - - - - - - - - -
  
  methods = 
    list(
      
      print_prms = function(){
        cat('\n--- Parameters for REEM `',name,'`\n')
        for(i in seq_along(prms)){
          cat(names(prms)[i],' = ', prms[[i]],'\n')
        }
        cat(' --------------------------\n')
      },
      
      simulate = function(deterministic){
        res = reem_simulate(prms, deterministic)  
        return(res)  
      },
      
      simulate_epi = function(deterministic){
        return(
          reem_simulate_epi(prms, deterministic)
        )
      },
      
      
      # = = = = = = = = = = 
      # ====   Fit   ====
      # = = = = = = = = = = 
      
      traj_dist_obs = function(use.cl, 
                               use.ww, 
                               err.type,
                               deterministic,
                               n.sim = 8, 
                               verbose = FALSE) {
        return(
          reem_traj_dist_obs(
            obj           = .self, 
            use.cl        = use.cl, 
            use.ww        = use.ww, 
            err.type      = err.type, 
            deterministic = deterministic, 
            n.sim         = n.sim, 
            verbose       = verbose
          )
        )
      },
      
      fit_abc = function(prm.abc,
                         prms.to.fit){
        res = reem_fit_abc(
          obj         = .self,
          prm.abc     = prm.abc,
          prms.to.fit = prms.to.fit)
        
        .self$is.fitted = TRUE
        .self$fit.obj   = res
        .self$fit.prm   = prm.abc
        return(res)
      },
      
      
      plot_fit = function(){
        return(reem_plot_fit(obj = .self))  
      },
      
      # = = = = = = = = = = 
      # ==== Forecasts ====
      # = = = = = = = = = = 
      
      forecast = function(prm.fcst, verbose = FALSE){
        
        if(!.self$is.fitted) 
          stop('Model cannot forecast because it is not fitted.')
        
        res = reem_forecast(obj      = .self, 
                            prm.fcst = prm.fcst,
                            verbose  = verbose) 
        
        .self$fcst.prm <- prm.fcst
        .self$fcst.obj <- res
        
        return(res)
      },
      
      plot_forecast = function(
                        date_breaks = '2 months',
                        date_labels = '%b`%y'){
        res = reem_plot_forecast(
          obj = .self,
          date_breaks = date_breaks,
          date_labels = date_labels
        )
       return(res) 
      },
      
      # = = = = = = = = = = 
      # ==== Backtest ====
      # = = = = = = = = = = 
      
      calc_scores = function( 
                      var, 
                      obs.new, 
                      density.n = 100, 
                      density.adjust = 0.4,
                      aggr.window = NULL)
        {
            res = reem_calc_scores(
              var            = var,
              obs.new        = obs.new,
              fcst           = .self$fcst.obj,
              density.n      = density.n,
              density.adjust = density.adjust,
              aggr.window    = aggr.window
            )
            return(res)
      },
      
      forecast_densities = function(
                      var, 
                      obs.new, 
                      density.n = 100, 
                      density.adjust = 0.4,
                      aggr.window)
      {
        res = reem_forecast_densities(
          var            = var,
          obs.new        = obs.new,
          fcst           = .self$fcst.obj,
          density.n      = density.n,
          density.adjust = density.adjust, 
          aggr.window    = aggr.window
          )
        return(res)
      }
    
    )
)




