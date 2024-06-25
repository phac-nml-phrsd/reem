
#' Class to represent a Renewal Equation Epidemic Model.
#'
#' @field name String for the name of the model. 
#' Just for convenience when working with multiple models,
#' this does not affect any calculation.
#' 
#' @field prms List of model parameters. (NEEDS TO BE CHECKED!)
#' \itemize{
#'   \item{\code{horizon}: }{horizon time of the simulation.}
#'   \item{\code{last.obs}: }{last observation time (must be < \code{horizon}).}
#'   \item{\code{B}: }{numerical vector representing the time 
#'   dependent multiplicative factor for transmission.}
#'   \item{\code{date.start}: }{start date of the epidemic.}
#'   \item{\code{R0}: }{basic reproduction number.}
#'   \item{\code{N}: }{population size.}
#'   \item{\code{alpha}: }{transmission heterogeneity (alpha=-Inf: homogeneous).}
#'   \item{\code{I.init}: }{initial incidence (overwritten in fit ABC).}
#'   \item{\code{lag}: }{aggregation window for clinical reports.}
#'   \item{\code{rho}: }{ mean reporting ratio.}
#'   \item{\code{g}: }{intrinsic generation interval distribution.}
#'   \item{\code{fec}: }{fecal shedding kinetics.}
#'   \item{\code{h.prop}: }{total proportion of infections that will be hospitalized.}
#'   \item{\code{h.lags}: }{time lags repartition weights between infection and hospital admission times.}
#'   \item{\code{kappa}: }{pathogen decay in wastewater.}
#'   \item{\code{psi}: }{plug flow simulation.}
#'   \item{\code{shed.mult}: }{deposited fecal shedding multiplier.}
#' }
#'   
#' @field is.fitted Logical. Is the model fitted to the observation
#' data attached (\code{obs.cl} and/or \code{obs.ww})? 
#' @field obs.cl Dataframe representing the clinical observations. 
#' Must have \code{date} and \code{obs} variables. 
#' @field obs.ha Dataframe representing the hospital admissions. 
#' Must have \code{date} and \code{obs} variables. 
#' @field obs.ww Dataframe representing the wastewater observations. 
#' Must have \code{date} and \code{obs} variables. 
#' @field fit.obj List containing the fitted object. 
#' @field fit.prm List that defines the parameters for the fitting algorithm.
#' 
#' @return An object of class \code{reem}.
#' 
#' @export
#' 
#' 
#'
#' 
setRefClass(
  Class = 'reem', 
  fields =  list(
    name      = "character",
    prms      = "list",
    obs.cl    = "data.frame",
    obs.ha    = "data.frame",
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
        print(paste0('--- Parameters for REEM `',name,'`'))
        for(i in seq_along(prms)){
          nam = names(prms)[i]
          val = prms[[i]]
          print(paste(nam,' = ', paste(val,collapse = ' ')))
        }
        print(' --------------------------')
      },
      
      simulate = function(deterministic){
        res = reem_simulate(prms, deterministic)  
        return(res)  
      },
      
      simulate_epi = function(deterministic){
        return(
          reem_simulate_epi(obj = .self, deterministic)
        )
      },

      
      # = = = = = = = = = = 
      # ====   Fit   ====
      # = = = = = = = = = = 
      
      traj_dist_obs = function(use.cl, 
                               use.ha, 
                               use.ww, 
                               err.type,
                               deterministic,
                               n.sim = 8, 
                               verbose = FALSE) {
        return(
          reem_traj_dist_obs(
            obj           = .self, 
            use.cl        = use.cl, 
            use.ha        = use.ha, 
            use.ww        = use.ww, 
            err.type      = err.type, 
            deterministic = deterministic, 
            n.sim         = n.sim, 
            verbose       = verbose
          )
        )
      },
      
      fit_abc = function(prm.abc,
                         prms.to.fit,
                         verbose = FALSE){
        res = reem_fit_abc(
          obj         = .self,
          prm.abc     = prm.abc,
          prms.to.fit = prms.to.fit,
          verbose     = verbose)
        
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
        
        .self$fcst.prm <- prm.fcst
        
        res = reem_forecast(obj      = .self, 
                            prm.fcst = prm.fcst,
                            verbose  = verbose) 
        
        .self$fcst.obj <- res
        
        return(res)
      },
      
      plot_forecast = function(
                        date_breaks = '2 months',
                        date_labels = '%b`%y', 
                        logscale = FALSE){
        res = reem_plot_forecast(
          obj         = .self,
          date_breaks = date_breaks,
          date_labels = date_labels,
          logscale    = logscale
        )
       return(res) 
      },
      
      
      plot_peak = function(var, 
                           date_labels = '%d %b`%y', 
                           logscale = FALSE){
        res = reem_plot_peak(
          var,
          obj         = .self,
          date_labels = date_labels,
          logscale    = logscale
        )
        return(res)
      },
      
      
      get_fcst_density = function(
              var, 
              date.future, 
              density.n = 100, 
              density.adjust = 0.4,
              aggr.window = NULL)
      {
        res = reem_get_fcst_density(
          var            = var,
          date.future    = date.future,
          fcst           = .self$fcst.obj,
          density.n      = density.n,
          density.adjust = density.adjust, 
          aggr.window    = aggr.window
        )
        return(res)
      },
      
      # ???DELETE BECAUSE DUPLICATES `get_fcst_density`???
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
      },
      
      forecast_peak = function(var){
        res = reem_forecast_peak(var = var, fcst = .self$fcst.obj)
        return(res)
      },
      
      proba_box = function(var, 
                           date.lower, 
                           date.upper,
                           val.lower, 
                           val.upper){
        res = reem_proba_box(var = var, 
                             date.lower = date.lower, 
                             date.upper = date.upper,
                             val.lower  = val.lower, 
                             val.upper  = val.upper,
                             fcst = .self$fcst.obj)
        return(res)
      },
      
      # = = = = = = = = = = 
      # ==== Backtest ====
      # = = = = = = = = = = 
     
      
      inside_CI = function(
                    var,
                    obs.new,
                    aggr.window,
                    ci.width)
      {
        res = reem_inside_CI(var, 
                             obs.new, 
                             ci.width, 
                             aggr.window,
                             fcst = .self$fcst.obj)
          return(res)
      },
      
       
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
      }
    )
)




