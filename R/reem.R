
#' Class to represent a Renewal Equation Epidemic Model.
#'
#' @field name character. 
#' @field prms list. 
#' @field is.fitted logical. 
#' 
#' @return
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
    is.fitted = "logical"
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
      
      traj_dist_obs = function(use.cl, 
                               use.ww, 
                               err.type,
                               deterministic,
                               n.sim = 10, 
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
          
        return(res)
      }
    )
  )




