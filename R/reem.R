
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
    is.fitted = "logical",
    fit.obj   = "list"
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
        return(res)
      },
      
      
      plot_fit = function(){
        ps = fit.obj$post.simulations %>% 
          dplyr::bind_rows() %>% 
          dplyr::group_by(date) %>%
          dplyr::summarise(Y.m = mean(Y),
                           Y.lo = min(Y),
                           Y.hi = max(Y),
                           Wr.m = mean(Wr),
                           Wr.lo = min(Wr),
                           Wr.hi = max(Wr))
        
        ggplot2::theme_set(ggplot2::theme_bw())
        
        # ---- Trajectories
        
        g.cl = ps %>% 
          ggplot2::ggplot(ggplot2::aes(x=date)) + 
          ggplot2::geom_line(ggplot2::aes(y = Y.m), color = 'chartreuse3')+
          ggplot2::geom_ribbon(ggplot2::aes(ymin=Y.lo, ymax=Y.hi), 
                               alpha=0.2, fill='chartreuse')+
          ggplot2::geom_point(data = obs.cl, ggplot2::aes(y=obs)) +
          labs(title = 'Fit to clinical data', y = 'cases')
        
        g.ww = ps %>%
          drop_na(starts_with('Wr')) %>%
          ggplot2::ggplot(ggplot2::aes(x=date)) + 
          ggplot2::geom_line(ggplot2::aes(y = Wr.m), color = 'chocolate3')+
          ggplot2::geom_ribbon(ggplot2::aes(ymin=Wr.lo, ymax=Wr.hi), 
                               alpha=0.2, fill='chocolate')+
          ggplot2::geom_point(data = obs.ww, ggplot2::aes(y=obs)) +
          labs(title = 'Fit to wastewater data', y = 'concentration')
        
        # ---- Posterior parameters
        
        # -- 1D density
        
        dp = rbind(mutate(fit.obj$all.distances, type = 'prior'),
                   mutate(fit.obj$post.prms, type = 'posterior')) %>%
          select(-starts_with('abc')) %>% 
          pivot_longer(cols = -type)
        
        g.post = dp %>% 
          ggplot2::ggplot(ggplot2::aes(x= value, 
                                       color = type, fill = type)) + 
          ggplot2::geom_density(alpha = 0.3)+
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
                ggplot2::theme(panel.grid = element_blank())+
                ggplot2::labs(x=nam[i], y=nam[j]) + 
                guides(fill = 'none')
              
              k = k+1
            }
          }
        }
        gpall2d = patchwork::wrap_plots(gp) +
          patchwork::plot_annotation(title = 'Posterior parameters 2D density')
        
        
        g.list = list(
          traj.cl = g.cl,
          traj.ww = g.ww,
          post.prms = g.post,
          post.prms.2d = gpall2d
        )
        g.all = patchwork::wrap_plots(g.list) 
        
        res = c(g.list, list(all = g.all))
        return(res)
      }
      
    )
)




