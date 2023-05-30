
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
        
        # Prepare dataframes for plotting
        
        ps = fit.obj$post.simulations 
        
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
          labs(title = 'Fit to clinical data', y = 'cases')
        # g.cl
        
        g.ww = ps.ww %>%
          drop_na(starts_with('Wr')) %>%
          ggplot2::ggplot(ggplot2::aes(x=date)) + 
          ggplot2::geom_line(ggplot2::aes(y = Wr.m), 
                             color = 'chocolate3')+
          ggplot2::geom_ribbon(ggplot2::aes(
            ymin = Wr.lo, 
            ymax = Wr.hi), 
            alpha=0.2, fill='chocolate')+
          ggplot2::geom_point(data = obs.ww, ggplot2::aes(y=obs)) +
          labs(title = 'Fit to wastewater data', y = 'concentration')
        # g.ww
        
        # ---- Posterior parameters
        
        # -- 1D density
        
        dp = rbind(mutate(fit.obj$all.distances, type = 'prior'),
                   mutate(fit.obj$post.prms, type = 'posterior')) %>%
          select(-starts_with('abc')) %>% 
          pivot_longer(cols = -type)
        
        col.pp =  c(posterior = 'indianred3', 
                    prior = 'gray70')
        
        g.post = dp %>% 
          ggplot2::ggplot(ggplot2::aes(x     = value, 
                                       color = type, 
                                       fill  = type)) + 
          ggplot2::geom_density(alpha = 0.3)+
          ggplot2::scale_color_manual(values = col.pp)+
          ggplot2::scale_fill_manual(values = col.pp)+
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
                ggplot2::theme(panel.grid = ggplot2::element_blank())+
                ggplot2::labs(x=nam[i], y=nam[j]) + 
                ggplot2::guides(fill = 'none')
              
              k = k+1
            }
          }
        }
        gpall2d = patchwork::wrap_plots(gp) +
          patchwork::plot_annotation(title = 'Posterior parameters 2D density')
        
        
        # -- Ordered ABC distances
        
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
          ggplot2::scale_color_manual(values = c('red2', 'gray80'))+
          ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) + 
          ggplot2::labs(
            title = 'ABC distances from data',
            x = 'ordered ABC iteration',
            y = 'distance'
          )
        
        
        g.list = list(
          traj.cl = g.cl,
          traj.ww = g.ww,
          post.prms = g.post,
          post.prms.2d = gpall2d,
          dist = g.dist
        )
        g.all = patchwork::wrap_plots(g.list) 
        
        res = c(g.list, list(all = g.all))
        return(res)
      },
      
      
      forecast = function(prm.fcst){
        
        if(!.self$is.fitted) 
          stop('Model cannot forecast because it is not fitted.')
        
        res = reem_forecast(obj = .self, prm.fcst = prm.fcst) 
        
        .self$fcst.prm <- prm.fcst
        .self$fcst.obj <- res
        
        return(res)
      },
      
      plot_forecast = function(){
        if(0){
          obs.cl = obj$obs.cl
          fcst.obj = obj$fcst.obj
          fcst.prm = obj$fcst.prm
        }
        
        col.fcst = 'steelblue2'
        col.fit  = 'tan2'
        
        sf = fcst.obj$summary.fcst %>% 
          filter(date >= fcst.prm$asof)
        
        obs.ww.before = filter(obs.ww, date <= fcst.prm$asof)
        obs.ww.after  = filter(obs.ww, date > fcst.prm$asof)
        
        fitsim.ww = obj$fit.obj$post.simulations %>% 
          bind_rows() %>% 
          group_by(date) %>% 
          summarize(m = mean(Wr), lo = min(Wr), hi = max(Wr)) %>% 
          drop_na(m)
        
        g = ggplot(data = drop_na(sf, Wr_mean),aes(x=date))+ 
          geom_line(data = fitsim.ww, aes(y=m), 
                    color = col.fit, linetype = 'dashed') + 
          geom_ribbon(data = fitsim.ww, aes(ymin=lo, ymax=hi), 
                      fill=col.fit, alpha = 0.1) + 
          geom_point(data = obs.ww.before, aes(y=obs))+ 
          geom_point(data = obs.ww.after, aes(y=obs), 
                     color='gray80')+ 
          geom_line( aes(y = Wr_mean), color= col.fcst, 
                     linetype = 'dotted') + 
          geom_ribbon(aes(ymin = Y_lo, ymax = Y_hi), 
                      alpha = 0.2, 
                      fill= col.fcst,
                      color= col.fcst) +
          geom_vline(xintercept = fcst.prm$asof, 
                     linetype = 'dashed', 
                     color = 'gray50') + 
          annotate(geom = 'text', y=1, x=fcst.prm$asof, 
                   label = fcst.prm$asof, size = 2) + 
          labs(title = 'Forecast wastewater concentration', 
               y = 'concentration')
        # g
         return(g)
        
      }
      
    )
)




