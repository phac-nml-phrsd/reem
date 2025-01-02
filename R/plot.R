#' Plot simulated epidemic
#'
#' @param simepi 
#'
#' @return A list of two ggplot objects plotting
#' the time series of the population groups and
#' pathogen concentration in wastewater.
#' 
#' @export
#'
plot_epi <- function(simepi) {
  
  # Populations
  sim.pop = simepi$sim |> 
    dplyr::select(date, S,I,A,Y,H) |>
    tidyr::pivot_longer(cols = -date)
  
  sim.pop$variable = NA
  sim.pop$variable[sim.pop$name == 'S'] = 'S: susceptible'
  sim.pop$variable[sim.pop$name == 'I'] = 'I: incidence'
  sim.pop$variable[sim.pop$name == 'A'] = 'A: aggregated incidence'
  sim.pop$variable[sim.pop$name == 'Y'] = 'Y: observed aggr. incidence'
  sim.pop$variable[sim.pop$name == 'H'] = 'H: daily hosp admissions'
  
  col.pop = c(`S: susceptible` = 'blue2', 
              `I: incidence` = 'red3', 
              `A: aggregated incidence` = 'gold', 
              `Y: observed aggr. incidence` = 'green4',
              `H: daily hosp admissions` = 'black')
  
  g.pop = sim.pop %>% 
    ggplot2::ggplot(ggplot2::aes(x=date, y = value, color = variable)) + 
    ggplot2::geom_line(linewidth = 1) + 
    ggplot2::scale_y_log10() + 
    ggplot2::scale_color_manual(values = col.pop) +
    ggplot2::theme_bw()+
    ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank() ) + 
    ggplot2::labs(title = 'Populations Simulated')
  # g.pop 
  
  # Wastewater
  sim.ww = simepi$sim |> 
    dplyr::select(date, dplyr::starts_with('W')) |>
    tidyr::pivot_longer(cols = dplyr::starts_with('W'))
  
  sim.ww$variable = NA
  sim.ww$variable[sim.ww$name == 'Wd'] = 'Wd: concentration deposited'
  sim.ww$variable[sim.ww$name == 'Wp'] = 'Wp: conc. present at sampling site'
  sim.ww$variable[sim.ww$name == 'Wr'] = 'Wr: conc. reportable (lab errors)'
  
  g.ww = sim.ww %>% 
    ggplot2::ggplot(ggplot2::aes(x=date, y = value, color = variable)) + 
    ggplot2::geom_line(data = dplyr::filter(sim.ww, name!='Wr'),
                       linewidth = 1) + 
    ggplot2::geom_point(data = dplyr::filter(sim.ww, name=='Wr')) + 
    ggplot2::theme_bw()+
    ggplot2::scale_color_manual(values = c('tan4', 'tan', 'tan2')) +
    ggplot2::labs(title = 'Pathogen concentration in wastewater')
  # g.ww
  
  # Observations 
  
  obs = rbind(
    dplyr::mutate(simepi$obs.cl, type = 'clinical reports'),
    dplyr::mutate(simepi$obs.ha, type = 'hosp. admissions'),
    dplyr::mutate(simepi$obs.ww, type = 'wastewater')
    )
  
  g.obs = obs |> 
    ggplot2::ggplot(ggplot2::aes(x=date, y=obs)) + 
    ggplot2::geom_step(color = 'grey') + 
    ggplot2::geom_point() +  
    ggplot2::theme_bw()+
    ggplot2::facet_wrap(~type, ncol = 1, scales = 'free_y') + 
    ggplot2::labs(title = 'Simulated observations', x='', y='value')
  
  g = list(
    populations = g.pop, 
    wastewater = g.ww,
    observations = g.obs)
  return(g)
}



#' Plot observed data associated with a `reem` object.
#'
#' @param obj A `reem` object.
#'
#' @return A ggplot object. 
#' @export
#' 
plot_obs <- function(obj) {
  tmp = list() 
  tmp[['cl']] = obj$obs.cl |> dplyr::mutate(type = 'case')
  tmp[['ha']] = obj$obs.ha |> dplyr::mutate(type = 'hosp. adm.')
  tmp[['ww']] = obj$obs.ww |> dplyr::mutate(type = 'ww')
  df = dplyr::bind_rows(tmp)  
  
  g = df |> ggplot(aes(x=date, y=obs)) + 
    geom_step(color = 'grey60') + 
    geom_point() + 
    facet_wrap(~type, scales = 'free_y', ncol = 1) + 
    labs(title = 'Observed data', y = 'value')
  return(g)
}

