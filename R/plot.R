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
#' @examples
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
              `Y: observed aggr. incidence` = 'gold4',
              `H: daily hosp admissions` = 'black')
  
  g.pop = sim.pop %>% 
    ggplot2::ggplot(ggplot2::aes(x=date, y = value, color = variable)) + 
    ggplot2::geom_line(linewidth = 1) + 
    ggplot2::scale_y_log10() + 
    ggplot2::scale_color_manual(values = col.pop) +
    ggplot2::theme_bw()+
    ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank() ) + 
    ggplot2::labs(title = 'Population')
  # g.pop 
  
  # Wastewater
  sim.ww = simepi$sim |> 
    dplyr::select(date, starts_with('W')) |>
    tidyr::pivot_longer(cols = starts_with('W'))
  
  sim.ww$variable = NA
  sim.ww$variable[sim.ww$name == 'Wd'] = 'Wd: concentration deposited'
  sim.ww$variable[sim.ww$name == 'Wp'] = 'Wd: concentration present at sampling site'
  sim.ww$variable[sim.ww$name == 'Wr'] = 'Wr: concentration reported at sampling site'
  
  g.ww = sim.ww %>% 
    ggplot2::ggplot(ggplot2::aes(x=date, y = value, color = variable)) + 
    ggplot2::geom_line(linewidth = 1) + 
    ggplot2::theme_bw()+
    ggplot2::labs(title = 'Pathogen concentration in wastewater')
  
  g = list(populations = g.pop, wastewater = g.ww)
  return(g)
}



#' Plot observed data associated with a `reem` object.
#'
#' @param obj A `reem` object.
#'
#' @return A ggplot object. 
#' @export
#'
#' @examples
#' 
plot_obs <- function(obj) {
  tmp = list() 
  tmp[['cl']] = obj$obs.cl |> mutate(type = 'case')
  tmp[['ha']] = obj$obs.ha |> mutate(type = 'hosp. adm.')
  tmp[['ww']] = obj$obs.ww |> mutate(type = 'ww')
  df = bind_rows(tmp)  
  
  g = df |> ggplot(aes(x=date, y=obs)) + 
    geom_step(color = 'grey60') + 
    geom_point() + 
    facet_wrap(~type, scales = 'free_y', ncol = 1) + 
    labs(title = 'Observed data', y = 'value')
  return(g)
}

