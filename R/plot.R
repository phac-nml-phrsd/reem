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
    dplyr::select(date, S,I,A,Y) |>
    tidyr::pivot_longer(cols = -date)
  
  sim.pop$variable = NA
  sim.pop$variable[sim.pop$name == 'S'] = 'S: susceptible'
  sim.pop$variable[sim.pop$name == 'I'] = 'I: incidence'
  sim.pop$variable[sim.pop$name == 'A'] = 'A: aggregated incidence'
  sim.pop$variable[sim.pop$name == 'Y'] = 'Y: observed aggr. incidence'
  
  g.pop = sim.pop %>% 
    ggplot2::ggplot(ggplot2::aes(x=date, y = value, color = variable)) + 
    ggplot2::geom_line() + 
    ggplot2::scale_y_log10() + 
    ggplot2::labs(title = 'Population')
  
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
    ggplot2::geom_line() + 
    ggplot2::labs(title = 'Pathogen concentration in wastewater')
  
  g = list(populations = g.pop, wastewater = g.ww)
  return(g)
}
