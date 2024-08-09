library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(lubridate)
library(stringr)
library(patchwork)

x = seq(0,1, by = 1e-3)

ff <- function(alpha, x) {
  y = x^exp(alpha)
  df = data.frame(x=x, y=y, alpha = as.factor(alpha))
  return(df)
}

as <- c(-3, -1, 0, 1, 3)

d = lapply(as, ff, x=x) |> bind_rows() 


g = ggplot(d, aes(x=x, y = y, color = alpha )) + 
  geom_line(linewidth = 3) + 
  theme(panel.grid.minor = element_blank()) + 
  labs(x = 'Susceptible fraction available', 
       y = 'Susceptible fraction exposed')+
  scale_color_brewer(palette = 'RdYlGn')
# g
pdf('alpha.pdf', width = 4, height = 3)
plot(g)
dev.off()