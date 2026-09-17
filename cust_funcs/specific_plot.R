library(ggplot2)

make_auto_effect_plot <- function(tip_cross, shock_type = c("Independent", "Correlated"),
                                  xmax = 5, var = "posaff_W", ncol = 1) {
  shock_type <- match.arg(shock_type)
  p <- tip_cross$Dynamics[[shock_type]][[1]]   # grabs whatever tipred is in there, by position not name
  
  df_auto <- p$data[p$data$row == var & p$data$col == var, ]
  df_auto$Effect <- droplevels(df_auto$Effect)
  
  p + df_auto +
    facet_wrap(~ Effect, ncol = ncol) +
    coord_cartesian(xlim = c(0, xmax))
}