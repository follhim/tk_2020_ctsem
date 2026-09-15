save_postpred_plots <- function(pp_plot, folder, width = 10, height = 7, dpi = 150) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  for (panel_name in names(pp_plot)) {
    p <- pp_plot[[panel_name]]
    if (!inherits(p, "ggplot")) next
    ggsave(
      filename = file.path(folder, paste0(panel_name, ".png")),
      plot = p, width = width, height = height, dpi = dpi
    )
    cat("saved:", panel_name, "\n")
  }
}