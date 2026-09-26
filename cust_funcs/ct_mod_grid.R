#' Moderation (simple-slopes) plots for every TIpred x DRIFT effect
#'
#' Loops plot_ribbon_multi() (mod_plot.R) over TIpreds and named
#' "predictor>outcome" effects, so any fit can get the full grid in one call.
#' Each plot shows the discrete-time effect at -1SD / median / +1SD of the
#' TIpred (other TIpreds held at 0 = their mean, since they are z-scored).
#'
#' @param fit      fitted ctsem model
#' @param latentNames latent names in model order (e.g. syn$latentNames)
#' @param effects  character vector "predictor>outcome"; bare "pa" = auto-effect
#' @param tipreds  TIpred names to loop over
#' @param labels   optional named vector: tipred name -> legend title
#' @param ...      passed to plot_ribbon_multi() (times, nsamp_use, ...)
#' @return nested list: out[[tipred]][[effect]] = plot_ribbon_multi() result
ct_mod_grid <- function(fit, latentNames, effects, tipreds,
                        labels = NULL, ...) {
  parse_eff <- function(s) {
    p <- trimws(strsplit(s, ">", fixed = TRUE)[[1]])
    if (length(p) == 1) p <- c(p, p)
    bad <- setdiff(p, latentNames)
    if (length(bad)) stop("Unknown latent(s) in '", s, "': ", paste(bad, collapse = ", "))
    list(predictor = p[1], outcome = p[2],
         row = match(p[2], latentNames),       # DRIFT row = outcome
         col = match(p[1], latentNames),       # DRIFT col = predictor
         label = if (p[1] == p[2]) paste0(p[1], " (auto)") else paste0(p[1], " -> ", p[2]))
  }
  effs <- lapply(effects, parse_eff)

  out <- list()
  for (tp in tipreds) {
    out[[tp]] <- list()
    for (e in effs) {
      message("== ", tp, " x ", e$label, " ==")
      res <- plot_ribbon_multi(
        fit, tipred_name = tp, row = e$row, col = e$col,
        plot_title   = e$label,
        y_lab        = "Effect",
        legend_title = if (!is.null(labels) && tp %in% names(labels)) labels[[tp]] else tp,
        ...
      )
      out[[tp]][[e$label]] <- res
    }
  }
  out
}
