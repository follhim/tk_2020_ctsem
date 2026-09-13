ct_effects_plot <- function(fit,
                            effects,
                            times = seq(0, 2, .1),
                            standardise = TRUE,
                            xlim = NULL,
                            ylim = NULL,
                            title = "Discrete-time effects over time",
                            legend_labels = NULL,
                            latentNames = NULL) {
  
  # ---- resolve latent names (VERIFY this field name against your fit object) ----
  if (is.null(latentNames)) {
    latentNames <- tryCatch(fit$ctstanmodelbase$latentNames, error = function(e) NULL)
    if (is.null(latentNames)) {
      stop("Could not auto-detect latentNames from `fit`. Run ",
           "str(fit, max.level = 1) to find the right field, or just pass ",
           "latentNames = c(...) explicitly (e.g. the vector from your ",
           "ct_lavaan() spec -- m$latentNames).")
    }
  }
  
  # ---- parse "A>B" (A predicts B), "A>auto" or bare "A" (auto-effect) ----
  parse_one <- function(spec) {
    spec <- trimws(spec)
    if (grepl(">", spec, fixed = TRUE)) {
      parts <- trimws(strsplit(spec, ">", fixed = TRUE)[[1]])
      if (length(parts) != 2) stop("Could not parse effect spec: '", spec, "'")
      predictor <- parts[1]
      outcome   <- if (tolower(parts[2]) == "auto") parts[1] else parts[2]
    } else {
      predictor <- spec; outcome <- spec
    }
    if (!predictor %in% latentNames) stop("Unknown latent '", predictor, "' in: '", spec, "'")
    if (!outcome   %in% latentNames) stop("Unknown latent '", outcome,   "' in: '", spec, "'")
    
    # DRIFT[row, col] = effect of column (predictor) on row (outcome) --
    # per the doc's own convention ("columns are predictor, rows are outcome")
    list(
      predictor = predictor, outcome = outcome,
      row = match(outcome, latentNames),
      col = match(predictor, latentNames),
      label = if (predictor == outcome) paste0(predictor, " (auto)")
      else paste0(predictor, " -> ", outcome),
      spec = spec
    )
  }
  
  parsed <- lapply(effects, parse_one)
  cr_indices <- do.call(rbind, lapply(parsed, function(p) c(p$row, p$col)))
  
  # ---- compute discrete-time parameters once ----
  x <- ctStanDiscretePars(fit, times = times, standardise = standardise, plot = FALSE)
  
  # ---- pull raw draws for exactly these effects ----
  # CHECK: run `str(out$dt)` / `unique(out$dt$Effect)` the first time you use
  # this, to confirm indices + ggcode=TRUE behaves as assumed here, and that
  # dt's row order lines up with `parsed` (used below to relabel Effect).
  out <- ctStanDiscreteParsPlot(x, indices = cr_indices, ggcode = TRUE)
  dt <- out$dt
  
  # Remap ctsem's own Effect labels to yours, by POSITION (order effects
  # were requested in), not by string-matching -- ctsem's auto-generated
  # names won't necessarily match your "A -> B" labels.
  effect_levels <- unique(dt$Effect)
  if (length(effect_levels) != length(parsed)) {
    warning("Number of effects in dt$Effect (", length(effect_levels),
            ") doesn't match number requested (", length(parsed), ") -- ",
            "inspect `out$dt` manually before trusting the relabeling below.")
  }
  relabel <- setNames(vapply(parsed, function(p) p$label, character(1)), effect_levels)
  dt$Effect <- relabel[as.character(dt$Effect)]
  
  # ---- summarize: median + 95% CrI at each timepoint, per effect ----
  summ <- dt %>%
    dplyr::filter(is.finite(value)) %>%
    dplyr::group_by(Effect, `Time interval`) %>%
    dplyr::summarize(
      median = median(value),
      lower  = quantile(value, .025),
      upper  = quantile(value, .975),
      .groups = "drop"
    )
  
  # ---- peak + CrI + time-of-peak, per effect ----
  peak_tbl <- summ %>%
    dplyr::group_by(Effect) %>%
    dplyr::slice_max(order_by = abs(median), n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::rename(peak_time = `Time interval`, peak_median = median,
                  peak_lower_2.5 = lower, peak_upper_97.5 = upper) %>%
    dplyr::arrange(dplyr::desc(abs(peak_median)))
  
  print(peak_tbl)
  
  # ---- one combined plot, full manual control ----
  p <- ggplot2::ggplot(summ, ggplot2::aes(x = `Time interval`, y = median,
                                          color = Effect, fill = Effect)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = .15, color = NA) +
    ggplot2::geom_line(ggplot2::aes(y = lower), linetype = "dotted", linewidth = .5) +
    ggplot2::geom_line(ggplot2::aes(y = upper), linetype = "dotted", linewidth = .5) +
    ggplot2::geom_line(linewidth = .8) +
    ggplot2::geom_hline(yintercept = 0, color = "grey50", linewidth = .3) +
    ggplot2::labs(title = title, x = "Time interval", y = "Effect size") +
    ggplot2::theme_bw()
  
  if (!is.null(legend_labels)) {
    p <- p + ggplot2::scale_color_discrete(name = "Effect", labels = legend_labels) +
      ggplot2::scale_fill_discrete(name = "Effect", labels = legend_labels)
  }
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
  }
  
  print(p)
  invisible(list(plot = p, summary = summ, peaks = peak_tbl))
}