plot_ribbon_multi <- function(fit, tipred_name = NULL, interaction = NULL,
                              levels_val = NULL,
                              times = seq(0, 5, by = 0.1),
                              row = 1, col = 1,
                              nsamp_use = 200,
                              colors = c("blue", "black", "red"),
                              plot_title = NULL,
                              x_lab = "Time interval",
                              y_lab = NULL,
                              legend_title = NULL,
                              line_width = 0.5,
                              bound_width = 0.3,
                              ribbon_alpha = 0.1) {

  ctmb <- ctsem:::.ctFitModelObject(fit)
  tipred_names <- ctmb$TIpredNames
  ntipred <- length(tipred_names)

  # ---- "interaction = 'tip_<tipred>_<param>'" shortcut ----
  # Lets you pass the exact row name printed by summary(fit)$tipreds (e.g.
  # "tip_zpos_auto_bin_psx_W") instead of looking up tipred_name/row/col by
  # hand. NOTE: this deliberately does NOT look anything up in ctmb$pars --
  # on a fitted ctStanFit object, .ctFitModelObject()'s pars$param column
  # has been overwritten with compiled Stan expressions (e.g. "state[5]"),
  # not the original "auto_bin_psx_W"-style label, so it can't be matched
  # against post-fit. Instead this relies on how ct_apply_crosslevel()
  # (lavaan_parser.R) always targets cross-level interactions: only
  # "auto_<latent>" (the DRIFT diagonal / self-effect cell, row == col ==
  # that latent's index) is supported by the DSL, so the latent name alone
  # is enough to recover row/col from ctmb$latentNames.
  if (!is.null(interaction)) {
    match_tp <- tipred_names[vapply(tipred_names, function(tp)
      startsWith(interaction, paste0("tip_", tp, "_")), logical(1))]
    if (length(match_tp) == 0) {
      stop("interaction '", interaction, "' doesn't start with 'tip_<tipred>_' for ",
           "any declared TIpred: ", paste(tipred_names, collapse = ", "))
    }
    # if one tipred name happens to be a prefix of another, prefer the
    # longer (more specific) match
    tipred_name <- match_tp[which.max(nchar(match_tp))]
    param_label <- sub(paste0("^tip_", tipred_name, "_"), "", interaction)

    if (!startsWith(param_label, "auto_")) {
      stop("interaction '", interaction, "' names DRIFT parameter '", param_label,
           "', which isn't an 'auto_<latent>' self-effect -- only within-latent ",
           "autoregression interactions (the ones ct_apply_crosslevel() can ",
           "target) are supported by this shortcut. Pass tipred_name/row/col ",
           "directly for anything else.")
    }
    latent_name <- sub("^auto_", "", param_label)
    idx <- match(latent_name, ctmb$latentNames)
    if (is.na(idx)) {
      stop("interaction '", interaction, "' -- latent '", latent_name, "' not ",
           "found in ctmb$latentNames: ", paste(ctmb$latentNames, collapse = ", "))
    }
    row <- idx
    col <- idx
  }

  if (is.null(tipred_name)) {
    stop("Provide either tipred_name (+ row/col), or interaction = 'tip_<tipred>_<param>'.")
  }

  tipred_index <- match(tipred_name, tipred_names)

  if (is.na(tipred_index)) {
    stop("tipred_name '", tipred_name, "' not found among: ",
         paste(tipred_names, collapse = ", "))
  }

  if (is.null(levels_val)) {
    tipred_raw <- fit$standata$tipredsdata[, tipred_index]
    m  <- mean(tipred_raw, na.rm = TRUE)
    s  <- sd(tipred_raw, na.rm = TRUE)
    md <- median(tipred_raw, na.rm = TRUE)
    levels_val <- c(m - s, md, m + s)
    message("Auto-detected levels (-1SD, median, +1SD) for ", tipred_name, ": ",
            paste(round(levels_val, 3), collapse = ", "))
  }

  rp_full <- ctsem:::.ctFitRawPosterior(fit)
  nsamp <- min(nsamp_use, nrow(rp_full))
  rp <- rp_full[sample(nrow(rp_full), nsamp), , drop = FALSE]

  if (is.null(legend_title)) legend_title <- tipred_name
  if (is.null(y_lab)) y_lab <- ctmb$latentNames[row]
  if (is.null(plot_title)) {
    plot_title <- paste0("Effect of ", tipred_name, " on ", ctmb$latentNames[row], " dynamics")
  }

  out <- vector("list", length(levels_val))

  for (li in seq_along(levels_val)) {
    tipreds <- rep(0, ntipred)
    tipreds[tipred_index] <- levels_val[li]

    message("Level ", li, "/", length(levels_val), " (", tipred_name, " = ", round(levels_val[li], 3), "):")

    draws_list <- pbapply::pblapply(seq_len(nsamp), function(s) {
      DRIFT <- suppressMessages(
        ctsem:::ctBackendParMatrices(fit, raw = rp[s, ], tipreds = tipreds, trim = FALSE)$DRIFT
      )
      sapply(times, function(t) expm::expm(DRIFT * t)[row, col])
    })
    draws <- do.call(rbind, draws_list)

    df <- data.frame(
      `Time interval` = times,
      med = apply(draws, 2, median),
      lo  = apply(draws, 2, quantile, 0.025),
      hi  = apply(draws, 2, quantile, 0.975),
      CovariateValue = factor(round(levels_val[li], 3),
                              levels = round(sort(levels_val), 3)),
      check.names = FALSE
    )
    out[[li]] <- df
  }

  rb <- do.call(rbind, out)

  p <- ggplot(rb, aes(x = `Time interval`, color = CovariateValue, fill = CovariateValue)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = ribbon_alpha, color = NA) +
    geom_line(aes(y = lo), linetype = "dotted", linewidth = bound_width) +
    geom_line(aes(y = hi), linetype = "dotted", linewidth = bound_width) +
    geom_line(aes(y = med), linewidth = line_width) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(title = plot_title, x = x_lab, y = y_lab,
         color = legend_title, fill = legend_title) +
    theme_minimal()

  list(plot = p, data = rb)
}
