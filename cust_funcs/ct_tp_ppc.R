#' Two-part posterior predictive check (occurrence + magnitude)
#'
#' ctPostPredPlots() compares simulated magnitudes at every beep where the
#' magnitude was OBSERVED, including simulated moments where the model's own
#' occurrence was 0 -- which biases the magnitude comparison low. This check
#' keeps a simulated magnitude only where the model's simulated occurrence is 1.
#'
#' Call after ctGenerateFromFit() and before setting fit$generated <- NULL.
#'
#' @param fit    fit with fit$generated$Y (samples x rows x manifests)
#' @param data   the data frame the model was fitted to (same rows)
#' @param bin    name of the binary occurrence manifest
#' @param mag    name of the magnitude manifest
#' @param breaks NULL for ordinal / integer magnitude (tabulate rounded values);
#'               numeric cut points for a continuous magnitude
#' @param labels optional labels for the `breaks` categories
#' @param time   time column (days) for the mean-by-day comparison
#' @return list(prop = category proportions, byday = mean magnitude by day,
#'              plot = ggplot of byday)
ct_tp_ppc <- function(fit, data, bin, mag, breaks = NULL, labels = NULL,
                      time = "time_day") {
  Y <- fit$generated$Y
  if (is.null(Y)) stop("fit$generated$Y is missing -- run ctGenerateFromFit() first.")
  if (dim(Y)[2] != nrow(data)) stop("generated rows (", dim(Y)[2],
                                    ") != nrow(data) (", nrow(data), ").")
  g_bin <- Y[, , bin]
  g_mag <- Y[, , mag]
  g_mag[is.na(g_bin) | g_bin == 0] <- NA          # model's own occurrence decides
  o_bin <- as.numeric(data[[bin]])
  o_mag <- as.numeric(data[[mag]])
  g_bin_obs <- g_bin[, !is.na(o_bin), drop = FALSE]   # occurrence where observed

  categorise <- if (is.null(breaks)) {
    function(x) round(x)
  } else {
    lab <- if (is.null(labels)) FALSE else labels
    function(x) cut(x, breaks = breaks, labels = lab)
  }
  ptab <- function(o, g, v, f) {
    fo <- f(o); fg <- f(as.vector(g))
    lev <- if (is.factor(fo)) levels(fo) else
      sort(unique(c(fo[!is.na(fo)], fg[!is.na(fg)])))
    ot <- prop.table(table(factor(fo, levels = lev)))
    gt <- prop.table(table(factor(fg, levels = lev)))
    data.frame(variable = v, category = as.character(lev),
               observed = as.numeric(ot), model = as.numeric(gt))
  }
  prop <- rbind(
    ptab(o_bin, g_bin_obs, bin, round),
    ptab(o_mag, g_mag, paste0(mag, " (given occurrence)"), categorise))
  rownames(prop) <- NULL

  day <- floor(as.numeric(data[[time]]))
  byday <- data.frame(day = day, observed = o_mag,
                      model = colMeans(g_mag, na.rm = TRUE))
  byday <- do.call(rbind, lapply(split(byday, byday$day), function(d)
    data.frame(day = d$day[1],
               observed = mean(d$observed, na.rm = TRUE),
               model = mean(d$model, na.rm = TRUE))))
  long <- rbind(data.frame(day = byday$day, source = "observed", mean = byday$observed),
                data.frame(day = byday$day, source = "model", mean = byday$model))
  plot <- ggplot2::ggplot(long, ggplot2::aes(day, mean, colour = source)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::labs(x = "Study day", y = paste0("Mean ", mag, " (given occurrence)"),
                  colour = NULL) +
    ggplot2::theme_bw()
  list(prop = prop, byday = byday, plot = plot)
}
