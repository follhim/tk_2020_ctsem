#' Discrete-time effect plots for any ctsem fit, by latent name
#'
#' Wraps ctDiscretePars() (ctsem >= 3.12). Effects are written
#' "predictor>outcome" (e.g. "sit>pa" = sit -> pa); a bare name ("pa") or
#' "pa>auto" is the auto-effect. Leave `effects` NULL to take every
#' auto/cross effect among the latents not listed in `exclude`.
#'
#' @param fit        a ctFit() / ctStanFit() object
#' @param effects    character vector of "predictor>outcome" specs, or NULL
#' @param type       "all", "auto" or "cross" -- which effects to keep
#' @param exclude    latent names to drop (e.g. trend / random-intercept factors)
#' @param drop       effects to leave out, as "predictor>outcome"; "a<>b" drops
#'                   both directions (e.g. "bin_psx<>con_psx" for cross-effects
#'                   the model does not estimate)
#' @param times      time intervals to compute (model time units)
#' @param facet      TRUE = auto and cross effects in separate panels
#' @param interval   width of the credible band
#' @param ...        passed to ctDiscretePars() (standardise, impulseType,
#'                   subjects, nsamples, ...)
#' @return invisibly, list(plot_auto, plot_cross, plot, summary, peaks):
#'   separate auto- and cross-effect plots (NULL when there are none), plus the
#'   old combined two-panel `plot`; the separate plots are also printed.
#'   summary = median + interval per effect x time; peaks = largest |median|
#'   of each cross-effect and the time interval where it occurs.
ct_dpars_plot <- function(fit, effects = NULL, type = c("all", "auto", "cross"),
                          exclude = NULL, drop = NULL, times = seq(0, 1, by = 0.01),
                          facet = TRUE, interval = 0.95, title = NULL, ...) {
  type <- match.arg(type)

  g  <- ctsem::ctDiscretePars(fit, times = times, plot = TRUE, ...)
  dt <- tryCatch(g@data, error = function(e) g$data)
  dt <- as.data.frame(dt)
  dt$Effect <- NULL                # ctsem's own label -- replaced by ours below
  dt$row <- as.character(dt$row)   # outcome
  dt$col <- as.character(dt$col)   # predictor
  latents <- unique(c(dt$row, dt$col))

  bad_ex <- setdiff(exclude, latents)
  if (length(bad_ex))
    stop("`exclude` takes latent (process) names, not parameter names -- unknown: ",
         paste(bad_ex, collapse = ", "), ". Latents in this fit: ",
         paste(latents, collapse = ", "))

  # ---- which (predictor, outcome) pairs ----
  if (is.null(effects)) {
    keep <- setdiff(latents, exclude)
    pairs <- expand.grid(predictor = keep, outcome = keep, stringsAsFactors = FALSE)
  } else {
    pairs <- do.call(rbind, lapply(trimws(effects), function(s) {
      p <- trimws(strsplit(s, ">", fixed = TRUE)[[1]])
      if (length(p) == 1 || tolower(p[2]) == "auto") p <- c(p[1], p[1])
      if (length(p) != 2) stop("Could not parse effect '", s, "' -- use 'predictor>outcome'.")
      bad <- setdiff(p, latents)
      if (length(bad)) stop("Unknown latent(s) in '", s, "': ", paste(bad, collapse = ", "),
                            ". Available: ", paste(latents, collapse = ", "))
      data.frame(predictor = p[1], outcome = p[2], stringsAsFactors = FALSE)
    }))
  }
  # ---- drop requested effects ("a>b", or "a<>b" for both directions) ----
  if (length(drop)) {
    dl <- do.call(rbind, lapply(trimws(drop), function(s) {
      both <- grepl("<>", s, fixed = TRUE)
      p <- trimws(strsplit(s, if (both) "<>" else ">", fixed = TRUE)[[1]])
      if (length(p) != 2) stop("Could not parse drop '", s, "' -- use 'a>b' or 'a<>b'.")
      bad <- setdiff(p, latents)
      if (length(bad)) stop("Unknown latent(s) in drop '", s, "': ", paste(bad, collapse = ", "))
      out <- data.frame(predictor = p[1], outcome = p[2], stringsAsFactors = FALSE)
      if (both) out <- rbind(out, data.frame(predictor = p[2], outcome = p[1], stringsAsFactors = FALSE))
      out
    }))
    key <- paste(pairs$predictor, pairs$outcome)
    pairs <- pairs[!key %in% paste(dl$predictor, dl$outcome), , drop = FALSE]
  }
  pairs$type <- ifelse(pairs$predictor == pairs$outcome, "Auto-effects", "Cross-effects")
  if (type == "auto")  pairs <- pairs[pairs$type == "Auto-effects", ]
  if (type == "cross") pairs <- pairs[pairs$type == "Cross-effects", ]
  if (!nrow(pairs)) stop("No effects left to plot after filtering.")
  pairs$Effect <- ifelse(pairs$type == "Auto-effects",
                         paste0(pairs$predictor, " (auto)"),
                         paste0(pairs$predictor, " -> ", pairs$outcome))

  d <- merge(dt, pairs, by.x = c("col", "row"), by.y = c("predictor", "outcome"))

  # ---- summarise over posterior samples ----
  a <- (1 - interval) / 2
  summ <- d %>%
    dplyr::filter(is.finite(value)) %>%
    dplyr::group_by(type, Effect, Subject, `Time interval`) %>%
    dplyr::summarise(median = stats::median(value),
                     lower  = stats::quantile(value, a),
                     upper  = stats::quantile(value, 1 - a), .groups = "drop") %>%
    dplyr::mutate(Effect = factor(Effect, levels = unique(pairs$Effect)))

  # peak of each cross-effect (auto-effects always peak at 1 at interval 0)
  peaks <- summ %>%
    dplyr::filter(type == "Cross-effects") %>%
    dplyr::group_by(type, Effect, Subject) %>%
    dplyr::slice_max(abs(median), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  # ---- plots: one for auto-effects, one for cross-effects ----
  multi_subj <- length(unique(summ$Subject)) > 1
  base_title <- if (is.null(title))
    tryCatch(g@labels$title, error = function(e) g$labels$title) else title
  make_plot <- function(sub, ttl) {
    if (!nrow(sub)) return(NULL)
    sub$Effect <- droplevels(sub$Effect)
    q <- ggplot2::ggplot(sub, ggplot2::aes(`Time interval`, median,
                                           colour = Effect, fill = Effect,
                                           group = interaction(Effect, Subject))) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = .15, colour = NA) +
      # dotted lines on the interval bounds, same colour as the effect
      ggplot2::geom_line(ggplot2::aes(y = lower), linetype = "dotted", linewidth = .6) +
      ggplot2::geom_line(ggplot2::aes(y = upper), linetype = "dotted", linewidth = .6) +
      ggplot2::geom_line(linewidth = .8) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey50", linewidth = .3) +
      ggplot2::labs(title = ttl, subtitle = base_title,
                    x = "Time interval", y = "Effect") +
      ggplot2::theme_bw()
    if (multi_subj) q <- q + ggplot2::aes(linetype = Subject)
    q
  }
  plot_auto  <- make_plot(summ[summ$type == "Auto-effects", ],  "Auto-effects")
  plot_cross <- make_plot(summ[summ$type == "Cross-effects", ], "Cross-effects")

  # combined version kept for backwards compatibility (`$plot`)
  p <- ggplot2::ggplot(summ, ggplot2::aes(`Time interval`, median,
                                          colour = Effect, fill = Effect,
                                          group = interaction(Effect, Subject))) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = .15, colour = NA) +
      # dotted lines on the interval bounds, same colour as the effect
      ggplot2::geom_line(ggplot2::aes(y = lower), linetype = "dotted", linewidth = .6) +
      ggplot2::geom_line(ggplot2::aes(y = upper), linetype = "dotted", linewidth = .6) +
    ggplot2::geom_line(linewidth = .8) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50", linewidth = .3) +
    ggplot2::labs(title = base_title, x = "Time interval", y = "Effect") +
    ggplot2::theme_bw()
  if (multi_subj) p <- p + ggplot2::aes(linetype = Subject)
  if (facet && length(unique(summ$type)) > 1)
    p <- p + ggplot2::facet_wrap(~type, scales = "free_y")

  if (!is.null(plot_auto))  print(plot_auto)
  if (!is.null(plot_cross)) print(plot_cross)
  invisible(list(plot_auto = plot_auto, plot_cross = plot_cross, plot = p,
                 summary = summ, peaks = peaks))
}
