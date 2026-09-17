# cust_funcs/ct_halflife.R
#
# Half-life for each within-person DRIFT-diagonal (autoregression)
# parameter of a fitted ctStanFit.
#
# For a stable continuous-time process, a deviation from equilibrium decays
# as exp(drift * t); half-life is the time at which that decay reaches
# 0.5 of its starting size:
#   0.5 = exp(drift * t)  =>  t = log(0.5) / drift = log(2) / -drift
# This is only meaningful when drift < 0 (a genuinely stable/mean-reverting
# process) -- a zero or positive diagonal never returns to equilibrium, and
# is left NA rather than silently printed as a negative "half-life".
#
# CI note: log(2)/-x is monotonic *decreasing* in x for x < 0, so it's
# valid to map the DRIFT parameter's 2.5%/97.5% bounds through it
# endpoint-by-endpoint to get the half-life's interval -- but because the
# map is decreasing, the more-negative (faster-decaying) drift bound
# becomes the *shorter* half-life bound, so the two swap order relative to
# the drift columns they came from.
#
# Usage inside a .qmd chunk (after source("cust_funcs/ct_halflife.R")):
#   ct_halflife(psx_twopart_fit, psx_twopart_syn$latentNames) %>% gt::gt()
# Every DRIFT-diagonal latent is returned by default; pass pattern = NULL
# to include between-person latents too (they're fixed at 0 here, so
# they'll just show NA).
ct_halflife <- function(fit, latentNames, pattern = "_W$", digits = 2) {
  pm <- summary(fit)$parmatrices
  pm <- pm[pm$matrix == "DRIFT" & pm$row == pm$col, , drop = FALSE]
  pm$latent <- latentNames[pm$row]
  if (!is.null(pattern)) pm <- pm[grepl(pattern, pm$latent), , drop = FALSE]

  pm %>%
    dplyr::transmute(
      latent,
      drift_mean  = Mean,
      drift_2.5   = `2.5%`,
      drift_97.5  = `97.5%`,
      halflife_days       = dplyr::if_else(drift_mean < 0, log(2) / -drift_mean, NA_real_),
      halflife_days_2.5   = dplyr::if_else(drift_97.5 < 0, log(2) / -drift_97.5, NA_real_),
      halflife_days_97.5  = dplyr::if_else(drift_2.5  < 0, log(2) / -drift_2.5,  NA_real_),
      halflife_hours      = halflife_days      * 24,
      halflife_hours_2.5  = halflife_days_2.5  * 24,
      halflife_hours_97.5 = halflife_days_97.5 * 24
    ) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~round(.x, digits)))
}
