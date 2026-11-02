# cust_funcs/ct_dyn_checks.R
#
# Helpers for checking whether a first-order CT model is adequate
# (oscillation / second-order dynamics) and for timescale ("relay") models.
#
#   expct_prep()        within-person z-scores in the format expct() wants
#   expct_tidy()        expct() output -> long tibble (est, lower, upper)
#   ct_implied_lagcor() model-implied lagged correlations of the OBSERVED
#                       variables, to overlay on the expct curves
#   ct_drift_eigen()    are the fitted dynamics oscillatory? (complex eigenvalues)
#   ct_osc2()           period / damping ratio of a second-order (level + velocity) fit
#   ct_peak_draws()     per-draw peak time + height of two cross-effect curves,
#                       and the difference between directions
#
# All model-based helpers work from posterior/uncertainty draws via ctExtract(),
# using the pop_ arrays: pop_DRIFT, pop_asymDIFFUSIONcov, pop_MANIFESTcov,
# each [draw, row, col]. Model time is assumed to be in DAYS; results are
# reported in hours.

# ---------------------------------------------------------------------------
# internal: draws of the population system matrices
.ct_draws <- function(fit, mats = c("DRIFT", "asymDIFFUSIONcov", "MANIFESTcov"),
                      nsamples = 300) {
  e <- ctsem::ctExtract(fit)
  out <- list()
  for (m in mats) {
    a <- e[[paste0("pop_", m)]]
    if (is.null(a))
      stop("ctExtract() has no pop_", m, ". Available: ",
           paste(grep("^pop_", names(e), value = TRUE), collapse = ", "))
    out[[m]] <- a
  }
  nd <- dim(out[[1]])[1]
  keep <- if (is.numeric(nsamples) && nsamples < nd) sort(sample.int(nd, nsamples)) else seq_len(nd)
  lapply(out, function(a) a[keep, , , drop = FALSE])
}

# internal: matrix exponential exp(A * t) for a vector of t, via eigen-decomposition
# (fast), falling back to expm::expm when A is (near) defective.
.mexp_seq <- function(A, times) {
  ev <- eigen(A)
  V  <- ev$vectors
  kv <- tryCatch(kappa(V), error = function(e) Inf)
  ok <- is.finite(kv) && kv < 1e8
  if (ok) {
    Vi <- solve(V)
    lapply(times, function(t) Re(V %*% (exp(ev$values * t) * Vi)))
  } else {
    lapply(times, function(t) expm::expm(A * t))
  }
}

# ---------------------------------------------------------------------------
#' Prepare data for expct(): id, time (hours since each person's first beep),
#' and within-person z-scores. People with < min_obs observations or zero
#' variance on any variable are dropped (their z-scores are undefined).
#' Optionally removes a time-of-day rhythm (two harmonics of a 24-h cycle)
#' so that daily rhythms are not mistaken for oscillating dynamics.
expct_prep <- function(data, vars, id = "id", time_hr = "time_hour",
                       hour_of_day = NULL, remove_tod = FALSE, min_obs = 10) {
  d <- data %>%
    dplyr::transmute(id = .data[[id]], time = as.numeric(.data[[time_hr]]),
                     dplyr::across(dplyr::all_of(vars), as.numeric),
                     hod = if (!is.null(hour_of_day)) as.numeric(.data[[hour_of_day]]) else NA_real_) %>%
    dplyr::group_by(id) %>%
    dplyr::arrange(time, .by_group = TRUE) %>%
    dplyr::mutate(time = time - min(time, na.rm = TRUE))

  if (remove_tod) {
    if (all(is.na(d$hod))) stop("remove_tod = TRUE needs `hour_of_day`.")
    d <- d %>% dplyr::ungroup() %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(vars), function(v) {
        r <- rep(NA_real_, length(v)); ok <- !is.na(v) & !is.na(hod)
        h <- hod[ok] * 2 * pi / 24
        r[ok] <- stats::residuals(stats::lm(v[ok] ~ sin(h) + cos(h) + sin(2 * h) + cos(2 * h)))
        r
      })) %>% dplyr::group_by(id)
  }

  d <- d %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(vars),
                                ~ (.x - mean(.x, na.rm = TRUE)) / stats::sd(.x, na.rm = TRUE))) %>%
    dplyr::filter(dplyr::n() >= min_obs) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(vars), ~ any(is.finite(.x)))) %>%
    dplyr::ungroup() %>%
    dplyr::select(id, time, dplyr::all_of(vars))

  # expct cannot take missing values in the analysed variables
  d <- d[stats::complete.cases(d), ]
  dropped <- dplyr::n_distinct(data[[id]]) - dplyr::n_distinct(d$id)
  if (dropped > 0) message(dropped, " people dropped (too few obs or zero within-person variance).")
  as.data.frame(d)
}

#' expct() output -> long tibble. Names in expct are "XtoY":
#' correlation between X at time t and Y at time t + lag.
expct_tidy <- function(out, Tpred) {
  purrr::imap_dfr(out$est, function(est, nm) {
    tibble::tibble(
      effect = nm,
      from   = sub("to.*$", "", nm),
      to     = sub("^.*to", "", nm),
      lag_hr = Tpred,
      est    = est,
      lower  = if (!is.null(out$lowCI))  out$lowCI[[nm]]  else NA_real_,
      upper  = if (!is.null(out$highCI)) out$highCI[[nm]] else NA_real_
    )
  }) %>%
    dplyr::mutate(type = ifelse(from == to, "auto", "cross"))
}

# ---------------------------------------------------------------------------
#' Model-implied lagged correlations between OBSERVED variables.
#'
#' cov(latent_i(t + d), latent_j(t)) = [expm(DRIFT d) %*% asymDIFFUSIONcov]_ij
#' Measurement error is uncorrelated over time, so it only enters the
#' denominators: var(obs_i) = asymDIFFUSIONcov_ii + MANIFESTcov_ii.
#' Person-level trend states have no diffusion, so their asymptotic variance is
#' ~0 and they drop out -- matching expct's within-person centring.
#'
#' Caveat: for a censored indicator the observed (floored) correlation is
#' attenuated relative to this latent-scale value. Compare SHAPES (decay,
#' sign changes, where it peaks), not exact heights.
#'
#' @param state_idx  indices of the dynamic latent states (e.g. c(1, 2))
#' @param manifest_idx indices of their indicators in MANIFESTcov (e.g. c(1, 2))
#' @param names      labels, in the same order (should match the expct names)
ct_implied_lagcor <- function(fit, state_idx, manifest_idx, names,
                              lags_hr = 1:48, nsamples = 200, interval = 0.95,
                              lambda = NULL) {
  # lambda: optional manifest x latent loading matrix (rows = the manifests in
  # `names` order, cols = ALL latents). Needed when an indicator loads on more
  # than one dynamic state (e.g. fast + slow). If NULL, each manifest is taken
  # to load only on state_idx[k].
  dr <- .ct_draws(fit, nsamples = nsamples)
  t_day <- lags_hr / 24
  n_lat <- dim(dr$DRIFT)[2]
  if (is.null(lambda)) {
    lambda <- matrix(0, length(state_idx), n_lat)
    lambda[cbind(seq_along(state_idx), state_idx)] <- 1
  }
  res <- vector("list", dim(dr$DRIFT)[1])
  for (s in seq_along(res)) {
    A  <- dr$DRIFT[s, , ];  S <- dr$asymDIFFUSIONcov[s, , ];  Th <- as.matrix(dr$MANIFESTcov[s, , ])
    v  <- diag(lambda %*% S %*% t(lambda)) + diag(Th)[manifest_idx]
    Ms <- .mexp_seq(A, t_day)
    res[[s]] <- purrr::map2_dfr(Ms, lags_hr, function(M, L) {
      C <- lambda %*% M %*% S %*% t(lambda)
      g <- expand.grid(i = seq_along(names), j = seq_along(names))
      tibble::tibble(from = names[g$j], to = names[g$i], lag_hr = L,
                     r = C[cbind(g$i, g$j)] / sqrt(v[g$i] * v[g$j]))
    }) %>% dplyr::mutate(draw = s)
  }
  a <- (1 - interval) / 2
  dplyr::bind_rows(res) %>%
    dplyr::group_by(from, to, lag_hr) %>%
    dplyr::summarise(implied = stats::median(r), implied_lower = stats::quantile(r, a),
                     implied_upper = stats::quantile(r, 1 - a), .groups = "drop") %>%
    dplyr::mutate(effect = paste0(from, "to", to))
}

# ---------------------------------------------------------------------------
#' Can the fitted (first-order, multivariate) dynamics oscillate?
#' Oscillation <=> complex eigenvalues of the dynamic block of DRIFT.
#' Period (hours) = 2*pi / |Im(lambda)| * 24.
ct_drift_eigen <- function(fit, state_idx, nsamples = 500) {
  dr <- .ct_draws(fit, mats = "DRIFT", nsamples = nsamples)
  ev <- purrr::map_dfr(seq_len(dim(dr$DRIFT)[1]), function(s) {
    l <- eigen(dr$DRIFT[s, state_idx, state_idx])$values
    tibble::tibble(draw = s, max_im = max(abs(Im(l))), max_re = max(Re(l)))
  })
  tibble::tibble(
    p_oscillatory   = mean(ev$max_im > 1e-8),   # share of draws with complex eigenvalues
    p_stable        = mean(ev$max_re < 0),      # all eigenvalues negative real part
    period_hr_median = if (any(ev$max_im > 1e-8))
      stats::median(2 * pi / ev$max_im[ev$max_im > 1e-8] * 24) else NA_real_
  )
}

#' Oscillation summary for a second-order (level + velocity) model.
#' d(level) = velocity ; d(velocity) = a * level + b * velocity
#' a = DRIFT[vel, level] (restoring force, < 0), b = DRIFT[vel, vel] (damping, < 0)
#' Oscillates when b^2 + 4a < 0 (equivalently damping ratio zeta < 1).
ct_osc2 <- function(fit, level_idx = 1, vel_idx = 2, nsamples = 1000, interval = 0.95) {
  dr <- .ct_draws(fit, mats = "DRIFT", nsamples = nsamples)
  a <- dr$DRIFT[, vel_idx, level_idx]
  b <- dr$DRIFT[, vel_idx, vel_idx]
  disc   <- b^2 + 4 * a
  zeta   <- -b / (2 * sqrt(-a))
  period <- ifelse(disc < 0, 2 * pi / (sqrt(-disc) / 2) * 24, NA_real_)   # hours
  q <- function(x) stats::quantile(x, c((1 - interval) / 2, .5, 1 - (1 - interval) / 2), na.rm = TRUE)
  tibble::tibble(
    quantity = c("restoring force a", "damping b", "damping ratio zeta",
                 "period (hours, oscillating draws only)"),
    lower  = c(q(a)[1], q(b)[1], q(zeta)[1], q(period)[1]),
    median = c(q(a)[2], q(b)[2], q(zeta)[2], q(period)[2]),
    upper  = c(q(a)[3], q(b)[3], q(zeta)[3], q(period)[3])
  ) %>%
    dplyr::bind_rows(tibble::tibble(quantity = "P(oscillating) = P(zeta < 1)",
                                    median = mean(disc < 0)))
}

# ---------------------------------------------------------------------------
#' Peak time and peak height of two cross-effect curves, per draw.
#'
#' Effect of state `a` on state `b` after interval t = expm(DRIFT t)[b, a].
#' This automatically includes paths through relay (timescale) states.
#' standardise = TRUE rescales by sqrt(asymDIFFUSIONcov[a,a] / [b,b]) so the
#' two directions are on a comparable (SD-unit) scale.
#'
#' @param from,to   state indices: direction 1 is from[1] -> to[1], direction 2
#'                  is from[2] -> to[2]
#' @param labels    two labels, e.g. c("dis -> ssit", "ssit -> dis")
#' @param times_hr  grid of intervals (hours) to search for the peak
#' @return list(draws, summary, difference)
ct_peak_draws <- function(fit, from, to, labels, times_hr = seq(0.05, 24, by = 0.05),
                          standardise = TRUE, nsamples = 300, interval = 0.95) {
  dr <- .ct_draws(fit, mats = c("DRIFT", "asymDIFFUSIONcov"), nsamples = nsamples)
  t_day <- times_hr / 24
  draws <- purrr::map_dfr(seq_len(dim(dr$DRIFT)[1]), function(s) {
    A <- dr$DRIFT[s, , ];  S <- dr$asymDIFFUSIONcov[s, , ]
    Ms <- .mexp_seq(A, t_day)
    purrr::map_dfr(1:2, function(k) {
      eff <- vapply(Ms, function(M) M[to[k], from[k]], numeric(1))
      if (standardise) eff <- eff * sqrt(S[from[k], from[k]] / S[to[k], to[k]])
      w <- which.max(abs(eff))
      tibble::tibble(draw = s, direction = labels[k],
                     peak_time_hr = times_hr[w], peak_height = eff[w])
    })
  })
  a <- (1 - interval) / 2
  summ <- draws %>%
    dplyr::group_by(direction) %>%
    dplyr::summarise(
      time_lower   = stats::quantile(peak_time_hr, a),
      time_median  = stats::median(peak_time_hr),
      time_upper   = stats::quantile(peak_time_hr, 1 - a),
      height_lower = stats::quantile(peak_height, a),
      height_median = stats::median(peak_height),
      height_upper = stats::quantile(peak_height, 1 - a),
      .groups = "drop")
  w <- tidyr::pivot_wider(draws, id_cols = draw, names_from = direction,
                          values_from = c(peak_time_hr, peak_height))
  dt <- w[[paste0("peak_time_hr_", labels[1])]] - w[[paste0("peak_time_hr_", labels[2])]]
  dh <- abs(w[[paste0("peak_height_", labels[1])]]) - abs(w[[paste0("peak_height_", labels[2])]])
  diffs <- tibble::tibble(
    contrast = c(paste0("peak time: ", labels[1], " minus ", labels[2], " (hours)"),
                 paste0("|peak height|: ", labels[1], " minus ", labels[2])),
    lower  = c(stats::quantile(dt, a), stats::quantile(dh, a)),
    median = c(stats::median(dt), stats::median(dh)),
    upper  = c(stats::quantile(dt, 1 - a), stats::quantile(dh, 1 - a)),
    p_greater_than_0 = c(mean(dt > 0), mean(dh > 0))
  )
  list(draws = draws, summary = summ, difference = diffs)
}

# ---------------------------------------------------------------------------
#' Compatibility patch for expct with newer mgcv.
#' expct calls predict(..., se.fit = "TRUE") -- a character string. Newer mgcv
#' versions test `se.fit || get.ci`, which errors on a string
#' ("invalid 'x' type in 'x || y'"). This wraps mgcv's predict.bam / predict.gam
#' so se.fit is coerced to logical. Session-only; safe to call repeatedly.
expct_patch_mgcv <- function() {
  ns <- asNamespace("mgcv")
  for (cls in c("bam", "gam")) {
    fname <- paste0("predict.", cls)
    orig  <- get(fname, envir = ns)
    if (isTRUE(attr(orig, "expct_patched"))) next
    fixed <- local({ o <- orig
      function(object, newdata, type = "link", se.fit = FALSE, ...)
        o(object, newdata, type = type, se.fit = isTRUE(as.logical(se.fit)), ...) })
    attr(fixed, "expct_patched") <- TRUE
    unlockBinding(fname, ns); assign(fname, fixed, envir = ns); lockBinding(fname, ns)
    assign(fname, fixed, envir = get(".__S3MethodsTable__.", envir = ns))
  }
  invisible(TRUE)
}
