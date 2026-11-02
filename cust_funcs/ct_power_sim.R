# cust_funcs/ct_power_sim.R
#
# Monte Carlo power for a peak standardized cross-lagged effect in a bivariate
# continuous-time model (psx, ssit), following the logic of Hecht & Zitzmann
# (2021): generate from a known CT model, fit the same model, and count the
# replications whose 95% interval excludes 0.
#
# Generating model (standardized, as in H&Z):
#   - within-person asymptotic covariance = I, so a raw discrete-time cross
#     effect IS the standardized one (their Eq. 10)
#   - DRIFT auto-effects taken from the empirical fit; the cross-effect
#     ssit -> psx solved so that its PEAK standardized discrete-time effect
#     equals `peak`; DIFFUSION then set to -(A + A') so the asymptotic
#     covariance stays I
#   - occasion-level measurement error (theta, in within-SD^2 units) and a
#     random intercept per variable (ICC), added on the observed scale
#   - ESM schedule: `days` x `beeps` stratified-random prompts, MCAR
#     missingness at `compliance`
#   - censored scenario: psx is floored at the quantile that puts `p_floor`
#     of observations at the floor (0.744 for esm_psx_sum)
#
# Index convention matches 03_ple_model.qmd: latent 1 = psx, latent 2 = ssit,
# DRIFT[1, 2] = ssit -> psx (label "ssit_psx").
#
# Requires: ctsem (juliaFit branch), expm

# ---------------------------------------------------------------------------
# population values

#' Auto-effects and standardized measurement-error variances from a fitted
#' psx/ssit model (the first two latents are taken as the dynamic pair).
ct_power_pars_from_fit <- function(fit, latents = 1:2, manifests = 1:2) {
  e   <- ctsem::ctExtract(fit)
  med <- function(a) apply(a, c(2, 3), stats::median)
  A   <- med(e$pop_DRIFT)[latents, latents]
  S   <- med(e$pop_asymDIFFUSIONcov)[latents, latents]
  Th  <- med(e$pop_MANIFESTcov)[manifests, manifests]
  list(
    auto  = diag(A),                         # per unit of model time (days)
    theta = diag(Th) / diag(S),              # error var / within process var
    reliability = diag(S) / (diag(S) + diag(Th)),
    drift_hat = A,
    std_cross_hat = A[1, 2] * sqrt(S[2, 2] / S[1, 1])
  )
}

# discrete-time effect of `from` on `to` after interval dt (standardized when
# the asymptotic covariance S is supplied)
ct_dt_effect <- function(A, dt, to = 1, from = 2, S = diag(nrow(A))) {
  vapply(dt, function(t) expm::expm(A * t)[to, from], numeric(1)) *
    sqrt(S[from, from] / S[to, to])
}

ct_peak <- function(A, to = 1, from = 2, S = diag(nrow(A)),
                    tmax = 5 / min(abs(diag(A)))) {   # peak sits near 1/|auto|
  o <- stats::optimize(function(t) -abs(ct_dt_effect(A, t, to, from, S)),
                       c(1e-4, tmax))
  list(dt = o$minimum, value = ct_dt_effect(A, o$minimum, to, from, S))
}

#' DRIFT matrix whose peak standardized ssit -> psx effect equals `peak`
#' (asymptotic covariance fixed at I).
ct_power_truth <- function(auto, peak, other_cross = 0,
                           tmax = 5 / min(abs(auto))) {
  build <- function(x) matrix(c(auto[1], x, other_cross, auto[2]), 2, byrow = TRUE)
  if (peak == 0) {
    A <- build(0)
  } else {
    upper <- 0.99 * sqrt(4 * auto[1] * auto[2]) - abs(other_cross)  # keeps Q PD
    f <- function(x) abs(ct_peak(build(x), tmax = tmax)$value) - abs(peak)
    if (f(upper) < 0) stop("A peak of ", peak, " is not reachable with these ",
                           "auto-effects while keeping asymptotic variance 1.")
    A <- build(sign(peak) * stats::uniroot(f, c(1e-8, upper), tol = 1e-10)$root)
  }
  Q <- -(A + t(A))
  if (any(eigen(Q, only.values = TRUE)$values <= 0)) stop("DIFFUSION not PD.")
  pk <- ct_peak(A, tmax = tmax)
  list(A = A, Q = Q, peak_dt = pk$dt, peak_value = pk$value)
}

#' True model taken from a fitted psx/ssit model: the full within-person
#' DRIFT (both auto- and both cross-effects), rescaled so each process has
#' within-person variance 1 (correlation R kept). The observed peak
#' standardized ssit -> psx effect is the target; `scale` shrinks/inflates the
#' ssit -> psx drift (e.g. 0.5 as a conservative sensitivity check).
ct_power_truth_from_fit <- function(fit, latents = 1:2, scale = 1) {
  e   <- ctsem::ctExtract(fit)
  med <- function(a) apply(a, c(2, 3), stats::median)
  A   <- med(e$pop_DRIFT)[latents, latents]
  S   <- med(e$pop_asymDIFFUSIONcov)[latents, latents]
  L   <- diag(sqrt(diag(S)))
  As  <- solve(L) %*% A %*% L               # dynamics in within-SD units
  R   <- stats::cov2cor(S)                  # within-person correlation
  As[1, 2] <- As[1, 2] * scale
  Q   <- -(As %*% R + R %*% t(As))          # keeps asymptotic cov = R
  Q   <- (Q + t(Q)) / 2
  if (any(eigen(Q, only.values = TRUE)$values <= 0))
    stop("Implied DIFFUSION is not positive definite (try scale closer to 1).")
  pk  <- ct_peak(As)                        # S = I: units already standardized
  list(A = As, Q = Q, R = R, peak_dt = pk$dt, peak_value = pk$value,
       raw_drift = A, asym_cov = S)
}

# ---------------------------------------------------------------------------
# design and data

# one person's prompt times (in days): `beeps` stratified-random prompts per
# day between start_hr and end_hr
ct_power_times <- function(days = 7, beeps = 8, start_hr = 10, end_hr = 22) {
  edges <- seq(start_hr, end_hr, length.out = beeps + 1)
  hrs <- unlist(lapply(seq_len(days) - 1, function(d)
    d * 24 + stats::runif(beeps, edges[-(beeps + 1)], edges[-1])))
  (hrs - hrs[1]) / 24
}

ct_power_generate <- function(N, truth, theta = c(0, 0), icc = c(.48, .38),
                              days = 7, beeps = 8, compliance = .66,
                              p_floor = NULL) {
  Tp <- days * beeps
  dtmat <- t(vapply(seq_len(N), function(i) diff(ct_power_times(days, beeps)),
                    numeric(Tp - 1)))

  gen <- ctsem::ctModel(
    type = "ct", Tpoints = Tp,
    manifestNames = c("y_psx", "y_ssit"), latentNames = c("psx", "ssit"),
    LAMBDA = diag(2),
    DRIFT = truth$A,
    DIFFUSION = t(chol(truth$Q)),          # Cholesky factors (ctGenerate docs)
    T0VAR = t(chol(if (is.null(truth$R)) diag(2) else truth$R)),  # stationary start
    T0MEANS = matrix(0, 2, 1), CINT = matrix(0, 2, 1),
    MANIFESTMEANS = matrix(0, 2, 1),
    MANIFESTVAR = diag(sqrt(pmax(theta, 1e-8)), 2)
  )
  # ctGenerate's arguments differ across ctsem versions (n vs n.subjects;
  # Tpoints/backend only in newer ones), so pass only what this version has
  fa <- names(formals(ctsem::ctGenerate))
  args <- list(ctmodelobj = gen, dtmat = dtmat)
  args[[if ("n" %in% fa) "n" else "n.subjects"]] <- N
  if ("Tpoints" %in% fa) args$Tpoints <- Tp
  if ("backend" %in% fa) args$backend <- "r"
  d <- as.data.frame(suppressMessages(do.call(ctsem::ctGenerate, args)))
  names(d)[1:2] <- c("id", "time")

  # random intercepts on the observed scale: ICC = v_b / (v_b + 1 + theta)
  v_b <- icc / (1 - icc) * (1 + theta)
  b <- cbind(stats::rnorm(N, 0, sqrt(v_b[1])), stats::rnorm(N, 0, sqrt(v_b[2])))
  d$y_psx  <- d$y_psx  + b[d$id, 1]
  d$y_ssit <- d$y_ssit + b[d$id, 2]

  # MCAR missed prompts (whole beep missing)
  d <- d[stats::runif(nrow(d)) < compliance, ]

  floor_val <- NA_real_
  if (!is.null(p_floor)) {
    floor_val <- stats::qnorm(p_floor, 0, sqrt(v_b[1] + 1 + theta[1]))
    d$y_psx <- pmax(d$y_psx, floor_val)
  }
  attr(d, "floor") <- floor_val
  d
}

# ---------------------------------------------------------------------------
# analysis model (mirrors the generating structure; censoring optional)

#' ri: how person differences in level are modelled
#'   "manifest" = random MANIFESTMEANS (random intercept on each indicator)
#'   "cint"     = random CINT + T0MEANS on each dynamic process (process mean)
#'   "trend"    = 03_ple_model.qmd structure: observed = process + trend
#'                factor, random T0MEANS and CINT on the trend factors
ct_power_model <- function(censored = FALSE, floor_val = -Inf,
                           ri = c("manifest", "cint", "trend")) {
  ri <- match.arg(ri)
  common <- list(
    type = "ct", id = "id", time = "time",
    manifestNames = c("y_psx", "y_ssit"),
    manifesttype = c(if (censored) 4 else 0, 0),
    censormin = c(if (censored) floor_val else -Inf, -Inf),
    censormax = c(Inf, Inf),
    MANIFESTVAR = matrix(c("mvar_psx", 0, 0, "mvar_ssit"), 2, 2)
  )
  if (ri != "trend") {
    m <- do.call(ctsem::ctModel, c(common, list(
      latentNames = c("psx", "ssit"),
      LAMBDA = diag(2),
      DRIFT = matrix(c("psx_psx",  "ssit_psx",
                       "psx_ssit", "ssit_ssit"), 2, 2, byrow = TRUE),
      DIFFUSION = matrix(c("diff_psx",      0,
                           "diff_psx_ssit", "diff_ssit"), 2, 2, byrow = TRUE),
      T0MEANS = if (ri == "cint") matrix(c("T0m_psx", "T0m_ssit"), 2, 1) else matrix(0, 2, 1),
      CINT    = if (ri == "cint") matrix(c("cint_psx", "cint_ssit"), 2, 1) else matrix(0, 2, 1),
      MANIFESTMEANS = if (ri == "manifest") matrix(c("mm_psx", "mm_ssit"), 2, 1) else matrix(0, 2, 1)
    )))
    re <- if (ri == "manifest") c("mm_psx", "mm_ssit") else
      c("T0m_psx", "T0m_ssit", "cint_psx", "cint_ssit")
  } else {
    m <- do.call(ctsem::ctModel, c(common, list(
      latentNames = c("psx", "ssit", "psx_trend", "ssit_trend"),
      LAMBDA = cbind(diag(2), diag(2)),
      DRIFT = matrix(c(
        "psx_psx" , "ssit_psx" , 0                , 0                 ,
        "psx_ssit", "ssit_ssit", 0                , 0                 ,
        0         , 0          , "drift_psx_trend", 0                 ,
        0         , 0          , 0                , "drift_ssit_trend"), 4, 4, byrow = TRUE),
      DIFFUSION = matrix(c(
        "diff_psx"     , 0          , 0, 0,
        "diff_psx_ssit", "diff_ssit", 0, 0,
        0              , 0          , 0, 0,
        0              , 0          , 0, 0), 4, 4, byrow = TRUE),
      T0VAR = matrix(c(
        "T0var_psx"     , 0           , 0, 0,
        "T0var_psx_ssit", "T0var_ssit", 0, 0,
        0               , 0           , 0, 0,
        0               , 0           , 0, 0), 4, 4, byrow = TRUE),
      T0MEANS = c(0, 0, "T0m_psx_trend", "T0m_ssit_trend"),
      CINT    = c(0, 0, "cint_psx_trend", "cint_ssit_trend"),
      MANIFESTMEANS = matrix(0, 2, 1)
    )))
    re <- c("T0m_psx_trend", "T0m_ssit_trend", "cint_psx_trend", "cint_ssit_trend")
  }
  m$pars$indvarying <- m$pars$param %in% re
  m
}

# per-draw peak standardized effect and its 95% interval (H&Z's criterion)
ct_power_peak_ci <- function(fit, nsamples = 200) {
  e <- ctsem::ctExtract(fit)
  A <- e$pop_DRIFT; S <- e$pop_asymDIFFUSIONcov
  idx <- seq_len(dim(A)[1])
  if (length(idx) > nsamples) idx <- sort(sample(idx, nsamples))
  pk <- vapply(idx, function(s)
    ct_peak(A[s, 1:2, 1:2], S = S[s, 1:2, 1:2])$value, numeric(1))
  c(peak_est = stats::median(pk),
    peak_lower = unname(stats::quantile(pk, .025)),
    peak_upper = unname(stats::quantile(pk, .975)))
}

# ---------------------------------------------------------------------------
# one replication

#' scenario: "continuous" = uncensored data, continuous model
#'           "censored"   = floored data, censored model
#'           "naive"      = floored data, continuous model (misspecified)
ct_power_rep <- function(N, truth, scenario = c("continuous", "censored", "naive"),
                         theta, icc, p_floor = .744, peak_ci = FALSE,
                         cores = 1, censor_eps = 0, fit_args = list(),
                         ri = "trend", ...) {
  scenario <- match.arg(scenario)
  t0 <- Sys.time()
  d <- suppressMessages(ct_power_generate(N, truth, theta = theta, icc = icc,
                         p_floor = if (scenario == "continuous") NULL else p_floor, ...))
  fl <- attr(d, "floor")
  mod <- suppressMessages(ct_power_model(censored = scenario == "censored",
                        floor_val = if (is.na(fl)) -Inf else fl + censor_eps,
                        ri = ri))

  fit <- tryCatch(
    suppressMessages(do.call(ctsem::ctFit, utils::modifyList(
      list(datalong = d, model = mod, backend = "julia", cores = cores,
           plot = FALSE, priors = TRUE, verbose = 0), fit_args))),
    error = function(e) e)

  out <- data.frame(N = N, scenario = scenario,
                    true_drift = truth$A[1, 2], true_peak = truth$peak_value,
                    prop_floor = if (is.na(fl)) 0 else mean(d$y_psx <= fl + 1e-12),
                    nobs = nrow(d), est = NA_real_, lower = NA_real_,
                    upper = NA_real_, converged = FALSE, error = NA_character_)
  if (inherits(fit, "error")) {
    out$error <- conditionMessage(fit)
  } else {
    pm <- summary(fit)$popmeans
    out$est   <- pm["ssit_psx", "mean"]
    out$lower <- pm["ssit_psx", "2.5%"]
    out$upper <- pm["ssit_psx", "97.5%"]
    out$converged <- all(is.finite(c(out$est, out$lower, out$upper)))
    if (peak_ci && out$converged) {
      pk <- tryCatch(ct_power_peak_ci(fit), error = function(e) NULL)
      if (!is.null(pk)) out <- cbind(out, t(pk))
    }
  }
  out$sig <- out$converged & (out$lower > 0 | out$upper < 0)
  out$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out
}

#' many replications of one design cell, saved as it goes
ct_power_cell <- function(N, scenario, nrep, truth, theta, icc, seed,
                          outfile = NULL, ...) {
  set.seed(seed)
  res <- list()
  for (r in seq_len(nrep)) {
    res[[r]] <- cbind(rep = r, ct_power_rep(N, truth, scenario,
                                            theta = theta, icc = icc, ...))
    if (!is.null(outfile)) saveRDS(do.call(rbind, res), outfile)
  }
  do.call(rbind, res)
}

#' power with Monte Carlo SE, per N x scenario
ct_power_summary <- function(res) {
  res |>
    dplyr::group_by(N, scenario) |>
    dplyr::summarise(
      reps = dplyr::n(),
      conv_rate = mean(converged),
      power = mean(sig[converged]),
      mc_se = sqrt(power * (1 - power) / sum(converged)),
      bias_drift = mean(est[converged]) - true_drift[1],
      coverage = mean((lower <= true_drift & upper >= true_drift)[converged]),
      median_sec = stats::median(seconds),
      .groups = "drop")
}
