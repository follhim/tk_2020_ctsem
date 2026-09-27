# How the 02 Positive Affect auto model (`pa_trend`) is specified

A plain-language walkthrough of `pa_trend_mod` in `02_pa_model.qmd`: what each piece of the
`ctModel()` call does, why it is set the way it is, and what it corresponds to in models you
already know (AR(1)/DSEM, RI-CLPM, latent growth curves).

---

## 1. The big picture in one sentence

Every positive-affect (PA) rating at a beep is treated as the sum of three things:

> **observed PA = slow-moving personal baseline (`pa_trend`) + fast-moving momentary state (`pa`) + beep-specific noise**

**Analogy - the ocean.** The *tide* is the baseline: it moves slowly and predictably over the
week. The *waves* are the momentary state: they get kicked up by events and die back down. The
*spray* is measurement noise: it is there at one beep and gone at the next.

**Comparisons**

| This model | What you would call it elsewhere |
|---|---|
| `pa_trend` (the tide) | The random intercept in an RI-CLPM, or the person mean in DSEM, except that here it is allowed to drift slowly over the week instead of being flat |
| `pa` (the waves) | The within-person, person-mean-centered part in DSEM or RI-CLPM, which has an AR(1) auto-effect |
| Measurement noise (the spray) | Residual or measurement error variance, like the within-level residual in a DSEM with measurement error |

This is the "trend factor" setup from the ctsem EAS example. It separates *where a person
usually is* (and whether that is changing across the week) from *how their mood bounces around
that level*. Those two things are kept apart so the auto-effect of the momentary state is not
contaminated by slow baseline changes. If a person's PA declines steadily across the week, a
model without a trend would mistake that slow slide for very high inertia.

---

## 2. Data going in

- **Standardization:** `esm_pa`, `zpos`, `zneg` and `zdis` are standardized with `scale()`, so
  effects are in SD units.
- **Time:** `time_day` is time since study start, measured in days. Continuous time means
  unequal gaps between beeps are fine. The model knows that 20 minutes and 5 hours are
  different amounts of time.
- **Rows:**
  - One row per person per time point; duplicated timestamps are dropped.
  - Beeps where PA is missing are dropped.
  - Nobody is dropped as a person. The Kalman filter uses whatever beeps each person has, with
    FIML-style handling of missingness.

---

## 3. Continuous time vs. the usual AR(1): the one idea to hold onto

In a discrete-time model (DSEM or Mplus), the auto-effect is a single number, for example
phi = .40, meaning "40% of this beep's deviation carries to the next beep." That assumes every
gap between beeps is the same length.

In continuous time, the auto-effect is a **drift** (a rate per day). The carry-over for any gap
of length Δt is

> **carry-over(Δt) = exp(drift × Δt)**

So the same drift gives more carry-over across a short gap and less across a long one.

- **Worked example:** with `pa_pa` = -3 (per day), the carry-over is
  - about .83 across a 1.5-hour gap (0.0625 day),
  - about .47 across 6 hours,
  - about .05 overnight (about 1 day).
- **More negative drift** means faster return to baseline: less inertia, quicker emotional
  recovery.
- **Drift near 0** means the state lingers: high inertia.
- **Half-life = ln(2) / |drift|.** This is how long it takes for half of a PA bump to fade.
  `ct_halflife()` reports it, and it is the most intuitive way to talk about `pa_pa`.

`ct_dpars_plot()` draws exactly this curve, carry-over as a function of time interval. So
instead of one phi you get the whole "how long does a PA bump last" curve.

---

## 4. Each part of the specification

### `manifestNames`, `latentNames`, `TIpredNames`

```r
manifestNames = c("esm_pa")                 # what you observed
latentNames   = c("pa", "pa_trend")         # the two hidden processes
TIpredNames   = c("zpos", "zneg", "zdis")   # person-level (time-invariant) predictors
```

- There is **one observed variable** and **two latent variables**. Each latent is a process
  that evolves over time.
- **TIpreds** are the schizotypy scores. They are measured once per person and used as
  moderators (see Section 6).

### `manifesttype`, `censormin`, `censormax`

```r
manifesttype = c(0); censormin = c(-Inf); censormax = c(Inf)
```

- `0` means continuous (Gaussian) measurement. PA is a roughly continuous, symmetric composite,
  so it doesn't need the censored treatment used for NA, PSX, etc. in 03–06.
- Censoring bounds of +/-Inf mean no censoring is active.

**Comparison:** ordinary continuous indicators in SEM.

### `LAMBDA`: how the latents produce the observation

```r
LAMBDA = cbind(diag(1), diag(1))   # = [ 1  1 ]
```

- **What it does:** it says `esm_pa = 1*pa + 1*pa_trend (+ noise)`. Both loadings are fixed to 1.
- **Why fixed:** with a single indicator, there is nothing to estimate a loading from. The
  observed score is simply split into two additive parts.
- **Analogy:** the tide plus the waves equals the water level you see on the pier.
- **Comparison:** in an RI-CLPM, the observed score = random intercept (loading 1) + within
  component (loading 1). Same idea here.

### `DRIFT`: how each process changes on its own

```r
DRIFT = [ pa_pa   0              ]
        [ 0       drift_pa_trend ]
```

- **`pa_pa`** is the **auto-effect of the momentary state**. It is the main parameter of
  interest: how quickly a PA deviation from someone's baseline fades. Its interpretation is
  covered in Section 3.
  - **Comparison:** phi in DSEM or the AR path in RI-CLPM, in continuous-time form.
- **`drift_pa_trend`** is the trend's own drift. It controls how the baseline curves over the
  week (see T0MEANS/CINT below).
  - The trend has no noise, so this drift is not an "inertia". It sets the shape of a smooth
    growth curve.
- **Off-diagonal 0s** mean the baseline and the momentary state don't push each other around.
  The momentary state is always a deviation from the baseline, not a cause of it.

### `DIFFUSION`: random shocks

```r
DIFFUSION = [ diff_pa  0 ]
            [ 0        0 ]
```

- **`diff_pa`** is how much unpredictable "stuff happens" input hits the momentary state over
  time: a good conversation, a bad text. Shocks enter, then fade according to `pa_pa`.
  - **Comparison:** the innovation or residual variance of the within-person AR process in
    DSEM.
  - **Analogy:** wind that creates new waves.
- **Trend diffusion is fixed to 0.** The baseline is not buffeted by random shocks. It follows
  a smooth path set by its start (T0MEANS), its target (CINT) and its drift. This is what keeps
  "baseline" and "momentary" distinguishable: one is smooth, the other is jumpy.

### `T0VAR`: variability at each person's first beep

```r
T0VAR = [ T0var_pa  0 ]
        [ 0         0 ]
```

- **`T0var_pa`** captures the fact that people's momentary state at their very first beep is
  somewhere around their baseline, not exactly on it. This is the variance of that first-beep
  deviation.
- **Trend T0VAR is 0** because individual differences in where the baseline starts are
  handled by making `T0m_pa_trend` a random effect (below). Having both would estimate the same
  thing twice.
- **Comparison:** the initial-condition variance in a state-space or Kalman model. DSEM
  handles this implicitly by conditioning on or discarding the first observation.

### `T0MEANS`: where things start

```r
T0MEANS = c(0, "T0m_pa_trend")
```

- **`pa` starts at 0.** The momentary state is a deviation from baseline, and on average people
  are not above or below their own baseline at the first beep.
- **`T0m_pa_trend` is free and random** (person-specific). It is each person's **baseline PA
  level at the start of the study**.
- **Comparison:** the random intercept in an RI-CLPM, or the intercept factor in a latent
  growth curve.

### `CINT`: where things are heading

```r
CINT = c(0, "cint_pa_trend")
```

- **`pa` CINT is 0.** The momentary state is always pulled back toward 0, i.e. back to the
  person's own baseline. This is what "regulation back to baseline" means in the model.
- **`cint_pa_trend` is free and random.** Together with `drift_pa_trend`, it sets where each
  person's baseline is heading:
  - **Long-run baseline (asymptote)** = -cint_pa_trend / drift_pa_trend
  - **Baseline path:** trend(t) = asymptote + (start - asymptote) * exp(drift_pa_trend * t)

  Because both the start (T0MEANS) and the target (via CINT) vary by person, each person's
  baseline can rise, fall or stay flat across the week, following a smooth exponential-shaped
  curve.
- **Comparison:** a latent growth curve with a random intercept (T0MEANS) and a random
  "where it ends up" component (CINT), using an exponential rather than linear shape. If
  someone's start and target are the same, their baseline is flat. That is exactly an RI-CLPM
  random intercept.

### `MANIFESTMEANS`

```r
MANIFESTMEANS = 0
```

- The measurement intercept is fixed at 0 because the trend already carries the level of PA.
- Freeing both a measurement intercept and the trend's level would ask the model to estimate
  the same mean twice, which is not identified.
- **Comparison:** in an RI-CLPM, you put the mean on either the observed intercept or the
  random intercept, never both.

### `MANIFESTVAR`: measurement or beep-specific noise

```r
MANIFESTVAR = "mvar_pa"
```

- This is variance that appears at a single beep and **does not carry over** to the next:
  - rating error,
  - momentary distraction while answering,
  - using the scale slightly differently at that beep.
- **How it differs from DIFFUSION:**
  - DIFFUSION is a real change in PA that lingers and fades.
  - MANIFESTVAR is noise that is gone by the next beep.
- **Analogy:** spray (MANIFESTVAR) vs. a wave (DIFFUSION).
- **Comparison:** DSEM with measurement error, versus plain DSEM. In plain DSEM this noise gets
  lumped into the innovation and biases phi toward 0. Separating it gives a cleaner auto-effect.

---

## 5. Random effects (`indvarying`)

```r
pa_trend_mod$pars$indvarying <- pa_trend_mod$pars$param %in% c("T0m_pa_trend", "cint_pa_trend")
```

- **Only two parameters differ from person to person:** where each person's baseline starts,
  and where it is heading. Their correlation is also estimated, so the model can ask whether
  people who start high tend to end high.
- **Everything else is a population value shared by everyone**, including the auto-effect
  `pa_pa`.
- **Comparison with DSEM:** in DSEM you would often make phi a random slope. Here, individual
  differences in inertia are modeled only through the schizotypy moderators (next section),
  not through a free random effect. This keeps the model fast and focused on the question of
  interest.

---

## 6. Schizotypy moderators (TIpreds)

```r
pa_trend_mod$pars[, eff] <- FALSE
pa_trend_mod$pars[pa_trend_mod$pars$matrix == "DRIFT" &
                  pa_trend_mod$pars$param %in% "pa_pa", eff] <- TRUE
```

- **What it does:** `zpos`, `zneg` and `zdis` are allowed to shift **only the auto-effect
  `pa_pa`**. The question is: *do people higher in positive, negative or disorganized
  schizotypy show more or less PA inertia?* In other words, does a PA bump fade more slowly or
  more quickly for them?
- **Reading the result:**
  - A **positive** `pa_pa` effect makes the drift less negative, so PA lingers longer (more
    inertia).
  - A **negative** effect means PA returns to baseline faster.
  - `ct_mod_grid()` turns this into carry-over curves at low, mean and high schizotypy.
- **Comparison:** a cross-level interaction in DSEM (level-2 predictor on a random AR slope),
  but without also estimating a free random slope.
- **What it does not do:** schizotypy is not set to predict the *baseline level* of PA,
  because the TIpred effects on `T0m_pa_trend` and `cint_pa_trend` are off. So this model does
  not test whether high-zneg people simply have lower PA overall.
  - To answer that, you would turn on `eff` for `T0m_pa_trend` and/or `cint_pa_trend`.

---

## 7. Estimation settings (`ctFit`)

```r
ctFit(datalong = pa_trend_data, model = pa_trend_mod,
      backend = "julia", priors = TRUE,
      optimcontrol = list(maxiter = 10000), cores = 10)
```

- **`backend = "julia"`:** the fast Kalman-filter engine. It compiles once per model shape per
  session.
- **`priors = TRUE`:** this makes the fit **MAP** (maximum a posteriori) rather than pure ML.
  Weak priors on all parameters keep variances and drifts away from impossible boundaries.
  With this much ESM data the priors barely move the estimates.
- **`maxiter`:** in practice the log-posterior stops changing long before 10,000 iterations.
  1,000 is typically enough; compare the estimates once to confirm.
- **`cores`:** ctsem's tuner picks the fastest thread count itself, usually about 2. Asking for
  more does not make a single fit faster.

---

## 8. Cheat sheet

| Piece | Setting | Plain meaning | Familiar equivalent |
|---|---|---|---|
| LAMBDA | [1 1] | Observed PA = baseline + momentary state | RI-CLPM decomposition |
| DRIFT `pa_pa` | free | How fast a PA bump fades (half-life) | AR(1) phi |
| DRIFT `drift_pa_trend` | free | Curvature of the week-long baseline path | Growth-curve shape |
| DIFFUSION `diff_pa` | free | Size of new "stuff happens" shocks | Innovation variance |
| DIFFUSION trend | 0 | Baseline is smooth, not shocked | - |
| T0VAR `T0var_pa` | free | First-beep deviation from baseline | Initial-state variance |
| T0MEANS `T0m_pa_trend` | free, random | Person's starting baseline | Random intercept |
| CINT `cint_pa_trend` | free, random | Where the person's baseline is heading | Random growth component |
| CINT `pa` | 0 | Momentary state is pulled back to baseline | Person-mean centering |
| MANIFESTMEANS | 0 | Level lives in the trend, not here | Identification constraint |
| MANIFESTVAR `mvar_pa` | free | One-beep noise that doesn't carry over | Measurement error |
| TIpreds | on `pa_pa` only | Does schizotypy change PA inertia? | Cross-level interaction on phi |

---

## 9. How the other 02 models extend this

- **`sit_trend`** adds positive situation (`esm24`, latent `psit`) as a second dynamic process
  with its own trend.
- **`str_trend`** adds stressful situation (`esm25`, latent `ssit`) in the same way.
- The only new ingredient in both is the **cross-effects** in DRIFT: `psit_pa` means the
  situation predicts later PA, and `pa_psit` means PA predicts later situations. These play the
  role of cross-lagged paths in an RI-CLPM, again in continuous time. Schizotypy moderates
  these DRIFT cells as well.
