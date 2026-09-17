# ctsem report

Backend: julia. Latents: 4. Manifests: 2. Covariates: 3.
Subjects: 203. Rows: 12723. Uncertainty draws: 1000.
Log likelihood: -7573.6. Free parameters: 31. AIC: 15209.2.
Written 2026-09-16 21:45 by ctReport().

## Files

### 01-summary.txt

Parameter estimates with uncertainty intervals.

*An interval excluding zero is the usual claim of an effect; check 03 first.*

### 02-parameter-matrices.txt

System matrices at each requested quantile.

*DRIFT diagonals negative means stable; off-diagonals are cross-effects per unit time.*

### 03-identification.txt

Whether the optimizer arrived, and whether the data pin the parameters.

*Read this first. A weak direction or a huge condition number invalidates the intervals in 01.*

### 04-loglik-profile.txt, 04-loglik-profile.pdf

Log-probability slice through each parameter against the quadratic its curvature implies.

*Off-peak: not converged. Ridged: pinned alone, free jointly. Non-quadratic or asymmetric: the interval is the wrong shape.*

20 of 31 parameters flagged: ridged, non-quadratic, asymmetric.

### 05-discrete-parameters.pdf

Regression coefficients between latents as a function of time interval.

*This is what the continuous-time DRIFT means for a given gap between measurements.*

### 06-network.pdf, 06-network-edges.csv

The temporal and contemporaneous networks at one time interval.

*An arrow is a directed effect over that interval; an undirected edge is shared innovation.*

### 08-residuals.pdf, 08-residuals.txt

Standardised prior residuals, and their autocorrelation.

*Autocorrelation outside the interval means structure the model has not taken up.*

### 09-covariance-check.csv, 09-covariance-check.pdf

Lagged covariances, observed against model-generated.

*Empirical outside the generated interval is misfit at that lag; a whole row or column off is a variable the model gets wrong.*

2 of 16 observed covariances fall outside the generated interval.

### 10-posterior-predictive.pdf

Distribution of generated data against observed.

*Observed away from the generated mass means the model cannot produce data like yours.*

### 11-covariate-effects.pdf

Predicted parameter values across the range of each covariate.

*A band that excludes a flat line is a covariate effect on that parameter.*

## Not produced

- **predictions** -- failed: second argument of / cannot be a "difftime" object

---

Each file comes from one ctsem function; the help page for that function
has the arguments and the detail. Regenerate with `ctReport(fit)`.
