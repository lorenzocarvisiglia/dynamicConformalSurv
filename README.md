# dynamicConformalSurv

Research R code for dynamic prediction intervals for survival times with longitudinal covariates.

The implementation accompanies the manuscript and PhD thesis chapter:

**Dynamic prediction intervals for survival times**

by Lorenzo Carvisiglia, Saverio Ranciati, and Mirko Signorelli.

## Scope

The repository provides a reusable, subject-level implementation of the proposed dynamic conformal method. It is script-based research code rather than an R package.

The method combines:

1. a working dynamic survival model based on Penalized Regression Calibration (PRC) through the `pencal` package;
2. subject-level bootstrap refitting of the full working model;
3. inverse probability of censoring weighting (IPCW) to sample observed post-landmark failures;
4. calibration on the survival-probability scale;
5. inversion of calibrated survival thresholds to obtain one-sided or two-sided prediction intervals.

Prediction is defined for subjects who are observed and event-free at a landmark time.

## Repository structure

```text
dynamicConformalSurv/
├── R/
│   ├── dynamic_conformal_pi.R
│   ├── coverage_one_run.R
│   └── load_dynamic_conformal.R
├── examples/
│   ├── example_toy_data.R
│   └── run_coverage_one_run.R
├── tests/
│   └── test_core_helpers.R
├── data/
│   └── README.md
├── .gitignore
├── LICENSE
└── README.md
```

The public implementation is intentionally independent of personal paths, cluster schedulers, and private datasets.

## Dependencies

The core method requires:

- `survival`
- `pencal`

A recent version of `pencal` is recommended.

## Input data

### Training data

`surv_train` must contain one row per subject, including:

- subject identifier;
- observed event or censoring time;
- event indicator;
- baseline covariates used in the working survival model.

`long_train` must contain one row per subject-visit, including:

- subject identifier;
- visit time;
- longitudinal markers.

The function constructs the landmark risk set internally and uses only longitudinal measurements observed up to the landmark.

### Subjects to be predicted

For genuine prediction, `surv_new` needs only the subject identifier and baseline covariates. Future event or censoring times are not required.

`long_new` should contain only information available at or before the landmark. The prediction time grid is determined exclusively from post-landmark event times in the training data, so future follow-up information from new subjects is never used to determine interval endpoints.

If an observed follow-up-time column is supplied in `surv_new`, it is used only to restrict the data to subjects known to be observed and event-free at the landmark.

## Basic usage

```r
source("R/load_dynamic_conformal.R")

fit <- dynamic_conformal_pi(
  surv_train = surv_train,
  long_train = long_train,
  surv_new = surv_new,
  long_new = long_new,
  landmark = 2,
  id_var = "id",
  time_var = "time",
  event_var = "event",
  long_time_var = "t.from.base",
  baseline_covariates = c("baseline.age"),
  longitudinal_markers = c("y1", "y2", "y3"),
  alpha = 0.10,
  B = 500,
  side = "two",
  lmm_fixefs = ~ t.from.base,
  lmm_ranefs = ~ t.from.base | id,
  seed = 123
)

fit$intervals
fit$cutoffs
fit$m_eff
```

If `lmm_fixefs` and `lmm_ranefs` are omitted, the same random-intercept/random-slope specification in `t.from.base` is used by default.

## Output

`dynamic_conformal_pi()` returns:

- `intervals`: prediction intervals on the original time scale;
- `cutoffs`: calibrated survival-probability thresholds;
- `calibration_scores`: successful bootstrap conformity scores;
- `m_eff`: number of successful bootstrap calibration replicates;
- `prc_fit`: fitted PRC objects;
- `landmark`, `alpha`, and `side`.

For two-sided intervals, the lower and upper endpoints are obtained by first crossing of the corresponding calibrated survival thresholds. If the fitted survival curve does not cross the threshold within the estimable training-event grid, the corresponding endpoint is returned as `Inf`.

## One-run simulation evaluation

`R/coverage_one_run.R` contains `coverage_one_run_prc()`, a helper for evaluating one simulated train/validation split. Unlike genuine prediction, this evaluation helper requires the true event time in the validation data so that empirical coverage can be computed.

See:

```text
examples/run_coverage_one_run.R
```

## Implementation details

The calibration distribution is formed among observed post-landmark failures with probability proportional to `1 / Ghat(T*)`, where `Ghat` is the Kaplan-Meier estimate of the censoring survival function in the landmark risk set.

Each bootstrap replicate resamples subjects with replacement, preserves bootstrap multiplicities by assigning new bootstrap subject identifiers, refits the complete PRC working model, draws one IPCW-weighted observed failure from the original landmark sample, and evaluates the refitted survival function at that subject's observed failure time.

The numerical implementation uses a small lower bound for estimated censoring survival probabilities to avoid division by zero. This is a numerical safeguard and does not replace the censoring-positivity assumption used in the theoretical results.

## Data and reproducibility

No ADNI participant-level data are distributed in this repository. ADNI data are subject to the access conditions of the Alzheimer's Disease Neuroimaging Initiative.

Large generated simulation datasets are also not stored in Git. The reusable method code does not depend on those datasets.

The current public release contains the core method and generic evaluation examples. Manuscript-specific archived cluster job scripts are intentionally not part of the reusable interface because they contain environment-specific execution details.

## Tests

Core helper functions can be checked without fitting a PRC model:

```bash
Rscript tests/test_core_helpers.R
```

## Citation

If you use this code, please cite the accompanying manuscript:

Carvisiglia, L., Ranciati, S., and Signorelli, M. (2026).  
*Dynamic prediction intervals for survival times*. arXiv:2609.10409.
