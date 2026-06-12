# dynamicConformalSurv

This repository contains R scripts implementing dynamic conformal prediction intervals for survival times with longitudinal covariates.

The code accompanies the manuscript:

**On the computation of prediction intervals for survival times in dynamic prediction problems**

## Overview

The method constructs prediction intervals for the event time of a subject who is event-free at a landmark time.

The algorithm combines:

1. a working dynamic survival model based on Penalized Regression Calibration;
2. bootstrap conformal calibration on the survival-probability scale;
3. inverse probability of censoring weighting to account for right-censored outcomes;
4. inversion of calibrated survival thresholds to obtain prediction interval endpoints.

## Repository structure

```text
dynamicConformalSurv/
├── R/
│   ├── load_dynamic_conformal.R
│   └── dynamic_conformal_pi.R
├── examples/
│   └── example_toy_data.R
└── data/
    └── README.md
```

## Required input data

Users should provide already-cleaned data.

The survival data should contain one row per subject, with at least:

- subject identifier;
- observed survival or censoring time;
- event indicator;
- baseline covariates.

The longitudinal data should contain one row per subject-visit, with at least:

- subject identifier;
- visit time;
- longitudinal markers.

## Dependencies

The scripts require:

- `survival`
- `pencal`

Please use a recent version of `pencal`. If survival prediction fails at time points beyond the maximum observed training time, update `pencal`.

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
  lmm_fixefs = ~ age,
  lmm_ranefs = ~ age | id,
  seed = 123
)

fit$intervals
fit$cutoffs
fit$m_eff
```

## Output

The main function returns a list with:

- `intervals`: prediction intervals on the original time scale;
- `cutoffs`: calibrated survival-probability thresholds;
- `calibration_scores`: bootstrap conformal scores;
- `m_eff`: number of successful bootstrap calibration replicates;
- `prc_fit`: fitted PRC model objects.

The interval table contains:

```text
id
landmark
alpha
side
lower
upper
```

For two-sided intervals, `lower` and `upper` are obtained by inverting the calibrated survival thresholds. If the fitted survival curve does not cross the required threshold within the prediction grid, the corresponding upper endpoint is returned as `Inf`.

## Example

A minimal toy example is provided in:

```text
examples/example_toy_data.R
```

The toy example is intended to show the required data structure and function call.

## Citation

If you use this code, please cite the accompanying manuscript:

Carvisiglia, L., Ranciati, S., and Signorelli, M.  
*On the computation of prediction intervals for survival times in dynamic prediction problems.*
