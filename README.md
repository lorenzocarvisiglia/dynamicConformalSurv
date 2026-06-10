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
│   ├── dynamic_conformal_pi.R
│   ├── working_prc_model.R
│   ├── conformal_calibration.R
│   └── survival_inversion.R
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

The scripts do not include data-cleaning functions. Dataset-specific preprocessing should be done by the user before calling the main function.

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
  long_time_var = "time_fup",
  baseline_covariates = c("age", "sex"),
  longitudinal_markers = c("marker1", "marker2"),
  alpha = 0.10,
  B = 500,
  side = "two",
  seed = 123
)

fit$intervals
fit$cutoffs
```

## Example

A minimal toy example is provided in:

```text
examples/example_toy_data.R
```

The toy example is intended to show the required data structure and function call.

## Data availability

This repository does not include ADNI data, simulation outputs, cluster scripts, or large `.RData` files.

## Dependencies

The scripts use standard R packages, including:

- `survival`
- `data.table`

The working dynamic survival model is based on Penalized Regression Calibration. Users should install the required PRC implementation before running the full method.

## Citation

If you use this code, please cite the accompanying manuscript:

Carvisiglia, L., Ranciati, S., and Signorelli, M.  
*On the computation of prediction intervals for survival times in dynamic prediction problems.*
