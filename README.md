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

The implementation is intentionally script-based. It is not an R package.

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
