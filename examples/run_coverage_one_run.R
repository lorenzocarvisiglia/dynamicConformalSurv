# Example: one-run coverage computation
#
# This script shows how to compute coverage metrics for a single
# train/validation split, without running an outer Monte Carlo loop.

source("R/load_dynamic_conformal.R")

one_run <- coverage_one_run_prc(
  surv_train = surv_train,
  long_train = long_train,
  surv_valid = surv_valid,
  long_valid = long_valid,
  landmark = 2,
  id_var = "id",
  time_var = "time",
  event_var = "event",
  long_time_var = "t.from.base",
  true_time_var = "T_true",
  baseline_covariates = c("baseline.age"),
  longitudinal_markers = c("y1", "y2", "y3"),
  alpha = 0.10,
  B = 500,
  lmm_fixefs = ~ t.from.base,
  lmm_ranefs = ~ t.from.base | id,
  seed = 123,
  verbose = TRUE
)

one_run$metrics
