# Main function for dynamic conformal prediction intervals
# This script implements the algorithm described in the paper.
#
# Users are expected to provide already-cleaned survival and longitudinal data.

dynamic_conformal_pi <- function(
  surv_train,
  long_train,
  surv_new,
  long_new,
  landmark,
  id_var,
  time_var,
  event_var,
  long_time_var,
  baseline_covariates,
  longitudinal_markers,
  alpha = 0.10,
  B = 500,
  side = c("two", "lower", "upper"),
  seed = NULL,
  verbose = TRUE
) {
  side <- match.arg(side)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (verbose) {
    message("Step 1: fitting the working PRC model")
  }

  fit <- fit_working_prc_model(
    surv_train = surv_train,
    long_train = long_train,
    landmark = landmark,
    id_var = id_var,
    time_var = time_var,
    event_var = event_var,
    long_time_var = long_time_var,
    baseline_covariates = baseline_covariates,
    longitudinal_markers = longitudinal_markers
  )

  if (verbose) {
    message("Step 2: computing conformal calibration scores")
  }

  calibration <- conformal_calibration(
    fit = fit,
    surv_train = surv_train,
    long_train = long_train,
    landmark = landmark,
    id_var = id_var,
    time_var = time_var,
    event_var = event_var,
    long_time_var = long_time_var,
    baseline_covariates = baseline_covariates,
    longitudinal_markers = longitudinal_markers,
    alpha = alpha,
    B = B,
    side = side,
    verbose = verbose
  )

  if (verbose) {
    message("Step 3: inverting calibrated survival thresholds")
  }

  intervals <- invert_survival_intervals(
    fit = fit,
    surv_new = surv_new,
    long_new = long_new,
    landmark = landmark,
    id_var = id_var,
    time_var = time_var,
    event_var = event_var,
    long_time_var = long_time_var,
    alpha = alpha,
    side = side,
    cutoffs = calibration$cutoffs
  )

  out <- list(
    intervals = intervals,
    cutoffs = calibration$cutoffs,
    calibration_scores = calibration$scores,
    working_model = fit,
    call = match.call()
  )

  class(out) <- "dynamic_conformal_pi"
  out
}
