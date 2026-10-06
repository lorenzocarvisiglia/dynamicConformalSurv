# Sequential helper for one simulation cell.

source("R/load_dynamic_conformal.R")
source("reproducibility/simulations/simulate_three_predictor.R")
source("reproducibility/simulations/simulate_twenty_predictor.R")

run_simulation_cell <- function(
  setting = c("three", "twenty"),
  scenario = c("weakNPH", "strongNPH"),
  n_train = 300L,
  n_validation = 10000L,
  target_censoring = 0.20,
  landmark = 2,
  alpha = 0.10,
  B = 500L,
  seed = 1L,
  verbose = TRUE
) {
  setting <- match.arg(setting)
  scenario <- match.arg(scenario)

  if (setting == "twenty" && scenario != "strongNPH") {
    stop("The 20-predictor design is defined only for strongNPH.")
  }

  if (setting == "three") {
    cmax <- calibrate_three_predictor_cmax(
      scenario = scenario,
      target_censoring = target_censoring,
      seed = 7000L + seed
    )
    train <- simulate_three_predictor_dataset(
      n_train, scenario, target_censoring, cmax, seed
    )
    valid <- simulate_three_predictor_dataset(
      n_validation, scenario, target_censoring, cmax, 1000000L + seed
    )
    markers <- c("y1", "y2", "y3")
  } else {
    cmax <- calibrate_twenty_predictor_cmax(
      target_censoring = target_censoring,
      seed = 7000L + seed
    )
    train <- simulate_twenty_predictor_dataset(
      n_train, target_censoring, cmax, seed
    )
    valid <- simulate_twenty_predictor_dataset(
      n_validation, target_censoring, cmax, 1000000L + seed
    )
    markers <- paste0("y", 1:20)
  }

  coverage_one_run_prc(
    surv_train = train$surv,
    long_train = train$long,
    surv_valid = valid$surv,
    long_valid = valid$long,
    landmark = landmark,
    id_var = "id",
    time_var = "time",
    event_var = "event",
    long_time_var = "t.from.base",
    true_time_var = "true.time",
    baseline_covariates = "baseline.age",
    longitudinal_markers = markers,
    alpha = alpha,
    B = B,
    lmm_fixefs = ~ age,
    lmm_ranefs = ~ age | id,
    seed = 900000L + seed,
    verbose = verbose
  )
}
