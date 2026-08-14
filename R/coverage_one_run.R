coverage_one_run_prc <- function(
  surv_train,
  long_train,
  surv_valid,
  long_valid,
  landmark,
  id_var,
  time_var,
  event_var,
  long_time_var,
  true_time_var,
  baseline_covariates,
  longitudinal_markers,
  alpha = 0.10,
  B = 500,
  lmm_fixefs = NULL,
  lmm_ranefs = NULL,
  penalty = "ridge",
  standardize = TRUE,
  seed = NULL,
  verbose = TRUE
) {
  if (!true_time_var %in% names(surv_valid)) {
    stop("surv_valid must contain the true event time column: ", true_time_var)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # 1. Fit conformal PRC once, using the existing algorithm.
  fit_two <- dynamic_conformal_pi(
    surv_train = surv_train,
    long_train = long_train,
    surv_new = surv_valid,
    long_new = long_valid,
    landmark = landmark,
    id_var = id_var,
    time_var = time_var,
    event_var = event_var,
    long_time_var = long_time_var,
    baseline_covariates = baseline_covariates,
    longitudinal_markers = longitudinal_markers,
    alpha = alpha,
    B = B,
    side = "two",
    lmm_fixefs = lmm_fixefs,
    lmm_ranefs = lmm_ranefs,
    penalty = penalty,
    standardize = standardize,
    n_cores = 1,
    seed = seed,
    verbose = verbose
  )

  # 2. Rebuild the landmarked datasets.
  train_lmk <- make_landmark_data(
    surv_data = surv_train,
    long_data = long_train,
    landmark = landmark,
    id_var = id_var,
    time_var = time_var,
    long_time_var = long_time_var
  )

  valid_lmk <- make_new_landmark_data(
    surv_data = surv_valid,
    long_data = long_valid,
    landmark = landmark,
    id_var = id_var,
    time_var = time_var,
    long_time_var = long_time_var
  )

  # 3. Predict survival curves for the validation subjects using the same PRC fit.
  time_grid <- make_prediction_grid(
    landmark = landmark,
    surv_train = train_lmk$surv,
    surv_new = valid_lmk$surv,
    time_var = time_var,
    event_var = event_var
  )

  pred_valid <- predict_prc_survival(
    prc_fit = fit_two$prc_fit,
    new_surv = valid_lmk$surv,
    new_long = valid_lmk$long,
    times = time_grid,
    id_var = id_var,
    baseline_covariates = baseline_covariates
  )

  scores <- fit_two$calibration_scores

  # 4. Conformal cutoffs for two-sided, lower one-sided, and upper one-sided PIs.
  conformal_cut_two <- compute_survival_cutoffs(
    scores = scores,
    alpha = alpha,
    side = "two"
  )

  conformal_cut_lower <- compute_survival_cutoffs(
    scores = scores,
    alpha = alpha,
    side = "lower"
  )

  conformal_cut_upper <- compute_survival_cutoffs(
    scores = scores,
    alpha = alpha,
    side = "upper"
  )

  # 5. Naive cutoffs from direct survival-function inversion.
  naive_cut_two <- list(
    lower_surv = alpha / 2,
    upper_surv = 1 - alpha / 2
  )

  naive_cut_lower <- list(
    lower_surv = NA_real_,
    upper_surv = 1 - alpha
  )

  naive_cut_upper <- list(
    lower_surv = alpha,
    upper_surv = NA_real_
  )

  # 6. Build conformal intervals.
  conformal_two <- invert_prediction_intervals(
    pred = pred_valid,
    landmark = landmark,
    alpha = alpha,
    side = "two",
    cutoffs = conformal_cut_two
  )

  conformal_lower <- invert_prediction_intervals(
    pred = pred_valid,
    landmark = landmark,
    alpha = alpha,
    side = "lower",
    cutoffs = conformal_cut_lower
  )

  conformal_upper <- invert_prediction_intervals(
    pred = pred_valid,
    landmark = landmark,
    alpha = alpha,
    side = "upper",
    cutoffs = conformal_cut_upper
  )

  # 7. Build naive intervals.
  naive_two <- invert_prediction_intervals(
    pred = pred_valid,
    landmark = landmark,
    alpha = alpha,
    side = "two",
    cutoffs = naive_cut_two
  )

  naive_lower <- invert_prediction_intervals(
    pred = pred_valid,
    landmark = landmark,
    alpha = alpha,
    side = "lower",
    cutoffs = naive_cut_lower
  )

  naive_upper <- invert_prediction_intervals(
    pred = pred_valid,
    landmark = landmark,
    alpha = alpha,
    side = "upper",
    cutoffs = naive_cut_upper
  )

  # 8. Largest observed post-landmark event time in the training risk set.
  eta_l <- max(
    train_lmk$surv[[time_var]][train_lmk$surv[[event_var]] == 1],
    na.rm = TRUE
  )

  # 9. Function to compute coverage metrics.
  compute_metrics <- function(two, lower, upper, method_name) {
    ids <- two$id

    truth <- valid_lmk$surv[
      match(ids, valid_lmk$surv[[id_var]]),
      true_time_var
    ]

    upper_trunc <- pmin(two$upper, eta_l)
    trunc_length <- upper_trunc - two$lower
    trunc_length[!is.finite(trunc_length)] <- NA_real_
    trunc_length[trunc_length < 0] <- NA_real_

    data.frame(
      method = method_name,
      landmark = landmark,
      n_valid = length(truth),
      eta_landmark = eta_l,
      right_cov = mean(truth >= lower$lower, na.rm = TRUE),
      left_cov = mean(truth <= upper$upper, na.rm = TRUE),
      total_cov = mean(truth >= two$lower & truth <= two$upper, na.rm = TRUE),
      trunc_total_cov = mean(truth >= two$lower & truth <= upper_trunc, na.rm = TRUE),
      avg_trunc_length = mean(trunc_length, na.rm = TRUE),
      inf_upper = mean(is.infinite(two$upper), na.rm = TRUE)
    )
  }

  metrics <- rbind(
    compute_metrics(
      two = naive_two,
      lower = naive_lower,
      upper = naive_upper,
      method_name = "Naive PRC"
    ),
    compute_metrics(
      two = conformal_two,
      lower = conformal_lower,
      upper = conformal_upper,
      method_name = "Conformal PRC"
    )
  )

  list(
    metrics = metrics,
    intervals = list(
      naive_two = naive_two,
      naive_lower = naive_lower,
      naive_upper = naive_upper,
      conformal_two = conformal_two,
      conformal_lower = conformal_lower,
      conformal_upper = conformal_upper
    ),
    cutoffs = list(
      naive_two = naive_cut_two,
      naive_lower = naive_cut_lower,
      naive_upper = naive_cut_upper,
      conformal_two = conformal_cut_two,
      conformal_lower = conformal_cut_lower,
      conformal_upper = conformal_cut_upper
    ),
    calibration_scores = scores,
    m_eff = fit_two$m_eff,
    prc_fit = fit_two$prc_fit
  )
}
