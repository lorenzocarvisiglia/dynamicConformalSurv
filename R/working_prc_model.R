# Working PRC model
#
# This file contains the working dynamic survival model used by the
# conformal procedure. In the paper this is Penalized Regression
# Calibration: longitudinal histories are summarized by LMM random effects,
# then a post-landmark Cox model is fitted.

fit_working_prc_model <- function(
  surv_train,
  long_train,
  landmark,
  id_var,
  time_var,
  event_var,
  long_time_var,
  baseline_covariates,
  longitudinal_markers
) {
  risk_ids <- surv_train[[id_var]][
    surv_train[[time_var]] > landmark
  ]

  surv_lm <- surv_train[
    surv_train[[id_var]] %in% risk_ids,
    ,
    drop = FALSE
  ]

  long_lm <- long_train[
    long_train[[id_var]] %in% risk_ids &
      long_train[[long_time_var]] <= landmark,
    ,
    drop = FALSE
  ]

  fit <- list(
    surv_train = surv_lm,
    long_train = long_lm,
    landmark = landmark,
    id_var = id_var,
    time_var = time_var,
    event_var = event_var,
    long_time_var = long_time_var,
    baseline_covariates = baseline_covariates,
    longitudinal_markers = longitudinal_markers,
    prc_fit = NULL
  )

  class(fit) <- "working_prc_model"
  fit
}


predict_survival_grid <- function(
  fit,
  surv_new,
  long_new,
  landmark,
  id_var,
  long_time_var,
  time_grid = NULL
) {
  ids <- surv_new[[id_var]]

  if (is.null(time_grid)) {
    event_times <- fit$surv_train[[fit$time_var]][
      fit$surv_train[[fit$event_var]] == 1 &
        fit$surv_train[[fit$time_var]] > landmark
    ]

    time_grid <- sort(unique(event_times))
  }

  if (length(time_grid) == 0) {
    stop("No post-landmark event times available for survival inversion.")
  }

  # Placeholder survival matrix.
  # This will be replaced by PRC-based survival prediction.
  surv_mat <- matrix(
    NA_real_,
    nrow = length(ids),
    ncol = length(time_grid)
  )

  colnames(surv_mat) <- paste0("t", seq_along(time_grid))
  rownames(surv_mat) <- ids

  list(
    ids = ids,
    time_grid = time_grid,
    surv = surv_mat
  )
}
