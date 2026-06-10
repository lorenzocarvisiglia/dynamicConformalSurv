# Conformal calibration
#
# This file computes IPCW bootstrap conformity scores on the survival
# probability scale and returns calibrated survival-probability cutoffs.

conformal_calibration <- function(
  fit,
  surv_train,
  long_train,
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
  verbose = TRUE
) {
  side <- match.arg(side)

  surv_lm <- surv_train[
    surv_train[[time_var]] > landmark,
    ,
    drop = FALSE
  ]

  if (sum(surv_lm[[event_var]] == 1) < 2) {
    stop("Fewer than two post-landmark events are available for calibration.")
  }

  weights <- estimate_ipcw_weights(
    surv_lm = surv_lm,
    time_var = time_var,
    event_var = event_var,
    landmark = landmark
  )

  event_rows <- which(surv_lm[[event_var]] == 1)
  event_prob <- weights[event_rows] / sum(weights[event_rows])

  scores <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    if (verbose && b %% 50 == 0) {
      message("Bootstrap calibration replicate ", b, " / ", B)
    }

    scores[b] <- tryCatch({
      boot_ids <- sample(
        surv_lm[[id_var]],
        size = nrow(surv_lm),
        replace = TRUE
      )

      surv_boot <- surv_train[
        surv_train[[id_var]] %in% boot_ids,
        ,
        drop = FALSE
      ]

      long_boot <- long_train[
        long_train[[id_var]] %in% boot_ids,
        ,
        drop = FALSE
      ]

      fit_b <- fit_working_prc_model(
        surv_train = surv_boot,
        long_train = long_boot,
        landmark = landmark,
        id_var = id_var,
        time_var = time_var,
        event_var = event_var,
        long_time_var = long_time_var,
        baseline_covariates = baseline_covariates,
        longitudinal_markers = longitudinal_markers
      )

      selected_row <- sample(event_rows, size = 1, prob = event_prob)

      surv_selected <- surv_lm[selected_row, , drop = FALSE]
      long_selected <- long_train[
        long_train[[id_var]] %in% surv_selected[[id_var]] &
          long_train[[long_time_var]] <= landmark,
        ,
        drop = FALSE
      ]

      pred <- predict_survival_grid(
        fit = fit_b,
        surv_new = surv_selected,
        long_new = long_selected,
        landmark = landmark,
        id_var = id_var,
        long_time_var = long_time_var,
        time_grid = surv_selected[[time_var]]
      )

      as.numeric(pred$surv[1, 1])
    }, error = function(e) {
      NA_real_
    })
  }

  scores <- scores[is.finite(scores)]

  if (length(scores) < 2) {
    stop("Fewer than two successful bootstrap calibration scores.")
  }

  cutoffs <- compute_survival_cutoffs(
    scores = scores,
    alpha = alpha,
    side = side
  )

  list(
    scores = scores,
    cutoffs = cutoffs,
    n_successful_bootstrap = length(scores)
  )
}


estimate_ipcw_weights <- function(
  surv_lm,
  time_var,
  event_var,
  landmark
) {
  censor_event <- 1 - surv_lm[[event_var]]

  km_fit <- survival::survfit(
    survival::Surv(surv_lm[[time_var]], censor_event) ~ 1
  )

  ghat <- summary(
    km_fit,
    times = surv_lm[[time_var]],
    extend = TRUE
  )$surv

  ghat[is.na(ghat)] <- min(ghat, na.rm = TRUE)
  ghat <- pmax(ghat, 1e-6)

  ifelse(
    surv_lm[[event_var]] == 1,
    1 / ghat,
    0
  )
}


compute_survival_cutoffs <- function(
  scores,
  alpha,
  side = c("two", "lower", "upper")
) {
  side <- match.arg(side)

  if (side == "two") {
    return(list(
      lower_surv = as.numeric(stats::quantile(scores, probs = alpha / 2, names = FALSE, type = 8)),
      upper_surv = as.numeric(stats::quantile(scores, probs = 1 - alpha / 2, names = FALSE, type = 8))
    ))
  }

  if (side == "lower") {
    return(list(
      lower_surv = NA_real_,
      upper_surv = as.numeric(stats::quantile(scores, probs = 1 - alpha, names = FALSE, type = 8))
    ))
  }

  list(
    lower_surv = as.numeric(stats::quantile(scores, probs = alpha, names = FALSE, type = 8)),
    upper_surv = NA_real_
  )
}
