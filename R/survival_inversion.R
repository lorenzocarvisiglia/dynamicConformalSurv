# Survival-function inversion utilities
#
# The function below converts calibrated survival-probability thresholds
# into prediction interval endpoints.

invert_survival_intervals <- function(
  fit,
  surv_new,
  long_new,
  landmark,
  id_var,
  time_var,
  event_var,
  long_time_var,
  alpha = 0.10,
  side = c("two", "lower", "upper"),
  cutoffs
) {
  side <- match.arg(side)

  pred <- predict_survival_grid(
    fit = fit,
    surv_new = surv_new,
    long_new = long_new,
    landmark = landmark,
    id_var = id_var,
    long_time_var = long_time_var
  )

  ids <- pred$ids
  time_grid <- pred$time_grid
  surv_mat <- pred$surv

  if (side == "two") {
    lower <- invert_one_threshold(
      surv_mat = surv_mat,
      time_grid = time_grid,
      threshold = cutoffs$upper_surv
    )

    upper <- invert_one_threshold(
      surv_mat = surv_mat,
      time_grid = time_grid,
      threshold = cutoffs$lower_surv
    )
  }

  if (side == "lower") {
    lower <- invert_one_threshold(
      surv_mat = surv_mat,
      time_grid = time_grid,
      threshold = cutoffs$upper_surv
    )

    upper <- rep(Inf, length(lower))
  }

  if (side == "upper") {
    lower <- rep(landmark, length(ids))

    upper <- invert_one_threshold(
      surv_mat = surv_mat,
      time_grid = time_grid,
      threshold = cutoffs$lower_surv
    )
  }

  data.frame(
    id = ids,
    landmark = landmark,
    side = side,
    alpha = alpha,
    lower = lower,
    upper = upper
  )
}


invert_one_threshold <- function(
  surv_mat,
  time_grid,
  threshold
) {
  out <- rep(Inf, nrow(surv_mat))

  for (i in seq_len(nrow(surv_mat))) {
    idx <- which(surv_mat[i, ] <= threshold)

    if (length(idx) > 0) {
      out[i] <- time_grid[min(idx)]
    }
  }

  out
}
