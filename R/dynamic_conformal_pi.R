# Dynamic conformal prediction intervals for survival times
#
# Script-based implementation based on Penalized Regression Calibration
# through the pencal package.
#
# Main function:
# dynamic_conformal_pi()
#
# Required data:
# - surv_train: one row per subject
# - long_train: one row per subject-visit
# - surv_new: one row per new subject
# - long_new: one row per new subject-visit
#
# The function returns prediction intervals on the original time scale.

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
  lmm_fixefs = NULL,
  lmm_ranefs = NULL,
  penalty = "ridge",
  standardize = TRUE,
  n_cores = 1,
  seed = NULL,
  verbose = TRUE
) {
  side <- match.arg(side)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  check_packages()
  check_dynamic_inputs(
    surv_data = surv_train,
    long_data = long_train,
    id_var = id_var,
    time_var = time_var,
    event_var = event_var,
    long_time_var = long_time_var,
    baseline_covariates = baseline_covariates,
    longitudinal_markers = longitudinal_markers
  )

  check_dynamic_inputs(
    surv_data = surv_new,
    long_data = long_new,
    id_var = id_var,
    time_var = time_var,
    event_var = event_var,
    long_time_var = long_time_var,
    baseline_covariates = baseline_covariates,
    longitudinal_markers = longitudinal_markers,
    require_event = FALSE
  )

  if (is.null(lmm_fixefs)) {
    lmm_fixefs <- stats::as.formula("~ t.from.base")
  }

  if (is.null(lmm_ranefs)) {
    lmm_ranefs <- stats::as.formula("~ t.from.base | id")
  }

  if (verbose) {
    message("Preparing landmark training data")
  }

  train_lmk <- make_landmark_data(
    surv_data = surv_train,
    long_data = long_train,
    landmark = landmark,
    id_var = id_var,
    time_var = time_var,
    long_time_var = long_time_var
  )

  surv_lmk <- train_lmk$surv
  long_lmk <- train_lmk$long

  if (nrow(surv_lmk) < 5) {
    stop("Too few subjects are at risk at the landmark.")
  }

  if (sum(surv_lmk[[event_var]] == 1, na.rm = TRUE) < 2) {
    stop("Fewer than two post-landmark events are available.")
  }

  if (verbose) {
    message("Fitting original PRC model")
  }

  original_fit <- fit_prc_landmark(
    surv_data = surv_lmk,
    long_data = long_lmk,
    id_var = id_var,
    time_var = time_var,
    event_var = event_var,
    long_time_var = long_time_var,
    baseline_covariates = baseline_covariates,
    longitudinal_markers = longitudinal_markers,
    lmm_fixefs = lmm_fixefs,
    lmm_ranefs = lmm_ranefs,
    penalty = penalty,
    standardize = standardize,
    n_cores = n_cores,
    verbose = verbose
  )

  if (verbose) {
    message("Computing IPCW failure distribution")
  }

  failure_dist <- make_ipcw_failure_distribution(
    surv_data = surv_lmk,
    time_var = time_var,
    event_var = event_var
  )

  if (verbose) {
    message("Running bootstrap conformal calibration")
  }

  scores <- bootstrap_conformal_scores(
    surv_data = surv_lmk,
    long_data = long_lmk,
    failure_dist = failure_dist,
    B = B,
    id_var = id_var,
    time_var = time_var,
    event_var = event_var,
    long_time_var = long_time_var,
    baseline_covariates = baseline_covariates,
    longitudinal_markers = longitudinal_markers,
    lmm_fixefs = lmm_fixefs,
    lmm_ranefs = lmm_ranefs,
    penalty = penalty,
    standardize = standardize,
    n_cores = n_cores,
    verbose = verbose
  )

  scores <- scores[is.finite(scores)]

  if (length(scores) < 2) {
    stop("Fewer than two successful bootstrap calibration scores.")
  }

  cutoffs <- compute_survival_cutoffs(
    scores = scores,
    alpha = alpha,
    side = side
  )

  if (verbose) {
    message("Preparing new landmark subjects")
  }

  new_lmk <- make_new_landmark_data(
    surv_data = surv_new,
    long_data = long_new,
    landmark = landmark,
    id_var = id_var,
    time_var = time_var,
    long_time_var = long_time_var
  )

  if (nrow(new_lmk$surv) == 0) {
    stop("No new subjects are available at the landmark.")
  }

  if (verbose) {
    message("Predicting survival curves for new subjects")
  }

  time_grid <- make_prediction_grid(
    landmark = landmark,
    surv_train = surv_lmk,
    surv_new = new_lmk$surv,
    time_var = time_var,
    event_var = event_var
  )

  pred_new <- predict_prc_survival(
    prc_fit = original_fit,
    new_surv = new_lmk$surv,
    new_long = new_lmk$long,
    times = time_grid,
    id_var = id_var,
    baseline_covariates = baseline_covariates
  )

  intervals <- invert_prediction_intervals(
    pred = pred_new,
    landmark = landmark,
    alpha = alpha,
    side = side,
    cutoffs = cutoffs
  )

  out <- list(
    intervals = intervals,
    cutoffs = cutoffs,
    calibration_scores = scores,
    m_eff = length(scores),
    prc_fit = original_fit,
    landmark = landmark,
    alpha = alpha,
    side = side,
    call = match.call()
  )

  class(out) <- "dynamic_conformal_pi"
  out
}


fit_prc_landmark <- function(
  surv_data,
  long_data,
  id_var,
  time_var,
  event_var,
  long_time_var,
  baseline_covariates,
  longitudinal_markers,
  lmm_fixefs,
  lmm_ranefs,
  penalty = "ridge",
  standardize = TRUE,
  n_cores = 1,
  verbose = FALSE
) {
  surv_p <- prepare_surv_for_pencal(
    surv_data = surv_data,
    id_var = id_var,
    time_var = time_var,
    event_var = event_var
  )

  long_p <- prepare_long_for_pencal(
    long_data = long_data,
    id_var = id_var,
    long_time_var = long_time_var
  )

  step1 <- pencal::fit_lmms(
    y.names = longitudinal_markers,
    fixefs = lmm_fixefs,
    ranefs = lmm_ranefs,
    t.from.base = t.from.base,
    long.data = long_p,
    surv.data = surv_p,
    n.boots = 0,
    n.cores = n_cores,
    verbose = verbose
  )

  step2 <- pencal::summarize_lmms(
    step1,
    n.cores = n_cores,
    verbose = verbose
  )

  baseline_formula <- stats::as.formula(
    paste("~", paste(baseline_covariates, collapse = " + "))
  )

  step3 <- pencal::fit_prclmm(
    object = step2,
    surv.data = surv_p,
    baseline.covs = baseline_formula,
    penalty = penalty,
    standardize = standardize,
    n.cores = n_cores,
    verbose = verbose
  )

  list(
    step1 = step1,
    step2 = step2,
    step3 = step3,
    id_var = id_var,
    time_var = time_var,
    event_var = event_var,
    long_time_var = long_time_var,
    baseline_covariates = baseline_covariates,
    longitudinal_markers = longitudinal_markers,
    lmm_fixefs = lmm_fixefs,
    lmm_ranefs = lmm_ranefs,
    penalty = penalty,
    standardize = standardize
  )
}


predict_prc_survival <- function(
  prc_fit,
  new_surv,
  new_long,
  times,
  id_var,
  baseline_covariates
) {
  new_base <- as.data.frame(
    new_surv[, c(id_var, baseline_covariates), drop = FALSE]
  )

  new_long <- as.data.frame(new_long)

  new_base <- prepare_new_base_for_pencal(
    new_base = new_base,
    id_var = id_var
  )

  new_long <- prepare_long_for_pencal(
    long_data = new_long,
    id_var = id_var,
    long_time_var = prc_fit$long_time_var
  )

  pred <- pencal::survpred_prclmm(
    step1 = prc_fit$step1,
    step2 = prc_fit$step2,
    step3 = prc_fit$step3,
    new.basecovs = new_base,
    new.longdata = new_long,
    times = times
  )

  pm <- as.data.frame(
    pred$predicted_survival,
    check.names = FALSE
  )

  surv_info <- extract_survival_matrix(pm)

  list(
    ids = pm$id,
    times = surv_info$times,
    surv = surv_info$surv
  )
}


bootstrap_conformal_scores <- function(
  surv_data,
  long_data,
  failure_dist,
  B,
  id_var,
  time_var,
  event_var,
  long_time_var,
  baseline_covariates,
  longitudinal_markers,
  lmm_fixefs,
  lmm_ranefs,
  penalty = "ridge",
  standardize = TRUE,
  n_cores = 1,
  verbose = TRUE
) {
  ids <- unique(surv_data[[id_var]])
  scores <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    if (verbose && (b == 1 || b %% 50 == 0 || b == B)) {
      message("Bootstrap replicate ", b, " / ", B)
    }

    scores[b] <- tryCatch({
      draw_ids <- sample(
        ids,
        size = length(ids),
        replace = TRUE
      )

      boot_data <- make_bootstrap_data(
        surv_data = surv_data,
        long_data = long_data,
        draw_ids = draw_ids,
        id_var = id_var
      )

      boot_fit <- fit_prc_landmark(
        surv_data = boot_data$surv,
        long_data = boot_data$long,
        id_var = id_var,
        time_var = time_var,
        event_var = event_var,
        long_time_var = long_time_var,
        baseline_covariates = baseline_covariates,
        longitudinal_markers = longitudinal_markers,
        lmm_fixefs = lmm_fixefs,
        lmm_ranefs = lmm_ranefs,
        penalty = penalty,
        standardize = standardize,
        n_cores = n_cores,
        verbose = FALSE
      )

      selected <- sample_failure_subject(
        failure_dist = failure_dist,
        long_data = long_data,
        id_var = id_var,
        time_var = time_var,
        long_time_var = long_time_var
      )

      pred <- predict_prc_survival(
        prc_fit = boot_fit,
        new_surv = selected$surv,
        new_long = selected$long,
        times = selected$time,
        id_var = id_var,
        baseline_covariates = baseline_covariates
      )

      as.numeric(pred$surv[1, 1])
    }, error = function(e) {
      if (verbose) {
        message("Bootstrap replicate ", b, " failed: ", conditionMessage(e))
      }
      NA_real_
    })
  }

  scores
}


make_landmark_data <- function(
  surv_data,
  long_data,
  landmark,
  id_var,
  time_var,
  long_time_var
) {
  surv_data <- as.data.frame(surv_data)
  long_data <- as.data.frame(long_data)

  risk_ids <- surv_data[[id_var]][surv_data[[time_var]] > landmark]

  surv_lmk <- surv_data[
    surv_data[[id_var]] %in% risk_ids,
    ,
    drop = FALSE
  ]

  long_lmk <- long_data[
    long_data[[id_var]] %in% risk_ids &
      long_data[[long_time_var]] <= landmark,
    ,
    drop = FALSE
  ]

  list(
    surv = surv_lmk,
    long = long_lmk
  )
}


make_new_landmark_data <- function(
  surv_data,
  long_data,
  landmark,
  id_var,
  time_var,
  long_time_var
) {
  surv_data <- as.data.frame(surv_data)
  long_data <- as.data.frame(long_data)

  if (time_var %in% names(surv_data)) {
    risk_ids <- surv_data[[id_var]][surv_data[[time_var]] > landmark]
  } else {
    risk_ids <- surv_data[[id_var]]
  }

  surv_lmk <- surv_data[
    surv_data[[id_var]] %in% risk_ids,
    ,
    drop = FALSE
  ]

  long_lmk <- long_data[
    long_data[[id_var]] %in% risk_ids &
      long_data[[long_time_var]] <= landmark,
    ,
    drop = FALSE
  ]

  list(
    surv = surv_lmk,
    long = long_lmk
  )
}


make_bootstrap_data <- function(
  surv_data,
  long_data,
  draw_ids,
  id_var
) {
  surv_list <- vector("list", length(draw_ids))
  long_list <- vector("list", length(draw_ids))

  for (k in seq_along(draw_ids)) {
    old_id <- draw_ids[k]
    new_id <- k

    srow <- surv_data[
      surv_data[[id_var]] == old_id,
      ,
      drop = FALSE
    ]

    srow[[id_var]] <- new_id
    surv_list[[k]] <- srow

    lrows <- long_data[
      long_data[[id_var]] == old_id,
      ,
      drop = FALSE
    ]

    lrows[[id_var]] <- new_id
    long_list[[k]] <- lrows
  }

  list(
    surv = as.data.frame(do.call(rbind, surv_list)),
    long = as.data.frame(do.call(rbind, long_list))
  )
}


make_ipcw_failure_distribution <- function(
  surv_data,
  time_var,
  event_var
) {
  failures <- surv_data[
    surv_data[[event_var]] == 1,
    ,
    drop = FALSE
  ]

  if (nrow(failures) < 2) {
    stop("Fewer than two failures are available.")
  }

  km_cens <- survival::survfit(
    survival::Surv(surv_data[[time_var]], 1 - surv_data[[event_var]]) ~ 1
  )

  ghat <- summary(
    km_cens,
    times = failures[[time_var]],
    extend = TRUE
  )$surv

  ghat[!is.finite(ghat)] <- NA_real_

  if (all(is.na(ghat))) {
    stop("Could not estimate censoring survival probabilities.")
  }

  ghat[is.na(ghat)] <- min(ghat, na.rm = TRUE)
  ghat <- pmax(ghat, 1e-6)

  weights <- 1 / ghat
  weights <- weights / sum(weights)

  list(
    failures = failures,
    weights = weights
  )
}


sample_failure_subject <- function(
  failure_dist,
  long_data,
  id_var,
  time_var,
  long_time_var
) {
  pick <- sample.int(
    nrow(failure_dist$failures),
    size = 1,
    prob = failure_dist$weights
  )

  srow <- failure_dist$failures[pick, , drop = FALSE]
  old_id <- srow[[id_var]]
  selected_time <- srow[[time_var]]

  srow[[id_var]] <- 1

  lrows <- long_data[
    long_data[[id_var]] == old_id,
    ,
    drop = FALSE
  ]

  lrows[[id_var]] <- 1

  list(
    time = selected_time,
    surv = srow,
    long = lrows
  )
}


make_prediction_grid <- function(
  landmark,
  surv_train,
  surv_new,
  time_var,
  event_var
) {
  grid <- sort(unique(c(
    landmark,
    surv_train[[time_var]][
      surv_train[[time_var]] > landmark &
        surv_train[[event_var]] == 1
    ],
    if (time_var %in% names(surv_new)) {
      surv_new[[time_var]][surv_new[[time_var]] > landmark]
    } else {
      numeric(0)
    }
  )))

  grid <- grid[is.finite(grid)]

  if (length(grid) == 0) {
    stop("The prediction time grid is empty.")
  }

  grid
}


compute_survival_cutoffs <- function(
  scores,
  alpha,
  side = c("two", "lower", "upper")
) {
  side <- match.arg(side)

  if (side == "two") {
    return(list(
      lower_surv = as.numeric(stats::quantile(
        scores,
        probs = alpha / 2,
        names = FALSE,
        na.rm = TRUE
      )),
      upper_surv = as.numeric(stats::quantile(
        scores,
        probs = 1 - alpha / 2,
        names = FALSE,
        na.rm = TRUE
      ))
    ))
  }

  if (side == "lower") {
    return(list(
      lower_surv = NA_real_,
      upper_surv = as.numeric(stats::quantile(
        scores,
        probs = 1 - alpha,
        names = FALSE,
        na.rm = TRUE
      ))
    ))
  }

  list(
    lower_surv = as.numeric(stats::quantile(
      scores,
      probs = alpha,
      names = FALSE,
      na.rm = TRUE
    )),
    upper_surv = NA_real_
  )
}


invert_prediction_intervals <- function(
  pred,
  landmark,
  alpha,
  side = c("two", "lower", "upper"),
  cutoffs
) {
  side <- match.arg(side)

  if (side == "two") {
    lower <- invert_first_crossing(
      surv_mat = pred$surv,
      times = pred$times,
      threshold = cutoffs$upper_surv
    )

    upper <- invert_first_crossing(
      surv_mat = pred$surv,
      times = pred$times,
      threshold = cutoffs$lower_surv
    )
  }

  if (side == "lower") {
    lower <- invert_first_crossing(
      surv_mat = pred$surv,
      times = pred$times,
      threshold = cutoffs$upper_surv
    )

    upper <- rep(Inf, length(lower))
  }

  if (side == "upper") {
    lower <- rep(landmark, nrow(pred$surv))

    upper <- invert_first_crossing(
      surv_mat = pred$surv,
      times = pred$times,
      threshold = cutoffs$lower_surv
    )
  }

  data.frame(
    id = pred$ids,
    landmark = landmark,
    alpha = alpha,
    side = side,
    lower = lower,
    upper = upper
  )
}


invert_first_crossing <- function(
  surv_mat,
  times,
  threshold
) {
  M <- surv_mat <= threshold
  has_crossing <- rowSums(M, na.rm = TRUE) > 0

  first_idx <- max.col(M, ties.method = "first")

  out <- rep(Inf, nrow(surv_mat))
  out[has_crossing] <- times[first_idx[has_crossing]]

  out
}


extract_survival_matrix <- function(pm) {
  if (!("id" %in% names(pm))) {
    names(pm)[1] <- "id"
  }

  surv_cols <- setdiff(names(pm), "id")

  times <- suppressWarnings(
    as.numeric(
      sub("^S\\((.*)\\)$", "\\1", surv_cols)
    )
  )

  keep <- is.finite(times)

  if (!any(keep)) {
    stop("Could not parse survival prediction time columns.")
  }

  surv_cols <- surv_cols[keep]
  times <- times[keep]

  ord <- order(times)

  list(
    times = times[ord],
    surv = as.matrix(pm[, surv_cols[ord], drop = FALSE])
  )
}


prepare_surv_for_pencal <- function(
  surv_data,
  id_var,
  time_var,
  event_var
) {
  out <- as.data.frame(surv_data)

  if (id_var != "id") {
    out$id <- out[[id_var]]
  }

  if (time_var != "time") {
    out$time <- out[[time_var]]
  }

  if (event_var != "event") {
    out$event <- out[[event_var]]
  }

  out
}


prepare_long_for_pencal <- function(
  long_data,
  id_var,
  long_time_var
) {
  out <- as.data.frame(long_data)

  if (id_var != "id") {
    out$id <- out[[id_var]]
  }

  if (long_time_var != "t.from.base") {
    out$t.from.base <- out[[long_time_var]]
  }

  out
}


prepare_new_base_for_pencal <- function(
  new_base,
  id_var
) {
  out <- as.data.frame(new_base)

  if (id_var != "id") {
    out$id <- out[[id_var]]
  }

  out
}


check_packages <- function() {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
  }

  if (!requireNamespace("pencal", quietly = TRUE)) {
    stop("Package 'pencal' is required.")
  }

  invisible(TRUE)
}


check_dynamic_inputs <- function(
  surv_data,
  long_data,
  id_var,
  time_var,
  event_var,
  long_time_var,
  baseline_covariates,
  longitudinal_markers,
  require_event = TRUE
) {
  surv_required <- c(id_var, time_var, baseline_covariates)

  if (require_event) {
    surv_required <- c(surv_required, event_var)
  } else {
    if (event_var %in% names(surv_data)) {
      surv_required <- c(surv_required, event_var)
    }
  }

  check_columns(
    data = surv_data,
    required = surv_required,
    data_name = "survival data"
  )

  check_columns(
    data = long_data,
    required = c(id_var, long_time_var, longitudinal_markers),
    data_name = "longitudinal data"
  )

  invisible(TRUE)
}


check_columns <- function(
  data,
  required,
  data_name
) {
  missing <- setdiff(required, names(data))

  if (length(missing) > 0) {
    stop(
      data_name,
      " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }

  invisible(TRUE)
}
