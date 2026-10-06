# Sequential ADNI cross-validation for the final Chapter 4 comparison.

suppressPackageStartupMessages({
  library(survival)
  library(data.table)
})

source("R/load_dynamic_conformal.R")

.prepare_adni_data <- function(data_file) {
  e <- new.env(parent = emptyenv())
  load(data_file, envir = e)

  required <- c(
    "surv", "long.transformed",
    "all.fixed.cov.names", "all.long.cov.names"
  )
  missing <- required[!vapply(required, exists, logical(1), envir = e)]
  if (length(missing)) {
    stop("ADNI workspace is missing: ", paste(missing, collapse = ", "))
  }

  surv <- as.data.frame(e$surv)
  long <- as.data.frame(e$long.transformed)

  names(surv) <- make.names(names(surv))
  names(long) <- make.names(names(long))

  base_names <- make.names(e$all.fixed.cov.names)
  long_names <- make.names(e$all.long.cov.names)

  surv <- surv[, c("id", "time.to.event", "event", base_names), drop = FALSE]
  names(surv)[names(surv) == "time.to.event"] <- "time"

  long <- long[, c("id", "time.fup", long_names), drop = FALSE]
  names(long)[names(long) == "time.fup"] <- "t.from.base"

  surv$time <- as.numeric(surv$time)
  surv$event <- as.integer(as.character(surv$event))
  long$t.from.base <- as.numeric(long$t.from.base)

  for (v in base_names) {
    if (is.character(surv[[v]])) surv[[v]] <- factor(surv[[v]])
  }
  for (v in long_names) {
    long[[v]] <- suppressWarnings(as.numeric(as.character(long[[v]])))
  }

  ids <- sort(unique(surv$id))
  map <- data.frame(original_id = ids, id_new = seq_along(ids))
  surv$id <- map$id_new[match(surv$id, map$original_id)]
  long$id <- map$id_new[match(long$id, map$original_id)]

  ok <- complete.cases(surv[, base_names, drop = FALSE])
  keep <- surv$id[ok]
  surv <- surv[surv$id %in% keep, , drop = FALSE]
  long <- long[long$id %in% keep, , drop = FALSE]

  list(
    surv = surv[order(surv$id), , drop = FALSE],
    long = long[order(long$id, long$t.from.base), , drop = FALSE],
    base_names = base_names,
    long_names = long_names
  )
}

.landmark_split <- function(surv, long, ids, landmark) {
  s <- surv[surv$id %in% ids & surv$time > landmark, , drop = FALSE]
  l <- long[long$id %in% s$id & long$t.from.base <= landmark, , drop = FALSE]
  keep <- intersect(s$id, unique(l$id))
  list(
    surv = s[s$id %in% keep, , drop = FALSE],
    long = l[l$id %in% keep, , drop = FALSE]
  )
}

.baseline_formula <- function(base_names) {
  if (!length(base_names)) {
    return(stats::as.formula("Surv(time, event) ~ 1"))
  }
  stats::reformulate(base_names, response = "Surv(time, event)")
}

.baseline_cox_calibration_scores <- function(surv_train, base_names, B, seed) {
  set.seed(seed)
  form <- .baseline_formula(base_names)

  fail <- surv_train[surv_train$event == 1, , drop = FALSE]
  km <- survfit(Surv(time, 1 - event) ~ 1, data = surv_train)
  G <- summary(km, times = fail$time, extend = TRUE)$surv
  G[!is.finite(G) | G <= 0] <- 1e-8
  w <- (1 / G) / sum(1 / G)

  ids <- surv_train$id
  out <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    draw <- sample(ids, length(ids), replace = TRUE)

    boot_list <- lapply(seq_along(draw), function(k) {
      z <- surv_train[surv_train$id == draw[k], , drop = FALSE]
      z$id <- k
      z
    })
    boot <- do.call(rbind, boot_list)

    fit <- try(coxph(form, data = boot, x = TRUE), silent = TRUE)
    if (inherits(fit, "try-error")) next

    pick <- sample.int(nrow(fail), 1L, prob = w)
    sf <- try(
      survfit(
        fit,
        newdata = fail[pick, , drop = FALSE],
        se.fit = FALSE,
        conf.int = FALSE
      ),
      silent = TRUE
    )
    if (inherits(sf, "try-error")) next

    ss <- summary(sf, times = fail$time[pick], extend = TRUE)$surv
    if (length(ss)) out[b] <- as.numeric(ss[1])
  }

  out[is.finite(out)]
}

.predict_baseline_cox <- function(surv_train, surv_test, base_names, times) {
  fit <- coxph(.baseline_formula(base_names), data = surv_train, x = TRUE)
  S <- matrix(NA_real_, nrow(surv_test), length(times))

  for (i in seq_len(nrow(surv_test))) {
    sf <- survfit(
      fit,
      newdata = surv_test[i, , drop = FALSE],
      se.fit = FALSE,
      conf.int = FALSE
    )
    S[i, ] <- summary(sf, times = times, extend = TRUE)$surv
  }

  list(ids = surv_test$id, times = times, surv = S)
}

.diag_rows <- function(method, fold, landmark, intervals, surv_test, eta) {
  truth <- surv_test[match(intervals$id, surv_test$id), , drop = FALSE]
  upper <- intervals$upper
  lower <- intervals$lower

  event_obs <- truth$event == 1
  upper_ok <- is.infinite(upper) | is.na(upper) | truth$time <= upper
  lower_trunc_ok <- is.finite(lower) & truth$time >= lower & truth$time <= eta

  upper_length <- upper - landmark
  upper_length[!is.finite(upper_length) | upper_length < 0] <- NA_real_

  data.frame(
    method = method,
    fold = fold,
    landmark = landmark,
    n_test = nrow(truth),
    events = sum(event_obs),
    upper_cov = if (any(event_obs)) mean(upper_ok[event_obs]) else NA_real_,
    lower_trunc_cov = if (any(event_obs)) mean(lower_trunc_ok[event_obs]) else NA_real_,
    upper_length = mean(upper_length, na.rm = TRUE),
    inf_upper = mean(is.na(upper) | is.infinite(upper)),
    eta = eta
  )
}

run_adni_cv <- function(
  data_file,
  landmarks = c(2, 4),
  K = 5L,
  B = 500L,
  alpha = 0.10,
  seed = 2026L,
  output_dir = NULL,
  verbose = TRUE
) {
  dat <- .prepare_adni_data(data_file)

  set.seed(seed)
  ids <- sort(unique(dat$surv$id))
  fold_assign <- data.frame(
    id = ids,
    fold = sample(rep(seq_len(K), length.out = length(ids)))
  )

  result_list <- list()
  k_out <- 1L

  for (landmark in landmarks) {
    for (fold in seq_len(K)) {
      if (verbose) message("ADNI: landmark=", landmark, ", fold=", fold, "/", K)

      train_ids <- fold_assign$id[fold_assign$fold != fold]
      test_ids <- fold_assign$id[fold_assign$fold == fold]

      tr <- .landmark_split(dat$surv, dat$long, train_ids, landmark)
      te <- .landmark_split(dat$surv, dat$long, test_ids, landmark)

      if (nrow(tr$surv) < 5 || nrow(te$surv) == 0) next
      eta <- max(tr$surv$time[tr$surv$event == 1], na.rm = TRUE)

      prc <- dynamic_conformal_pi(
        surv_train = tr$surv,
        long_train = tr$long,
        surv_new = te$surv,
        long_new = te$long,
        landmark = landmark,
        id_var = "id",
        time_var = "time",
        event_var = "event",
        long_time_var = "t.from.base",
        baseline_covariates = dat$base_names,
        longitudinal_markers = dat$long_names,
        alpha = alpha,
        B = B,
        side = "two",
        lmm_fixefs = ~ t.from.base,
        lmm_ranefs = ~ t.from.base | id,
        seed = seed + 1000L * fold + as.integer(100 * landmark),
        verbose = FALSE
      )

      grid <- make_prediction_grid(
        landmark = landmark,
        surv_train = tr$surv,
        time_var = "time",
        event_var = "event"
      )

      prc_pred <- predict_prc_survival(
        prc_fit = prc$prc_fit,
        new_surv = te$surv,
        new_long = te$long,
        times = grid,
        id_var = "id",
        baseline_covariates = dat$base_names
      )

      conf_lower <- invert_prediction_intervals(
        prc_pred, landmark, alpha, "lower",
        compute_survival_cutoffs(prc$calibration_scores, alpha, "lower")
      )
      conf_upper <- invert_prediction_intervals(
        prc_pred, landmark, alpha, "upper",
        compute_survival_cutoffs(prc$calibration_scores, alpha, "upper")
      )
      conf_int <- data.frame(
        id = conf_lower$id,
        lower = conf_lower$lower,
        upper = conf_upper$upper
      )

      naive_lower <- invert_prediction_intervals(
        prc_pred, landmark, alpha, "lower",
        list(lower_surv = NA_real_, upper_surv = 1 - alpha)
      )
      naive_upper <- invert_prediction_intervals(
        prc_pred, landmark, alpha, "upper",
        list(lower_surv = alpha, upper_surv = NA_real_)
      )
      naive_int <- data.frame(
        id = naive_lower$id,
        lower = naive_lower$lower,
        upper = naive_upper$upper
      )

      base_scores <- .baseline_cox_calibration_scores(
        tr$surv, dat$base_names, B,
        seed + 50000L + 1000L * fold + as.integer(100 * landmark)
      )
      base_pred <- .predict_baseline_cox(
        tr$surv, te$surv, dat$base_names, grid
      )
      base_lower <- invert_prediction_intervals(
        base_pred, landmark, alpha, "lower",
        compute_survival_cutoffs(base_scores, alpha, "lower")
      )
      base_upper <- invert_prediction_intervals(
        base_pred, landmark, alpha, "upper",
        compute_survival_cutoffs(base_scores, alpha, "upper")
      )
      base_int <- data.frame(
        id = base_lower$id,
        lower = base_lower$lower,
        upper = base_upper$upper
      )

      result_list[[k_out]] <- rbind(
        .diag_rows("Naive PRC", fold, landmark, naive_int, te$surv, eta),
        .diag_rows("Conformal PRC", fold, landmark, conf_int, te$surv, eta),
        .diag_rows("Conformal baseline Cox", fold, landmark, base_int, te$surv, eta)
      )
      k_out <- k_out + 1L
    }
  }

  fold_results <- do.call(rbind, result_list)

  summary_results <- aggregate(
    cbind(
      n_test, events, upper_cov, lower_trunc_cov,
      upper_length, inf_upper
    ) ~ method + landmark,
    data = fold_results,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(
      fold_results,
      file.path(output_dir, "adni_cv_fold_results.csv")
    )
    data.table::fwrite(
      summary_results,
      file.path(output_dir, "adni_cv_summary.csv")
    )
  }

  list(
    fold_results = fold_results,
    summary = summary_results,
    fold_assignment = fold_assign
  )
}
