# Serial reproduction of the ADNI cross-validation analysis.
#
# Requires the processed ADNI .RData object used in the thesis. Participant-level
# data are not distributed with this repository.

suppressPackageStartupMessages({
  library(survival)
})

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1L]
}

repo_root <- normalizePath(get_arg("--repo", "."), mustWork = TRUE)
data_file <- get_arg("--data", Sys.getenv("ADNI_DATA_FILE", ""))
out_dir <- get_arg("--output", file.path("results", "adni_cv"))
B <- as.integer(get_arg("--B", "500"))
K <- as.integer(get_arg("--K", "5"))
seed <- as.integer(get_arg("--seed", "2026"))
alpha <- 0.10
landmarks <- c(2, 3, 4)

if (!nzchar(data_file) || !file.exists(data_file)) {
  stop("Provide the processed ADNI .RData file with --data or ADNI_DATA_FILE.")
}
if (B < 2L) stop("--B must be at least 2")
if (K < 2L) stop("--K must be at least 2")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(repo_root, "R", "load_dynamic_conformal.R"))
load(data_file)

# Expected objects: surv, long.transformed, all.fixed.cov.names, all.long.cov.names.
surv_data <- as.data.frame(surv)
long_data <- as.data.frame(long.transformed)
names(surv_data) <- make.names(names(surv_data))
names(long_data) <- make.names(names(long_data))
base_names <- make.names(all.fixed.cov.names)
long_names <- make.names(all.long.cov.names)

surv_data <- surv_data[, c("id", "time.to.event", "event", base_names), drop = FALSE]
names(surv_data)[names(surv_data) == "time.to.event"] <- "time"
long_data <- long_data[, c("id", "time.fup", long_names), drop = FALSE]
names(long_data)[names(long_data) == "time.fup"] <- "t.from.base"

surv_data$time <- as.numeric(surv_data$time)
surv_data$event <- as.integer(as.character(surv_data$event))
long_data$t.from.base <- as.numeric(long_data$t.from.base)
for (v in base_names) if (is.character(surv_data[[v]])) surv_data[[v]] <- factor(surv_data[[v]])
for (v in long_names) long_data[[v]] <- suppressWarnings(as.numeric(as.character(long_data[[v]])))

old_ids <- sort(unique(surv_data$id))
id_map <- data.frame(original_id = old_ids, id = seq_along(old_ids))
surv_data$id <- id_map$id[match(surv_data$id, id_map$original_id)]
long_data$id <- id_map$id[match(long_data$id, id_map$original_id)]

keep_ids <- surv_data$id[complete.cases(surv_data[, base_names, drop = FALSE])]
surv_data <- surv_data[surv_data$id %in% keep_ids, , drop = FALSE]
long_data <- long_data[long_data$id %in% keep_ids, , drop = FALSE]
surv_data <- surv_data[order(surv_data$id), , drop = FALSE]
long_data <- long_data[order(long_data$id, long_data$t.from.base), , drop = FALSE]

set.seed(seed)
ids_all <- sort(unique(surv_data$id))
fold_id <- sample(rep(seq_len(K), length.out = length(ids_all)))
fold_map <- data.frame(id = ids_all, fold = fold_id)

landmark_subset <- function(ids, landmark) {
  s <- surv_data[surv_data$id %in% ids & surv_data$time > landmark, , drop = FALSE]
  l <- long_data[long_data$id %in% s$id & long_data$t.from.base <= landmark, , drop = FALSE]
  keep <- intersect(s$id, unique(l$id))
  list(
    surv = s[s$id %in% keep, , drop = FALSE],
    long = l[l$id %in% keep, , drop = FALSE]
  )
}

prediction_from_fit <- function(fit, tr_s, te_s, te_l, landmark) {
  grid <- make_prediction_grid(landmark, tr_s, "time", "event")
  predict_prc_survival(
    prc_fit = fit,
    new_surv = te_s,
    new_long = te_l,
    times = grid,
    id_var = "id",
    baseline_covariates = base_names
  )
}

interval_rows <- function(method, fold, landmark, te_s, lower, upper, eta) {
  ids <- lower$id
  data.frame(
    fold = fold,
    landmark = landmark,
    method = method,
    id = ids,
    obs_time = te_s$time[match(ids, te_s$id)],
    event = te_s$event[match(ids, te_s$id)],
    lower_endpoint = lower$lower,
    upper_endpoint = upper$upper,
    eta = eta,
    stringsAsFactors = FALSE
  )
}

# Baseline Cox helpers. The same folds and landmark risk sets are used as for PRC.
make_baseline_design <- function(train_raw, test_raw) {
  tr <- train_raw
  te <- test_raw

  for (v in base_names) {
    if (is.numeric(tr[[v]]) || is.integer(tr[[v]])) {
      med <- stats::median(tr[[v]], na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      tr[[v]][is.na(tr[[v]])] <- med
      te[[v]][is.na(te[[v]])] <- med
    } else {
      a <- as.character(tr[[v]])
      b <- as.character(te[[v]])
      tab <- sort(table(a), decreasing = TRUE)
      mode_value <- if (length(tab)) names(tab)[1] else "missing"
      a[is.na(a) | a == ""] <- mode_value
      b[is.na(b) | b == ""] <- mode_value
      lev <- sort(unique(c(a, b)))
      tr[[v]] <- factor(a, levels = lev)
      te[[v]] <- factor(b, levels = lev)
    }
  }

  tr$.set <- "train"
  te$.set <- "test"
  both <- rbind(tr, te)
  form <- stats::as.formula(paste("~", paste(sprintf("`%s`", base_names), collapse = " + ")))
  mm <- stats::model.matrix(form, data = both)
  if ("(Intercept)" %in% colnames(mm)) mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]

  xnames <- paste0("x", seq_len(ncol(mm)))
  colnames(mm) <- xnames
  d <- cbind(both[, c("id", "time", "event", ".set")], as.data.frame(mm))
  tr_d <- d[d$.set == "train", , drop = FALSE]
  te_d <- d[d$.set == "test", , drop = FALSE]
  tr_d$.set <- NULL
  te_d$.set <- NULL

  keep <- vapply(xnames, function(v) all(is.finite(tr_d[[v]])) && stats::sd(tr_d[[v]]) > 0, logical(1))
  xnames <- xnames[keep]
  list(train = tr_d[, c("id", "time", "event", xnames), drop = FALSE], test = te_d[, c("id", "time", "event", xnames), drop = FALSE], xnames = xnames)
}

fit_baseline_cox <- function(d, xnames) {
  form <- if (length(xnames)) {
    stats::as.formula(paste("survival::Surv(time, event) ~", paste(xnames, collapse = " + ")))
  } else {
    survival::Surv(time, event) ~ 1
  }
  fit <- survival::coxph(form, data = d, ties = "breslow", singular.ok = TRUE, x = TRUE, model = TRUE)
  if (length(stats::coef(fit))) fit$coefficients[is.na(fit$coefficients)] <- 0
  list(fit = fit, xnames = xnames)
}

baseline_survival_matrix <- function(obj, newdata, times, landmark) {
  bh <- survival::basehaz(obj$fit, centered = FALSE)
  H <- if (nrow(bh)) stats::approx(bh$time, bh$hazard, xout = times, method = "constant", f = 0, rule = 2, yleft = 0)$y else rep(0, length(times))
  H[times <= landmark] <- 0
  lp <- if (length(obj$xnames)) as.numeric(as.matrix(newdata[, obj$xnames, drop = FALSE]) %*% stats::coef(obj$fit)[obj$xnames]) else rep(0, nrow(newdata))
  exp(-outer(exp(lp), H, "*"))
}

baseline_survival_at <- function(obj, newdata, time, landmark) {
  bh <- survival::basehaz(obj$fit, centered = FALSE)
  H <- if (nrow(bh)) stats::approx(bh$time, bh$hazard, xout = time, method = "constant", f = 0, rule = 2, yleft = 0)$y else rep(0, length(time))
  H[time <= landmark] <- 0
  lp <- if (length(obj$xnames)) as.numeric(as.matrix(newdata[, obj$xnames, drop = FALSE]) %*% stats::coef(obj$fit)[obj$xnames]) else rep(0, nrow(newdata))
  exp(-exp(lp) * H)
}

run_baseline <- function(tr_s, te_s, landmark, fold) {
  design <- make_baseline_design(tr_s, te_s)
  tr <- design$train
  te <- design$test
  fit0 <- fit_baseline_cox(tr, design$xnames)
  eta <- max(tr$time[tr$event == 1], na.rm = TRUE)
  grid <- sort(unique(c(landmark, tr$time[tr$event == 1])))

  km <- survival::survfit(survival::Surv(time, 1 - event) ~ 1, data = tr)
  failures <- tr[tr$event == 1, , drop = FALSE]
  G <- summary(km, times = failures$time, extend = TRUE)$surv
  G[!is.finite(G) | G <= 0] <- 1e-8
  prob <- (1 / G) / sum(1 / G)

  set.seed(seed + 100000L + 1000L * fold + as.integer(100 * landmark))
  boot_draws <- replicate(B, sample(seq_len(nrow(tr)), nrow(tr), replace = TRUE), simplify = FALSE)
  fail_picks <- sample.int(nrow(failures), B, replace = TRUE, prob = prob)
  scores <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    scores[b] <- tryCatch({
      fit_b <- fit_baseline_cox(tr[boot_draws[[b]], , drop = FALSE], design$xnames)
      cal <- failures[fail_picks[b], , drop = FALSE]
      baseline_survival_at(fit_b, cal, cal$time, landmark)[1]
    }, error = function(e) NA_real_)
  }
  scores <- scores[is.finite(scores)]
  if (length(scores) < 2L) stop("Too few baseline-Cox calibration scores")

  S <- baseline_survival_matrix(fit0, te, grid, landmark)
  pred <- list(ids = te$id, times = grid, surv = S)
  lower <- invert_prediction_intervals(pred, landmark, alpha, "lower", compute_survival_cutoffs(scores, alpha, "lower"))
  upper <- invert_prediction_intervals(pred, landmark, alpha, "upper", compute_survival_cutoffs(scores, alpha, "upper"))
  interval_rows("Conformal baseline Cox", fold, landmark, te_s, lower, upper, eta)
}

rows <- list()
k <- 1L

for (landmark in landmarks) {
  for (fold in seq_len(K)) {
    message("landmark=", landmark, " fold=", fold)
    train_ids <- fold_map$id[fold_map$fold != fold]
    test_ids <- fold_map$id[fold_map$fold == fold]
    tr <- landmark_subset(train_ids, landmark)
    te <- landmark_subset(test_ids, landmark)
    if (!nrow(tr$surv) || !nrow(te$surv)) stop("Empty landmark train/test set")

    fit <- dynamic_conformal_pi(
      surv_train = tr$surv,
      long_train = tr$long,
      surv_new = te$surv,
      long_new = te$long,
      landmark = landmark,
      id_var = "id",
      time_var = "time",
      event_var = "event",
      long_time_var = "t.from.base",
      baseline_covariates = base_names,
      longitudinal_markers = long_names,
      alpha = alpha,
      B = B,
      side = "two",
      lmm_fixefs = ~ t.from.base,
      lmm_ranefs = ~ t.from.base | id,
      n_cores = 1,
      seed = seed + 100000L * fold + 1000L * as.integer(round(100 * landmark)) + B,
      verbose = FALSE
    )

    pred <- prediction_from_fit(fit$prc_fit, tr$surv, te$surv, te$long, landmark)
    eta <- max(tr$surv$time[tr$surv$event == 1], na.rm = TRUE)

    naive_lower <- invert_prediction_intervals(pred, landmark, alpha, "lower", list(lower_surv = NA_real_, upper_surv = 1 - alpha))
    naive_upper <- invert_prediction_intervals(pred, landmark, alpha, "upper", list(lower_surv = alpha, upper_surv = NA_real_))
    conf_lower <- invert_prediction_intervals(pred, landmark, alpha, "lower", compute_survival_cutoffs(fit$calibration_scores, alpha, "lower"))
    conf_upper <- invert_prediction_intervals(pred, landmark, alpha, "upper", compute_survival_cutoffs(fit$calibration_scores, alpha, "upper"))

    rows[[k]] <- interval_rows("Naive PRC", fold, landmark, te$surv, naive_lower, naive_upper, eta); k <- k + 1L
    rows[[k]] <- interval_rows("Conformal PRC", fold, landmark, te$surv, conf_lower, conf_upper, eta); k <- k + 1L
    rows[[k]] <- run_baseline(tr$surv, te$surv, landmark, fold); k <- k + 1L
  }
}

intervals <- do.call(rbind, rows)
utils::write.csv(intervals, file.path(out_dir, "adni_cv_intervals.csv"), row.names = FALSE)

summarise_group <- function(d) {
  ev <- d$event == 1
  upper_inf <- is.na(d$upper_endpoint) | is.infinite(d$upper_endpoint)
  upper_len <- d$upper_endpoint - d$landmark
  upper_len[!is.finite(upper_len) | upper_len < 0] <- NA_real_
  upper_trunc <- ifelse(upper_inf, d$eta, pmin(d$upper_endpoint, d$eta)) - d$landmark
  upper_trunc[!is.finite(upper_trunc) | upper_trunc < 0] <- NA_real_

  data.frame(
    n_test = nrow(d),
    n_events = sum(ev),
    upper_cov = mean(is.na(d$upper_endpoint[ev]) | is.infinite(d$upper_endpoint[ev]) | d$obs_time[ev] <= d$upper_endpoint[ev], na.rm = TRUE),
    upper_length_finite = mean(upper_len, na.rm = TRUE),
    upper_length_truncated = mean(upper_trunc, na.rm = TRUE),
    inf_upper = mean(upper_inf),
    lower_trunc_cov = mean(
      is.finite(d$lower_endpoint[ev]) & d$obs_time[ev] >= d$lower_endpoint[ev] & d$obs_time[ev] <= d$eta[ev],
      na.rm = TRUE
    )
  )
}

groups <- split(intervals, interaction(intervals$landmark, intervals$method, drop = TRUE))
summary <- do.call(rbind, lapply(groups, function(d) {
  ans <- summarise_group(d)
  ans$landmark <- d$landmark[1]
  ans$method <- d$method[1]
  ans
}))
rownames(summary) <- NULL
summary$method_order <- match(summary$method, c("Naive PRC", "Conformal PRC", "Conformal baseline Cox"))
summary <- summary[order(summary$landmark, summary$method_order), c("landmark", "method", "n_test", "n_events", "upper_cov", "upper_length_finite", "upper_length_truncated", "inf_upper", "lower_trunc_cov")]
utils::write.csv(summary, file.path(out_dir, "adni_cv_summary.csv"), row.names = FALSE)
print(summary)

cat("\nupper_length_finite is the historical table convention (finite endpoints only).\n")
cat("upper_length_truncated truncates infinite endpoints at the largest training event time.\n")
