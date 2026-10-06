# Alpha-sensitivity analysis for the 20-marker strong-NPH setting.
# The manuscript analysis used N=1000, censoring target 35%, landmarks 2 and 4,
# 1000 Monte Carlo replications, and alpha in {0.05, 0.10, 0.15, 0.20}.

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1L]
}

repo_root <- normalizePath(get_arg("--repo", "."), mustWork = TRUE)
source(file.path(repo_root, "R", "load_dynamic_conformal.R"))
source(file.path(repo_root, "reproduction", "simulations", "dgp_3marker.R"))
source(file.path(repo_root, "reproduction", "simulations", "dgp_20marker.R"))

N <- as.integer(get_arg("--N", "1000"))
target <- as.numeric(get_arg("--censor", "0.35"))
n_rep <- as.integer(get_arg("--reps", "1"))
B <- as.integer(get_arg("--B", "20"))
validation_n <- as.integer(get_arg("--validation-n", "10000"))
seed <- as.integer(get_arg("--seed", "2026"))
out_file <- get_arg("--output", file.path("results", "alpha_sensitivity_20marker.csv"))

alphas <- c(0.05, 0.10, 0.15, 0.20)
landmarks <- c(2, 4)
markers <- paste0("y", 1:20)

cmap <- sim20_censor_map(targets = target, seed = seed + 10L)
cmax <- cmap[[as.character(target)]]
val_lat <- sim20_latent(validation_n, seed = seed + 700000L)
val <- sim20_observe(val_lat, target, cmax, seed = seed + 800000L)

rows <- list()
k <- 1L

for (rep_id in seq_len(n_rep)) {
  message("replicate ", rep_id, " / ", n_rep)
  lat <- sim20_latent(N, seed = seed + 100000L + rep_id)
  train <- sim20_observe(lat, target, cmax, seed = seed + 200000L + rep_id)

  for (landmark in landmarks) {
    # Fit/calibrate once at alpha=0.10; calibration scores themselves do not depend on alpha.
    base_fit <- dynamic_conformal_pi(
      surv_train = train$surv.data,
      long_train = train$long.data,
      surv_new = val$surv.data,
      long_new = val$long.data,
      landmark = landmark,
      id_var = "id",
      time_var = "time",
      event_var = "event",
      long_time_var = "t.from.base",
      baseline_covariates = "baseline.age",
      longitudinal_markers = markers,
      alpha = 0.10,
      B = B,
      side = "two",
      lmm_fixefs = ~ age,
      lmm_ranefs = ~ age | id,
      seed = seed + 900000L + rep_id + 100L * landmark,
      verbose = FALSE
    )

    train_lmk <- make_landmark_data(train$surv.data, train$long.data, landmark, "id", "time", "t.from.base")
    val_lmk <- make_new_landmark_data(val$surv.data, val$long.data, landmark, "id", "time", "t.from.base")
    grid <- make_prediction_grid(landmark, train_lmk$surv, "time", "event")
    pred <- predict_prc_survival(base_fit$prc_fit, val_lmk$surv, val_lmk$long, grid, "id", "baseline.age")
    truth <- val_lmk$surv$true.time[match(pred$ids, val_lmk$surv$id)]
    eta <- max(train_lmk$surv$time[train_lmk$surv$event == 1], na.rm = TRUE)

    for (alpha in alphas) {
      cut_conf <- compute_survival_cutoffs(base_fit$calibration_scores, alpha, "two")
      cut_naive <- list(lower_surv = alpha / 2, upper_surv = 1 - alpha / 2)

      for (method in c("Naive PRC", "Conformal PRC")) {
        cut <- if (method == "Naive PRC") cut_naive else cut_conf
        int <- invert_prediction_intervals(pred, landmark, alpha, "two", cut)
        upper_trunc <- pmin(int$upper, eta)
        len <- pmax(upper_trunc - int$lower, 0)

        rows[[k]] <- data.frame(
          rep = rep_id,
          landmark = landmark,
          alpha = alpha,
          nominal_coverage = 1 - alpha,
          method = method,
          mean_truncated_length = mean(len[is.finite(len)], na.rm = TRUE)
        )
        k <- k + 1L
      }
    }
  }
}

out <- do.call(rbind, rows)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(out, out_file, row.names = FALSE)
print(aggregate(mean_truncated_length ~ landmark + alpha + nominal_coverage + method, out, mean, na.rm = TRUE))
message("Saved: ", out_file)
