# Serial Monte Carlo runner for the simulation study.
#
# Examples:
#   Rscript reproduction/simulations/run_simulation.R --setting 3 --scenario weakNPH --N 300 --censor 0.20 --landmark 2 --reps 10 --B 50
#   Rscript reproduction/simulations/run_simulation.R --setting 20 --N 300 --censor 0.20 --landmark 2 --reps 10 --B 50
#
# Manuscript settings use reps=1000, B=500, validation-n=10000.

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

setting <- as.integer(get_arg("--setting", "3"))
scenario <- get_arg("--scenario", if (setting == 3L) "weakNPH" else "strongNPH")
N <- as.integer(get_arg("--N", "300"))
target <- as.numeric(get_arg("--censor", "0.20"))
landmark <- as.numeric(get_arg("--landmark", "2"))
n_rep <- as.integer(get_arg("--reps", "1"))
B <- as.integer(get_arg("--B", "20"))
alpha <- as.numeric(get_arg("--alpha", "0.10"))
validation_n <- as.integer(get_arg("--validation-n", "10000"))
seed <- as.integer(get_arg("--seed", "2026"))
out_file <- get_arg("--output", file.path("results", sprintf("sim_setting%d_%s_N%d_cens%02d_lmk%s.csv", setting, scenario, N, round(100 * target), landmark)))

if (!(setting %in% c(3L, 20L))) stop("--setting must be 3 or 20")
if (!(target %in% c(0.20, 0.35, 0.50))) stop("--censor must be 0.20, 0.35, or 0.50")
if (!(landmark %in% c(2, 4))) stop("--landmark must be 2 or 4 for the main study")
if (setting == 20L && scenario != "strongNPH") stop("The 20-marker study uses strongNPH only")

if (setting == 3L) {
  cmap <- sim3_censor_map(targets = target, seed = seed + 10L)
  cmax <- cmap[[scenario]][[as.character(target)]]
  val_lat <- sim3_latent(validation_n, scenario = scenario, seed = seed + 700000L)
  val <- sim3_observe(val_lat, target, cmax, seed = seed + 710000L)
  markers <- c("y1", "y2", "y3")
} else {
  cmap <- sim20_censor_map(targets = target, seed = seed + 10L)
  cmax <- cmap[[as.character(target)]]
  val_lat <- sim20_latent(validation_n, seed = seed + 700000L)
  val <- sim20_observe(val_lat, target, cmax, seed = seed + 800000L)
  markers <- paste0("y", 1:20)
}

rows <- vector("list", n_rep * 2L)
k <- 1L

for (rep_id in seq_len(n_rep)) {
  message("replicate ", rep_id, " / ", n_rep)

  if (setting == 3L) {
    lat <- sim3_latent(N, scenario = scenario, seed = seed + 100000L + rep_id)
    train <- sim3_observe(lat, target, cmax, seed = seed + 200000L + rep_id)
  } else {
    lat <- sim20_latent(N, seed = seed + 100000L + rep_id)
    train <- sim20_observe(lat, target, cmax, seed = seed + 200000L + rep_id)
  }

  ans <- coverage_one_run_prc(
    surv_train = train$surv.data,
    long_train = train$long.data,
    surv_valid = val$surv.data,
    long_valid = val$long.data,
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
    penalty = "ridge",
    standardize = TRUE,
    seed = seed + 900000L + rep_id,
    verbose = FALSE
  )

  m <- ans$metrics
  m$rep <- rep_id
  m$setting <- setting
  m$scenario <- scenario
  m$N <- N
  m$censor_target <- target
  m$B <- B
  m$validation_n <- validation_n
  m$realized_train_censoring <- mean(train$surv.data$event == 0)
  m$realized_validation_censoring <- mean(val$surv.data$event == 0)

  rows[[k]] <- m[1, ]; k <- k + 1L
  rows[[k]] <- m[2, ]; k <- k + 1L
}

out <- do.call(rbind, rows)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(out, out_file, row.names = FALSE)

summary <- aggregate(
  cbind(right_cov, left_cov, total_cov, trunc_total_cov, avg_trunc_length, inf_upper) ~ method + setting + scenario + N + censor_target + landmark,
  data = out,
  FUN = mean,
  na.rm = TRUE
)

print(summary)
message("Saved per-replicate results to: ", out_file)
