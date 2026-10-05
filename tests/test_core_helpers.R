# Lightweight checks for core helper functions.
#
# These tests do not fit PRC models and therefore do not require pencal.

source("R/dynamic_conformal_pi.R")

# Prediction grid must depend only on the training data.
train <- data.frame(
  id = 1:4,
  time = c(3, 5, 7, 4),
  event = c(1, 0, 1, 1)
)

grid <- make_prediction_grid(
  landmark = 2,
  surv_train = train,
  time_var = "time",
  event_var = "event"
)

stopifnot(identical(grid, c(2, 3, 4, 7)))
stopifnot(!("surv_new" %in% names(formals(make_prediction_grid))))

# Subject-level bootstrap must preserve multiplicities.
surv <- data.frame(
  id = 1:3,
  time = c(3, 4, 5),
  event = c(1, 1, 0)
)

long <- data.frame(
  id = c(1, 1, 2, 2, 3, 3),
  t.from.base = rep(c(0, 1), 3),
  y1 = seq_len(6)
)

boot <- make_bootstrap_data(
  surv_data = surv,
  long_data = long,
  draw_ids = c(1, 1, 3),
  id_var = "id"
)

stopifnot(nrow(boot$surv) == 3)
stopifnot(identical(boot$surv$id, 1:3))
stopifnot(sum(boot$long$id == 1) == 2)
stopifnot(sum(boot$long$id == 2) == 2)
stopifnot(sum(boot$long$id == 3) == 2)

# Survival cutoffs must be ordered for a two-sided interval.
scores <- seq(0.01, 0.99, length.out = 99)
cuts <- compute_survival_cutoffs(scores, alpha = 0.10, side = "two")
stopifnot(cuts$lower_surv < cuts$upper_surv)

# First-crossing inversion.
surv_mat <- rbind(
  c(1.00, 0.90, 0.70, 0.40),
  c(1.00, 0.95, 0.92, 0.91)
)
times <- c(2, 3, 4, 5)

inv <- invert_first_crossing(
  surv_mat = surv_mat,
  times = times,
  threshold = 0.80
)

stopifnot(inv[1] == 4)
stopifnot(is.infinite(inv[2]))

cat("All core helper checks passed.\n")
