# Toy example for dynamicConformalSurv
#
# This file shows the expected data structure and a genuine prediction call.
# The numerical results are illustrative and are not intended to reproduce
# the manuscript simulations.

set.seed(1)

source("R/load_dynamic_conformal.R")

n <- 100
landmark <- 2

surv_data <- data.frame(
  id = 1:n,
  time = rexp(n, rate = 0.12) + 1,
  event = rbinom(n, size = 1, prob = 0.65),
  age = rnorm(n, mean = 65, sd = 8),
  sex = rbinom(n, size = 1, prob = 0.5)
)

long_data <- do.call(
  rbind,
  lapply(1:n, function(i) {
    visit_times <- sort(c(0, runif(4, 0, min(5, surv_data$time[i]))))

    data.frame(
      id = i,
      time_fup = visit_times,
      marker1 = 10 + 0.3 * visit_times + rnorm(length(visit_times), 0, 1),
      marker2 = 20 - 0.2 * visit_times + rnorm(length(visit_times), 0, 1)
    )
  })
)

train_ids <- sample(surv_data$id, size = 70)
new_ids <- setdiff(surv_data$id, train_ids)

surv_train <- surv_data[surv_data$id %in% train_ids, ]
long_train <- long_data[long_data$id %in% train_ids, ]

# In a real prediction problem, eligibility at the landmark is known at the
# prediction time. Future outcomes are therefore not passed to the prediction
# function.
new_risk_ids <- surv_data$id[
  surv_data$id %in% new_ids & surv_data$time > landmark
]

surv_new <- surv_data[
  surv_data$id %in% new_risk_ids,
  c("id", "age", "sex"),
  drop = FALSE
]

long_new <- long_data[
  long_data$id %in% new_risk_ids &
    long_data$time_fup <= landmark,
  ,
  drop = FALSE
]

fit <- dynamic_conformal_pi(
  surv_train = surv_train,
  long_train = long_train,
  surv_new = surv_new,
  long_new = long_new,
  landmark = landmark,
  id_var = "id",
  time_var = "time",
  event_var = "event",
  long_time_var = "time_fup",
  baseline_covariates = c("age", "sex"),
  longitudinal_markers = c("marker1", "marker2"),
  alpha = 0.10,
  B = 20,
  side = "two",
  seed = 123,
  verbose = TRUE
)

head(fit$intervals)
fit$cutoffs
