# Final 3-predictor simulation data-generating mechanism.

if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required.")
}

three_predictor_parameters <- function() {
  list(
    visit_min = 4L,
    visit_max = 9L,
    visit_horizon = 5,
    age_mean = 60,
    age_sd = 10,
    age_min = 40,
    age_max = 90,
    re_sd = c(1.00, 0.20, 6.00, 0.80, 3.50, 0.55),
    re_rho_within = -0.20,
    long_intercept = c(10, 100, 50),
    long_slope = c(-0.30, 1.80, -1.20),
    long_error_sd = c(0.8, 5.0, 2.5),
    censor_min = 1,
    scenarios = list(
      weakNPH = list(
        mu0 = log(5), sigma_T = 0.50, beta_age = 0.18,
        beta_10 = 0.20, beta_11 = 0.18, beta_20 = 0.22, beta_31 = -0.20
      ),
      strongNPH = list(
        mu0 = log(5), sigma_T = 0.40, beta_age = 0.30,
        beta_10 = 0.40, beta_11 = 0.36, beta_20 = 0.42, beta_31 = -0.38
      )
    )
  )
}

.make_three_re_cov <- function(par) {
  R <- diag(6)
  R[1, 2] <- R[2, 1] <- par$re_rho_within
  R[3, 4] <- R[4, 3] <- par$re_rho_within
  R[5, 6] <- R[6, 5] <- par$re_rho_within
  D <- diag(par$re_sd)
  D %*% R %*% D
}

.draw_three_visits <- function(par) {
  m <- sample(par$visit_min:par$visit_max, 1L)
  c(0, sort(runif(m - 1L, 0, par$visit_horizon)))
}

simulate_three_predictor_latent <- function(
  n,
  scenario = c("weakNPH", "strongNPH"),
  seed = NULL
) {
  scenario <- match.arg(scenario)
  if (!is.null(seed)) set.seed(seed)

  par <- three_predictor_parameters()
  sc <- par$scenarios[[scenario]]
  id <- seq_len(n)

  baseline_age <- pmin(
    pmax(rnorm(n, par$age_mean, par$age_sd), par$age_min),
    par$age_max
  )
  age_std <- (baseline_age - par$age_mean) / par$age_sd

  b <- MASS::mvrnorm(
    n = n,
    mu = rep(0, 6),
    Sigma = .make_three_re_cov(par)
  )
  colnames(b) <- c("b10", "b11", "b20", "b21", "b30", "b31")

  eta <- sc$beta_age * age_std +
    sc$beta_10 * b[, "b10"] / par$re_sd[1] +
    sc$beta_11 * b[, "b11"] / par$re_sd[2] +
    sc$beta_20 * b[, "b20"] / par$re_sd[3] +
    sc$beta_31 * b[, "b31"] / par$re_sd[6]

  true_time <- 1 + exp(sc$mu0 + eta + sc$sigma_T * rnorm(n))

  long_list <- vector("list", n)
  for (i in seq_len(n)) {
    tt <- .draw_three_visits(par)

    y1 <- par$long_intercept[1] + par$long_slope[1] * tt +
      b[i, "b10"] + b[i, "b11"] * tt +
      rnorm(length(tt), 0, par$long_error_sd[1])

    y2 <- par$long_intercept[2] + par$long_slope[2] * tt +
      b[i, "b20"] + b[i, "b21"] * tt +
      rnorm(length(tt), 0, par$long_error_sd[2])

    y3 <- par$long_intercept[3] + par$long_slope[3] * tt +
      b[i, "b30"] + b[i, "b31"] * tt +
      rnorm(length(tt), 0, par$long_error_sd[3])

    long_list[[i]] <- data.frame(
      id = i,
      t.from.base = tt,
      age = baseline_age[i] + tt,
      y1 = y1,
      y2 = y2,
      y3 = y3
    )
  }

  list(
    baseline = data.frame(id = id, baseline.age = baseline_age),
    long_full = do.call(rbind, long_list),
    true_time = true_time,
    random_effects = data.frame(id = id, b)
  )
}

.censoring_probability <- function(T, cmax, cmin = 1) {
  if (cmax <= cmin) stop("cmax must be larger than cmin.")
  ifelse(T <= cmin, 0, ifelse(T >= cmax, 1, (T - cmin) / (cmax - cmin)))
}

calibrate_three_predictor_cmax <- function(
  scenario = c("weakNPH", "strongNPH"),
  target_censoring,
  n_calibration = 200000L,
  seed = 5001L
) {
  scenario <- match.arg(scenario)
  latent <- simulate_three_predictor_latent(
    n = n_calibration,
    scenario = scenario,
    seed = seed
  )

  f <- function(cmax) {
    mean(.censoring_probability(latent$true_time, cmax, 1)) - target_censoring
  }

  lower <- 1 + 1e-6
  upper <- max(stats::quantile(latent$true_time, 0.999), 2)
  while (f(upper) > 0) upper <- upper * 1.5
  stats::uniroot(f, c(lower, upper), tol = 1e-6)$root
}

simulate_three_predictor_dataset <- function(
  n,
  scenario = c("weakNPH", "strongNPH"),
  target_censoring = 0.20,
  cmax = NULL,
  seed = 1L,
  calibration_seed = 5001L,
  n_calibration = 200000L
) {
  scenario <- match.arg(scenario)
  par <- three_predictor_parameters()

  latent <- simulate_three_predictor_latent(n, scenario, seed)

  if (is.null(cmax)) {
    cmax <- calibrate_three_predictor_cmax(
      scenario, target_censoring, n_calibration, calibration_seed
    )
  }

  set.seed(seed + 100000L)
  censor_time <- runif(n, min = par$censor_min, max = cmax)
  obs_time <- pmin(latent$true_time, censor_time)
  event <- as.integer(latent$true_time <= censor_time)

  surv <- data.frame(
    id = latent$baseline$id,
    baseline.age = latent$baseline$baseline.age,
    time = obs_time,
    event = event,
    true.time = latent$true_time
  )

  long <- latent$long_full[
    latent$long_full$t.from.base <=
      obs_time[match(latent$long_full$id, surv$id)],
    ,
    drop = FALSE
  ]

  list(
    surv = surv,
    long = long,
    random_effects = latent$random_effects,
    meta = list(
      scenario = scenario,
      n = n,
      target_censoring = target_censoring,
      cmax = cmax,
      realized_censoring = mean(event == 0)
    )
  )
}
