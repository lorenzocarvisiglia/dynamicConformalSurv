# Final 20-predictor simulation data-generating mechanism.

if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required.")
}

twenty_predictor_parameters <- function() {
  beta_t_4_20 <- c(
    0.12, 0.08, -0.06, -0.09, -0.12, 0.10, 0.07, 0.03, -0.04,
    0.13, -0.10, 0.05, -0.07, 0.11, -0.05, 0.09, -0.02
  )

  active_names <- c(
    as.vector(rbind(paste0("b0_y", 1:5), paste0("b1_y", 1:5))),
    paste0("b0_y", 6:10),
    paste0("b1_y", 11:15)
  )

  active_sd <- c(
    0.85, 0.24, 0.80, 0.20, 0.75, 0.22, 0.70, 0.18, 0.65, 0.16,
    seq(0.60, 0.80, length.out = 5),
    seq(0.14, 0.24, length.out = 5)
  )

  list(
    visit_min = 4L, visit_max = 9L, visit_horizon = 5,
    age_mean = 60, age_sd = 10, age_min = 40, age_max = 85,
    beta0 = c(0.80, -0.10, 0.40, seq(0.65, -0.45, length.out = 17)),
    beta_t = c(0.16, 0.10, -0.06, beta_t_4_20),
    beta_age = c(0.20, -0.10, 0.08, seq(0.16, -0.10, length.out = 17)),
    err_sd = c(0.35, 0.30, 0.35, seq(0.30, 0.45, length.out = 17)),
    active_names = active_names,
    active_sd = active_sd,
    re_rho = 0.35,
    mu0 = 1.65, sigma_T = 0.90, gamma_age = 0.22,
    gamma_b10 = 0.52, gamma_b11 = 0.36,
    gamma_b20 = -0.42, gamma_b31 = 0.40,
    censor_min = 1
  )
}

.draw_twenty_visits <- function(par) {
  m <- sample(par$visit_min:par$visit_max, 1L)
  c(0, sort(runif(m - 1L, 0, par$visit_horizon)))
}

.make_twenty_re_cov <- function(par) {
  p <- length(par$active_sd)
  R <- outer(seq_len(p), seq_len(p), function(i, j) par$re_rho^abs(i - j))
  D <- diag(par$active_sd)
  D %*% R %*% D
}

simulate_twenty_predictor_latent <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  par <- twenty_predictor_parameters()

  id <- seq_len(n)
  baseline_age <- pmin(
    pmax(rnorm(n, par$age_mean, par$age_sd), par$age_min),
    par$age_max
  )
  age_std <- (baseline_age - par$age_mean) / par$age_sd

  b_active <- MASS::mvrnorm(
    n = n,
    mu = rep(0, length(par$active_names)),
    Sigma = .make_twenty_re_cov(par)
  )
  colnames(b_active) <- par$active_names

  b0 <- matrix(0, n, 20)
  b1 <- matrix(0, n, 20)

  for (j in 1:5) {
    b0[, j] <- b_active[, paste0("b0_y", j)]
    b1[, j] <- b_active[, paste0("b1_y", j)]
  }
  for (j in 6:10) b0[, j] <- b_active[, paste0("b0_y", j)]
  for (j in 11:15) b1[, j] <- b_active[, paste0("b1_y", j)]

  mu_T <- par$mu0 +
    par$gamma_age * age_std +
    par$gamma_b10 * b0[, 1] +
    par$gamma_b11 * b1[, 1] +
    par$gamma_b20 * b0[, 2] +
    par$gamma_b31 * b1[, 3]

  true_time <- 1 + exp(mu_T + par$sigma_T * rnorm(n))

  long_list <- vector("list", n)
  for (i in seq_len(n)) {
    tt <- .draw_twenty_visits(par)
    out <- data.frame(
      id = i,
      t.from.base = tt,
      age = baseline_age[i] + tt
    )

    for (j in 1:20) {
      out[[paste0("y", j)]] <-
        par$beta0[j] + par$beta_t[j] * tt + par$beta_age[j] * age_std[i] +
        b0[i, j] + b1[i, j] * tt +
        rnorm(length(tt), 0, par$err_sd[j])
    }
    long_list[[i]] <- out
  }

  list(
    baseline = data.frame(id = id, baseline.age = baseline_age),
    long_full = do.call(rbind, long_list),
    true_time = true_time,
    random_effects = data.frame(id = id, b_active, check.names = FALSE)
  )
}

.censor_prob_20 <- function(T, cmax) {
  p <- ifelse(T < cmax, (T - 1) / (cmax - 1), 1)
  pmin(pmax(p, 0), 1)
}

calibrate_twenty_predictor_cmax <- function(
  target_censoring,
  n_calibration = 200000L,
  seed = 5001L
) {
  latent <- simulate_twenty_predictor_latent(n_calibration, seed)
  f <- function(cmax) mean(.censor_prob_20(latent$true_time, cmax)) - target_censoring

  lower <- 1 + 1e-6
  upper <- max(stats::quantile(latent$true_time, 0.999), 2)
  while (f(upper) > 0) upper <- upper * 1.5
  stats::uniroot(f, c(lower, upper), tol = 1e-6)$root
}

simulate_twenty_predictor_dataset <- function(
  n,
  target_censoring = 0.20,
  cmax = NULL,
  seed = 1L,
  calibration_seed = 5001L,
  n_calibration = 200000L
) {
  latent <- simulate_twenty_predictor_latent(n, seed)

  if (is.null(cmax)) {
    cmax <- calibrate_twenty_predictor_cmax(
      target_censoring, n_calibration, calibration_seed
    )
  }

  set.seed(seed + 100000L)
  censor_time <- runif(n, 1, cmax)
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
      scenario = "strongNPH",
      n = n,
      target_censoring = target_censoring,
      cmax = cmax,
      realized_censoring = mean(event == 0)
    )
  )
}
