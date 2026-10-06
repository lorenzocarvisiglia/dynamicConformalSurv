# Data-generating mechanism for the 20-marker strong-NPH simulation.
#
# This is a user-independent version of the final generator used for Scenario C.

if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required.")
}

sim20_parameters <- function() {
  list(
    landmarks_descriptive = c(2, 2.5, 3, 4, 4.5),
    landmarks_methods = c(2, 4),
    censor_targets = c(0.20, 0.35, 0.50),
    n_values = c(300L, 1000L, 1500L),
    n_validation = 10000L,
    n_markers = 20L,
    visit_min = 4L,
    visit_max = 9L,
    visit_horizon = 5,
    age_mean = 60,
    age_sd = 10,
    age_min = 40,
    age_max = 85,
    censor_min = 1,
    mu0 = 1.65,
    sigma = 0.90,
    gamma_age = 0.22,
    gamma_b10 = 0.52,
    gamma_b11 = 0.36,
    gamma_b20 = -0.42,
    gamma_b31 = 0.40,
    re_rho = 0.35
  )
}

sim20_fixed_parameters <- function() {
  beta0 <- rep(NA_real_, 20)
  beta_t <- rep(NA_real_, 20)
  beta_age <- rep(NA_real_, 20)
  err_sd <- rep(NA_real_, 20)

  beta0[1:3] <- c(0.80, -0.10, 0.40)
  beta_t[1:3] <- c(0.16, 0.10, -0.06)
  beta_age[1:3] <- c(0.20, -0.10, 0.08)
  err_sd[1:3] <- c(0.35, 0.30, 0.35)

  beta0[4:20] <- seq(0.65, -0.45, length.out = 17)
  beta_t[4:20] <- c(
    0.12, 0.08, -0.06, -0.09, -0.12, 0.10, 0.07,
    0.03, -0.04, 0.13, -0.10, 0.05, -0.07, 0.11,
    -0.05, 0.09, -0.02
  )
  beta_age[4:20] <- seq(0.16, -0.10, length.out = 17)
  err_sd[4:20] <- seq(0.30, 0.45, length.out = 17)

  list(beta0 = beta0, beta_t = beta_t, beta_age = beta_age, err_sd = err_sd)
}

sim20_re_spec <- function() {
  has_intercept <- rep(FALSE, 20)
  has_slope <- rep(FALSE, 20)
  has_intercept[1:5] <- TRUE
  has_slope[1:5] <- TRUE
  has_intercept[6:10] <- TRUE
  has_slope[11:15] <- TRUE

  sd_intercept <- rep(0, 20)
  sd_slope <- rep(0, 20)
  sd_intercept[1:3] <- c(0.85, 0.80, 0.75)
  sd_slope[1:3] <- c(0.24, 0.20, 0.22)
  sd_intercept[4:5] <- c(0.70, 0.65)
  sd_slope[4:5] <- c(0.18, 0.16)
  sd_intercept[6:10] <- seq(0.60, 0.80, length.out = 5)
  sd_slope[11:15] <- seq(0.14, 0.24, length.out = 5)

  list(
    has_intercept = has_intercept,
    has_slope = has_slope,
    sd_intercept = sd_intercept,
    sd_slope = sd_slope
  )
}

sim20_latent <- function(n, seed = NULL, keep_long = TRUE) {
  par <- sim20_parameters()
  fixed <- sim20_fixed_parameters()
  spec <- sim20_re_spec()
  if (!is.null(seed)) set.seed(seed)

  id <- seq_len(n)
  baseline.age <- pmin(
    pmax(rnorm(n, par$age_mean, par$age_sd), par$age_min),
    par$age_max
  )
  age_std <- (baseline.age - par$age_mean) / par$age_sd

  active_names <- character(0)
  active_sd <- numeric(0)
  active_marker <- integer(0)
  active_type <- character(0)

  for (j in seq_len(par$n_markers)) {
    if (spec$has_intercept[j]) {
      active_names <- c(active_names, paste0("b0_y", j))
      active_sd <- c(active_sd, spec$sd_intercept[j])
      active_marker <- c(active_marker, j)
      active_type <- c(active_type, "intercept")
    }
    if (spec$has_slope[j]) {
      active_names <- c(active_names, paste0("b1_y", j))
      active_sd <- c(active_sd, spec$sd_slope[j])
      active_marker <- c(active_marker, j)
      active_type <- c(active_type, "slope")
    }
  }

  R <- outer(seq_along(active_sd), seq_along(active_sd), function(i, j) par$re_rho^abs(i - j))
  D <- diag(active_sd)
  b_active <- MASS::mvrnorm(n, mu = rep(0, length(active_sd)), Sigma = D %*% R %*% D)
  colnames(b_active) <- active_names

  b0 <- matrix(0, nrow = n, ncol = 20)
  b1 <- matrix(0, nrow = n, ncol = 20)
  colnames(b0) <- paste0("b0_y", 1:20)
  colnames(b1) <- paste0("b1_y", 1:20)

  for (k in seq_along(active_names)) {
    j <- active_marker[k]
    if (active_type[k] == "intercept") b0[, j] <- b_active[, k]
    if (active_type[k] == "slope") b1[, j] <- b_active[, k]
  }

  mu <- par$mu0 +
    par$gamma_age * age_std +
    par$gamma_b10 * b0[, 1] +
    par$gamma_b11 * b1[, 1] +
    par$gamma_b20 * b0[, 2] +
    par$gamma_b31 * b1[, 3]

  true.time <- 1 + exp(mu + par$sigma * rnorm(n))

  basecov.data <- data.frame(id = id, baseline.age = baseline.age)
  if (!keep_long) {
    return(list(id = id, baseline.age = baseline.age, true.time = true.time, basecov.data = basecov.data))
  }

  long_list <- vector("list", n)
  for (i in seq_len(n)) {
    m_i <- sample(par$visit_min:par$visit_max, 1L)
    tt <- c(0, sort(runif(m_i - 1L, 0, par$visit_horizon)))
    Y <- matrix(NA_real_, nrow = m_i, ncol = 20)

    for (j in 1:20) {
      Y[, j] <- fixed$beta0[j] + fixed$beta_t[j] * tt + fixed$beta_age[j] * age_std[i] +
        b0[i, j] + b1[i, j] * tt + rnorm(m_i, 0, fixed$err_sd[j])
    }
    colnames(Y) <- paste0("y", 1:20)

    long_list[[i]] <- data.frame(
      id = id[i],
      t.from.base = tt,
      age = baseline.age[i] + tt,
      Y,
      check.names = FALSE
    )
  }

  list(
    id = id,
    baseline.age = baseline.age,
    true.time = true.time,
    basecov.data = basecov.data,
    long.full = do.call(rbind, long_list),
    random.effects = data.frame(id = id, b0, b1, check.names = FALSE)
  )
}

if (!exists("sim_censor_probability", mode = "function")) {
  sim_censor_probability <- function(true_time, cmax, cmin = 1) {
    if (cmax <= cmin) stop("cmax must be greater than cmin")
    p <- ifelse(
      true_time <= cmin,
      0,
      ifelse(true_time >= cmax, 1, (true_time - cmin) / (cmax - cmin))
    )
    pmin(pmax(p, 0), 1)
  }
}

if (!exists("sim_calibrate_cmax", mode = "function")) {
  sim_calibrate_cmax <- function(true_time, target, cmin = 1) {
    f <- function(cmax) mean(sim_censor_probability(true_time, cmax, cmin)) - target
    lower <- cmin + 1e-6
    upper <- max(as.numeric(stats::quantile(true_time, 0.999)), cmin + 1)
    while (f(upper) > 0) upper <- upper * 1.5
    stats::uniroot(f, interval = c(lower, upper), tol = 1e-7)$root
  }
}

sim20_censor_map <- function(targets = sim20_parameters()$censor_targets, pilot_n = 200000L, seed = 5002L) {
  latent <- sim20_latent(pilot_n, seed = seed, keep_long = FALSE)
  setNames(
    vapply(targets, function(target) sim_calibrate_cmax(latent$true.time, target), numeric(1)),
    as.character(targets)
  )
}

sim20_observe <- function(latent, target_censoring, cmax, seed = NULL) {
  par <- sim20_parameters()
  if (!is.null(seed)) set.seed(seed)

  C <- runif(length(latent$true.time), min = par$censor_min, max = cmax)
  time <- pmin(latent$true.time, C)
  event <- as.integer(latent$true.time <= C)

  surv.data <- data.frame(
    id = latent$id,
    baseline.age = latent$baseline.age,
    time = time,
    event = event,
    true.time = latent$true.time,
    censor.time = C
  )

  list(
    surv.data = surv.data,
    long.data = latent$long.full,
    basecov.data = latent$basecov.data,
    meta = list(
      scenario = "strongNPH",
      target_censoring = target_censoring,
      cmax = cmax,
      realized_censoring = mean(event == 0),
      n_markers = 20L,
      survival_depends_on = c("b0_y1", "b1_y1", "b0_y2", "b1_y3")
    )
  )
}
