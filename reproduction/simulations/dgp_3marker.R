# Data-generating mechanism for the 3-marker simulation study.
#
# This public implementation follows the final thesis specification and the
# metadata stored with the simulation datasets used for the reported results.

if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required.")
}

sim3_parameters <- function() {
  list(
    landmarks_descriptive = c(2, 2.5, 3, 4, 4.5),
    landmarks_methods = c(2, 4),
    censor_targets = c(0.20, 0.35, 0.50),
    n_values = c(300L, 1000L, 1500L),
    n_validation = 10000L,
    visit_min = 4L,
    visit_max = 9L,
    visit_horizon = 5,
    age_mean = 60,
    age_sd = 10,
    age_min = 40,
    age_max = 90,
    censor_min = 1,
    re_sd = c(1.0, 0.20, 6.0, 0.80, 3.5, 0.55),
    corr_within = -0.20,
    scenarios = list(
      weakNPH = list(
        sigma = 0.50,
        beta_age = 0.18,
        beta_y1_b0 = 0.20,
        beta_y1_b1 = 0.18,
        beta_y2_b0 = 0.22,
        beta_y3_b1 = -0.20
      ),
      strongNPH = list(
        sigma = 0.40,
        beta_age = 0.30,
        beta_y1_b0 = 0.40,
        beta_y1_b1 = 0.36,
        beta_y2_b0 = 0.42,
        beta_y3_b1 = -0.38
      )
    )
  )
}

sim3_re_covariance <- function(par = sim3_parameters()) {
  R <- diag(6)
  for (idx in list(1:2, 3:4, 5:6)) {
    R[idx[1], idx[2]] <- par$corr_within
    R[idx[2], idx[1]] <- par$corr_within
  }
  D <- diag(par$re_sd)
  D %*% R %*% D
}

sim3_latent <- function(n, scenario = c("weakNPH", "strongNPH"), seed = NULL) {
  scenario <- match.arg(scenario)
  par <- sim3_parameters()
  sc <- par$scenarios[[scenario]]

  if (!is.null(seed)) set.seed(seed)

  id <- seq_len(n)
  baseline.age <- pmin(
    pmax(rnorm(n, par$age_mean, par$age_sd), par$age_min),
    par$age_max
  )

  b <- MASS::mvrnorm(
    n = n,
    mu = rep(0, 6),
    Sigma = sim3_re_covariance(par)
  )
  colnames(b) <- c("y1_b0", "y1_b1", "y2_b0", "y2_b1", "y3_b0", "y3_b1")

  eta <- sc$beta_age * (baseline.age - 60) / 10 +
    sc$beta_y1_b0 * b[, "y1_b0"] / par$re_sd[1] +
    sc$beta_y1_b1 * b[, "y1_b1"] / par$re_sd[2] +
    sc$beta_y2_b0 * b[, "y2_b0"] / par$re_sd[3] +
    sc$beta_y3_b1 * b[, "y3_b1"] / par$re_sd[6]

  true.time <- 1 + exp(log(5) + eta + sc$sigma * rnorm(n))

  long_list <- vector("list", n)
  for (i in seq_len(n)) {
    m_i <- sample(par$visit_min:par$visit_max, 1L)
    tt <- c(0, sort(runif(m_i - 1L, 0, par$visit_horizon)))

    long_list[[i]] <- data.frame(
      id = id[i],
      baseline.age = baseline.age[i],
      t.from.base = tt,
      age = baseline.age[i] + tt,
      y1 = 10 - 0.30 * tt + b[i, "y1_b0"] + b[i, "y1_b1"] * tt + rnorm(m_i, 0, 0.8),
      y2 = 100 + 1.80 * tt + b[i, "y2_b0"] + b[i, "y2_b1"] * tt + rnorm(m_i, 0, 5.0),
      y3 = 50 - 1.20 * tt + b[i, "y3_b0"] + b[i, "y3_b1"] * tt + rnorm(m_i, 0, 2.5),
      stringsAsFactors = FALSE
    )
  }

  list(
    id = id,
    baseline.age = baseline.age,
    random.effects = b,
    true.time = true.time,
    long.full = do.call(rbind, long_list),
    scenario = scenario
  )
}

sim_censor_probability <- function(true_time, cmax, cmin = 1) {
  if (cmax <= cmin) stop("cmax must be greater than cmin")
  p <- ifelse(
    true_time <= cmin,
    0,
    ifelse(true_time >= cmax, 1, (true_time - cmin) / (cmax - cmin))
  )
  pmin(pmax(p, 0), 1)
}

sim_calibrate_cmax <- function(true_time, target, cmin = 1) {
  f <- function(cmax) mean(sim_censor_probability(true_time, cmax, cmin)) - target
  lower <- cmin + 1e-6
  upper <- max(as.numeric(stats::quantile(true_time, 0.999)), cmin + 1)
  while (f(upper) > 0) upper <- upper * 1.5
  stats::uniroot(f, interval = c(lower, upper), tol = 1e-7)$root
}

sim3_censor_map <- function(
  targets = sim3_parameters()$censor_targets,
  pilot_n = 200000L,
  seed = 5010L
) {
  out <- list()
  for (scenario in names(sim3_parameters()$scenarios)) {
    latent <- sim3_latent(pilot_n, scenario = scenario, seed = seed + match(scenario, names(sim3_parameters()$scenarios)))
    out[[scenario]] <- setNames(
      vapply(targets, function(target) sim_calibrate_cmax(latent$true.time, target), numeric(1)),
      as.character(targets)
    )
  }
  out
}

sim3_observe <- function(latent, target_censoring, cmax, seed = NULL) {
  par <- sim3_parameters()
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
    censor.time = C,
    stringsAsFactors = FALSE
  )

  long.data <- latent$long.full[
    latent$long.full$t.from.base <= time[match(latent$long.full$id, latent$id)],
    ,
    drop = FALSE
  ]

  list(
    surv.data = surv.data,
    long.data = long.data,
    basecov.data = surv.data[, c("id", "baseline.age"), drop = FALSE],
    meta = list(
      scenario = latent$scenario,
      target_censoring = target_censoring,
      cmax = cmax,
      realized_censoring = mean(event == 0)
    )
  )
}
