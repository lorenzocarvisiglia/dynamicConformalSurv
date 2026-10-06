# Sequential Monte Carlo driver for the simulation study.
#
# This script is intentionally simple and portable. It does not use a cluster,
# task array, manifest, or parallel backend.

source("reproducibility/simulations/run_simulation_cell.R")

run_simulation_study <- function(
  setting = c("three", "twenty"),
  scenarios = NULL,
  n_values = c(300L, 1000L, 1500L),
  censoring_values = c(0.20, 0.35, 0.50),
  landmarks = c(2, 4),
  n_rep = 1000L,
  n_validation = 10000L,
  alpha_values = 0.10,
  B = 500L,
  seed_offset = 20260000L,
  output_file = NULL,
  verbose = TRUE
) {
  setting <- match.arg(setting)

  if (is.null(scenarios)) {
    scenarios <- if (setting == "three") {
      c("weakNPH", "strongNPH")
    } else {
      "strongNPH"
    }
  }

  grid <- expand.grid(
    scenario = scenarios,
    N = n_values,
    censoring = censoring_values,
    landmark = landmarks,
    alpha = alpha_values,
    rep = seq_len(n_rep),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  rows <- vector("list", nrow(grid))

  # The original study used one censoring calibration and one large validation
  # sample for each scenario/censoring design, reused across Monte Carlo
  # training replications.
  cache <- new.env(parent = emptyenv())

  get_design_objects <- function(scenario, censoring) {
    key <- paste(setting, scenario, censoring, sep = "_")
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }

    design_seed <- seed_offset + match(scenario, unique(scenarios)) * 10000L +
      as.integer(round(1000 * censoring))

    if (setting == "three") {
      cmax <- calibrate_three_predictor_cmax(
        scenario = scenario,
        target_censoring = censoring,
        seed = design_seed
      )
      validation <- simulate_three_predictor_dataset(
        n = n_validation,
        scenario = scenario,
        target_censoring = censoring,
        cmax = cmax,
        seed = design_seed + 500000L
      )
    } else {
      cmax <- calibrate_twenty_predictor_cmax(
        target_censoring = censoring,
        seed = design_seed
      )
      validation <- simulate_twenty_predictor_dataset(
        n = n_validation,
        target_censoring = censoring,
        cmax = cmax,
        seed = design_seed + 500000L
      )
    }

    out <- list(cmax = cmax, validation = validation)
    assign(key, out, envir = cache)
    out
  }

  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]

    if (verbose) {
      message(
        "Simulation ", i, "/", nrow(grid),
        ": setting=", setting,
        ", scenario=", g$scenario,
        ", N=", g$N,
        ", censoring=", g$censoring,
        ", landmark=", g$landmark,
        ", alpha=", g$alpha,
        ", rep=", g$rep
      )
    }

    design <- get_design_objects(g$scenario, as.numeric(g$censoring))

    ans <- try(
      run_simulation_cell(
        setting = setting,
        scenario = g$scenario,
        n_train = as.integer(g$N),
        n_validation = as.integer(n_validation),
        target_censoring = as.numeric(g$censoring),
        landmark = as.numeric(g$landmark),
        alpha = as.numeric(g$alpha),
        B = as.integer(B),
        seed = as.integer(seed_offset + g$rep),
        cmax = design$cmax,
        validation_data = design$validation,
        verbose = FALSE
      ),
      silent = TRUE
    )

    if (inherits(ans, "try-error")) {
      rows[[i]] <- data.frame(
        setting = setting,
        scenario = g$scenario,
        N = g$N,
        censoring = g$censoring,
        landmark = g$landmark,
        alpha = g$alpha,
        rep = g$rep,
        method = NA_character_,
        error = as.character(ans)
      )
      next
    }

    z <- ans$metrics
    z$setting <- setting
    z$scenario <- g$scenario
    z$N <- g$N
    z$censoring <- g$censoring
    z$alpha <- g$alpha
    z$rep <- g$rep
    z$error <- NA_character_

    rows[[i]] <- z
  }

  out <- do.call(rbind, rows)

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(out, output_file, row.names = FALSE)
  }

  out
}

# Main 3-predictor design:
# res3 <- run_simulation_study(
#   setting = "three",
#   n_rep = 1000,
#   alpha_values = 0.10,
#   B = 500,
#   output_file = "results/simulation_three_predictor.csv"
# )
#
# Main 20-predictor design:
# res20 <- run_simulation_study(
#   setting = "twenty",
#   scenarios = "strongNPH",
#   n_rep = 1000,
#   alpha_values = 0.10,
#   B = 500,
#   output_file = "results/simulation_twenty_predictor.csv"
# )
#
# Alpha sensitivity used in the 20-predictor study:
# res_alpha <- run_simulation_study(
#   setting = "twenty",
#   scenarios = "strongNPH",
#   n_values = 1000,
#   censoring_values = 0.35,
#   landmarks = c(2, 4),
#   n_rep = 1000,
#   alpha_values = c(0.05, 0.10, 0.15, 0.20),
#   B = 500,
#   output_file = "results/alpha_sensitivity.csv"
# )
