# ADNI cross-validation

The ADNI participant-level data are not distributed in this repository.

The script run_adni_cv.R expects a prepared RData file containing the same objects used in the original analysis:

- surv: one row per subject, including id, time.to.event, event, and the baseline covariates;
- long.transformed: one row per subject-visit, including id, time.fup, and the longitudinal markers;
- all.fixed.cov.names: character vector naming the baseline covariates;
- all.long.cov.names: character vector naming the longitudinal markers.

The public script is sequential and does not use SLURM, task files, or cluster-specific paths.

Example:

source("reproducibility/adni/run_adni_cv.R")

results <- run_adni_cv(
  data_file = "path/to/adni2023.RData",
  landmarks = c(2, 4),
  K = 5,
  B = 500,
  alpha = 0.10,
  output_dir = "results/adni"
)

The script reproduces the final comparison used in the chapter: Naive PRC inversion, Conformal PRC, and Conformal baseline Cox.

The observed-event coverage quantities are descriptive diagnostics and do not estimate population-level predictive coverage under censoring.
