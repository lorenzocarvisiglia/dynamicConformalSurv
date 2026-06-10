# Data

This folder is reserved for small example datasets.

The scripts in this repository do not include the ADNI data or the simulation data used in the paper.

Users should provide already-cleaned data in two data frames:

1. a survival data frame, with one row per subject;
2. a longitudinal data frame, with one row per subject-visit.

The survival data frame should contain at least:

- subject identifier;
- observed survival or censoring time;
- event indicator;
- baseline covariates.

The longitudinal data frame should contain at least:

- subject identifier;
- visit time;
- longitudinal markers observed up to follow-up.

See `examples/example_toy_data.R` for a minimal example.
