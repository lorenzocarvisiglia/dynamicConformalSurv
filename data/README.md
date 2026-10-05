# Data

No participant-level study data or large simulation datasets are stored in this repository.

For the reusable method, users provide:

1. a training survival data frame with one row per subject;
2. a training longitudinal data frame with one row per subject-visit;
3. baseline covariates for the subjects to be predicted;
4. longitudinal measurements observed no later than the chosen landmark.

Training survival data must include the observed event/censoring time and event indicator. Future outcomes are not required for subjects being predicted.

The ADNI data used in the accompanying analysis are not redistributed here because access is governed by the Alzheimer's Disease Neuroimaging Initiative.

See `examples/example_toy_data.R` for the expected structure.
