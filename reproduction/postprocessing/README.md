# Post-processing checks for archived results

These scripts are provided for re-summarising already-computed result files from the original cluster analyses. They do not refit PRC models.

- `recompute_3marker_left_coverage.R` corrects the one-sided upper-endpoint coverage summary for the 3-marker simulations. An infinite upper endpoint represents an unbounded interval and therefore covers every finite event time.
- `reaggregate_adni_existing.R` keeps only conformal rows from the baseline-Cox files and reports both the historical finite-only upper-length summary and the alternative length obtained by truncating infinite upper endpoints at the largest observed training event time.

The scripts accept the original analysis root as a command-line argument and contain no personal paths.
