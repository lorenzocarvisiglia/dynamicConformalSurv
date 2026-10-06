# ADNI analysis reproduction

The ADNI participant-level data are not distributed in this repository. Access must be obtained separately from the Alzheimer's Disease Neuroimaging Initiative.

`run_adni_cv.R` is a serial, user-independent version of the final cross-validation analysis. It compares naive PRC, conformal PRC, and conformal baseline Cox using the same subject-level five-fold split at landmarks 2, 3, and 4. The default calibration size is `B = 500`.

The script expects the processed `adni2023.RData` object used in the analysis. It must contain:

- `surv`, with subject ID, `time.to.event`, `event`, and the baseline covariates;
- `long.transformed`, with subject ID, `time.fup`, and the longitudinal markers;
- `all.fixed.cov.names`, the names of the baseline covariates;
- `all.long.cov.names`, the names of the longitudinal markers.

Run:

```bash
Rscript reproduction/adni/run_adni_cv.R \
  --data /path/to/adni2023.RData \
  --output results/adni_cv \
  --B 500
```

The output reports both `upper_length_finite`, the historical table convention that averages upper one-sided length over finite upper endpoints, and `upper_length_truncated`, where infinite endpoints are truncated at the largest observed training event time.

The public script also avoids an old development aggregation error in which naive and conformal baseline-Cox rows could be pooled under the same label. Only conformal baseline Cox is used in the three-method comparison.
