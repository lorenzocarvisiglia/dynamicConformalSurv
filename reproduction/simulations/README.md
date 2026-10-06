# Simulation reproduction

These scripts provide serial, user-independent versions of the simulation designs used in Chapter 4. They do not contain cluster paths, SLURM job arrays, or parallel execution.

## Main study

The main simulation settings are:

- 3-marker weak NPH and strong NPH scenarios;
- 20-marker strong NPH scenario;
- training sample sizes 300, 1000, and 1500;
- censoring targets 20%, 35%, and 50% before landmarking;
- landmarks 2 and 4 for interval evaluation;
- validation sample size 10,000;
- 1,000 Monte Carlo replications;
- 500 bootstrap calibration replicates;
- nominal coverage 90% in the main tables.

For a small local run:

```bash
Rscript reproduction/simulations/run_simulation.R \
  --setting 3 --scenario weakNPH --N 300 --censor 0.20 \
  --landmark 2 --reps 2 --B 20 --validation-n 1000
```

For the full manuscript settings, set `--reps 1000 --B 500 --validation-n 10000` and run each design cell serially.

`dgp_3marker.R` follows the final thesis equations and the metadata stored with the datasets used in the study. `dgp_20marker.R` is a cleaned version of the final Scenario C generator.

## Alpha sensitivity

```bash
Rscript reproduction/simulations/run_alpha_sensitivity.R \
  --N 1000 --censor 0.35 --reps 1000 --B 500 --validation-n 10000
```

The alpha-sensitivity analysis uses `alpha = 0.05, 0.10, 0.15, 0.20` and landmarks 2 and 4.

## Note on historical post-processing

The reusable code treats an infinite upper endpoint as an unbounded upper interval, so any finite event time is covered by that one-sided interval. This is the convention implemented in `R/coverage_one_run.R` and should be used for new analyses.
