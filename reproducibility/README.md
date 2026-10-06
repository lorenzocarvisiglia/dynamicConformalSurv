# Reproducibility scripts

This directory contains portable, sequential R scripts corresponding to the simulation and ADNI analyses reported with the dynamic conformal survival method.

The original computations were run on a computing cluster. The public scripts here are portable reconstructions of the final analyses: cluster-specific paths, task arrays, manifests, and parallel execution have been removed. The statistical specifications are preserved. For the 3-predictor generator, the final saved datasets and the thesis specification are the reference because a later helper script retained on the cluster had already been modified after those datasets were generated. The 20-predictor generator was reconstructed from the final June 2026 generator and checked against the saved validation data.

## Structure

reproducibility/
- simulations/simulate_three_predictor.R
- simulations/simulate_twenty_predictor.R
- simulations/run_simulation_cell.R
- simulations/run_simulation_study.R
- adni/README.md
- adni/run_adni_cv.R

## Simulation study

The 3-predictor data-generating mechanism is the final design reported in the thesis appendix. It uses baseline age truncated to 40--90 years, 4--9 visits over five years, three longitudinal markers with marker-specific random intercepts and slopes, within-marker random intercept/slope correlation -0.20, shifted-lognormal event times, weak and strong non-proportional-hazards settings, and independent uniform censoring.

The 20-predictor generator is a separate strong-NPH design. It uses baseline age truncated to 40--85 years, 20 longitudinal markers with heterogeneous random-effects structures, and survival depending directly on selected random effects from the first three markers.

The main simulation settings reported in the thesis are training sample sizes 300, 1000, and 1500; target censoring 20%, 35%, and 50%; landmarks 2 and 4; validation sample size 10,000; 1,000 Monte Carlo replications; nominal coverage 90%; and 500 bootstrap calibration replicates in the final analyses.

A single cell can be run with run_simulation_cell(). The full Monte Carlo grid, including the 20-predictor alpha-sensitivity settings 0.05, 0.10, 0.15, and 0.20, can be run with run_simulation_study(). Running the full design sequentially is computationally intensive.

## ADNI analysis

The ADNI participant-level data are not distributed. The public ADNI script expects the same prepared R workspace used in the analysis, with objects described in adni/README.md.

The cross-validation script runs sequentially and does not require a scheduler.

## Reproducibility note

Because the public repository intentionally excludes private data, large generated simulation datasets, and cluster execution files, the scripts are designed to regenerate the statistical analyses from inputs rather than reproduce the original filesystem layout.
