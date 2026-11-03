# Placebo atlas: simulation study

Code for the registered simulation study (OSF registration [osf.io/yv6hf](https://osf.io/yv6hf)).

The study checks how well the main model estimates the share of labelled side effects that are
not more common on the medicine than on placebo. It uses simulated data only. No real data are stored here.

- `04_code/R/model_census.stan`: the model (Stan).
- `04_code/R/sim_functions.R`, `sim_study.R`: data generation, fitting and saving.
- `04_code/R/pilot_templates.R`: table sizes and rules taken from the pilot products.
- `05_results/simulation/calibration.csv`: settings for the 32 scenarios.
- `05_results/simulation/fits/`: one file per fitted data set (32 scenarios x 200 data sets).

The fits run on GitHub Actions (`.github/workflows/simulation.yml`). Each fit is seeded,
so any single fit can be reproduced with `Rscript sim_study.R`.

Authors: Marwan Alqumbaey and Yousef M. Nasher, Faculty of Medicine and Health Sciences,
Sana'a University, Yemen.
