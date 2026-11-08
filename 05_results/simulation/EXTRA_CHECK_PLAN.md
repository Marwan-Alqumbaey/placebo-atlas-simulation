# Extra simulation check: theta for fully observed rows

Written 8 November 2026, after the registered simulation was complete and before any of the extra fits were run. No census outcome data had been extracted or analysed.

## Why we added it

The registered simulation (OSF yv6hf) reported bias, root mean squared error and coverage of theta and theta_conv. Coverage of theta was below 90% in 8 of 32 scenarios. Under the registered rule, theta for fully observed rows therefore becomes the primary analysis. The registered simulation did not measure the performance of that estimate. This check measures it, so that the primary estimate has its own evidence. It is an addition to the registered plan and will be reported as such.

## What stays the same

The model file (model_census.stan, unchanged; see CODE_SHA256_v2.txt), the data-generating function (read from sim_study.R, not copied), the 32 scenarios, the calibration, the seeds (scenario x 100000 + replication), the sampler settings (1 chain, 500 warm-up and 1,000 sampling iterations) and 200 data sets per scenario. Each data set is therefore the same as in the registered simulation.

## Estimands

- theta_obs: in each data set, the share of fully observed rows (type 1) whose true excess E is below 1%. Its estimate is the posterior of the same share computed from the full model fit.
- theta_cens: the same for censored rows (type 3), reported for information.
- theta and theta_conv are kept again to check that the refits reproduce the registered results.

## Performance measures

Bias, root mean squared error and coverage of the 95% credible interval for theta_obs and theta_cens, by scenario, with Monte Carlo standard errors. Convergence is summarised as in the registered simulation (R-hat, effective sample size, divergent transitions).

## How we will read the result

- If coverage of theta_obs is at least 90% in every scenario, theta for fully observed rows is supported as the primary estimate.
- If it is below 90% in any scenario, we keep theta for fully observed rows as the primary estimate (as registered), report this as a limitation, and give extra weight to the registered observed analysis of fully observed rows (Newcombe intervals, no model) in the interpretation.
- If the refits do not reproduce the registered theta (posterior means differing by more than 0.01 on average), we report this and investigate before using either set.
