# Extra simulation check (added 2026-11-08, after the registered simulation was complete).
# Why: the registered simulation found coverage of theta below 90% in 8 of 32 scenarios, so by the
#   registered rule "theta for fully observed rows" becomes the primary estimate. The registered
#   simulation did not measure that estimate. This check measures it with the same data sets.
# What is unchanged: the model (model_census.stan, same file and hash), the data-generating code
#   (simulate_census_rules is read from sim_study.R, not copied), seeds, scenarios, calibration,
#   sampler settings and 200 data sets per scenario.
# What is added: from each fit we keep the posterior of
#   theta_obs  = share of fully observed rows (type 1) with E < 1%
#   theta_cens = share of censored rows (type 3) with E < 1%
#   and their true values in the same data set. theta and theta_conv are kept again as a
#   reproducibility check against 05_results/simulation/fits.
# Plan written before any of these fits: 05_results/simulation/EXTRA_CHECK_PLAN.md
# Usage: Rscript sim_obs.R <worker_id> <n_workers> <minutes_budget>

source("sim_functions.R"); source("init_fun.R"); source("pilot_templates.R")
suppressMessages(library(rstan)); rstan_options(auto_write = TRUE)

src <- readLines("sim_study.R")                       # use the frozen data-generating function as it is
a <- grep("^simulate_census_rules <- function", src); b <- grep("^# Calibrate", src) - 1
eval(parse(text = src[a:b]))
J_PER_DATASET <- 150; N_REPS <- 200

out_dir <- "../../05_results/simulation"
fit_dir <- file.path(out_dir, "fits_obs"); dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
cal <- read.csv(file.path(out_dir, "calibration.csv"))

args <- commandArgs(trailingOnly = TRUE)
worker <- as.integer(args[1]); n_workers <- as.integer(args[2]); budget_min <- as.numeric(args[3])
mod <- stan_model("model_census.stan")
t_start <- Sys.time()
jobs <- expand.grid(rep = 1:N_REPS, sid = cal$id)
jobs <- jobs[order(jobs$rep, jobs$sid), ]
done <- file.exists(file.path(fit_dir, sprintf("s%02d_r%03d.csv", jobs$sid, jobs$rep)))
jobs <- jobs[!done, ]
jobs <- jobs[seq(worker, nrow(jobs), by = n_workers), ]
jobs <- jobs[!is.na(jobs$rep), ]
ci <- function(v) c(mean(v), quantile(v, .025), quantile(v, .975))
for (q in seq_len(nrow(jobs))) {
  if (as.numeric(difftime(Sys.time(), t_start, units = "mins")) > budget_min) break
  sid <- jobs$sid[q]; rep <- jobs$rep[q]
  f <- file.path(fit_dir, sprintf("s%02d_r%03d.csv", sid, rep))
  if (file.exists(f)) next
  s <- cal[cal$id == sid, ]; seed <- sid * 100000 + rep
  x <- simulate_census_rules(J_PER_DATASET, s$alpha, s$beta, s$sdv, s$rho, pilot_templates, s$p_rule, s$w10,
                             seed = seed, rows_per_product = as.integer(pilot_rows_per_product), cut_hi = s$cut_hi)
  sd <- make_standata(x); t0 <- Sys.time()
  fit <- suppressWarnings(sampling(mod, data = sd, chains = 1, iter = 1500, warmup = 500, refresh = 0,
                                   seed = seed, init = list(make_init(sd)())))
  E <- rstan::extract(fit, "E")$E                      # draws x rows
  th <- rstan::extract(fit, "theta")$theta; tc <- rstan::extract(fit, "theta_conv")$theta_conv
  o <- x$rows$type == 1; cz <- x$rows$type == 3
  to <- rowMeans(E[, o, drop = FALSE] < 0.01)
  tz <- if (any(cz)) rowMeans(E[, cz, drop = FALSE] < 0.01) else rep(NA_real_, nrow(E))
  sm <- summary(fit, pars = c("alpha", "beta", "sd_prod", "sd_pt", "sd_row", "sd_c", "sd_s", "theta"))$summary
  sp <- get_sampler_params(fit, inc_warmup = FALSE)[[1]]
  r <- function(v) unname(round(v, 6))
  res <- data.frame(sid = sid, rep = rep, seed = seed, rows = nrow(x$rows), n_obs = sum(o), n_cens = sum(cz),
    obs_true = mean(x$rows$E_true[o] < 0.01), obs_mean = r(ci(to)[1]), obs_lo = r(ci(to)[2]), obs_hi = r(ci(to)[3]),
    cens_true = if (any(cz)) mean(x$rows$E_true[cz] < 0.01) else NA, cens_mean = r(ci(tz)[1]), cens_lo = r(ci(tz)[2]), cens_hi = r(ci(tz)[3]),
    theta_true = x$theta_true, theta_mean = r(mean(th)), theta_lo = r(quantile(th, .025)), theta_hi = r(quantile(th, .975)),
    conv_mean = r(mean(tc)),
    max_rhat = max(sm[, "Rhat"], na.rm = TRUE), min_ess = min(sm[, "n_eff"], na.rm = TRUE),
    divergent = sum(sp[, "divergent__"]), minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  write.csv(res, paste0(f, ".tmp"), row.names = FALSE); file.rename(paste0(f, ".tmp"), f)
  cat(format(Sys.time()), basename(f), "done\n")
}
cat("worker", worker, "stopped at", format(Sys.time()), "\n")
