# Registered simulation study (OSF yv6hf; protocol 10.1; ADEMP).
# Aims: bias, RMSE and 95% coverage of theta and theta_conv from the main model.
# Data-generating mechanism: census-like data sets; arm sizes and table rules from the pilot
#   products; side effects per product resampled from pilot label counts; label selection
#   (row kept only if observed drug rate >= 1%). Factors: alpha {-3.5,-2.5}, beta {0.3,0.8},
#   all SDs {0.5,1.0}, corr(e_r,w_r) {0,-0.5}, censored share {20%,40%}: 32 scenarios.
# Estimands: theta and theta_conv as realised in each data set (finite population).
# Methods: the registered Stan model with the label-selection term (decision log 73-74; rstan;
#   1 chain, 500 warm-up + 1000 sampling per fit).
# Performance: bias, RMSE, coverage of the 95% credible interval, with Monte Carlo SEs.
# Usage: Rscript sim_study.R <worker_id> <n_workers> <minutes_budget>

source("sim_functions.R"); source("init_fun.R"); source("pilot_templates.R")
suppressMessages(library(rstan)); rstan_options(auto_write = TRUE)

J_PER_DATASET <- 150     # lower end of the expected census; a harder test than the full census
N_REPS <- 200
out_dir <- "../../05_results/simulation"; dir.create(file.path(out_dir, "fits"), recursive = TRUE, showWarnings = FALSE)

scen <- expand.grid(alpha = c(-3.5, -2.5), beta = c(0.3, 0.8), sdv = c(0.5, 1.0),
                    rho = c(0, -0.5), cens = c(0.20, 0.40))
scen$id <- seq_len(nrow(scen))

# Rules: each product's table gets a rule with probability p_rule; the rule is drawn from the
# pilot rules (either-arm >=5%: Leqembi, Xgeva, Onerji; either-arm >10%: Beyonttra), with
# weight w10 on the >10% rule. Arm sizes always come from the six pilot templates.
sim_one <- function(s, seed, p_rule, w10) {
  tm <- pilot_templates
  set.seed(seed)
  tm$rule <- "none"
  sizes <- tm
  x <- simulate_census_rules(J = J_PER_DATASET, alpha = s$alpha, beta = s$beta, sd_all = s$sdv,
                             rho_row = s$rho, sizes = sizes, p_rule = p_rule, w10 = w10, seed = seed,
                             rows_per_product = as.integer(pilot_rows_per_product))
  x
}

simulate_census_rules <- function(J, alpha, beta, sd_all, rho_row, sizes, p_rule, w10, seed, rows_per_product,
                                  K = 300, S = 20, cut_hi = 10) {
  set.seed(seed)
  term_w <- 1 / (1:K)^0.8; term_soc <- sample(1:S, K, replace = TRUE)
  mk_pair <- function(n, rho) { z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n); cbind(sd_all * z1, sd_all * z2) }
  re_prod <- mk_pair(J, 0); re_pt <- mk_pair(K, 0); c_s <- rnorm(S, 0, sd_all); s_s <- rnorm(S, 0, sd_all)
  # product-level draws
  m <- sample(rows_per_product, J, replace = TRUE)
  tpl <- sample(nrow(sizes), J, replace = TRUE)
  NPj <- sizes$NP[tpl]; NDj <- sizes$ND[tpl]
  ruled <- runif(J) < p_rule
  is10 <- runif(J) < w10
  cut <- ifelse(is10, cut_hi, 5); ge <- !is10
  kPj <- ifelse(ge, ceiling(cut * NPj / 100) - 1, floor(cut * NPj / 100))
  kDj <- ifelse(ge, ceiling(cut * NDj / 100) - 1, floor(cut * NDj / 100))
  # row-level draws (all candidate rows at once)
  prod <- rep(1:J, m); n <- length(prod)
  pt <- unlist(lapply(m, function(mm) sample(1:K, mm, prob = term_w)))
  soc <- term_soc[pt]
  ew <- mk_pair(n, rho_row)
  mu <- alpha + re_prod[prod, 1] + re_pt[pt, 1] + c_s[soc] + ew[, 1]
  de <- beta + re_prod[prod, 2] + re_pt[pt, 2] + s_s[soc] + ew[, 2]
  pP <- plogis(mu); pD <- plogis(mu + de)
  NP <- NPj[prod]; ND <- NDj[prod]
  yP <- rbinom(n, NP, pP); yD <- rbinom(n, ND, pD)
  keep <- yD / ND >= 0.01                                   # label selection
  shown <- !ruled[prod] | yD > kDj[prod] | yP > kPj[prod]
  u2 <- runif(n) < 0.05
  type <- ifelse(shown, ifelse(u2, 2, 1), 3)
  idx <- which(keep); R <- length(idx)
  rows <- data.frame(r = seq_len(R), prod = prod[idx], pt = pt[idx], soc = soc[idx], type = type[idx], NP = NP[idx],
                     yP = ifelse(type[idx] == 3, kPj[prod[idx]], yP[idx]), yP_hi = yP[idx],
                     E_true = (pD - pP)[idx], pD_true = pD[idx])
  arms <- data.frame(row = seq_len(R), ND = ND[idx], yD = ifelse(type[idx] == 3, kDj[prod[idx]], yD[idx]), yD_hi = yD[idx])
  i2 <- rows$type == 2
  rows$yP[i2] <- pmax(0, rows$yP_hi[i2] - 1); rows$yP_hi[i2] <- pmin(rows$NP[i2], rows$yP_hi[i2] + 1)
  arms$yD_hi[i2] <- pmin(arms$ND[i2], arms$yD_hi[i2] + 1); arms$yD[i2] <- pmax(0, arms$yD_hi[i2] - 2)
  for (v in c("prod", "pt", "soc")) rows[[v]] <- as.integer(factor(rows[[v]]))
  list(rows = rows, arms = arms, theta_true = mean(rows$E_true < 0.01),
       theta_mis_true = mean(rows$pD_true < 0.01), censored_share = mean(rows$type == 3))
}

# Calibrate (p_rule, w10) per scenario so the expected censored share matches the target.
calibrate <- function(s) {
  share <- function(p, w, ch = 10) mean(sapply(1:20, function(i)
    simulate_census_rules(J_PER_DATASET, s$alpha, s$beta, s$sdv, s$rho, pilot_templates, p, w,
                          seed = 9000 + i, rows_per_product = as.integer(pilot_rows_per_product), cut_hi = ch)$censored_share))
  w10 <- 0.25; cut_hi <- 10                           # pilot mix: 1 of 4 ruled tables uses ">10%"
  if (share(1, w10) < s$cens) {                       # step 1: raise the share of tables with the strict pilot rule
    for (w in seq(0.25, 1, by = 0.05)) { w10 <- w; if (share(1, w) >= s$cens) break }
  }
  if (share(1, w10, cut_hi) < s$cens) {               # step 2 (stress test): raise the strict cutoff above 10%
    for (ch in 11:30) { cut_hi <- ch; if (share(1, 1, ch) >= s$cens) break }
  }
  lo <- 0; hi <- 1
  for (it in 1:8) { mid <- (lo + hi) / 2; if (share(mid, w10, cut_hi) < s$cens) lo <- mid else hi <- mid }
  c(p_rule = hi, w10 = w10, cut_hi = cut_hi, achieved = share(hi, w10, cut_hi))
}

calib_file <- file.path(out_dir, "calibration.csv")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) && args[1] == "calibrate") {
  cal <- t(sapply(seq_len(nrow(scen)), function(i) calibrate(scen[i, ])))
  write.csv(cbind(scen, cal), calib_file, row.names = FALSE); print(cbind(scen, round(cal, 3))); quit(save = "no")
}

worker <- as.integer(args[1]); n_workers <- as.integer(args[2]); budget_min <- as.numeric(args[3])
# Job scheduling only (added 2026-11-02 so a second computer can share the work; the method is unchanged).
# "fwd" (default) works through replications 1, 2, 3, ...; "rev" works through each block of 50
# replications from its end (50, 49, ..., then 100, 99, ...). Both skip fits that already exist,
# so two computers meet in the middle of each block.
direction <- if (length(args) >= 4) args[4] else "fwd"
if (direction == "rev") {
  lg <- file(file.path(out_dir, sprintf("pc_worker%02d.log", worker)), open = "at"); sink(lg); sink(lg, type = "message")
}
cal <- read.csv(calib_file)
mod <- stan_model("model_census.stan")
t_start <- Sys.time()
jobs <- expand.grid(rep = 1:N_REPS, sid = cal$id)
jobs <- if (direction == "rev") jobs[order(ceiling(jobs$rep / 50), -jobs$rep, jobs$sid), ] else
  jobs[order(jobs$rep, jobs$sid), ]                  # interleave scenarios so partial results are balanced
# "pending" (added 2026-11-03 for many machines run in rounds): drop finished fits first, then split.
if (direction == "pending") {
  done <- file.exists(file.path(out_dir, "fits", sprintf("s%02d_r%03d.csv", jobs$sid, jobs$rep)))
  jobs <- jobs[!done, ]
}
jobs <- jobs[seq(worker, nrow(jobs), by = n_workers), ]
jobs <- jobs[!is.na(jobs$rep), ]
for (q in seq_len(nrow(jobs))) {
  if (as.numeric(difftime(Sys.time(), t_start, units = "mins")) > budget_min) break
  sid <- jobs$sid[q]; rep <- jobs$rep[q]
  f <- file.path(out_dir, "fits", sprintf("s%02d_r%03d.csv", sid, rep))
  if (file.exists(f)) next
  s <- cal[cal$id == sid, ]; seed <- sid * 100000 + rep
  x <- simulate_census_rules(J_PER_DATASET, s$alpha, s$beta, s$sdv, s$rho, pilot_templates, s$p_rule, s$w10,
                             seed = seed, rows_per_product = as.integer(pilot_rows_per_product), cut_hi = s$cut_hi)
  sd <- make_standata(x); t0 <- Sys.time()
  fit <- suppressWarnings(sampling(mod, data = sd, chains = 1, iter = 1500, warmup = 500, refresh = 0,
                                   seed = seed, init = list(make_init(sd)())))
  th <- rstan::extract(fit, "theta")$theta; tc <- rstan::extract(fit, "theta_conv")$theta_conv
  sm <- summary(fit, pars = c("alpha", "beta", "sd_prod", "sd_pt", "sd_row", "sd_c", "sd_s", "theta"))$summary
  sp <- get_sampler_params(fit, inc_warmup = FALSE)[[1]]
  theta_conv_true <- x$theta_true - x$theta_mis_true
  res <- data.frame(sid = sid, rep = rep, seed = seed, rows = nrow(x$rows), censored = x$censored_share,
    theta_true = x$theta_true, theta_mean = mean(th), theta_lo = quantile(th, .025), theta_hi = quantile(th, .975),
    conv_true = theta_conv_true, conv_mean = mean(tc), conv_lo = quantile(tc, .025), conv_hi = quantile(tc, .975),
    max_rhat = max(sm[, "Rhat"], na.rm = TRUE), min_ess = min(sm[, "n_eff"], na.rm = TRUE),
    divergent = sum(sp[, "divergent__"]), minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  write.csv(res, paste0(f, ".tmp"), row.names = FALSE); file.rename(paste0(f, ".tmp"), f)
  cat(format(Sys.time()), basename(f), "done\n")
}
cat("worker", worker, "stopped at", format(Sys.time()), "\n")
