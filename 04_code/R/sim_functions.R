# Simulation of census-like data sets for the simulation study (protocol 10.1; ADEMP).
# Each data set mimics: products with several labelled side effects, crossed MedDRA terms,
# organ classes, one or two dose arms sharing one placebo arm, table inclusion rules that
# hide rows, and label selection (a row enters only if the observed drug rate is >= 1%).

simulate_census <- function(J = 150, alpha = -3, beta = 0.5, sd_all = 0.5, rho_row = 0,
                            p_rule_none = 0.6, K = 300, S = 20, seed = 1,
                            sizes = NULL, rows_per_product = NULL) {
  set.seed(seed)
  # term frequencies (Zipf-like) and term -> organ class map
  term_w <- 1 / (1:K)^0.8
  term_soc <- sample(1:S, K, replace = TRUE)
  # random effects
  mk_pair <- function(n, rho) {
    z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
    cbind(sd_all * z1, sd_all * z2)
  }
  re_prod <- mk_pair(J, 0)
  re_pt <- mk_pair(K, 0)
  c_s <- rnorm(S, 0, sd_all); s_s <- rnorm(S, 0, sd_all)

  rows <- list(); arms <- list(); r <- 0
  for (j in 1:J) {
    m <- if (is.null(rows_per_product)) 1 + rpois(1, 5) else sample(rows_per_product, 1)
    if (!is.null(sizes)) {
      # resample a pilot table structure: arm sizes and its row rule
      sz <- sizes[sample(nrow(sizes), 1), ]
      NP <- sz$NP; NDv <- sz$ND; n_arm <- 1
      rule <- if (sz$rule == "none" || runif(1) < p_rule_none) "none" else sz$rule
      cut <- if (is.na(sz$cut)) 5 else sz$cut
      ge <- if (is.na(sz$ge)) TRUE else sz$ge
    } else {
      ge <- TRUE
      n_arm <- if (runif(1) < 0.2) 2 else 1
      base <- round(exp(runif(1, log(80), log(1000))))
      NDv <- rep(base, n_arm)
      NP <- max(30, round(base * sample(c(1, 0.5), 1)))
      rule <- if (runif(1) < p_rule_none) "none" else
        sample(c("either", "drug", "drug_gt_pbo"), 1)
      cut <- sample(c(2, 5, 10), 1)
    }
    terms <- sample(1:K, m, prob = term_w)
    for (k in terms) {
      ew <- mk_pair(1, rho_row)
      s <- term_soc[k]
      mu <- alpha + re_prod[j, 1] + re_pt[k, 1] + c_s[s] + ew[1]
      de <- beta + re_prod[j, 2] + re_pt[k, 2] + s_s[s] + ew[2]
      pP <- plogis(mu); pD <- plogis(mu + de)
      yP <- rbinom(1, NP, pP); yD <- rbinom(n_arm, NDv, pD)
      # label selection: listed only if crude drug rate >= 1%
      if (sum(yD) / sum(NDv) < 0.01) next
      kP <- if (ge) ceiling(cut * NP / 100) - 1 else floor(cut * NP / 100)
      kD <- if (ge) ceiling(cut * NDv / 100) - 1 else floor(cut * NDv / 100)
      shown <- switch(rule,
        none = TRUE,
        either = any(yD > kD) || yP > kP,
        drug = any(yD > kD),
        drug_gt_pbo = any(yD > kD & yD / NDv > yP / NP))
      r <- r + 1
      if (shown) {
        type <- if (runif(1) < 0.05) 2 else 1
      } else {
        type <- switch(rule, either = 3, drug = 4, drug_gt_pbo = 5)
      }
      rows[[r]] <- data.frame(r = r, prod = j, pt = k, soc = s, type = type, NP = NP,
                              yP = if (type == 3) kP else yP,
                              yP_hi = yP, E_true = pD - pP, pD_true = pD)
      for (i in seq_along(NDv)) {
        arms[[length(arms) + 1]] <- data.frame(row = r, ND = NDv[i],
          yD = if (type %in% 3:5) kD[i] else yD[i], yD_hi = yD[i])
      }
    }
  }
  rows <- do.call(rbind, rows); arms <- do.call(rbind, arms)
  # interval rows: widen by +-1 around the truth (ambiguous back-calculation)
  i2 <- rows$type == 2
  rows$yP[i2] <- pmax(0, rows$yP_hi[i2] - 1); rows$yP_hi[i2] <- pmin(rows$NP[i2], rows$yP_hi[i2] + 1)
  a2 <- arms$row %in% rows$r[i2]
  arms$yD_hi[a2] <- pmin(arms$ND[a2], arms$yD_hi[a2] + 1); arms$yD[a2] <- pmax(0, arms$yD_hi[a2] - 2)
  # re-index products and terms that actually appear
  rows$prod <- as.integer(factor(rows$prod)); rows$pt <- as.integer(factor(rows$pt))
  rows$soc <- as.integer(factor(rows$soc))
  list(rows = rows, arms = arms,
       theta_true = mean(rows$E_true < 0.01),
       theta_mis_true = mean(rows$pD_true < 0.01),
       censored_share = mean(rows$type %in% 3:5))
}

# Build the Stan data list, including index sets and type-5 support vectors.
make_standata <- function(sim, prior_only = 0, ge = FALSE, sel = TRUE) {
  rows <- sim$rows; arms <- sim$arms
  R <- nrow(rows); A <- nrow(arms)
  # label selection (row listed only if pooled drug count >= 1% of pooled N). Applied to observed rows
  # and to hidden rows with one drug arm; kmin = 0 means no adjustment.
  Ntot <- as.numeric(tapply(arms$ND, factor(arms$row, levels = 1:R), sum))
  narm <- as.integer(table(factor(arms$row, levels = 1:R)))
  kmin <- if (sel) as.integer(ceiling(Ntot / 100 - 1e-9)) else integer(R)
  kmin[narm > 1 & rows$type != 1] <- 0L
  if (sel) {   # observed rows must meet the rule; otherwise they were selected on other data
    obs_ok <- rows$type != 1 | as.numeric(tapply(arms$yD, factor(arms$row, levels = 1:R), sum)) >= kmin
    kmin[!obs_ok] <- 0L
  }
  a2 <- which(rows$type[arms$row] == 2 & narm[arms$row] == 1)
  low2 <- arms$row[a2][arms$yD_hi[a2] < kmin[arms$row[a2]]]; kmin[low2] <- 0L   # interval wholly below 1%
  kmin[rows$type == 2 & narm > 1] <- 0L
  arms$yD[a2] <- pmax(arms$yD[a2], kmin[arms$row[a2]])
  if (any(arms$yD_hi[a2] < arms$yD[a2])) stop("type-2 interval below the selection threshold")
  rows1 <- which(rows$type == 1); arms1 <- which(rows$type[arms$row] == 1)
  rows5 <- which(rows$type == 5)
  lchP <- numeric(0); P_start <- integer(R)
  lchD <- numeric(0); D_start <- integer(A)
  n_arm5 <- integer(R); arm5 <- matrix(0L, R, 2)
  thr <- integer(0); t5_start <- integer(R)
  for (r in rows5) {
    a <- which(arms$row == r); NP <- rows$NP[r]
    if (length(a) > 2) stop("more than 2 dose arms in a type-5 row")
    n_arm5[r] <- length(a); arm5[r, seq_along(a)] <- a
    P_start[r] <- length(lchP) + 1; lchP <- c(lchP, lchoose(NP, 0:NP))
    for (i in a) { D_start[i] <- length(lchD) + 1; lchD <- c(lchD, lchoose(arms$ND[i], 0:arms$ND[i])) }
    t5_start[r] <- length(thr) + 1
    for (y in 0:NP) for (i in a) {
      lim <- if (ge) ceiling(y * arms$ND[i] / NP) - 1 else floor(y * arms$ND[i] / NP)
      thr <- c(thr, min(arms$ND[i], max(arms$yD[i], lim)))
    }
  }
  # hidden arms ("units"): placebo arms of type-3 rows and drug arms of type-3/4 rows
  hp <- which(rows$type == 3); hd <- which(rows$type[arms$row] %in% 3:4)
  uN <- c(rows$NP[hp], arms$ND[hd]); uk <- c(rows$yP[hp], arms$yD[hd]); ulin <- c(hp, R + arms$row[hd])
  ulo <- c(integer(length(hp)), kmin[arms$row[hd]])          # hidden drug arms: counts kmin..k
  if (any(ulo > uk)) stop("hidden row whose cutoff is below the selection threshold")
  ulen <- uk - ulo + 1; ustart <- if (length(ulen)) cumsum(c(1, head(ulen, -1))) else integer(0)
  c_y <- as.numeric(unlist(lapply(seq_along(uk), function(i) ulo[i]:uk[i]))); c_N <- as.numeric(rep(uN, ulen))
  # selection: rows with kmin <= 4 exact (quadrature in Stan); kmin >= 5 by a normal-CDF mixture
  sr <- which(kmin > 0); sa <- sr[kmin[sr] <= 4]; sb <- sr[kmin[sr] >= 5]
  a_lc <- outer(Ntot[sa], 1:3, lchoose); a_use <- outer(kmin[sa], 1:3, ">") * 1
  a_lc[a_use == 0] <- 0
  Q <- 12; jm <- matrix(0, Q, Q); for (k in 1:(Q - 1)) jm[k, k + 1] <- jm[k + 1, k] <- sqrt(k / 2)
  ev <- eigen(jm, symmetric = TRUE); gh_x <- ev$values; gh_lw <- log(ev$vectors[1, ]^2)   # weights / sqrt(pi)
  mx <- lapply(sb, function(r) sel_mix(Ntot[r], kmin[r]))
  b_la <- do.call(rbind, lapply(mx, function(z) log(z$a))); b_c <- do.call(rbind, lapply(mx, function(z) z$c))
  b_t2 <- do.call(rbind, lapply(mx, function(z) z$t^2))
  if (!length(sb)) b_la <- b_c <- b_t2 <- matrix(0, 0, 3)
  if (!length(sa)) a_lc <- a_use <- matrix(0, 0, 3)
  list(U = length(uN), u_start = as.array(as.integer(ustart)), u_len = as.array(as.integer(ulen)),
       TC = length(c_y), c_lch = as.array(lchoose(c_N, c_y)), c_y = as.array(as.numeric(c_y)), c_N = as.array(as.numeric(c_N)),
       c_lin = as.array(as.integer(rep(ulin, ulen))),
       R = R, J = max(rows$prod), K = max(rows$pt), S = max(rows$soc),
       prod = rows$prod, pt = rows$pt, soc = rows$soc, type = rows$type,
       NP = rows$NP, yP = rows$yP, yP_hi = rows$yP_hi,
       A = A, arm_row = arms$row, ND = arms$ND, yD = arms$yD, yD_hi = arms$yD_hi,
       R1 = length(rows1), rows1 = as.array(rows1), A1 = length(arms1), arms1 = as.array(arms1),
       R5 = length(rows5), rows5 = as.array(rows5),
       LP = length(lchP), lchP = as.array(lchP), P_start = P_start,
       LD = length(lchD), lchD = as.array(lchD), D_start = D_start,
       n_arm5 = n_arm5, arm5 = arm5,
       T5 = length(thr), thr = as.array(thr), t5_start = t5_start,
       kmin = kmin, SA = length(sa), a_row = as.array(sa), a_N = as.array(as.numeric(Ntot[sa])), a_lc = a_lc, a_use = a_use,
       Q = Q, gh_x = gh_x, gh_lw = gh_lw, SB = length(sb), b_row = as.array(sb), b_la = b_la, b_c = b_c, b_t2 = b_t2,
       prior_only = prior_only)
}

# P(Binomial(N, plogis(eta)) >= k) ~ sum_h a_h Phi((eta - c_h) / t_h), fitted on the log scale where
# log P > -14. Cached by (N, k); the global random-number state is left untouched.
.sel_cache <- new.env()
sel_mix <- function(N, k, H = 3) {
  key <- paste(N, k); if (!is.null(.sel_cache[[key]])) return(.sel_cache[[key]])
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv()) else NULL
  eta <- seq(qlogis(1e-6), qlogis(0.6), length.out = 400)
  ex <- pbinom(k - 1, N, plogis(eta), lower.tail = FALSE, log.p = TRUE); keep <- ex > -14; eta <- eta[keep]; ex <- ex[keep]
  f <- function(par) { a <- exp(par[1:H]) / sum(exp(par[1:H])); cc <- par[H + 1:H]; tt <- exp(par[2 * H + 1:H])
    sum((log(Reduce(`+`, lapply(1:H, function(i) a[i] * pnorm((eta - cc[i]) / tt[i])))) - ex)^2) }
  best <- NULL
  for (st in 1:5) { set.seed(st)
    p0 <- c(rep(0, H), qlogis(k / N) + rnorm(H, 0, 0.5), log(rep(0.3, H)) + rnorm(H, 0, 0.3))
    o <- optim(p0, f, method = "BFGS", control = list(maxit = 2000, reltol = 1e-14))
    if (is.null(best) || o$value < best$value) best <- o }
  if (!is.null(old)) assign(".Random.seed", old, globalenv())
  res <- list(a = exp(best$par[1:H]) / sum(exp(best$par[1:H])), c = best$par[H + 1:H], t = exp(best$par[2 * H + 1:H]))
  .sel_cache[[key]] <- res; res
}
