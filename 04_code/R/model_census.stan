// Hierarchical binomial model for drug vs placebo side-effect rates, with censoring
// implied by table inclusion rules (protocol section 10.1; OSF registration osf.io/yv6hf).
//
// Row r = one labelled side effect of one product. Each row has one placebo arm and
// one or more drug (dose) arms sharing one drug rate.
// Row types:
//   1 fully observed counts
//   2 interval-censored counts (partial matches, ambiguous back-calculation)
//   3 hidden by a cutoff in either arm
//   4 hidden by a cutoff in the drug arm
//   5 hidden by "cutoff and drug > placebo" (or ">=", encoded in thr); shown if ANY dose arm passes
// Label selection: a side effect is listed as common only if its observed drug rate is at least 1%.
// Rows enter the data only if selected, so the population model for the row effects is truncated:
// for each row with kmin > 0, target += -log P(selected | m_mu, m_delta, Sigma_row), where
// P(selected) = integral of P(pooled drug count >= kmin | eta) over eta ~ N(m_mu + m_delta, s_eta^2),
// Rows with kmin <= 4: exact binomial tail, integrated by Gauss-Hermite quadrature.
// Rows with kmin >= 5: P(count >= kmin | eta) is approximated by a 3-part mixture of normal CDFs in eta
// (fitted in R for each pooled N), which makes the integral over eta closed form.
// Hidden single-arm rows also know their count is >= kmin.
// kmin = 0 switches the adjustment off for that row (the earlier model).
// The likelihood is exact. For type 5, P(not shown) = sum_yP Bin(yP) * prod_i F_i(thr_i(yP)),
// computed from full probability vectors (closed-form binomial pmf with precomputed log choose).

data {
  int<lower=1> R;
  int<lower=1> J;
  int<lower=1> K;
  int<lower=1> S;
  array[R] int<lower=1, upper=J> prod;
  array[R] int<lower=1, upper=K> pt;
  array[R] int<lower=1, upper=S> soc;
  array[R] int<lower=1, upper=5> type;

  array[R] int<lower=1> NP;
  array[R] int<lower=0> yP;      // type 1 count; type 2 lower; type 3 k_P
  array[R] int<lower=0> yP_hi;   // type 2 upper

  int<lower=1> A;
  array[A] int<lower=1, upper=R> arm_row;
  array[A] int<lower=1> ND;
  array[A] int<lower=0> yD;      // type 1 count; type 2 lower; types 3-5 k_D
  array[A] int<lower=0> yD_hi;   // type 2 upper

  // index sets
  int<lower=0> R1;  array[R1] int<lower=1, upper=R> rows1;   // type-1 rows
  int<lower=0> A1;  array[A1] int<lower=1, upper=A> arms1;   // drug arms of type-1 rows
  int<lower=0> R5;  array[R5] int<lower=1, upper=R> rows5;   // type-5 rows

  // type-5 support: log choose vectors and thresholds, stored flat
  int<lower=0> LP; vector[LP] lchP; array[R] int<lower=0> P_start;   // placebo 0..NP
  int<lower=0> LD; vector[LD] lchD; array[A] int<lower=0> D_start;   // drug arm 0..ND
  array[R] int<lower=0> n_arm5;                                       // arms per type-5 row
  array[R, 2] int<lower=0> arm5;                                      // up to 2 arms per row
  int<lower=0> T5; array[T5] int<lower=0> thr; array[R] int<lower=0> t5_start;

  // types 3-4: every hidden arm is a "unit" u with log CDF = log sum_{y=0}^{k_u} Bin(y; N_u, p_u).
  // All terms are stored flat; lin_u indexes p_u in [mu; eta_D_row] (1..R placebo, R+1..2R drug).
  int<lower=0> U; array[U] int<lower=1> u_start; array[U] int<lower=1> u_len;
  int<lower=0> TC; vector[TC] c_lch; vector[TC] c_y; vector[TC] c_N; array[TC] int<lower=1> c_lin;

  // label selection: kmin per row (0 = none); denominator units with terms y = 0..kmin-1 on pooled N
  array[R] int<lower=0> kmin;
  int<lower=0> SA; array[SA] int<lower=1, upper=R> a_row; vector[SA] a_N;     // kmin <= 4
  matrix[SA, 3] a_lc; matrix[SA, 3] a_use;              // log choose(N, j) and 1 if j < kmin, j = 1..3
  int<lower=1> Q; vector[Q] gh_x; vector[Q] gh_lw;      // Gauss-Hermite nodes and log(weight / sqrt(pi))
  int<lower=0> SB; array[SB] int<lower=1, upper=R> b_row;                       // kmin >= 5
  matrix[SB, 3] b_la; matrix[SB, 3] b_c; matrix[SB, 3] b_t2;                   // log weight, centre, scale^2

  int<lower=0, upper=1> prior_only;
}

parameters {
  real alpha;
  real beta;
  matrix[2, J] z_prod;  vector<lower=0>[2] sd_prod; cholesky_factor_corr[2] L_prod;
  matrix[2, K] z_pt;    vector<lower=0>[2] sd_pt;   cholesky_factor_corr[2] L_pt;
  vector[R] mu;         vector[R] delta;           // row-level, centred parameterisation
  vector<lower=0>[2] sd_row;  cholesky_factor_corr[2] L_row;
  vector[S] z_c; vector[S] z_s; real<lower=0> sd_c; real<lower=0> sd_s;
}

transformed parameters {
  vector[R] m_mu;      // expected placebo log-odds without the row effect e_r
  vector[R] m_delta;   // expected log-odds ratio without the row effect w_r
  {
    matrix[2, J] re_prod = diag_pre_multiply(sd_prod, L_prod) * z_prod;
    matrix[2, K] re_pt = diag_pre_multiply(sd_pt, L_pt) * z_pt;
    m_mu = alpha + re_prod[1, prod]' + re_pt[1, pt]' + sd_c * z_c[soc];
    m_delta = beta + re_prod[2, prod]' + re_pt[2, pt]' + sd_s * z_s[soc];
  }
}

model {
  alpha ~ normal(-3, 1.5);
  beta ~ normal(0, 1);
  to_vector(z_prod) ~ std_normal();
  to_vector(z_pt) ~ std_normal();
  {
    // (e_r, w_r) = (mu_r - m_mu_r, delta_r - m_delta_r) ~ bivariate normal(0, Sigma_row)
    matrix[2, 2] L_S = diag_pre_multiply(sd_row, L_row);
    array[R] vector[2] res;
    for (r in 1:R) res[r] = [mu[r] - m_mu[r], delta[r] - m_delta[r]]';
    res ~ multi_normal_cholesky(rep_vector(0, 2), L_S);
  }
  z_c ~ std_normal(); z_s ~ std_normal();
  sd_prod ~ normal(0, 1); sd_pt ~ normal(0, 1); sd_row ~ normal(0, 1);
  sd_c ~ normal(0, 1); sd_s ~ normal(0, 1);
  L_prod ~ lkj_corr_cholesky(2); L_pt ~ lkj_corr_cholesky(2); L_row ~ lkj_corr_cholesky(2);

  if (!prior_only) {
    vector[R] eta_D = mu + delta;
    // type 1, vectorised
    if (R1 > 0) target += binomial_logit_lpmf(yP[rows1] | NP[rows1], mu[rows1]);
    if (A1 > 0) target += binomial_logit_lpmf(yD[arms1] | ND[arms1], eta_D[arm_row[arms1]]);

    for (r in 1:R) {
      if (type[r] == 2) {
        vector[yP_hi[r] - yP[r] + 1] lp;
        for (y in yP[r]:yP_hi[r]) lp[y - yP[r] + 1] = binomial_logit_lpmf(y | NP[r], mu[r]);
        target += log_sum_exp(lp);
      }
    }
    for (a in 1:A) {
      int r = arm_row[a];
      if (type[r] == 2) {
        vector[yD_hi[a] - yD[a] + 1] lp;
        for (y in yD[a]:yD_hi[a]) lp[y - yD[a] + 1] = binomial_logit_lpmf(y | ND[a], eta_D[r]);
        target += log_sum_exp(lp);
      }
    }
    if (U > 0) {
      vector[2 * R] lin = append_row(mu, eta_D);
      vector[2 * R] lp = log_inv_logit(lin);
      vector[2 * R] l1mp = log1m_inv_logit(lin);
      vector[TC] terms = c_lch + c_y .* lp[c_lin] + (c_N - c_y) .* l1mp[c_lin];
      for (u in 1:U) target += log_sum_exp(segment(terms, u_start[u], u_len[u]));
    }
    for (q in 1:R5) {
      int r = rows5[q];
      int nP = NP[r];
      vector[nP + 1] yseq = linspaced_vector(nP + 1, 0, nP);
      vector[nP + 1] terms = segment(lchP, P_start[r], nP + 1)
                             + yseq * log_inv_logit(mu[r])
                             + (nP - yseq) * log1m_inv_logit(mu[r]);
      for (i in 1:n_arm5[r]) {
        int a = arm5[r, i];
        int nD = ND[a];
        vector[nD + 1] ys = linspaced_vector(nD + 1, 0, nD);
        vector[nD + 1] lpD = segment(lchD, D_start[a], nD + 1)
                             + ys * log_inv_logit(eta_D[r])
                             + (nD - ys) * log1m_inv_logit(eta_D[r]);
        int k0 = kmin[r];                       // 0 unless the selection adjustment applies (single arm)
        vector[nD + 1] lF;                      // lF[t] = log P(k0 <= Y <= t - 1)
        lF[k0 + 1] = lpD[k0 + 1];
        for (t in (k0 + 2):(nD + 1)) lF[t] = log_sum_exp(lF[t - 1], lpD[t]);
        for (y in 0:nP) {
          int th = thr[t5_start[r] + y * n_arm5[r] + i - 1];
          terms[y + 1] += (th >= k0) ? lF[th + 1] : negative_infinity();
        }
      }
      target += log_sum_exp(terms);
    }
    if (SA + SB > 0) {
      real s2 = square(sd_row[1]) + square(sd_row[2]) + 2 * sd_row[1] * sd_row[2] * L_row[2, 1];
      vector[R] m_eta = m_mu + m_delta;
      if (SA > 0) {
        vector[SA] me = m_eta[a_row];
        vector[SA] Ps = rep_vector(0, SA);
        real se = sqrt(s2);
        for (q in 1:Q) {
          vector[SA] eq = me + sqrt2() * se * gh_x[q];
          vector[SA] lp = log_inv_logit(eq);
          vector[SA] l1p = log1m_inv_logit(eq);
          vector[SA] low = exp(a_N .* l1p);
          for (j in 1:3) low += a_use[:, j] .* exp(a_lc[:, j] + j * lp + (a_N - j) .* l1p);
          Ps += exp(gh_lw[q]) * (1 - low);
        }
        target += -sum(log(Ps));
      }
      for (u in 1:SB) {
        real m = m_eta[b_row[u]];
        vector[3] lt;
        for (h in 1:3) lt[h] = b_la[u, h] + std_normal_lcdf((m - b_c[u, h]) / sqrt(b_t2[u, h] + s2));
        target += -log_sum_exp(lt);
      }
    }
  }
}

generated quantities {
  vector[R] pD_lat = inv_logit(mu + delta);
  vector[R] E = pD_lat - inv_logit(mu);
  real theta = 0;
  real theta_mis = 0;
  real theta_conv;
  for (r in 1:R) {
    theta += (E[r] < 0.01);
    theta_mis += (pD_lat[r] < 0.01);
  }
  theta /= R;
  theta_mis /= R;
  theta_conv = theta - theta_mis;
}
