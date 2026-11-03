# Starting values for the sampler: each side effect starts at its observed rates
# (with a small continuity adjustment). Hidden rows start just below their table cutoff.
make_init <- function(sd) {
  R <- sd$R
  pP <- (sd$yP + 0.5) / (sd$NP + 1)
  yD_row <- tapply(sd$yD, sd$arm_row, sum)[as.character(1:R)]
  ND_row <- tapply(sd$ND, sd$arm_row, sum)[as.character(1:R)]
  pD <- (as.numeric(yD_row) + 0.5) / (as.numeric(ND_row) + 1)
  hid <- sd$type %in% 3:5
  pP[hid] <- pmax(0.002, pmin(pP[hid], 0.5) * 0.6)
  pD[hid] <- pmax(0.002, pmin(pD[hid], 0.5) * 0.6)
  mu0 <- qlogis(pP); de0 <- qlogis(pD) - mu0
  function() list(
    alpha = median(mu0), beta = median(de0),
    z_prod = matrix(0, 2, sd$J), z_pt = matrix(0, 2, sd$K),
    sd_prod = c(0.3, 0.3), sd_pt = c(0.3, 0.3), sd_row = c(0.5, 0.5),
    L_prod = diag(2), L_pt = diag(2), L_row = diag(2),
    z_c = rep(0, sd$S), z_s = rep(0, sd$S), sd_c = 0.3, sd_s = 0.3,
    mu = as.numeric(mu0), delta = as.numeric(de0))
}
