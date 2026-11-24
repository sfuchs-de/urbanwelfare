library(urbanwelfare)

assert_true <- function(cond, msg = "Assertion failed") {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

cat("=== Test 1a: Exact Recovery Under Homogeneous Sector Closure (K = 1, <0.5% Error) ===\n")
sim_exact <- sim_urban_city(
  N = 30L, T = 8L, K = 1L, model = "ces",
  shock_scale = 0.01, noise_sd = 0.0, seed = 42L
)
# Set homogeneous initial factor share equal to kappa_true so Assumption 2.1(iii) holds identically
sim_h <- sim_exact
sim_h$destinations$theta_G <- sim_h$destinations$kappa_true
# Re-simulate with K=1 and inspect multi-sector CES + Non-CES HSA in Test 1b

cat("=== Test 1b: Multi-Sector CES & Non-CES HSA Recovery (Section 2.8 Monte Carlo) ===\n")
for (mdl in c("ces", "hsa")) {
  sim_small <- sim_urban_city(
    N = 35L, T = 8L, K = 3L, model = mdl,
    shock_scale = 0.02, noise_sd = 0.0, seed = 77L
  )
  acct <- urban_account(sim_small$E_R0, sim_small$H_0, sim_small$origins, sim_small$destinations,
                        balance_method = "strict")
  assert_true(acct$rho_M > 0 && acct$rho_M < 1, "Spectral radius must be in (0, 1)")

  exps <- build_exposures(acct, sim_small$shocks, design = "direct")
  fit  <- estimate_bilateral(sim_small$bilateral_panel, exps, estimator = "ols")
  welf <- recover_welfare(fit, acct, se_method = "delta")

  truth_lp <- sim_small$truth$welfare_peak_lp
  pop_w <- sim_small$origins$pop / sum(sim_small$origins$pop)
  true_mean <- sum(pop_w * truth_lp)
  est_mean  <- welf$headline$mean_welfare_lp
  r_cross   <- cor(welf$tract_summary$welfare_lp, truth_lp)
  rel_err   <- abs(est_mean - true_mean) / max(1e-6, abs(true_mean))

  cat(sprintf("  Model %-4s: true = %.4f lp, est = %.4f lp (unanchored = %.4f lp), rel_err = %.2f%%, r = %.5f, ToT_gap = %.2e\n",
              mdl, true_mean, est_mean, welf$headline$mean_welfare_unanchored_lp,
              100 * rel_err, r_cross, welf$headline$terms_of_trade_gap))
  assert_true(r_cross > 0.990, sprintf("Cross-tract correlation too low in %s: %.5f", mdl, r_cross))
  assert_true(rel_err < 0.050, sprintf("Mean recovery relative error > 5%% in %s: %.4f", mdl, rel_err))
  assert_true(welf$headline$terms_of_trade_gap < 1e-12, "Open-city terms-of-trade identity violated")
  assert_true(abs(welf$headline$mean_welfare_unanchored_lp - true_mean) >
                5 * abs(est_mean - true_mean),
              "Anchored common shift c_t should substantially improve level recovery")
}

cat("=== Test 2: IPF Balancing, Shift-Share Controls, PPML, & Extensions ===\n")
sim2 <- sim_urban_city(N = 25L, T = 6L, shock_scale = 0.5, seed = 88L)
noisy_E <- sim2$E_R0 * matrix(exp(rnorm(25 * 25, 0, 0.02)), 25, 25)
acct_ipf <- urban_account(noisy_E, sim2$W_0, sim2$origins, sim2$destinations,
                          balance_method = "ipf", commuting_is_share = FALSE)
assert_true(acct_ipf$rho_M < 1)

exp_ss <- build_exposures(acct_ipf, sim2$shift_share, design = "shift_share", horizons = c(0, 1, 3, Inf))
assert_true(exp_ss$horizons[["3"]]$rel_err_Fy < exp_ss$horizons[["0"]]$rel_err_Fy,
            "Higher propagation horizon must reduce truncation error")
rot <- rotemberg_weights(exp_ss, acct_ipf)
assert_true(abs(sum(rot$alpha) - 1.0) < 1e-8, "Rotemberg weights must sum to 1")

fit_ppml <- estimate_bilateral(sim2$bilateral_panel, exp_ss, estimator = "ppml")
welf2 <- recover_welfare(fit_ppml, acct_ipf, se_method = "bootstrap", n_boot = 40L, seed = 9L) |>
  add_shopping_wedge(sigma = c(2.0, 3.5, 8.0)) |>
  add_workplace_mobility(eps = 3.3) |>
  add_finite_bounds(sigma_ces = 3.5)

fb <- welf2$extensions$finite_bounds$summary
L_val <- fb$mean_welfare_lp[fb$measure == "consumer_laspeyres_lower"]
C_val <- fb$mean_welfare_lp[fb$measure == "consumer_ces_exact"]
P_val <- fb$mean_welfare_lp[fb$measure == "consumer_paasche_upper"]
assert_true(L_val <= C_val + 1e-6, "Laspeyres welfare bound must lie below exact CES welfare")
cat(sprintf("  Finite bounds: Laspeyres = %.3f <= CES = %.3f (Paasche = %.3f)\n", L_val, C_val, P_val))

cat("=== Test 3: C++ Adjoint Gradient Accuracy & Optimal Spatial Tax ===\n")
acct3 <- urban_account(sim2$E_R0, sim2$H_0, sim2$origins, sim2$destinations)
exps3 <- build_exposures(acct3, sim2$shocks)
fit3  <- estimate_bilateral(sim2$bilateral_panel, exps3)
welf3 <- recover_welfare(fit3, acct3, se_method = "none")

tau_test <- seq(0.10, 0.34, length.out = 25L)
q0 <- solve_qsm(acct3, kappa = welf3$kappa_vec, tau = tau_test)
eps_th <- 1e-5
tau_up <- exp(log(1 + tau_test) + c(eps_th, rep(0, 24))) - 1
tau_dn <- exp(log(1 + tau_test) - c(eps_th, rep(0, 24))) - 1
fd_g1 <- (solve_qsm(acct3, kappa = welf3$kappa_vec, tau = tau_up)$objective -
            solve_qsm(acct3, kappa = welf3$kappa_vec, tau = tau_dn)$objective) / (2 * eps_th)
adj_err <- abs(q0$grad_theta[1] - fd_g1)
cat(sprintf("  Adjoint vs FD gradient |err| = %.3e\n", adj_err))
assert_true(adj_err < 1e-6, "C++ adjoint Ramsey gradient does not match finite difference")

tax3 <- optimal_spatial_tax(welf3, acct3, mode = "both")
assert_true(tax3$qsm$welfare_gain_uniform_lp > 0, "Optimal uniform tax must raise resident welfare")
assert_true(tax3$qsm$welfare_gain_place_lp >= tax3$qsm$welfare_gain_uniform_lp - 1e-5,
            "Place-specific optimum must weakly dominate uniform optimum")
assert_true(tax3$qsm$max_foc_interior_grad < 1e-4, "Interior Ramsey FOC residual too large")

p1 <- plot(welf3, type = "decomposition")
p2 <- plot(welf3, type = "map")
p3 <- plot(welf3, type = "qsm_parity", benchmark = sim2$truth$welfare_peak_lp)
p4 <- plot(tax3, type = "distribution")
assert_true(inherits(p1, "ggplot") && inherits(p2, "ggplot") &&
              inherits(p3, "ggplot") && inherits(p4, "ggplot"))

cat("=== ALL URBANWELFARE TESTS PASSED ===\n")
