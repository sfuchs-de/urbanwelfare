#' Simulate a Synthetic Multi-Neighborhood City with Known Ground-Truth Welfare
#'
#' Generates a complete synthetic open-city economy following the Section 2.8
#' Monte Carlo data generating process of Allen et al. (2026), supporting CES,
#' non-CES Homothetic Single Aggregator (HSA + capacity technology), and
#' heterogeneous-elasticity stress-test worlds, and solves for exact nonlinear
#' (or first-order) equilibrium bilateral spending, income, prices, and
#' money-metric resident welfare.
#'
#' @param N Number of residential/commercial neighborhoods (default \code{80}).
#' @param T Number of time periods / months (default \code{12}).
#' @param K Number of commercial activity sector groups (default \code{4}).
#' @param n_visitor_groups Number of external visitor origin groups for the
#'   shift-share shock structure (default \code{6}).
#' @param model Character: \code{"ces"} (default), \code{"hsa"} (non-CES HSA
#'   demand + capacity production), \code{"iia_het"} (origin-varying \eqn{\sigma_n}),
#'   or \code{"kappa_het"} (within-group \eqn{\kappa_i^S} dispersion).
#' @param sigma Mean resident shopping substitution elasticity (default \code{3.4555}).
#' @param nu Non-housing share of resident expenditure (default \code{0.68}).
#' @param tourism_share Ratio of aggregate base-month visitor spending to
#'   aggregate local resident spending (default \code{0.10}).
#' @param shock_scale Multiplier on seasonal + group shock amplitude (default \code{1.0}).
#' @param noise_sd Log-normal measurement noise SD added to observed bilateral
#'   spending (default \code{0.0}).
#' @param seed Integer random seed (default \code{42}).
#'
#' @return A named list with the exact schema of \code{\link{sim_city_160}}.
#' @export
#' @examples
#' sim <- sim_urban_city(N = 30, T = 6, seed = 123)
#' acct <- urban_account(sim$E_R0, sim$H_0, sim$origins, sim$destinations)
#' exps <- build_exposures(acct, sim$shocks)
#' fit  <- estimate_bilateral(sim$bilateral_panel, exps)
#' welf <- recover_welfare(fit, acct, se_method = "none")
#' cor(welf$tract_summary$welfare_lp, sim$truth$welfare_peak_lp)
sim_urban_city <- function(N = 80L,
                           T = 12L,
                           K = 4L,
                           n_visitor_groups = 6L,
                           model = c("ces", "hsa", "iia_het", "kappa_het"),
                           sigma = 3.4555,
                           nu = 0.68,
                           tourism_share = 0.10,
                           shock_scale = 1.0,
                           noise_sd = 0.0,
                           seed = 42L) {
  model <- match.arg(model)
  set.seed(as.integer(seed))
  N <- as.integer(N)
  T <- as.integer(T)
  K <- as.integer(K)
  G_vis <- as.integer(n_visitor_groups)

  # 1. Spatial coordinates & monocentric commercial attractiveness
  xy <- matrix(stats::runif(2L * N), nrow = N, ncol = 2L)
  dist_mat <- as.matrix(stats::dist(xy))
  centre_d <- sqrt((xy[, 1] - 0.5)^2 + (xy[, 2] - 0.5)^2)

  pop <- exp(stats::rnorm(N, 0, 0.5))
  pop <- pop / mean(pop)
  inc_pc <- exp(stats::rnorm(N, 0, 0.35))
  v_guess <- pop * inc_pc

  attract <- exp(stats::rnorm(N, 0, 0.5)) * (1.0 + 2.0 * exp(-centre_d^2 / 0.05))
  local_w <- sweep(exp(-3.0 * dist_mat), 2, attract, "*")
  s_out_NH <- stats::runif(N, 0.08, 0.20)
  S_NH0 <- sweep(local_w, 1, rowSums(local_w), "/")
  S_NH0 <- sweep(S_NH0, 1, 1.0 - s_out_NH, "*")
  res_guess <- as.numeric(crossprod(S_NH0, nu * v_guess))

  # 2. Visitor hotspots & seasonal shift-share demand
  vis_weight <- exp(-centre_d / 0.15) + 0.05
  hotspots <- matrix(stats::runif(2L * G_vis, 0.2, 0.8), nrow = G_vis, ncol = 2L)
  pi_vis <- matrix(0.0, nrow = N, ncol = G_vis)
  for (g in seq_len(G_vis)) {
    d_h <- sqrt((xy[, 1] - hotspots[g, 1])^2 + (xy[, 2] - hotspots[g, 2])^2)
    w_g <- vis_weight * exp(-d_h / 0.20)
    pi_vis[, g] <- w_g / sum(w_g)
  }
  grp_raw <- stats::rgamma(G_vis, shape = 3.0, rate = 1.0)
  grp_size <- grp_raw / sum(grp_raw)
  base_tourism <- tourism_share * sum(res_guess)

  months <- seq(0L, T - 1L)
  season <- 0.5 * exp(-((months %% 12) - 6.0)^2 / 4.0) - 0.5 * exp(-6.0^2 / 4.0)
  shocks_g <- matrix(stats::rnorm(G_vis * T, 0, 0.15), nrow = G_vis, ncol = T)
  shocks_g[, 1L] <- 0.0
  log_shift_gt <- shock_scale * sweep(shocks_g, 2, season, "+")
  spend_gt <- sweep(exp(log_shift_gt), 1, base_tourism * grp_size, "*")
  ET_mat <- pi_vis %*% spend_gt # N x T

  # Exposure shares pi_{ig} = initial share of destination i's visitor spending from group g
  ET_i0 <- ET_mat[, 1L]
  shares_ig <- sweep(sweep(pi_vis, 2, base_tourism * grp_size, "*"), 1, ET_i0, "/")
  expected_gt <- matrix(rep(shock_scale * season, each = G_vis), nrow = G_vis, ncol = T)

  # 3. Technology, commuting & balanced initial equilibrium at unit prices
  A_fixed0 <- stats::runif(N, 0.10, 0.40) * res_guess
  sector_group <- sample.int(K, N, replace = TRUE)
  kappa_group <- stats::runif(K, 0.32, 0.58)
  kappa <- kappa_group[sector_group]
  if (model == "kappa_het") {
    kappa <- pmin(0.80, pmax(0.20, kappa + stats::runif(N, -0.12, 0.12)))
  }

  mu <- if (model == "hsa") {
    eps_s <- 1.0 / kappa - 1.0
    eps_s * log(1.0 + 1.0 / eps_s)
  } else {
    stats::runif(N, 0.35, 0.60)
  }
  theta_G <- 1.0 - mu

  commute <- sweep(exp(-4.0 * dist_mat), 1, pop, "*")
  commute <- sweep(commute, 2, colSums(commute), "/")
  res_f_share <- stats::runif(N, 0.65, 0.92)
  H_0 <- sweep(commute, 2, res_f_share, "*")
  v_out <- 0.45 * v_guess

  # Solve initial linear account closure: v0 = H_0 diag(1 - mu) y0 + v_out, y0 = nu S_NH0' v0 + ET_0 + A_fixed0
  B_mat <- sweep(H_0, 2, 1.0 - mu, "*")
  v0 <- as.numeric(solve(diag(N) - B_mat %*% (nu * t(S_NH0)), B_mat %*% (ET_i0 + A_fixed0) + v_out))
  y0 <- as.numeric(crossprod(S_NH0, nu * v0)) + ET_i0 + A_fixed0
  G0 <- (1.0 - mu) * y0
  E_R0 <- sweep(S_NH0, 1, nu * v0, "*")
  W_0 <- sweep(H_0, 2, G0, "*")

  sigma_n <- if (model == "iia_het") {
    sigma * exp(stats::rnorm(N, -0.5 * 0.25^2, 0.25))
  } else {
    rep(sigma, N)
  }

  # 4. Solve nonlinear monthly equilibria t = 1..T (month 1 is base equilibrium)
  logv_true <- matrix(0.0, nrow = N, ncol = T)
  logp_true <- matrix(0.0, nrow = N, ncol = T)
  logy_true <- matrix(0.0, nrow = N, ncol = T)
  welf_true <- matrix(0.0, nrow = N, ncol = T)
  logE_arr <- array(0.0, dim = c(N, N, T))
  for (t_idx in seq_len(T)) {
    if (t_idx == 1L) {
      logE_arr[, , 1L] <- log(pmax(1e-18, E_R0))
      next
    }
    if (model == "ces") {
      sol <- cpp_solve_qsm(
        y0 = y0, mu = mu, kappa = kappa, A_fixed = A_fixed0,
        ET_shock = ET_mat[, t_idx], v_out = v_out, nu = rep(nu, N),
        s_out0 = s_out_NH, S_NH0 = S_NH0, H_alloc = H_0,
        pop = pop, is_res = rep(1.0, N), w_welf = pop / sum(pop),
        theta = rep(0.0, N), use_tax_demand = FALSE,
        sigma_R = sigma, sigma_T = 7.0, eta_out = 2.5,
        logp_init = logp_true[, t_idx - 1L], compute_jacobian = FALSE
      )
      lp <- as.numeric(sol$logp)
      vt <- as.numeric(sol$v)
      yt <- as.numeric(sol$y)
      logP_n <- as.numeric(sol$logP_NH)
      exp_p <- exp((1.0 - sigma) * lp)
      denom <- s_out_NH + as.numeric(S_NH0 %*% exp_p)
      S_t <- sweep(sweep(S_NH0, 2, exp_p, "*"), 1, denom, "/")
    } else {
      sol_r <- .solve_general_dgp(model, N, y0, mu, kappa, A_fixed0, ET_mat[, t_idx],
                                  v_out, v0, nu, s_out_NH, S_NH0, H_0, sigma_n,
                                  logp_true[, t_idx - 1L])
      lp <- sol_r$logp
      vt <- sol_r$v
      yt <- sol_r$y
      logP_n <- sol_r$logP_NH
      S_t <- sol_r$S_NH
    }
    logp_true[, t_idx] <- lp
    logv_true[, t_idx] <- log(vt / v0)
    logy_true[, t_idx] <- log(yt / y0)
    welf_true[, t_idx] <- log(vt / v0) - nu * logP_n
    E_Rt <- sweep(S_t, 1, nu * vt, "*")
    noise <- if (noise_sd > 0) matrix(stats::rnorm(N * N, 0, noise_sd), N, N) else 0.0
    logE_arr[, , t_idx] <- log(pmax(1e-18, E_Rt)) + noise
  }

  ids <- sprintf("T%03d", seq_len(N))
  t_labels <- sprintf("M%02d", seq_len(T))
  dimnames(E_R0) <- dimnames(W_0) <- dimnames(H_0) <- list(ids, ids)
  rownames(shares_ig) <- ids
  colnames(shares_ig) <- paste0("G", seq_len(G_vis))
  dimnames(log_shift_gt) <- dimnames(expected_gt) <- list(colnames(shares_ig), t_labels)

  origins_df <- data.frame(
    id = ids, x = xy[, 1], y = xy[, 2], pop = pop,
    v0 = v0, v_out = v_out, s_H = rep(1.0 - nu, N),
    s_out_NH = s_out_NH, is_city_resident = rep(TRUE, N),
    stringsAsFactors = FALSE
  )
  destinations_df <- data.frame(
    id = ids, x = xy[, 1], y = xy[, 2], y0 = y0, G0 = G0,
    E_ext0 = ET_i0, A_fixed0 = A_fixed0, sector_group = sector_group,
    theta_G = theta_G, theta_L = rep(0.65, N), kappa_true = kappa,
    stringsAsFactors = FALSE
  )

  h_mat <- log(sweep(ET_mat, 1, ET_i0, "/"))
  d_mat <- sweep(h_mat, 1, ET_i0 / y0, "*")
  shocks_df <- data.frame(
    dest = rep(ids, times = T),
    time = rep(t_labels, each = N),
    h = as.vector(h_mat),
    d = as.vector(d_mat),
    E_ext = as.vector(ET_mat),
    stringsAsFactors = FALSE
  )

  grid_od <- expand.grid(origin = ids, dest = ids, time = t_labels,
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid_od$log_spending <- as.vector(logE_arr)
  grid_od$spending <- exp(grid_od$log_spending)

  peak_idx <- which.max(colMeans(d_mat))
  dimnames(logv_true) <- dimnames(logp_true) <- dimnames(logy_true) <- dimnames(welf_true) <- list(ids, t_labels)

  list(
    E_R0 = E_R0,
    W_0 = W_0,
    H_0 = H_0,
    origins = origins_df,
    destinations = destinations_df,
    shocks = shocks_df,
    shift_share = list(
      shares = shares_ig,
      shifters = log_shift_gt,
      shifter_expectations = expected_gt
    ),
    bilateral_panel = grid_od,
    truth = list(
      model = model,
      peak_time = t_labels[peak_idx],
      peak_index = peak_idx,
      logv = logv_true,
      logp = logp_true,
      logy = logy_true,
      welfare = welf_true,
      welfare_peak_lp = 100 * welf_true[, peak_idx]
    )
  )
}

.solve_general_dgp <- function(model, N, y0, mu, kappa, A_fixed0, ET_t,
                               v_out, v0, nu, s_out_NH, S_NH0, H_0, sigma_n, lp_init) {
  hsa_sig <- 1.8
  a_all <- cbind(s_out_NH, S_NH0)

  eval_shares <- function(lp) {
    lp_all <- c(0.0, lp)
    if (model == "hsa") {
      beta_n <- sigma_n - hsa_sig
      logA <- rep(0.0, N)
      for (it in 1:50) {
        logz <- outer(rep(1.0, N), lp_all) - logA
        z <- exp(logz)
        S <- a_all * exp((1.0 - hsa_sig) * logz - sweep(z - 1.0, 1, beta_n, "*"))
        f <- rowSums(S) - 1.0
        df <- rowSums(S * (-(1.0 - hsa_sig) + sweep(z, 1, beta_n, "*")))
        step <- f / df
        logA <- logA - step
        if (max(abs(step)) < 1e-13) break
      }
      logz <- outer(rep(1.0, N), lp_all) - logA
      S <- a_all * exp((1.0 - hsa_sig) * logz - sweep(exp(logz) - 1.0, 1, beta_n, "*"))
      return(S)
    }
    w <- a_all * exp(outer(1.0 - sigma_n, lp_all))
    sweep(w, 1, rowSums(w), "/")
  }

  eval_supply <- function(lp) {
    p <- exp(lp)
    if (model == "hsa") {
      m_cap <- (1.0 / kappa - 1.0) * y0
      Q_cap <- y0 + m_cap
      M_in <- m_cap * log(pmax(1e-8, p * Q_cap / m_cap))
      y <- p * Q_cap - m_cap
      return(list(y = y, G = y - M_in))
    }
    s_prod <- (1.0 / kappa - 1.0) * (1.0 - mu) / mu
    sh_m <- pmin(0.99999, mu * exp((s_prod - 1.0) * lp))
    q_ratio <- ((1.0 - mu) / (1.0 - sh_m))^(s_prod / (s_prod - 1.0))
    y <- y0 * p * q_ratio
    list(y = y, G = (1.0 - sh_m) * y)
  }

  excess_fn <- function(lp) {
    sup <- eval_supply(lp)
    v <- as.numeric(H_0 %*% sup$G) + v_out
    S_all <- eval_shares(lp)
    S_loc <- S_all[, -1L, drop = FALSE]
    (sup$y - as.numeric(crossprod(S_loc, nu * v)) - ET_t - A_fixed0) / y0
  }

  # Damped Newton in R for general non-CES stress models
  lp <- lp_init
  f <- excess_fn(lp)
  eps_fd <- 1e-6
  for (it in 1:30) {
    if (max(abs(f)) < 1e-10) break
    J <- matrix(0.0, N, N)
    for (j in seq_len(N)) {
      lp_u <- lp; lp_u[j] <- lp_u[j] + eps_fd
      lp_d <- lp; lp_d[j] <- lp_d[j] - eps_fd
      J[, j] <- (excess_fn(lp_u) - excess_fn(lp_d)) / (2.0 * eps_fd)
    }
    step <- -solve(J, f)
    lp <- lp + step
    f <- excess_fn(lp)
  }

  sup <- eval_supply(lp)
  v <- as.numeric(H_0 %*% sup$G) + v_out
  S_all <- eval_shares(lp)
  if (model == "hsa") {
    # 8-point Gauss-Legendre path integral of homothetic price index
    nodes <- c(-0.96028986, -0.79666648, -0.52553241, -0.18343464,
               0.18343464, 0.52553241, 0.79666648, 0.96028986)
    wts <- c(0.10122854, 0.22238103, 0.31370665, 0.36268378,
             0.36268378, 0.31370665, 0.22238103, 0.10122854)
    logP_NH <- rep(0.0, N)
    for (q in seq_along(nodes)) {
      Sq <- eval_shares(0.5 * (nodes[q] + 1.0) * lp)[, -1L, drop = FALSE]
      logP_NH <- logP_NH + 0.5 * wts[q] * as.numeric(Sq %*% lp)
    }
  } else {
    lp_all <- c(0.0, lp)
    logP_NH <- log(rowSums(a_all * exp(outer(1.0 - sigma_n, lp_all)))) / (1.0 - sigma_n)
  }
  list(logp = lp, y = sup$y, v = v, logP_NH = logP_NH, S_NH = S_all[, -1L, drop = FALSE])
}
