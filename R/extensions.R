#' Add CES Shopping Amenity/Congestion Residual Valuation (Section 5.2.1)
#'
#' Evaluates Redding-Weinstein (2020) CES perceived destination shopping cost
#' changes \eqn{\widehat{c}^{CES}_{it} = \widehat{\Delta}^D_{it} / (1 - \sigma)},
#' isolates the shopping residual \eqn{\widehat{\zeta}_{it} = \widehat{c}^{CES}_{it} - \widehat{p}^{prod}_{it}},
#' and computes augmented resident welfare across one or more candidate
#' substitution elasticities \code{sigma}.
#'
#' @param welfare An \code{\link{urban_welfare}} object.
#' @param sigma Numeric vector of substitution elasticities \eqn{\sigma \ne 1}
#'   (default \code{c(2.0, 3.456, 8.235)}).
#' @return The input \code{urban_welfare} object with \code{extensions$shopping_wedge}
#'   populated.
#' @export
add_shopping_wedge <- function(welfare, sigma = c(2.0, 3.456, 8.235)) {
  if (!inherits(welfare, "urban_welfare")) stop("welfare must be an 'urban_welfare' object.")
  if (any(abs(sigma - 1.0) < 1e-6)) stop("sigma cannot equal 1 in CES share inversion.")

  acct <- welfare$account
  p_idx <- welfare$headline$period_index
  is_res <- acct$origins$is_city_resident
  pop_w <- acct$origins$pop[is_res] / sum(acct$origins$pop[is_res])

  Delta_D <- welfare$responses$Delta_D_hat[, p_idx]
  p_prod <- welfare$anchored$dlog_p[, p_idx]
  u_base <- welfare$anchored$welfare[is_res, p_idx]
  S_res <- acct$S_tot[is_res, , drop = FALSE]

  rows <- lapply(sigma, function(sig) {
    c_ces <- Delta_D / (1.0 - sig)
    zeta <- c_ces - p_prod
    u_aug <- u_base - as.numeric(S_res %*% zeta)
    data.frame(
      sigma = sig,
      mean_producer_price_lp = 100 * mean(p_prod),
      mean_ces_cost_lp = 100 * mean(c_ces),
      mean_residual_zeta_lp = 100 * mean(zeta),
      mean_welfare_aug_lp = 100 * sum(pop_w * u_aug),
      p10_welfare_aug_lp = 100 * .weighted_quantile(u_aug, pop_w, 0.10),
      p90_welfare_aug_lp = 100 * .weighted_quantile(u_aug, pop_w, 0.90)
    )
  })
  summary_df <- do.call(rbind, rows)
  welfare$extensions$shopping_wedge <- list(summary = summary_df)
  welfare
}

#' Add Workplace Mobility Le Chatelier Price Rules (Section 5.2.2)
#'
#' Compares fixed-workplace factor-share living costs and resident welfare
#' against finite Frechet workplace reallocation (\code{eps}) and perfect
#' city-wide labor mobility, holding estimated income and sales responses fixed.
#'
#' @param welfare An \code{\link{urban_welfare}} object.
#' @param eps Positive Frechet workplace sorting elasticity \eqn{\varepsilon}
#'   (default \code{3.3}).
#' @param include_perfect Logical (default \code{TRUE}): also compute the
#'   \eqn{\varepsilon \to \infty} city-wide wage equalization benchmark.
#' @return The input \code{urban_welfare} object with
#'   \code{extensions$workplace_mobility} populated.
#' @export
add_workplace_mobility <- function(welfare, eps = 3.3, include_perfect = TRUE) {
  if (!inherits(welfare, "urban_welfare")) stop("welfare must be an 'urban_welfare' object.")
  if (eps <= 0) stop("eps must be strictly positive.")

  acct <- welfare$account
  p_idx <- welfare$headline$period_index
  is_res <- acct$origins$is_city_resident
  pop_w <- acct$origins$pop[is_res] / sum(acct$origins$pop[is_res])
  Nd <- acct$N_dest

  y0 <- acct$destinations$y0
  theta_G <- acct$destinations$theta_G
  theta_L <- acct$destinations$theta_L
  a_i <- theta_L * theta_G
  b_i <- (1.0 - theta_L) * theta_G

  dlog_y <- welfare$anchored$dlog_y[, p_idx]
  dlog_v <- welfare$anchored$dlog_v[is_res, p_idx]
  S_res <- acct$S_tot[is_res, , drop = FALSE]

  # 1. Fixed-workplace factor-share benchmark: dlog p = (a_i + b_i) dlog y
  p_fixed <- (a_i + b_i) * dlog_y
  cost_fixed <- as.numeric(S_res %*% p_fixed)
  welf_fixed <- dlog_v - cost_fixed

  # 2. Finite Frechet workplace sorting (eps): solve linearized workplace labor market clearing
  payroll0 <- sweep(acct$H_m, 2, a_i * y0, "*")
  W0_orig <- rowSums(payroll0)
  pi0 <- sweep(payroll0, 1, ifelse(W0_orig > 0, W0_orig, 1.0), "/")
  L_res_i <- colSums(payroll0)
  L_nonres_i <- pmax(0.0, a_i * y0 - L_res_i)

  # Linearized labor supply response d L_i / d log w_j:
  # d L_i = eps * L_res_i * dlog w_i + (1/eps - eps) * sum_n W0_n pi0_{ni} sum_j pi0_{nj} dlog w_j + L_nonres_i dlog w_i
  cross_mat <- crossprod(pi0, sweep(pi0, 1, W0_orig, "*")) # Nd x Nd: sum_n W0_n pi0_{ni} pi0_{nj}
  Jac_L <- diag(eps * L_res_i + L_nonres_i, Nd) + (1.0 / eps - eps) * cross_mat
  rhs_L <- a_i * y0 * dlog_y
  dlog_w_fin <- as.numeric(solve(Jac_L, rhs_L))
  p_fin <- a_i * dlog_w_fin + b_i * dlog_y
  cost_fin <- as.numeric(S_res %*% p_fin)
  welf_fin <- dlog_v - cost_fin

  out_rows <- list(
    data.frame(regime = "fixed_factor_share", eps = 0.0,
               mean_living_cost_lp = 100 * sum(pop_w * cost_fixed),
               mean_welfare_lp = 100 * sum(pop_w * welf_fixed)),
    data.frame(regime = "frechet_mobility", eps = eps,
               mean_living_cost_lp = 100 * sum(pop_w * cost_fin),
               mean_welfare_lp = 100 * sum(pop_w * welf_fin))
  )

  if (include_perfect) {
    dlog_w_perf <- sum(a_i * y0 * dlog_y) / sum(a_i * y0)
    p_perf <- a_i * dlog_w_perf + b_i * dlog_y
    cost_perf <- as.numeric(S_res %*% p_perf)
    welf_perf <- dlog_v - cost_perf
    out_rows[[3L]] <- data.frame(
      regime = "perfect_mobility", eps = Inf,
      mean_living_cost_lp = 100 * sum(pop_w * cost_perf),
      mean_welfare_lp = 100 * sum(pop_w * welf_perf)
    )
  }

  welfare$extensions$workplace_mobility <- list(
    summary = do.call(rbind, out_rows),
    dlog_p_frechet = p_fin
  )
  welfare
}

#' Add Consumer and Producer Finite-Change Welfare Bounds (Section 5.2.2)
#'
#' Computes exact finite-change consumer Laspeyres, Paasche, and CES welfare
#' measures (Diewert 1976) alongside C++ bounded-variable least-squares (BVLS)
#' producer-side concave-technology projection bounds ensuring global producer
#' feasibility \eqn{q_i \in [\min(0, y_{i1}/y_{i0}-1), \max(0, y_{i1}/y_{i0}-1)]}.
#'
#' @param welfare An \code{\link{urban_welfare}} object.
#' @param sigma_ces Numeric substitution elasticity for the exact finite-change
#'   CES consumer benchmark (default \code{3.456}).
#' @return The input \code{urban_welfare} object with
#'   \code{extensions$finite_bounds} populated.
#' @export
add_finite_bounds <- function(welfare, sigma_ces = 3.456) {
  if (!inherits(welfare, "urban_welfare")) stop("welfare must be an 'urban_welfare' object.")
  acct <- welfare$account
  p_idx <- welfare$headline$period_index
  is_res <- acct$origins$is_city_resident
  pop_w <- acct$origins$pop[is_res] / sum(acct$origins$pop[is_res])

  dlog_v <- welfare$anchored$dlog_v[is_res, p_idx]
  dlog_p <- welfare$anchored$dlog_p[, p_idx]
  dlog_y <- welfare$anchored$dlog_y[, p_idx]
  Delta_D <- welfare$responses$Delta_D_hat[, p_idx]

  nu_res <- acct$origins$nu[is_res]
  S_NH_res <- acct$S_NH[is_res, , drop = FALSE]
  s_out_res <- acct$s_out_NH[is_res]
  S_tot_res <- acct$S_tot[is_res, , drop = FALSE]

  # 1. Consumer finite-change bounds at exp(dlog_p)
  p_ratio <- exp(dlog_p)
  s_rem0 <- 1.0 - rowSums(S_tot_res) # housing + outside good share at fixed unit price
  L_index <- s_rem0 + as.numeric(S_tot_res %*% p_ratio)
  w_laspeyres <- dlog_v - log(L_index)

  # Post-shock conditional nonhousing shares from estimated Delta_D
  raw_post <- sweep(S_NH_res, 2, exp(Delta_D), "*")
  den_post <- s_out_res + rowSums(raw_post)
  S_tot_post <- sweep(raw_post, 1, nu_res / den_post, "*")
  s_rem1 <- 1.0 - rowSums(S_tot_post)
  P_inv_index <- s_rem1 + as.numeric(S_tot_post %*% (1.0 / p_ratio))
  w_paasche <- dlog_v + log(P_inv_index)

  # Exact CES log price index at exp(dlog_p)
  ces_inner <- s_out_res + as.numeric(S_NH_res %*% exp((1.0 - sigma_ces) * dlog_p))
  logP_ces <- log(ces_inner) / (1.0 - sigma_ces)
  w_ces_exact <- dlog_v - nu_res * logP_ces

  # 2. Producer concave-technology feasible finite-change projection (BVLS in C++)
  v0_res <- acct$origins$v0[is_res]
  y0 <- acct$destinations$y0
  H_res <- acct$H_m[is_res, , drop = FALSE]
  A_mat <- sweep(sweep(H_res, 2, y0, "*"), 1, v0_res, "/")
  y_change <- exp(dlog_y) - 1.0
  lo_q <- pmin(0.0, y_change)
  hi_q <- pmax(0.0, y_change)
  target_v <- exp(dlog_v) - 1.0

  w_sqrt <- sqrt(pop_w)
  Aw <- sweep(A_mat, 1, w_sqrt, "*")
  bw <- w_sqrt * target_v
  q_proj <- cpp_bvls_box(Aw, bw, lo_q, hi_q, max_iter = 3000L, tol = 1e-11)
  v_ratio_proj <- pmax(1e-6, 1.0 + as.numeric(A_mat %*% q_proj))
  # Under concave production with dG_i = y0_i q_i and y1_i = y0_i (1 + y_change_i),
  # finite price ratio p1_i / p0_i satisfies 1 + q_i / (1 + y_change_i) <= p1_i / p0_i <= 1 / (1 - q_i / 1) when bounded
  # Evaluate matched Laspeyres/Paasche producer-feasible interval:
  L_prod <- s_rem0 + as.numeric(S_tot_res %*% (1.0 + q_proj))
  P_prod_inv <- s_rem1 + as.numeric(S_tot_post %*% (1.0 - q_proj / (1.0 + y_change)))
  w_prod_lo <- log(v_ratio_proj) - log(pmax(1e-8, L_prod))
  w_prod_hi <- log(v_ratio_proj) + log(pmax(1e-8, P_prod_inv))

  summary_df <- data.frame(
    measure = c("first_order_log_linear", "consumer_laspeyres_lower",
                "consumer_ces_exact", "consumer_paasche_upper",
                "producer_concave_lower", "producer_concave_upper"),
    mean_welfare_lp = 100 * c(
       welfare$headline$mean_welfare_lp / 100,
      sum(pop_w * w_laspeyres),
      sum(pop_w * w_ces_exact),
      sum(pop_w * w_paasche),
      sum(pop_w * pmin(w_prod_lo, w_prod_hi)),
      sum(pop_w * pmax(w_prod_lo, w_prod_hi))
    ),
    stringsAsFactors = FALSE
  )

  welfare$extensions$finite_bounds <- list(
    summary = summary_df,
    producer_projection_rmse_pp = 100 * sqrt(sum(pop_w * ((v_ratio_proj - 1.0) - target_v)^2))
  )
  welfare
}
