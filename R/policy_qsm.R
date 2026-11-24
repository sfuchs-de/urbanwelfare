#' Solve the Calibrated Nonlinear Open-City Quantitative Spatial Model
#'
#' Solves the two-sector open-city nonlinear QSM calibrated to \code{account}
#' and marginal price shares \code{kappa} for an arbitrary external demand
#' vector \code{ET_shock} or destination ad valorem tax schedule \code{tau}
#' using the C++ damped-Newton and homotopy solver.
#'
#' @param account An \code{\link{urban_account}} object.
#' @param kappa Numeric vector of length \code{N_dest} of marginal price shares
#'   \eqn{\kappa_i^S \in (0, 1)} (defaults to \code{account$destinations$theta_G}).
#' @param ET_shock Optional numeric vector of destination external demand levels
#'   \eqn{E_i^T} (defaults to initial \code{account$destinations$E_ext0}).
#' @param tau Optional numeric vector of destination ad valorem visitor taxes
#'   \eqn{\tau_i \ge 0} (default all zeros).
#' @param sigma_R Resident CES substitution elasticity across destinations (default \code{3.456}).
#' @param sigma_T Visitor CES substitution elasticity across destinations (default \code{7.0}).
#' @param eta_out Visitor outer city-level demand elasticity (default \code{2.5}).
#' @param compute_jacobian Logical: also compute the general-equilibrium tax-base
#'   elasticity matrix \eqn{\mathcal{E}} and conditional welfare derivatives (default \code{FALSE}).
#' @return A list of nonlinear QSM equilibrium prices, sales, incomes, tax
#'   revenue, and exact resident log-utility changes.
#' @export
solve_qsm <- function(account,
                      kappa = NULL,
                      ET_shock = NULL,
                      tau = NULL,
                      sigma_R = 3.456,
                      sigma_T = 7.0,
                      eta_out = 2.5,
                      compute_jacobian = FALSE) {
  if (!inherits(account, "urban_account")) stop("account must be an 'urban_account' object.")
  Nd <- account$N_dest
  if (is.null(kappa)) kappa <- account$destinations$theta_G
  kappa <- pmin(0.99, pmax(0.02, as.numeric(kappa)))
  if (is.null(ET_shock)) ET_shock <- account$destinations$E_ext0
  use_tax <- !is.null(tau)
  theta <- if (use_tax) log(1.0 + pmax(-0.95, as.numeric(tau))) else rep(0.0, Nd)

  is_res_num <- as.numeric(account$origins$is_city_resident)
  pop <- account$origins$pop
  w_welf <- (is_res_num * pop) / sum(is_res_num * pop)
  mu_vec <- 1.0 - account$destinations$theta_G

  res <- cpp_solve_qsm(
    y0 = account$destinations$y0,
    mu = mu_vec,
    kappa = kappa,
    A_fixed = account$destinations$A_fixed0,
    ET_shock = as.numeric(ET_shock),
    v_out = account$origins$v_out,
    nu = account$origins$nu,
    s_out0 = account$s_out_NH,
    S_NH0 = account$S_NH,
    H_alloc = account$H_m,
    pop = pop,
    is_res = is_res_num,
    w_welf = w_welf,
    theta = theta,
    use_tax_demand = use_tax,
    sigma_R = sigma_R,
    sigma_T = sigma_T,
    eta_out = eta_out,
    logp_init = rep(0.0, Nd),
    compute_jacobian = compute_jacobian
  )
  is_res <- account$origins$is_city_resident
  pop_w <- pop[is_res] / sum(pop[is_res])
  res$welfare_lp <- 100 * res$logp * 0 # placeholder replaced below
  res$resident_welfare_lp <- 100 * (log(res$v[is_res] / res$v0[is_res]) - account$origins$nu[is_res] * res$logP_NH[is_res])
  res$mean_welfare_lp <- sum(pop_w * res$resident_welfare_lp)
  res
}

#' Compute Empirical and Nonlinear QSM Optimal Spatial Taxes (Section 6 & Appendix B)
#'
#' Computes optimal visitor/external demand taxes rebated lump-sum to city
#' residents under both the first-order empirical sufficient-statistics index
#' and the exact C++ adjoint nonlinear QSM Ramsey planner (both uniform city-wide
#' tax \eqn{\tau_{\mathrm{unif}}^\star} and place-specific schedule
#' \eqn{\boldsymbol{\tau}^\star \ge \mathbf{0}}), along with the four-channel
#' spatial Ramsey decomposition (equation 6.2).
#'
#' @param welfare An \code{\link{urban_welfare}} object.
#' @param account An \code{\link{urban_account}} object (defaults to \code{welfare$account}).
#' @param demand_elasticities Named list with \code{sigma_T} (visitor across-tract
#'   substitution elasticity, default \code{7.0}) and \code{eta_out} (outer
#'   city-level demand elasticity, default \code{2.5}).
#' @param sigma_R Resident across-tract substitution elasticity for the QSM
#'   benchmark (default \code{3.456}).
#' @param mode Character: \code{"both"} (default), \code{"empirical"}, or
#'   \code{"nonlinear_qsm"}.
#' @param max_tax Upper bound on ad valorem tax rate \eqn{\tau_i} during
#'   numerical optimization (default \code{3.0}, i.e. 300\%).
#'
#' @aliases urban_tax
#' @return An S3 object of class \code{"urban_tax"}.
#' @export
#' @examples
#' data(sim_city_160)
#' acct <- urban_account(sim_city_160$E_R0, sim_city_160$H_0,
#'                       sim_city_160$origins, sim_city_160$destinations)
#' exps <- build_exposures(acct, sim_city_160$shocks)
#' fit  <- estimate_bilateral(sim_city_160$bilateral_panel, exps)
#' welf <- recover_welfare(fit, acct, se_method = "none")
#' tax  <- optimal_spatial_tax(welf, acct, mode = "both")
#' summary(tax)
optimal_spatial_tax <- function(welfare,
                                account = welfare$account,
                                demand_elasticities = list(sigma_T = 7.0, eta_out = 2.5),
                                sigma_R = 3.456,
                                mode = c("both", "empirical", "nonlinear_qsm"),
                                max_tax = 3.0) {
  if (!inherits(welfare, "urban_welfare")) stop("welfare must be an 'urban_welfare' object.")
  mode <- match.arg(mode)
  sigma_T <- demand_elasticities$sigma_T
  eta_out <- demand_elasticities$eta_out
  if (sigma_T <= 1 || eta_out <= 1) {
    stop("sigma_T and eta_out must exceed 1 for interior optimal visitor taxes.")
  }

  Nd <- account$N_dest
  is_res <- account$origins$is_city_resident
  pop <- account$origins$pop
  Omega_pop <- sum(pop[is_res])
  w_welf <- (as.numeric(is_res) * pop) / Omega_pop
  kappa_vec <- welfare$kappa_vec
  ET0 <- account$destinations$E_ext0
  ET0_tot <- sum(ET0)
  if (ET0_tot <= 0) stop("Total external/visitor spending E_ext0 is zero; cannot compute visitor tax.")
  alpha0 <- ET0 / ET0_tot
  y0 <- account$destinations$y0
  v0_res <- account$origins$v0[is_res]

  # 1. Empirical first-order Ramsey tax schedule around initial equilibrium
  emp_schedule <- NULL
  if (mode %in% c("both", "empirical")) {
    # First-order welfare derivative wrt destination log visitor spending e_i = dE^T_i / E^{T,0}_i:
    # shock d_i = tau0_i * de_i => network sales response (I - M)^(-1) diag(tau0)
    # Using anchored producer envelope: dEV_n = sum_j (H_{nj}^m y0_j - E^{R,0}_{nj}) kappa_j dlog y_j
    net_dest_wedge <- as.numeric(crossprod(w_welf[is_res] / v0_res,
                                           sweep(account$H_m[is_res, , drop = FALSE], 2, y0 * kappa_vec, "*") -
                                             sweep(account$E_R0[is_res, , drop = FALSE], 2, kappa_vec, "*")))
    # Propagate through (I - M)^(-1)
    b_T_emp <- ET0 * as.numeric(crossprod(account$I_minus_M_inv, net_dest_wedge / y0)) * y0 / pmax(1e-12, y0)
    # Wait: dlog y / de_j = (I - M)^(-1)_{.j} * (ET0_j / y0_j)
    dlogy_de <- sweep(account$I_minus_M_inv, 2, ET0 / y0, "*")
    b_T_emp <- as.numeric(crossprod(dlogy_de, net_dest_wedge))

    # Marginal social value of lump-sum transfer T around tau=0:
    # Rebate raises resident income by pop_n dT, which induces local spending nu_n S_NH dT
    d_reb <- as.numeric(account$Omega_v[, is_res, drop = FALSE] %*% (pop[is_res] * account$origins$nu[is_res] * (1.0 - account$s_out_NH[is_res]))) * 0
    rhs_reb <- as.numeric(crossprod(account$S_NH[is_res, , drop = FALSE], pop[is_res] * account$origins$nu[is_res])) / y0
    dlogy_dT <- as.numeric(account$I_minus_M_inv %*% rhs_reb)
    mu_T_emp <- sum(w_welf[is_res] * pop[is_res] / v0_res) + sum(net_dest_wedge * dlogy_dT)

    # Partial + price-feedback tax-base elasticity matrix E_emp at tau=0:
    B_T <- (1.0 - sigma_T) * diag(Nd) + (sigma_T - 1.0 - eta_out) * outer(rep(1.0, Nd), alpha0)
    # First-order dlog E^T = B_T dlog p + (B_T - I) dtheta, with dlog p = diag(kappa) (I - M)^(-1) diag(tau0) dlog E^T
    P_fb <- sweep(sweep(account$I_minus_M_inv, 1, kappa_vec, "*"), 2, ET0 / y0, "*")
    E_emp <- -solve(diag(Nd) - B_T %*% P_fb, B_T - diag(Nd))

    # Solve linear Ramsey FOC: (I + diag(ET0)^(-1) E_emp' diag(ET0) - I ...):
    # From eq. (B.3): (1 + tau_j) ET0_j - sum_i E_{ij} tau_i ET0_i = (Omega / mu_T) sum_i E_{ij} b_{T,i}
    # => [(E_emp' - I) diag(ET0)] tau = ET0 - (Omega / mu_T) E_emp' b_T_emp
    A_foc <- sweep(t(E_emp) - diag(Nd), 2, ET0, "*")
    rhs_foc <- ET0 - (Omega_pop / mu_T_emp) * as.numeric(crossprod(E_emp, b_T_emp))
    active <- ET0 > 1e-10 * ET0_tot
    tau_emp <- rep(0.0, Nd)
    if (any(active)) {
      tau_raw <- solve(A_foc[active, active, drop = FALSE], rhs_foc[active])
      tau_emp[active] <- pmin(max_tax, pmax(0.0, as.numeric(tau_raw)))
    }
    e_unif_emp <- as.numeric(crossprod(alpha0, E_emp %*% rep(1.0, Nd)))
    b_tot_emp <- sum(as.numeric(crossprod(E_emp, b_T_emp)))
    tau_unif_emp <- max(0.0, (1.0 - (Omega_pop / (mu_T_emp * ET0_tot)) * b_tot_emp) / max(0.1, e_unif_emp - 1.0))
    emp_schedule <- list(
      tau_place = tau_emp,
      tau_uniform = tau_unif_emp,
      base_weighted_mean = sum(alpha0 * tau_emp)
    )
  }

  # 2. Exact Nonlinear QSM + Adjoint Ramsey Optimal Spatial Tax
  qsm_out <- NULL
  if (mode %in% c("both", "nonlinear_qsm")) {
    mu_vec <- 1.0 - account$destinations$theta_G
    is_res_num <- as.numeric(is_res)
    theta_max <- log(1.0 + max_tax)
    warm_lp <- rep(0.0, Nd)

    eval_qsm_theta <- function(theta_vec, need_jac = FALSE) {
      res <- cpp_solve_qsm(
        y0 = y0, mu = mu_vec, kappa = kappa_vec,
        A_fixed = account$destinations$A_fixed0, ET_shock = ET0,
        v_out = account$origins$v_out, nu = account$origins$nu,
        s_out0 = account$s_out_NH, S_NH0 = account$S_NH, H_alloc = account$H_m,
        pop = pop, is_res = is_res_num, w_welf = w_welf,
        theta = theta_vec, use_tax_demand = TRUE,
        sigma_R = sigma_R, sigma_T = sigma_T, eta_out = eta_out,
        logp_init = warm_lp, compute_jacobian = need_jac
      )
      warm_lp <<- as.numeric(res$logp)
      res
    }

    # 2a. Optimal Uniform Tax (1D L-BFGS-B using sum of adjoint gradient)
    last_th_u <- NULL
    last_res_u <- NULL
    get_u <- function(th_scalar) {
      if (is.null(last_th_u) || abs(th_scalar - last_th_u) > 1e-15) {
        last_th_u <<- th_scalar
        last_res_u <<- eval_qsm_theta(rep(th_scalar, Nd), FALSE)
      }
      last_res_u
    }
    opt_unif <- stats::optim(
      par = log(1.0 + 1.0 / max(0.5, eta_out - 1.0)) * 0.6,
      fn = function(th_s) -100.0 * get_u(th_s)$objective,
      gr = function(th_s) -100.0 * sum(get_u(th_s)$grad_theta),
      method = "L-BFGS-B", lower = 0.0, upper = theta_max,
      control = list(factr = 1e6, maxit = 60L)
    )
    th_unif_star <- opt_unif$par
    tau_unif_star <- exp(th_unif_star) - 1.0
    res_unif <- eval_qsm_theta(rep(th_unif_star, Nd), FALSE)

    # 2b. Optimal Place-Specific Schedule (N_dest-dimensional L-BFGS-B with exact C++ adjoint gradient)
    last_th_p <- NULL
    last_res_p <- NULL
    get_p <- function(th_vec) {
      if (is.null(last_th_p) || max(abs(th_vec - last_th_p)) > 1e-15) {
        last_th_p <<- th_vec
        last_res_p <<- eval_qsm_theta(th_vec, FALSE)
      }
      last_res_p
    }
    opt_place <- stats::optim(
      par = rep(th_unif_star, Nd),
      fn = function(th_v) -100.0 * get_p(th_v)$objective,
      gr = function(th_v) -100.0 * as.numeric(get_p(th_v)$grad_theta),
      method = "L-BFGS-B", lower = rep(0.0, Nd), upper = rep(theta_max, Nd),
      control = list(factr = 1e6, pgtol = 1e-7, maxit = 120L)
    )
    th_place_star <- opt_place$par
    tau_place_star <- exp(th_place_star) - 1.0
    res_place <- eval_qsm_theta(th_place_star, TRUE)

    # 2c. Exact 4-Channel Spatial Ramsey Decomposition at optimum (eq. 6.2 / Appendix B)
    E_ge <- res_place$E_ge
    ET_star <- as.numeric(res_place$ET)
    b_T <- as.numeric(res_place$b_T)
    mu_T <- res_place$mu_T_ge
    diag_E <- diag(E_ge)
    denom_j <- pmax(0.1, diag_E - 1.0) * pmax(1e-12, ET_star)

    ch1_markup <- 1.0 / pmax(0.1, diag_E - 1.0)
    ch2_own_inc <- -(Omega_pop / mu_T) * (diag_E * b_T) / denom_j
    E_off <- E_ge
    diag(E_off) <- 0.0
    ch3_cross_inc <- -(Omega_pop / mu_T) * as.numeric(crossprod(E_off, b_T)) / denom_j
    ch4_fiscal_spill <- as.numeric(crossprod(E_off, tau_place_star * ET_star)) / denom_j

    decomp_df <- data.frame(
      dest = account$dest_ids,
      E_T_base = ET0,
      E_T_opt = ET_star,
      tau_opt = tau_place_star,
      own_markup = ch1_markup,
      own_incidence = ch2_own_inc,
      cross_incidence = ch3_cross_inc,
      fiscal_spillover = ch4_fiscal_spill,
      foc_unconstrained = ch1_markup + ch2_own_inc + ch3_cross_inc + ch4_fiscal_spill,
      stringsAsFactors = FALSE
    )
    if (all(c("x", "y") %in% names(account$destinations))) {
      decomp_df$x <- account$destinations$x
      decomp_df$y <- account$destinations$y
    }

    w_star <- ET_star / sum(ET_star)
    qsm_out <- list(
      tau_uniform = tau_unif_star,
      welfare_gain_uniform_lp = 100.0 * res_unif$objective,
      revenue_uniform = res_unif$revenue,
      tau_place = tau_place_star,
      base_weighted_mean_tau = sum(w_star * tau_place_star),
      p10_tau = .weighted_quantile(tau_place_star, w_star, 0.10),
      p90_tau = .weighted_quantile(tau_place_star, w_star, 0.90),
      welfare_gain_place_lp = 100.0 * res_place$objective,
      uniform_share_of_place_gain = res_unif$objective / max(1e-12, res_place$objective),
      revenue_place = res_place$revenue,
      max_foc_interior_grad = max(abs(res_place$grad_theta[th_place_star > 1e-5 & th_place_star < theta_max - 1e-5]), 0),
      decomposition = decomp_df
    )
  }

  structure(
    list(
      mode = mode,
      demand_elasticities = demand_elasticities,
      sigma_R = sigma_R,
      empirical = emp_schedule,
      qsm = qsm_out,
      account = account
    ),
    class = "urban_tax"
  )
}

#' @method print urban_tax
#' @export
print.urban_tax <- function(x, ...) {
  cat(sprintf("<urban_tax: mode = '%s', sigma_T = %.2f, eta_out = %.2f>\n",
              x$mode, x$demand_elasticities$sigma_T, x$demand_elasticities$eta_out))
  if (!is.null(x$qsm)) {
    q <- x$qsm
    cat(sprintf("  Nonlinear QSM Optimal Uniform Tax     : %6.2f%%  (Welfare gain: +%.3f log pts)\n",
                100 * q$tau_uniform, q$welfare_gain_uniform_lp))
    cat(sprintf("  Nonlinear QSM Place-Specific Mean Tax : %6.2f%%  [P10: %.1f%%, P90: %.1f%%]  (Welfare gain: +%.3f log pts)\n",
                100 * q$base_weighted_mean_tau, 100 * q$p10_tau, 100 * q$p90_tau, q$welfare_gain_place_lp))
    cat(sprintf("  Uniform share of place-specific gain  : %6.2f%%\n",
                100 * q$uniform_share_of_place_gain))
  }
  if (!is.null(x$empirical)) {
    cat(sprintf("  Empirical Linear Ramsey Mean Tax      : %6.2f%%  (Uniform: %.2f%%)\n",
                100 * x$empirical$base_weighted_mean, 100 * x$empirical$tau_uniform))
  }
  invisible(x)
}

#' @method summary urban_tax
#' @export
summary.urban_tax <- function(object, ...) {
  cat("=== Optimal Spatial Tax Summary ===\n")
  print.urban_tax(object)
  if (!is.null(object$qsm)) {
    d <- object$qsm$decomposition
    w <- d$E_T_opt / sum(d$E_T_opt)
    cat("\nBase-Weighted 4-Channel Spatial Ramsey Decomposition (Ad Valorem %):\n")
    cat(sprintf("  1. Own inverse-elasticity markup : %+7.2f%%\n", 100 * sum(w * d$own_markup)))
    cat(sprintf("  2. Own resident incidence        : %+7.2f%%\n", 100 * sum(w * d$own_incidence)))
    cat(sprintf("  3. Cross-tract resident incidence: %+7.2f%%\n", 100 * sum(w * d$cross_incidence)))
    cat(sprintf("  4. Cross-base fiscal spillovers  : %+7.2f%%\n", 100 * sum(w * d$fiscal_spillover)))
  }
  invisible(object)
}
