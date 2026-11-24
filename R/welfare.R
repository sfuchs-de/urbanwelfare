#' Recover Producer Price Shares, Common Income Level, and Resident Welfare (Steps 2 & 3)
#'
#' Translates estimated bilateral spending slopes into neighborhood nominal
#' income growth, destination sales and producer price responses, the open-city
#' common income level \eqn{\widehat{c}_t} (Proposition 2.5), and first-order
#' money-metric resident welfare \eqn{\mathrm{d}EV_{nt} / v_n^0}. Also verifies
#' the open-city aggregate terms-of-trade identity and computes two-step
#' bootstrap or delta-method confidence intervals.
#'
#' @param fit An \code{\link{urban_bilateral}} object from \code{\link{estimate_bilateral}}.
#' @param account An \code{\link{urban_account}} object.
#' @param price_shares Either \code{"estimate_by_group"} (default, estimates
#'   group-specific marginal price shares \eqn{\kappa_g^S} via two-way FE
#'   spatial factor-income matching, eq. 2.7), \code{"estimate_common"} (single
#'   city-wide \eqn{\kappa^S}), or a numeric vector of length \code{N_dest} (or
#'   named by \code{sector_group}) supplying \eqn{\kappa_i^S \in (0, 1]} directly.
#' @param level_anchor Character: \code{"open_city_factor_balance"} (default,
#'   Proposition 2.5's \eqn{\widehat{c}_t = -\psi_t / (1 - \vartheta)}) or
#'   \code{"zero"} (cross-sectional differences at \eqn{c_t = 0}).
#' @param period Integer or character selecting the focal reporting period for
#'   cross-sectional summaries and plots (defaults to the period with largest
#'   mean external shock \eqn{\bar{d}_t}).
#' @param kappa_periods Optional integer vector of column indices of periods used
#'   when estimating \eqn{\kappa_g^S} (defaults to all periods).
#' @param se_method Character: \code{"bootstrap"} (default, parametric redraws of
#'   \eqn{\widehat{\boldsymbol{\beta}}} and \eqn{\widehat{\boldsymbol{\kappa}}^S}),
#'   \code{"delta"}, or \code{"none"}.
#' @param n_boot Number of parametric bootstrap draws (default 300).
#' @param seed Optional integer random seed for bootstrap reproducibility.
#'
#' @aliases urban_welfare
#' @return An S3 object of class \code{"urban_welfare"}.
#' @export
#' @examples
#' data(sim_city_160)
#' acct <- urban_account(sim_city_160$E_R0, sim_city_160$H_0,
#'                       sim_city_160$origins, sim_city_160$destinations)
#' exps <- build_exposures(acct, sim_city_160$shocks)
#' fit  <- estimate_bilateral(sim_city_160$bilateral_panel, exps)
#' welf <- recover_welfare(fit, acct, seed = 1)
#' summary(welf)
recover_welfare <- function(fit,
                            account,
                            price_shares = "estimate_by_group",
                            level_anchor = c("open_city_factor_balance", "zero"),
                            period = NULL,
                            kappa_periods = NULL,
                            se_method = c("bootstrap", "delta", "none"),
                            n_boot = 300L,
                            seed = 42L) {
  if (!inherits(fit, "urban_bilateral")) stop("fit must be an 'urban_bilateral' object.")
  if (!inherits(account, "urban_account")) stop("account must be an 'urban_account' object.")
  level_anchor <- match.arg(level_anchor)
  se_method <- match.arg(se_method)
  exps <- fit$exposures
  Nt <- exps$N_time

  if (is.null(period)) {
    p_idx <- which.max(colMeans(exps$d))
  } else if (is.character(period)) {
    p_idx <- match(period, exps$time_ids)
    if (is.na(p_idx)) stop("period not found in exposures$time_ids.")
  } else {
    p_idx <- as.integer(period)
  }
  if (is.null(kappa_periods)) kappa_periods <- seq_len(Nt)

  # Step 2a: Relative origin income & destination sales at c_t = 0
  resp <- .compute_relative_responses(fit$coefficients, account, exps)

  # Step 2b: Estimate or assign destination price shares kappa^S_i
  kap_fit <- .resolve_kappa(price_shares, account, resp$dlog_v_hat, resp$dlog_y_hat, kappa_periods)
  kappa_vec <- kap_fit$kappa_vec

  # Step 3: Common level c_t (Proposition 2.5) and welfare evaluation
  eval_out <- .evaluate_anchored_welfare(account, resp, kappa_vec, level_anchor)

  is_res <- account$origins$is_city_resident
  pop_res <- account$origins$pop[is_res]
  pop_weights <- pop_res / sum(pop_res)

  # Peak/focal period tract-level summary data.frame for city residents
  w_p <- 100 * eval_out$welfare[is_res, p_idx]
  inc_p <- 100 * eval_out$dlog_v[is_res, p_idx]
  cost_p <- 100 * eval_out$living_cost[is_res, p_idx]
  w_zero_p <- 100 * eval_out$welfare_zero[is_res, p_idx]

  tract_df <- data.frame(
    id = account$orig_ids[is_res],
    pop = pop_res,
    v0 = account$origins$v0[is_res],
    income_lp = inc_p,
    living_cost_lp = cost_p,
    welfare_lp = w_p,
    welfare_unanchored_lp = w_zero_p,
    stringsAsFactors = FALSE
  )
  if (all(c("x", "y") %in% names(account$origins))) {
    tract_df$x <- account$origins$x[is_res]
    tract_df$y <- account$origins$y[is_res]
  }

  headline <- list(
    period = exps$time_ids[p_idx],
    period_index = p_idx,
    mean_income_lp = sum(pop_weights * inc_p),
    mean_living_cost_lp = sum(pop_weights * cost_p),
    mean_welfare_lp = sum(pop_weights * w_p),
    mean_welfare_unanchored_lp = sum(pop_weights * w_zero_p),
    common_shift_c_lp = 100 * eval_out$c_t[p_idx],
    vartheta = eval_out$vartheta,
    p10_welfare_lp = .weighted_quantile(w_p, pop_weights, 0.10),
    p50_welfare_lp = .weighted_quantile(w_p, pop_weights, 0.50),
    p90_welfare_lp = .weighted_quantile(w_p, pop_weights, 0.90),
    terms_of_trade_gap = eval_out$tot_check$max_abs_gap
  )

  # Confidence intervals for mean resident welfare in focal period
  ci <- c(lower = NA_real_, upper = NA_real_, se = NA_real_)
  if (se_method == "bootstrap" && n_boot > 1L) {
    if (!is.null(seed)) set.seed(as.integer(seed))
    b_draws <- .rmvnorm_safe(n_boot, fit$coefficients, fit$vcov)
    k_noise <- if (!is.null(kap_fit$vcov)) {
      .rmvnorm_safe(n_boot, rep(0.0, nrow(kap_fit$vcov)), kap_fit$vcov)
    } else {
      matrix(0.0, nrow = n_boot, ncol = length(kap_fit$coef))
    }
    boot_means <- numeric(n_boot)
    for (b in seq_len(n_boot)) {
      r_b <- .compute_relative_responses(b_draws[b, ], account, exps)
      if (kap_fit$estimated) {
        kf_b <- .resolve_kappa(price_shares, account, r_b$dlog_v_hat, r_b$dlog_y_hat, kappa_periods, compute_vcov = FALSE)
        c_new <- kf_b$coef + k_noise[b, ]
        kv_b <- .map_group_kappa(c_new, account, kap_fit$mode)
      } else {
        kv_b <- kappa_vec
      }
      ev_b <- .evaluate_anchored_welfare(account, r_b, kv_b, level_anchor)
      boot_means[b] <- 100 * sum(pop_weights * ev_b$welfare[is_res, p_idx])
    }
    qs <- stats::quantile(boot_means, probs = c(0.025, 0.975), names = FALSE, na.rm = TRUE)
    ci <- c(lower = qs[1], upper = qs[2], se = stats::sd(boot_means, na.rm = TRUE))
  } else if (se_method == "delta") {
    # Numerical gradient wrt 4 bilateral slopes holding kappa_vec fixed
    base_val <- headline$mean_welfare_lp
    grad <- numeric(4L)
    eps_d <- 1e-5
    for (k in 1:4) {
      b_up <- fit$coefficients
      b_up[k] <- b_up[k] + eps_d
      r_u <- .compute_relative_responses(b_up, account, exps)
      ev_u <- .evaluate_anchored_welfare(account, r_u, kappa_vec, level_anchor)
      val_u <- 100 * sum(pop_weights * ev_u$welfare[is_res, p_idx])
      grad[k] <- (val_u - base_val) / eps_d
    }
    se_val <- sqrt(max(0.0, as.numeric(t(grad) %*% fit$vcov %*% grad)))
    ci <- c(lower = base_val - 1.959964 * se_val, upper = base_val + 1.959964 * se_val, se = se_val)
  }
  headline$ci_lower_lp <- unname(ci["lower"])
  headline$ci_upper_lp <- unname(ci["upper"])
  headline$se_lp <- unname(ci["se"])

  # Decile summary table
  deciles <- .build_decile_table(tract_df, pop_weights)

  structure(
    list(
      fit = fit,
      account = account,
      level_anchor = level_anchor,
      kappa_fit = kap_fit,
      kappa_vec = kappa_vec,
      responses = resp,
      anchored = eval_out,
      headline = headline,
      tract_summary = tract_df,
      deciles = deciles,
      extensions = list()
    ),
    class = "urban_welfare"
  )
}

.compute_relative_responses <- function(b4, account, exps) {
  a_O <- unname(b4[1])
  b_v <- unname(b4[2])
  a_D <- unname(b4[3])
  b_D <- unname(b4[4])

  # Relative origin income growth: d log v_hat = a_O h_orig + b_v Fv
  dlog_v_hat <- a_O * exps$h_orig + b_v * exps$Fv
  Delta_D_hat <- a_D * exps$h + b_D * exps$Fy

  # Exact analytical sum_n Omega^v_{in} d log E^R_{nit}:
  # d log E^R_{nit} = dlog_v_hat[n, t] + Delta_D_hat[i, t] - bar_Delta[n, t]
  bar_Delta <- account$S_NH %*% Delta_D_hat
  orig_shift <- dlog_v_hat - bar_Delta
  Omega_row_sum <- rowSums(account$Omega_v)
  resident_sales_contrib <- (account$Omega_v %*% orig_shift) + sweep(Delta_D_hat, 1, Omega_row_sum, "*")
  dlog_y_hat <- exps$d + resident_sales_contrib

  list(
    dlog_v_hat = dlog_v_hat,
    dlog_y_hat = dlog_y_hat,
    Delta_D_hat = Delta_D_hat
  )
}

.resolve_kappa <- function(price_shares, account, dlog_v_hat, dlog_y_hat, kappa_periods, compute_vcov = TRUE) {
  Nd <- account$N_dest
  is_res <- account$origins$is_city_resident
  v0_res <- account$origins$v0[is_res]
  H_res <- account$H_m[is_res, , drop = FALSE]
  y0 <- account$destinations$y0
  groups <- account$destinations$sector_group
  grp_levels <- account$group_levels

  if (is.numeric(price_shares)) {
    if (length(price_shares) == Nd) {
      kv <- as.numeric(price_shares)
    } else if (!is.null(names(price_shares)) && all(as.character(grp_levels) %in% names(price_shares))) {
      kv <- as.numeric(price_shares[as.character(groups)])
    } else if (length(price_shares) == 1L) {
      kv <- rep(as.numeric(price_shares), Nd)
    } else {
      stop("Numeric price_shares must have length 1, length N_dest, or names matching sector_group.")
    }
    return(list(estimated = FALSE, mode = "user", coef = price_shares, vcov = NULL, kappa_vec = kv))
  }

  mode <- match.arg(price_shares, c("estimate_by_group", "estimate_common"))
  Y_sub <- dlog_v_hat[is_res, kappa_periods, drop = FALSE]
  Yd <- .demean_twoway_mat(Y_sub)
  No_res <- nrow(Y_sub)
  Nt_sub <- ncol(Y_sub)

  use_groups <- if (mode == "estimate_common") list(all = rep(TRUE, Nd)) else {
    setNames(lapply(grp_levels, function(g) groups == g), as.character(grp_levels))
  }
  K <- length(use_groups)
  Xd_list <- vector("list", K)
  X_mat <- matrix(0.0, nrow = No_res * Nt_sub, ncol = K)
  for (k in seq_len(K)) {
    mask <- as.numeric(use_groups[[k]])
    sales_g <- sweep(dlog_y_hat[, kappa_periods, drop = FALSE], 1, y0 * mask, "*")
    X_raw <- sweep(H_res %*% sales_g, 1, v0_res, "/")
    Xd_k <- .demean_twoway_mat(X_raw)
    Xd_list[[k]] <- Xd_k
    X_mat[, k] <- as.vector(Xd_k)
  }

  y_vec <- as.vector(Yd)
  XX <- crossprod(X_mat)
  cf <- as.numeric(solve(XX, crossprod(X_mat, y_vec)))
  names(cf) <- names(use_groups)
  kv <- .map_group_kappa(cf, account, mode)

  vc <- NULL
  if (compute_vcov) {
    res_mat <- matrix(y_vec - X_mat %*% cf, nrow = No_res, ncol = Nt_sub)
    om_scores <- matrix(0.0, nrow = No_res, ncol = K)
    tm_scores <- matrix(0.0, nrow = Nt_sub, ncol = K)
    cell_meat <- matrix(0.0, nrow = K, ncol = K)
    for (k in seq_len(K)) {
      sc_k <- Xd_list[[k]] * res_mat
      om_scores[, k] <- rowSums(sc_k)
      tm_scores[, k] <- colSums(sc_k)
    }
    for (k1 in seq_len(K)) {
      for (k2 in seq_len(K)) {
        cell_meat[k1, k2] <- sum((Xd_list[[k1]] * res_mat) * (Xd_list[[k2]] * res_mat))
      }
    }
    meat <- crossprod(om_scores) + crossprod(tm_scores) - cell_meat
    G_min <- max(2L, min(No_res, Nt_sub))
    n_tot <- length(y_vec)
    bread <- solve(XX)
    vc <- ((n_tot - 1) / max(1, n_tot - K)) * (G_min / (G_min - 1)) * (bread %*% meat %*% bread)
    eig <- eigen(0.5 * (vc + t(vc)), symmetric = TRUE)
    vc <- eig$vectors %*% (pmax(eig$values, 0.0) * t(eig$vectors))
    dimnames(vc) <- list(names(cf), names(cf))
  }

  list(estimated = TRUE, mode = mode, coef = cf, vcov = vc, kappa_vec = kv)
}

.map_group_kappa <- function(cf, account, mode) {
  if (mode == "estimate_common") {
    kv <- rep(unname(cf[1]), account$N_dest)
  } else {
    kv <- as.numeric(cf[as.character(account$destinations$sector_group)])
  }
  pmin(0.999, pmax(0.01, kv))
}

.evaluate_anchored_welfare <- function(account, resp, kappa_vec, level_anchor) {
  is_res <- account$origins$is_city_resident
  v0 <- account$origins$v0
  v0_res <- v0[is_res]
  V_R0 <- sum(v0_res)
  y0 <- account$destinations$y0
  H_res <- account$H_m[is_res, , drop = FALSE]
  H_col_res <- colSums(H_res)

  # Omega_i^v summed over city residents receiving c_t
  Omega_res_i <- rowSums(account$Omega_v[, is_res, drop = FALSE])
  vartheta <- sum(H_col_res * y0 * kappa_vec * Omega_res_i) / V_R0

  # Unanchored price & welfare
  dlog_p_zero <- sweep(resp$dlog_y_hat, 1, kappa_vec, "*")
  living_cost_zero <- account$S_tot %*% dlog_p_zero
  welfare_zero <- resp$dlog_v_hat - living_cost_zero

  if (level_anchor == "zero") {
    c_t <- rep(0.0, ncol(resp$dlog_v_hat))
  } else {
    if (vartheta >= 1.0) {
      stop(sprintf("Open-city recycling share vartheta = %.4f >= 1; common shift c_t is not uniquely anchored.", vartheta))
    }
    lhs_zero <- colSums(sweep(resp$dlog_v_hat[is_res, , drop = FALSE], 1, v0_res, "*"))
    rhs_zero <- as.numeric(crossprod(H_col_res * y0 * kappa_vec, resp$dlog_y_hat))
    psi_t <- (lhs_zero - rhs_zero) / V_R0
    c_t <- -psi_t / (1.0 - vartheta)
  }

  # Add common shift c_t to city residents
  dlog_v <- resp$dlog_v_hat
  dlog_v[is_res, ] <- sweep(dlog_v[is_res, , drop = FALSE], 2, c_t, "+")

  dlog_y <- resp$dlog_y_hat + outer(Omega_res_i, c_t)
  dlog_p <- sweep(dlog_y, 1, kappa_vec, "*")
  living_cost <- account$S_tot %*% dlog_p
  welfare <- dlog_v - living_cost

  # Verify open-city aggregate terms-of-trade identity (eq. 2.14) when anchored
  # Sum_{n in R} dEV_nt = Sum_i E^{ext,tot}_i dlog p_it - Sum_i (1 - H_{res,i}) y0_i dlog p_it
  # where E^{ext,tot}_i = y0_i - sum_{n in R} E_{ni}^{R,0}
  E_nonres_i <- y0 - colSums(account$E_R0[is_res, , drop = FALSE])
  sum_dEV <- colSums(sweep(welfare[is_res, , drop = FALSE], 1, v0_res, "*"))
  tot_rhs <- as.numeric(crossprod(E_nonres_i - (1.0 - H_col_res) * y0, dlog_p))
  tot_gap <- if (level_anchor == "open_city_factor_balance") max(abs(sum_dEV - tot_rhs)) / V_R0 else NA_real_

  list(
    vartheta = vartheta,
    c_t = c_t,
    dlog_v = dlog_v,
    dlog_y = dlog_y,
    dlog_p = dlog_p,
    living_cost = living_cost,
    welfare = welfare,
    welfare_zero = welfare_zero,
    tot_check = list(sum_dEV = sum_dEV, tot_rhs = tot_rhs, max_abs_gap = tot_gap)
  )
}

.weighted_quantile <- function(x, w, probs) {
  ord <- order(x)
  xo <- x[ord]
  wo <- w[ord] / sum(w)
  cw <- cumsum(wo)
  sapply(probs, function(p) xo[which(cw >= p)[1]])
}

.build_decile_table <- function(tract_df, weights) {
  cuts <- .weighted_quantile(tract_df$welfare_lp, weights, seq(0, 1, by = 0.1))
  cuts[1] <- -Inf
  cuts[length(cuts)] <- Inf
  cuts_u <- unique(cuts)
  if (length(cuts_u) < 3L) {
    return(data.frame(decile = 1L, income_lp = sum(weights * tract_df$income_lp),
                      living_cost_lp = sum(weights * tract_df$living_cost_lp),
                      welfare_lp = sum(weights * tract_df$welfare_lp)))
  }
  bins <- cut(tract_df$welfare_lp, breaks = cuts_u, include.lowest = TRUE, labels = FALSE)
  out <- lapply(sort(unique(bins)), function(d) {
    idx <- which(bins == d)
    wd <- weights[idx] / sum(weights[idx])
    data.frame(
      decile = d,
      n_tracts = length(idx),
      income_lp = sum(wd * tract_df$income_lp[idx]),
      living_cost_lp = sum(wd * tract_df$living_cost_lp[idx]),
      welfare_lp = sum(wd * tract_df$welfare_lp[idx])
    )
  })
  do.call(rbind, out)
}

.rmvnorm_safe <- function(n, mu, Sigma) {
  k <- length(mu)
  eig <- eigen(0.5 * (Sigma + t(Sigma)), symmetric = TRUE)
  vals <- pmax(eig$values, 0.0)
  R <- eig$vectors %*% (sqrt(vals) * t(eig$vectors))
  Z <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)
   sweep(Z %*% R, 2, mu, "+")
}

#' @export
as.data.frame.urban_welfare <- function(x, ...) {
  x$tract_summary
}

#' @method print urban_welfare
#' @export
print.urban_welfare <- function(x, ...) {
  h <- x$headline
  cat(sprintf("<urban_welfare: period = '%s', anchor = '%s'>\n", h$period, x$level_anchor))
  cat(sprintf("  Mean Resident Welfare : %6.2f log pts", h$mean_welfare_lp))
  if (is.finite(h$ci_lower_lp)) {
    cat(sprintf("  (95%% CI: [%.2f, %.2f])", h$ci_lower_lp, h$ci_upper_lp))
  }
  cat("\n")
  cat(sprintf("  Income vs Living Cost : %6.2f income  -  %6.2f living cost  (c_t = %.2f, vartheta = %.3f)\n",
              h$mean_income_lp, h$mean_living_cost_lp, h$common_shift_c_lp, h$vartheta))
  cat(sprintf("  Tract P10 / P50 / P90 : %6.2f / %6.2f / %6.2f log pts\n",
              h$p10_welfare_lp, h$p50_welfare_lp, h$p90_welfare_lp))
  invisible(x)
}

#' @method summary urban_welfare
#' @export
summary.urban_welfare <- function(object, ...) {
  cat("=== Resident Welfare Incidence Summary ===\n")
  print.urban_welfare(object)
  cat("\nEstimated / Assigned Price Shares (kappa^S):\n")
  print(round(object$kappa_fit$coef, 4))
  if (length(object$extensions) > 0L) {
    cat("\nAttached Extensions:\n")
    for (nm in names(object$extensions)) {
      cat("  -", nm, "\n")
    }
  }
  invisible(object)
}
