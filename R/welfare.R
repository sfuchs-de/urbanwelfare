#' Recover Producer Price Shares, Common Income Level, and Resident Welfare (Steps 2 & 3)
#'
#' Translates estimated bilateral spending slopes into neighborhood nominal
#' income growth, destination sales and producer price responses, the open-city
#' common income level \eqn{\widehat{c}_t} (Proposition 2.5), and first-order
#' money-metric resident welfare \eqn{\mathrm{d}EV_{nt} / v_n^0}. Supports
#' unbalanced weighted Step-2 matching samples (\code{kappa_support}),
#' non-homothetic marginal propensities to spend (\code{mpc_nonhousing}),
#' upstream outside-sector retention (\code{psi_O}), and two-period contrasts
#' (\code{period} vs. \code{base_period}).
#'
#' @param fit An \code{\link{urban_bilateral}} object from
#'   \code{\link{estimate_bilateral}} or \code{\link{bilateral_from_estimates}}.
#' @param account An \code{\link{urban_account}} object.
#' @param price_shares Either \code{"estimate_by_group"} (default, estimates
#'   group-specific marginal price shares \eqn{\kappa_g^S} via two-way FE
#'   spatial factor-income matching, eq. 2.7), \code{"estimate_common"} (single
#'   city-wide \eqn{\kappa^S}), or a numeric vector of length \code{N_dest} (or
#'   named by \code{sector_group}) supplying \eqn{\kappa_i^S \in (0, 1]} directly.
#' @param level_anchor Character: \code{"open_city_factor_balance"} (default,
#'   Proposition 2.5's \eqn{\widehat{c}_t = -\psi_t / (1 - \vartheta)}) or
#'   \code{"zero"} (cross-sectional differences at \eqn{c_t = 0}).
#' @param period Integer or character selecting the focal reporting period
#'   (defaults to the period with largest mean external shock \eqn{\bar{d}_t}).
#' @param base_period Optional integer or character selecting a reference period
#'   to difference against (when \code{NULL}, responses are reported as stored,
#'   e.g. relative to the already-normalized base month).
#' @param kappa_periods Optional integer vector of column indices of periods used
#'   when estimating \eqn{\kappa_g^S} on a balanced panel (defaults to all periods).
#' @param kappa_support Optional data frame with columns \code{origin} (or \code{id})
#'   and \code{time} (or \code{t}), plus optional \code{weight} (or
#'   \code{estimation_weight}) and additional control columns, specifying an
#'   exact unbalanced WLS/FE support for Step-2 \eqn{\kappa_g^S} estimation.
#' @param mpc_nonhousing Optional scalar or length-\code{N_orig} vector giving
#'   the marginal propensity to spend on non-housing goods \eqn{\mu_n = \mathrm{d}x_n / \mathrm{d}v_n}
#'   (defaults to \code{account$origins$nu}, i.e. \eqn{\mu_n = \nu_n}).
#' @param psi_O Numeric scalar in \eqn{[0, 1)} for upstream outside-sector factor
#'   retention per euro of flexible input expansion (default \code{0.0}).
#' @param se_method Character: \code{"bootstrap"} (default, parametric redraws of
#'   \eqn{\widehat{\boldsymbol{\beta}}} and \eqn{\widehat{\boldsymbol{\kappa}}^S}),
#'   \code{"delta"}, or \code{"none"}.
#' @param n_boot Number of parametric bootstrap draws (default 300).
#' @param seed Optional integer random seed for bootstrap reproducibility.
#'
#' @aliases urban_welfare
#' @return An S3 object of class \code{"urban_welfare"}.
#' @export
recover_welfare <- function(fit,
                            account,
                            price_shares = "estimate_by_group",
                            level_anchor = c("open_city_factor_balance", "zero"),
                            period = NULL,
                            base_period = NULL,
                            kappa_periods = NULL,
                            kappa_support = NULL,
                            mpc_nonhousing = NULL,
                            psi_O = 0.0,
                            se_method = c("bootstrap", "delta", "none"),
                            n_boot = 300L,
                            seed = 42L) {
  if (!inherits(fit, "urban_bilateral")) stop("fit must be an 'urban_bilateral' object.")
  if (!inherits(account, "urban_account")) stop("account must be an 'urban_account' object.")
  level_anchor <- match.arg(level_anchor)
  se_method <- match.arg(se_method)
  exps <- fit$exposures
  Nt <- exps$N_time

  p_idx <- .resolve_period_idx(period, exps, default_peak = TRUE)
  b_idx <- if (!is.null(base_period)) .resolve_period_idx(base_period, exps, default_peak = FALSE) else NULL
  if (is.null(kappa_periods)) kappa_periods <- seq_len(Nt)

  nu_vec <- account$origins$nu
  mu_vec <- if (is.null(mpc_nonhousing)) nu_vec else {
    if (length(mpc_nonhousing) == 1L) rep(as.numeric(mpc_nonhousing), account$N_orig) else as.numeric(mpc_nonhousing)
  }

  # Step 2a: Relative origin income & destination sales at c_t = 0
  resp <- .compute_relative_responses(fit$coefficients, account, exps, nu_vec, mu_vec)

  # Step 2b: Estimate or assign destination price shares kappa^S_i
  kap_fit <- .resolve_kappa(price_shares, account, resp$dlog_v_hat, resp$dlog_y_hat,
                            kappa_periods, kappa_support, exps$time_ids)
  kappa_vec <- kap_fit$kappa_vec

  # Step 3: Common level c_t (Proposition 2.5) and welfare evaluation
  eval_out <- .evaluate_anchored_welfare(account, resp, kappa_vec, level_anchor, nu_vec, mu_vec, psi_O)

  is_res <- account$origins$is_city_resident
  pop_res <- account$origins$pop[is_res]
  pop_weights <- pop_res / sum(pop_res)

  extract_contrast <- function(ev, col_p, col_b) {
    if (is.null(col_b)) {
      list(
        w = 100 * ev$welfare[is_res, col_p],
        inc = 100 * ev$dlog_v[is_res, col_p],
        cost = 100 * ev$living_cost[is_res, col_p],
        w_zero = 100 * ev$welfare_zero[is_res, col_p],
        c_lp = 100 * ev$c_t[col_p]
      )
    } else {
      list(
        w = 100 * (ev$welfare[is_res, col_p] - ev$welfare[is_res, col_b]),
        inc = 100 * (ev$dlog_v[is_res, col_p] - ev$dlog_v[is_res, col_b]),
        cost = 100 * (ev$living_cost[is_res, col_p] - ev$living_cost[is_res, col_b]),
        w_zero = 100 * (ev$welfare_zero[is_res, col_p] - ev$welfare_zero[is_res, col_b]),
        c_lp = 100 * (ev$c_t[col_p] - ev$c_t[col_b])
      )
    }
  }

  ct <- extract_contrast(eval_out, p_idx, b_idx)
  tract_df <- data.frame(
    id = account$orig_ids[is_res],
    pop = pop_res,
    v0 = account$origins$v0[is_res],
    income_lp = ct$inc,
    living_cost_lp = ct$cost,
    welfare_lp = ct$w,
    welfare_unanchored_lp = ct$w_zero,
    stringsAsFactors = FALSE
  )
  if (all(c("x", "y") %in% names(account$origins))) {
    tract_df$x <- account$origins$x[is_res]
    tract_df$y <- account$origins$y[is_res]
  }

  period_label <- if (is.null(b_idx)) exps$time_ids[p_idx] else paste0(exps$time_ids[p_idx], " - ", exps$time_ids[b_idx])
  headline <- list(
    period = period_label,
    period_index = p_idx,
    base_period_index = b_idx,
    mean_income_lp = sum(pop_weights * ct$inc),
    mean_living_cost_lp = sum(pop_weights * ct$cost),
    mean_welfare_lp = sum(pop_weights * ct$w),
    mean_welfare_unanchored_lp = sum(pop_weights * ct$w_zero),
    common_shift_c_lp = ct$c_lp,
    vartheta = eval_out$vartheta,
    p10_welfare_lp = .weighted_quantile(ct$w, pop_weights, 0.10),
    p50_welfare_lp = .weighted_quantile(ct$w, pop_weights, 0.50),
    p90_welfare_lp = .weighted_quantile(ct$w, pop_weights, 0.90),
    terms_of_trade_gap = eval_out$tot_check$max_abs_gap
  )

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
      r_b <- .compute_relative_responses(b_draws[b, ], account, exps, nu_vec, mu_vec)
      if (kap_fit$estimated) {
        kf_b <- .resolve_kappa(price_shares, account, r_b$dlog_v_hat, r_b$dlog_y_hat,
                               kappa_periods, kappa_support, exps$time_ids, compute_vcov = FALSE)
        c_new <- kf_b$coef + k_noise[b, ]
        kv_b <- .map_group_kappa(c_new, account, kap_fit$mode, clamp = FALSE)
      } else {
        kv_b <- kappa_vec
      }
      ev_b <- .evaluate_anchored_welfare(account, r_b, kv_b, level_anchor, nu_vec, mu_vec, psi_O)
      ct_b <- extract_contrast(ev_b, p_idx, b_idx)
      boot_means[b] <- sum(pop_weights * ct_b$w)
    }
    qs <- stats::quantile(boot_means, probs = c(0.025, 0.975), names = FALSE, na.rm = TRUE)
    ci <- c(lower = qs[1], upper = qs[2], se = stats::sd(boot_means, na.rm = TRUE))
  } else if (se_method == "delta") {
    base_val <- headline$mean_welfare_lp
    grad <- numeric(4L)
    eps_d <- 1e-5
    for (k in 1:4) {
      b_up <- fit$coefficients
      b_up[k] <- b_up[k] + eps_d
      r_u <- .compute_relative_responses(b_up, account, exps, nu_vec, mu_vec)
      ev_u <- .evaluate_anchored_welfare(account, r_u, kappa_vec, level_anchor, nu_vec, mu_vec, psi_O)
      ct_u <- extract_contrast(ev_u, p_idx, b_idx)
      grad[k] <- (sum(pop_weights * ct_u$w) - base_val) / eps_d
    }
    se_val <- sqrt(max(0.0, as.numeric(t(grad) %*% fit$vcov %*% grad)))
    ci <- c(lower = base_val - 1.959964 * se_val, upper = base_val + 1.959964 * se_val, se = se_val)
  }
  headline$ci_lower_lp <- unname(ci["lower"])
  headline$ci_upper_lp <- unname(ci["upper"])
  headline$se_lp <- unname(ci["se"])

  deciles <- .build_decile_table(tract_df, pop_weights)

  structure(
    list(
      fit = fit,
      account = account,
      level_anchor = level_anchor,
      mpc_nonhousing = mu_vec,
      psi_O = psi_O,
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

.resolve_period_idx <- function(period, exps, default_peak = TRUE) {
  if (is.null(period)) {
    return(if (default_peak) which.max(colMeans(exps$d)) else 1L)
  }
  if (is.character(period)) {
    idx <- match(period, exps$time_ids)
    if (is.na(idx)) stop(sprintf("Period '%s' not found in exposures$time_ids.", period))
    return(idx)
  }
  as.integer(period)
}

.compute_relative_responses <- function(b4, account, exps, nu_vec, mu_vec) {
  a_O <- unname(b4[1])
  b_v <- unname(b4[2])
  a_D <- unname(b4[3])
  b_D <- unname(b4[4])

  # Non-housing spending shift dlog x_nt = a_O h_orig + b_v Fv
  # Note: non-city metro origins do not receive city origin direct shock a_O h_orig
  is_res <- account$origins$is_city_resident
  h_orig_eff <- exps$h_orig
  if (any(!is_res)) h_orig_eff[!is_res, ] <- 0.0
  dlog_x_hat <- a_O * h_orig_eff + b_v * exps$Fv
  dlog_v_hat <- sweep(dlog_x_hat, 1, nu_vec / mu_vec, "*")
  Delta_D_hat <- a_D * exps$h + b_D * exps$Fy

  bar_Delta <- account$S_NH %*% Delta_D_hat
  orig_shift <- dlog_x_hat - bar_Delta
  Omega_row_sum <- rowSums(account$Omega_v)
  resident_sales_contrib <- (account$Omega_v %*% orig_shift) + sweep(Delta_D_hat, 1, Omega_row_sum, "*")
  dlog_y_hat <- exps$d + resident_sales_contrib

  list(
    dlog_x_hat = dlog_x_hat,
    dlog_v_hat = dlog_v_hat,
    dlog_y_hat = dlog_y_hat,
    Delta_D_hat = Delta_D_hat
  )
}

.resolve_kappa <- function(price_shares, account, dlog_v_hat, dlog_y_hat,
                           kappa_periods, kappa_support = NULL, time_ids = NULL,
                           compute_vcov = TRUE) {
  Nd <- account$N_dest
  is_res <- account$origins$is_city_resident
  res_ids <- account$orig_ids[is_res]
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
  use_groups <- if (mode == "estimate_common") list(all = rep(TRUE, Nd)) else {
    setNames(lapply(grp_levels, function(g) groups == g), as.character(grp_levels))
  }
  K <- length(use_groups)

  # Unbalanced weighted WLS/FE support when kappa_support is supplied
  if (!is.null(kappa_support)) {
    sup <- as.data.frame(kappa_support)
    id_col <- if ("origin" %in% names(sup)) "origin" else if ("id" %in% names(sup)) "id" else stop("kappa_support needs 'origin' or 'id'.")
    tm_col <- if ("time" %in% names(sup)) "time" else if ("t" %in% names(sup)) "t" else stop("kappa_support needs 'time' or 't'.")
    wt_col <- if ("weight" %in% names(sup)) "weight" else if ("estimation_weight" %in% names(sup)) "estimation_weight" else NULL
    ctrl_cols <- setdiff(names(sup), c(id_col, tm_col, wt_col))

    r_ix <- match(as.character(sup[[id_col]]), res_ids)
    t_ix <- match(as.character(sup[[tm_col]]), as.character(time_ids))
    keep <- !is.na(r_ix) & !is.na(t_ix)
    sup <- sup[keep, , drop = FALSE]
    r_ix <- r_ix[keep]
    t_ix <- t_ix[keep]

    Y_res_all <- dlog_v_hat[is_res, , drop = FALSE]
    reg_df <- data.frame(
      origin_fe = res_ids[r_ix],
      time_fe = as.character(time_ids)[t_ix],
      y_dep = Y_res_all[cbind(r_ix, t_ix)],
      w_est = if (!is.null(wt_col)) as.numeric(sup[[wt_col]]) else 1.0,
      stringsAsFactors = FALSE
    )
    g_col_names <- paste0("Xg_", seq_len(K))
    for (k in seq_len(K)) {
      mask <- as.numeric(use_groups[[k]])
      sales_g <- sweep(dlog_y_hat, 1, y0 * mask, "*")
      X_full <- sweep(H_res %*% sales_g, 1, v0_res, "/")
      reg_df[[g_col_names[k]]] <- X_full[cbind(r_ix, t_ix)]
    }
    for (cc in ctrl_cols) reg_df[[cc]] <- as.numeric(sup[[cc]])

    rhs_all <- c(g_col_names, if (length(ctrl_cols) > 0L) paste0("`", ctrl_cols, "`") else character(0))
    fml_k <- as.formula(sprintf("y_dep ~ %s | origin_fe + time_fe", paste(rhs_all, collapse = " + ")))
    est_k <- fixest::feols(fml_k, data = reg_df, weights = ~ w_est,
                           vcov = if (compute_vcov) ~ origin_fe + time_fe else "iid", notes = FALSE)
    cf <- as.numeric(coef(est_k)[g_col_names])
    names(cf) <- names(use_groups)
    vc <- if (compute_vcov) {
      V_sub <- as.matrix(vcov(est_k))[g_col_names, g_col_names, drop = FALSE]
      dimnames(V_sub) <- list(names(use_groups), names(use_groups))
      V_sub
    } else NULL
    kv <- .map_group_kappa(cf, account, mode, clamp = FALSE)
    return(list(estimated = TRUE, mode = mode, coef = cf, vcov = vc, kappa_vec = kv))
  }

  # Fast analytical balanced two-way FE path
  Y_sub <- dlog_v_hat[is_res, kappa_periods, drop = FALSE]
  Yd <- .demean_twoway_mat(Y_sub)
  No_res <- nrow(Y_sub)
  Nt_sub <- ncol(Y_sub)

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

.map_group_kappa <- function(cf, account, mode, clamp = TRUE) {
  if (mode == "estimate_common") {
    kv <- rep(unname(cf[1]), account$N_dest)
  } else {
    kv <- as.numeric(cf[as.character(account$destinations$sector_group)])
  }
  if (clamp) pmin(0.999, pmax(0.01, kv)) else kv
}

.evaluate_anchored_welfare <- function(account, resp, kappa_vec, level_anchor, nu_vec, mu_vec, psi_O = 0.0) {
  is_res <- account$origins$is_city_resident
  v0 <- account$origins$v0
  v0_res <- v0[is_res]
  V_R0 <- sum(v0_res)
  y0 <- account$destinations$y0
  H_res <- account$H_m[is_res, , drop = FALSE]
  H_col_res <- colSums(H_res)

  mu_over_nu_res <- mu_vec[is_res] / nu_vec[is_res]
  Omega_res_weighted <- as.numeric(account$Omega_v[, is_res, drop = FALSE] %*% mu_over_nu_res)
  eff_retention_i <- H_col_res * kappa_vec + psi_O * (1.0 - kappa_vec)
  vartheta <- sum(y0 * eff_retention_i * Omega_res_weighted) / V_R0

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
    rhs_zero <- as.numeric(crossprod(y0 * eff_retention_i, resp$dlog_y_hat))
    psi_t <- (lhs_zero - rhs_zero) / V_R0
    c_t <- -psi_t / (1.0 - vartheta)
  }

  dlog_v <- resp$dlog_v_hat
  dlog_v[is_res, ] <- sweep(dlog_v[is_res, , drop = FALSE], 2, c_t, "+")

  dlog_y <- resp$dlog_y_hat + outer(Omega_res_weighted, c_t)
  dlog_p <- sweep(dlog_y, 1, kappa_vec, "*")
  living_cost <- account$S_tot %*% dlog_p
  welfare <- dlog_v - living_cost

  E_nonres_i <- y0 - colSums(account$E_R0[is_res, , drop = FALSE])
  sum_dEV <- colSums(sweep(welfare[is_res, , drop = FALSE], 1, v0_res, "*"))
  dm_mat <- sweep(dlog_y, 1, y0, "*") - sweep(dlog_p, 1, y0, "*")
  tot_rhs <- as.numeric(crossprod(E_nonres_i - (1.0 - H_col_res) * y0, dlog_p)) + psi_O * colSums(dm_mat)
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
