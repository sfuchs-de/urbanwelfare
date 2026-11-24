#' Estimate the Bilateral Spending Regression (Step 1)
#'
#' Estimates Proposition 2.4's bilateral spending specification wrapping
#' \code{fixest::feols} (OLS or 2SLS IV) or \code{fixest::fepois} (PPML) with
#' origin-destination pair and time fixed effects and two-way origin-time and
#' destination-time clustered standard errors. Uses C++-accelerated basket
#' centering for destination exposures \eqn{\widetilde{h}_{nit}} and
#' \eqn{\widetilde{F}^y_{nit}}, and automatically includes Borusyak-Hull
#' expected network controls when present in \code{exposures}.
#'
#' @param panel Data frame with columns \code{origin}, \code{dest}, \code{time},
#'   and either \code{log_spending} or positive \code{spending}.
#' @param exposures An \code{\link{urban_exposures}} object.
#' @param estimator Character string: \code{"ols"} (default), \code{"ppml"}, or
#'   \code{"iv"} (2SLS when \code{exposures} contains shift-share instruments).
#' @param cluster Character vector of length 2 specifying multi-way cluster
#'   dimensions (default \code{c("origin^time", "dest^time")}).
#' @param include_direct Logical (default \code{TRUE}): whether to include the
#'   direct origin and destination terms \eqn{(h_{nt}, \widetilde{h}_{nit})}
#'   alongside the network exposures \eqn{(F^v_{nt}, \widetilde{F}^y_{nit})}.
#' @param use_expected_controls Logical (default \code{TRUE}): include
#'   Borusyak-Hull expected exposures when present in \code{exposures}.
#' @param extra_controls Optional character vector of additional covariate column
#'   names already present in \code{panel} to include on the right-hand side.
#' @param weights Optional character string naming a non-negative observation
#'   weight column in \code{panel}.
#'
#' @aliases urban_bilateral
#' @return An S3 object of class \code{"urban_bilateral"} containing the
#'   fitted \code{fixest} model, structural slope vector
#'   \eqn{(a_O, b_v, a_D, b_D)}, clustered covariance matrix, and metadata.
#' @export
#' @examples
#' data(sim_city_160)
#' acct <- urban_account(sim_city_160$E_R0, sim_city_160$H_0,
#'                       sim_city_160$origins, sim_city_160$destinations)
#' exps <- build_exposures(acct, sim_city_160$shocks)
#' fit <- estimate_bilateral(sim_city_160$bilateral_panel, exps, estimator = "ols")
#' summary(fit)
estimate_bilateral <- function(panel,
                               exposures,
                               estimator = c("ols", "ppml", "iv"),
                               cluster = c("origin^time", "dest^time"),
                               include_direct = TRUE,
                               use_expected_controls = TRUE,
                               extra_controls = character(0),
                               weights = NULL) {
  if (!inherits(exposures, "urban_exposures")) stop("exposures must be an 'urban_exposures' object.")
  estimator <- match.arg(estimator)
  df <- as.data.frame(panel)
  if (!all(c("origin", "dest", "time") %in% names(df))) {
    stop("panel must contain columns 'origin', 'dest', and 'time'.")
  }
  if (length(extra_controls) > 0L && !all(extra_controls %in% names(df))) {
    stop("All extra_controls columns must exist in panel.")
  }

  o_idx <- match(as.character(df$origin), exposures$orig_ids)
  d_idx <- match(as.character(df$dest), exposures$dest_ids)
  t_idx <- match(as.character(df$time), exposures$time_ids)
  if (anyNA(o_idx) || anyNA(d_idx) || anyNA(t_idx)) {
    stop("Unmatched origin, dest, or time IDs between panel and exposures.")
  }

  # Construct regressors using C++ basket-centering
  df$h_orig <- exposures$h_orig[cbind(o_idx, t_idx)]
  df$Fv_orig <- exposures$Fv[cbind(o_idx, t_idx)]
  if (isTRUE(exposures$center_dest)) {
    df$h_dest_tilde <- cpp_basket_center(exposures$h, exposures$S_NH, o_idx, d_idx, t_idx)
    df$Fy_dest_tilde <- cpp_basket_center(exposures$Fy, exposures$S_NH, o_idx, d_idx, t_idx)
  } else {
    df$h_dest_tilde <- exposures$h[cbind(d_idx, t_idx)]
    df$Fy_dest_tilde <- exposures$Fy[cbind(d_idx, t_idx)]
  }

  ctrl_terms <- character(0)
  if (isTRUE(use_expected_controls) && !is.null(exposures$expected_controls)) {
    ec <- exposures$expected_controls
    df$Fv_orig_exp <- ec$Fv_exp[cbind(o_idx, t_idx)]
    df$Fy_dest_tilde_exp <- cpp_basket_center(ec$Fy_exp, exposures$S_NH, o_idx, d_idx, t_idx)
    ctrl_terms <- c("Fv_orig_exp", "Fy_dest_tilde_exp")
    if (include_direct) {
      df$h_orig_exp <- ec$h_orig_exp[cbind(o_idx, t_idx)]
      df$h_dest_tilde_exp <- cpp_basket_center(ec$h_exp, exposures$S_NH, o_idx, d_idx, t_idx)
      ctrl_terms <- c("h_orig_exp", "Fv_orig_exp", "h_dest_tilde_exp", "Fy_dest_tilde_exp")
    }
  }
  if (length(extra_controls) > 0L) {
    ctrl_terms <- c(ctrl_terms, paste0("`", extra_controls, "`"))
  }

  main_terms <- if (include_direct) {
    c("h_orig", "Fv_orig", "h_dest_tilde", "Fy_dest_tilde")
  } else {
    c("Fv_orig", "Fy_dest_tilde")
  }
  rhs_ctrl <- if (length(ctrl_terms) > 0L) paste("+", paste(ctrl_terms, collapse = " + ")) else ""
  vcov_fml <- as.formula(paste("~", paste(cluster, collapse = " + ")))
  w_fml <- if (!is.null(weights)) as.formula(paste0("~ `", weights, "`")) else NULL

  if (estimator == "ppml") {
    if (!("spending" %in% names(df))) {
      if ("log_spending" %in% names(df)) df$spending <- exp(df$log_spending)
      else stop("PPML requires column 'spending' (or 'log_spending') in panel.")
    }
    fml <- as.formula(sprintf("spending ~ %s %s | origin^dest + time",
                              paste(main_terms, collapse = " + "), rhs_ctrl))
    est <- fixest::fepois(fml, data = df, vcov = vcov_fml, weights = w_fml, notes = FALSE)
  } else if (estimator == "iv") {
    if (is.null(exposures$iv_instruments)) {
      stop("estimator = 'iv' requires shift-share instruments in exposures (pass realized_h in build_exposures).")
    }
    if (!("log_spending" %in% names(df))) {
      if ("spending" %in% names(df)) df$log_spending <- log(df$spending)
      else stop("OLS/IV requires column 'log_spending' or positive 'spending' in panel.")
    }
    iv <- exposures$iv_instruments
    df$h_orig_iv <- iv$h_orig_iv[cbind(o_idx, t_idx)]
    df$Fv_orig_iv <- iv$Fv_iv[cbind(o_idx, t_idx)]
    df$h_dest_tilde_iv <- cpp_basket_center(iv$h_iv, exposures$S_NH, o_idx, d_idx, t_idx)
    df$Fy_dest_tilde_iv <- cpp_basket_center(iv$Fy_iv, exposures$S_NH, o_idx, d_idx, t_idx)
    iv_terms <- paste0(main_terms, "_iv")
    rhs_1 <- if (length(ctrl_terms) > 0L) paste(ctrl_terms, collapse = " + ") else "1"
    fml <- as.formula(sprintf("log_spending ~ %s | origin^dest + time | %s ~ %s",
                              rhs_1,
                              paste(main_terms, collapse = " + "),
                              paste(iv_terms, collapse = " + ")))
    est <- fixest::feols(fml, data = df, vcov = vcov_fml, weights = w_fml, notes = FALSE)
  } else {
    if (!("log_spending" %in% names(df))) {
      if ("spending" %in% names(df)) df$log_spending <- log(df$spending)
      else stop("OLS requires column 'log_spending' or positive 'spending' in panel.")
    }
    fml <- as.formula(sprintf("log_spending ~ %s %s | origin^dest + time",
                              paste(main_terms, collapse = " + "), rhs_ctrl))
    est <- fixest::feols(fml, data = df, vcov = vcov_fml, weights = w_fml, notes = FALSE)
  }

  cf <- coef(est)
  vc <- as.matrix(vcov(est))
  canon_names <- c("h_orig", "Fv_orig", "h_dest_tilde", "Fy_dest_tilde")
  fit_names <- if (estimator == "iv") paste0("fit_", canon_names) else canon_names
  b4 <- setNames(rep(0.0, 4L), c("a_O", "b_v", "a_D", "b_D"))
  V4 <- matrix(0.0, 4L, 4L, dimnames = list(names(b4), names(b4)))
  present <- which(fit_names %in% names(cf))
  b4[present] <- cf[fit_names[present]]
  V4[present, present] <- vc[fit_names[present], fit_names[present], drop = FALSE]

  structure(
    list(
      fixest_model = est,
      estimator = estimator,
      include_direct = include_direct,
      coefficients = b4,
      vcov = V4,
      n_obs = stats::nobs(est),
      exposures = exposures
    ),
    class = "urban_bilateral"
  )
}

#' Construct an \code{urban_bilateral} Object from Pre-Estimated Slopes
#'
#' Enables secure-enclave workflows where Step 1's multi-million-row bilateral
#' panel regression is estimated inside a restricted enclave and only the four
#' estimated slopes \eqn{(a_O, b_v, a_D, b_D)} and their \eqn{4 \times 4}
#' covariance matrix are exported for downstream welfare and policy analysis.
#'
#' @param coefficients Numeric vector of length 4 giving \eqn{(a_O, b_v, a_D, b_D)}.
#' @param exposures An \code{\link{urban_exposures}} object.
#' @param vcov Optional \eqn{4 \times 4} covariance matrix (defaults to zeros).
#' @param n_obs Optional integer observation count for reporting.
#' @return An S3 object of class \code{"urban_bilateral"}.
#' @export
bilateral_from_estimates <- function(coefficients, exposures, vcov = NULL, n_obs = NA_integer_) {
  if (!inherits(exposures, "urban_exposures")) stop("exposures must be an 'urban_exposures' object.")
  if (length(coefficients) != 4L) stop("coefficients must have length 4: c(a_O, b_v, a_D, b_D).")
  b4 <- setNames(as.numeric(coefficients), c("a_O", "b_v", "a_D", "b_D"))
  if (is.null(vcov)) {
    V4 <- matrix(0.0, 4L, 4L, dimnames = list(names(b4), names(b4)))
  } else {
    V4 <- as.matrix(vcov)
    if (!all(dim(V4) == c(4L, 4L))) stop("vcov must be 4 x 4.")
    dimnames(V4) <- list(names(b4), names(b4))
  }
  structure(
    list(
      fixest_model = NULL,
      estimator = "imported",
      include_direct = TRUE,
      coefficients = b4,
      vcov = V4,
      n_obs = as.integer(n_obs),
      exposures = exposures
    ),
    class = "urban_bilateral"
  )
}

#' Compute Rotemberg Sensitivity Weights for Shift-Share Network Designs
#'
#' Computes Aggregated Rotemberg (1983) / Goldsmith-Pinkham, Sorkin, and Swift
#' (2020) instrument relevance weights across shock origin groups \eqn{g = 1,\dots,G}
#' for the network income exposure \eqn{F^v_{nt}}.
#'
#' @param exposures An \code{\link{urban_exposures}} object created with
#'   \code{design = "shift_share"}.
#' @param account An \code{\link{urban_account}} object.
#' @return A data frame of shock groups sorted by absolute Rotemberg weight,
#'   along with effective shock count \eqn{1 / \sum_g \alpha_g^2}.
#' @export
rotemberg_weights <- function(exposures, account) {
  if (exposures$design != "shift_share" || is.null(exposures$shift_share)) {
    stop("rotemberg_weights requires an urban_exposures object with design = 'shift_share'.")
  }
  pi_ig <- exposures$shift_share$shares
  g_gt <- exposures$shift_share$shifters
  tau0 <- account$destinations$tau0
  G <- ncol(pi_ig)

  # Decompose Fv = sum_g Fv^{(g)}, where Fv^{(g)} = Gamma_v (I - M)^(-1) diag(tau0) pi_{.g} g_{g.}
  Op <- account$Gamma_v %*% account$I_minus_M_inv
  Fv_dm <- .demean_twoway_mat(exposures$Fv)
  denom <- sum(Fv_dm * Fv_dm)
  alpha <- numeric(G)
  for (g in seq_len(G)) {
    s_g <- Op %*% (tau0 * pi_ig[, g])
    Fv_g <- s_g %*% matrix(g_gt[g, ], nrow = 1L)
    Fv_g_dm <- .demean_twoway_mat(Fv_g)
    alpha[g] <- sum(Fv_g_dm * Fv_dm) / denom
  }
  grp_names <- if (!is.null(colnames(pi_ig))) colnames(pi_ig) else paste0("group_", seq_len(G))
  df <- data.frame(
    group = grp_names,
    alpha = alpha,
    abs_alpha_share = abs(alpha) / sum(abs(alpha)),
    stringsAsFactors = FALSE
  )
  df <- df[order(-abs(df$alpha)), ]
  rownames(df) <- NULL
  attr(df, "effective_shocks_hhi_inv") <- 1.0 / sum(alpha^2)
  df
}

#' Leave-One-Location-Out Jackknife Influence Diagnostics
#'
#' Re-estimates \code{\link{estimate_bilateral}} dropping each sampled origin
#' (or destination) in turn to check that no single hub drives the network
#' slope estimates \eqn{(b_v, b_D)}.
#'
#' @param panel Bilateral spending panel passed to \code{\link{estimate_bilateral}}.
#' @param exposures An \code{\link{urban_exposures}} object.
#' @param by Character: \code{"dest"} (default) or \code{"origin"}.
#' @param max_units Maximum number of highest-exposure units to test (default 15).
#' @param ... Additional arguments passed to \code{\link{estimate_bilateral}}.
#' @return A data frame of leave-one-out slope estimates.
#' @export
loo_influence <- function(panel, exposures, by = c("dest", "origin"), max_units = 15L, ...) {
  by <- match.arg(by)
  scores <- if (by == "dest") rowMeans(abs(exposures$d)) else rowMeans(abs(exposures$Fv))
  ids <- if (by == "dest") exposures$dest_ids else exposures$orig_ids
  ord <- head(order(-scores), min(length(ids), as.integer(max_units)))
  target_ids <- ids[ord]

  rows <- lapply(target_ids, function(uid) {
    sub_df <- panel[as.character(panel[[by]]) != uid, , drop = FALSE]
    fit_sub <- estimate_bilateral(sub_df, exposures, ...)
    cf <- fit_sub$coefficients
    data.frame(dropped_unit = uid, a_O = cf["a_O"], b_v = cf["b_v"],
               a_D = cf["a_D"], b_D = cf["b_D"], stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

.demean_twoway_mat <- function(X) {
  rm <- rowMeans(X)
  cm <- colMeans(X)
  gm <- mean(X)
  sweep(sweep(X, 1, rm, "-"), 2, cm, "-") + gm
}

#' @method print urban_bilateral
#' @export
print.urban_bilateral <- function(x, ...) {
  cat(sprintf("<urban_bilateral: estimator = '%s', N = %s>\n",
              x$estimator, format(x$n_obs, big.mark = ",")))
  se <- sqrt(pmax(0, diag(x$vcov)))
  tab <- cbind(Estimate = x$coefficients, `Std. Error` = se,
               `t value` = ifelse(se > 0, x$coefficients / se, NA_real_))
  print(round(tab, 4))
  invisible(x)
}

#' @method summary urban_bilateral
#' @export
summary.urban_bilateral <- function(object, ...) {
  print.urban_bilateral(object)
  invisible(object)
}
