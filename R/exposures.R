#' Construct Direct and Network-Propagated Shock Exposures
#'
#' Computes direct destination shocks \eqn{d_{it} = \tau_i^0 h_{it}}, infinite-horizon
#' network sales exposures \eqn{\mathbf{F}_t^y = (I - M)^{-1}\mathbf{d}_t}, and
#' origin income exposures \eqn{\mathbf{F}_t^v = \Gamma_v \mathbf{F}_t^y}
#' (Proposition 2.3a). Supports both direct destination shock panels and
#' shift-share designs with optional Borusyak-Hull (2023) expected-exposure
#' recentering propagated through \eqn{(I - M)^{-1}}.
#'
#' @param account An \code{\link{urban_account}} object.
#' @param shock_panel For \code{design = "direct"}, either a data frame with
#'   columns \code{dest}, \code{time}, and \code{h} (proportional log shift in
#'   external demand) or \code{d} (sales-scaled shock \eqn{\mathrm{d}E_i^T / y_i^0}),
#'   or an \eqn{N_I \times T} numeric matrix of \code{h}. For
#'   \code{design = "shift_share"}, a named list with \code{shares}
#'   (\eqn{N_I \times G} matrix of destination group shares \eqn{\pi_{ig}}) and
#'   \code{shifters} (\eqn{G \times T} matrix of group log shifts \eqn{g_{gt}}),
#'   plus optional \code{shifter_expectations} (\eqn{G \times T} matrix of
#'   expected shifts for Borusyak-Hull network controls) and optional
#'   \code{realized_h} (\eqn{N_I \times T} realized log shifts for 2SLS IV).
#' @param design Character string: \code{"direct"} or \code{"shift_share"}.
#' @param center_dest Logical (default \code{TRUE}): whether to pre-enable
#'   origin-basket centering \eqn{\widetilde{x}_{nit} = x_{it} - \sum_j s_{nj}^{NH} x_{jt}}.
#' @param horizons Numeric vector of propagation horizons to compute alongside
#'   the exact \code{Inf} inverse (e.g. \code{c(0, 1, 2, Inf)}).
#'
#' @aliases urban_exposures
#' @return An S3 object of class \code{"urban_exposures"}.
#' @export
#' @examples
#' data(sim_city_160)
#' acct <- urban_account(sim_city_160$E_R0, sim_city_160$H_0,
#'                       sim_city_160$origins, sim_city_160$destinations)
#' exp_dir <- build_exposures(acct, sim_city_160$shocks, design = "direct")
#' summary(exp_dir)
build_exposures <- function(account,
                            shock_panel,
                            design = c("direct", "shift_share"),
                            center_dest = TRUE,
                            horizons = c(0, 1, 2, Inf)) {
  if (!inherits(account, "urban_account")) stop("account must be an 'urban_account' object.")
  design <- match.arg(design)
  Nd <- account$N_dest
  No <- account$N_orig
  tau0 <- account$destinations$tau0

  expected_controls <- NULL
  iv_instruments <- NULL
  ss_meta <- NULL

  if (design == "direct") {
    parsed <- .parse_direct_shocks(shock_panel, account$dest_ids, tau0)
    h_mat <- parsed$h_mat
    d_mat <- parsed$d_mat
    time_ids <- parsed$time_ids
  } else {
    if (!is.list(shock_panel) || is.null(shock_panel$shares) || is.null(shock_panel$shifters)) {
      stop("For design = 'shift_share', shock_panel must be a list with 'shares' (Nd x G) and 'shifters' (G x T).")
    }
    pi_ig <- as.matrix(shock_panel$shares)
    g_gt <- as.matrix(shock_panel$shifters)
    if (nrow(pi_ig) != Nd || ncol(pi_ig) != nrow(g_gt)) {
      stop("Conformability error between shock_panel$shares (Nd x G) and shock_panel$shifters (G x T).")
    }
    time_ids <- if (!is.null(colnames(g_gt))) colnames(g_gt) else as.character(seq_len(ncol(g_gt)))
    h_ss <- pi_ig %*% g_gt
    d_ss <- sweep(h_ss, 1, tau0, "*")
    ss_meta <- list(shares = pi_ig, shifters = g_gt)

    if (!is.null(shock_panel$shifter_expectations)) {
      mu_gt <- as.matrix(shock_panel$shifter_expectations)
      h_exp <- pi_ig %*% mu_gt
      d_exp <- sweep(h_exp, 1, tau0, "*")
      Fy_exp <- account$I_minus_M_inv %*% d_exp
      Fv_exp <- account$Gamma_v %*% Fy_exp
      expected_controls <- list(
        h_exp = h_exp,
        d_exp = d_exp,
        Fy_exp = Fy_exp,
        Fv_exp = Fv_exp,
        h_orig_exp = .map_dest_to_orig(h_exp, account)
      )
    }

    if (!is.null(shock_panel$realized_h)) {
      # Realized shocks are endogenous; shift-share exposures serve as 2SLS instruments
      h_mat <- as.matrix(shock_panel$realized_h)
      d_mat <- sweep(h_mat, 1, tau0, "*")
      Fy_iv <- account$I_minus_M_inv %*% d_ss
      Fv_iv <- account$Gamma_v %*% Fy_iv
      iv_instruments <- list(
        h_iv = h_ss,
        d_iv = d_ss,
        Fy_iv = Fy_iv,
        Fv_iv = Fv_iv,
        h_orig_iv = .map_dest_to_orig(h_ss, account)
      )
    } else {
      h_mat <- h_ss
      d_mat <- d_ss
    }
  }

  rownames(h_mat) <- account$dest_ids
  colnames(h_mat) <- time_ids
  rownames(d_mat) <- account$dest_ids
  colnames(d_mat) <- time_ids

  # Exact infinite-horizon propagation: Fy = (I - M)^(-1) d, Fv = Gamma_v Fy
  Fy <- account$I_minus_M_inv %*% d_mat
  Fv <- account$Gamma_v %*% Fy
  h_orig <- .map_dest_to_orig(h_mat, account)
  rownames(Fv) <- account$orig_ids
  colnames(Fv) <- time_ids

  # Optional finite-round approximations sum_{r=0}^H M^r d_t
  fin_horizons <- sort(unique(horizons[is.finite(horizons) & horizons >= 0]))
  horizon_list <- list()
  if (length(fin_horizons) > 0L) {
    max_H <- max(fin_horizons)
    term <- d_mat
    accum <- d_mat
    if (0 %in% fin_horizons) {
      horizon_list[["0"]] <- list(Fy = accum, Fv = account$Gamma_v %*% accum,
                                  rel_err_Fy = max(abs(accum - Fy)) / pmax(1e-12, max(abs(Fy))))
    }
    if (max_H >= 1L) {
      for (r in seq_len(max_H)) {
        term <- account$M %*% term
        accum <- accum + term
        if (r %in% fin_horizons) {
          horizon_list[[as.character(r)]] <- list(
            Fy = accum,
            Fv = account$Gamma_v %*% accum,
            rel_err_Fy = max(abs(accum - Fy)) / pmax(1e-12, max(abs(Fy)))
          )
        }
      }
    }
  }

  structure(
    list(
      design = design,
      center_dest = center_dest,
      time_ids = time_ids,
      N_time = length(time_ids),
      orig_ids = account$orig_ids,
      dest_ids = account$dest_ids,
      S_NH = account$S_NH,
      h = h_mat,
      d = d_mat,
      h_orig = h_orig,
      Fy = Fy,
      Fv = Fv,
      horizons = horizon_list,
      expected_controls = expected_controls,
      iv_instruments = iv_instruments,
      shift_share = ss_meta
    ),
    class = "urban_exposures"
  )
}

.parse_direct_shocks <- function(shock_panel, dest_ids, tau0) {
  Nd <- length(dest_ids)
  if (is.matrix(shock_panel)) {
    if (nrow(shock_panel) != Nd) stop("Matrix shock_panel must have N_dest rows.")
    h_mat <- shock_panel
    d_mat <- sweep(h_mat, 1, tau0, "*")
    t_ids <- if (!is.null(colnames(shock_panel))) colnames(shock_panel) else as.character(seq_len(ncol(shock_panel)))
    return(list(h_mat = h_mat, d_mat = d_mat, time_ids = t_ids))
  }
  df <- as.data.frame(shock_panel)
  if (!all(c("dest", "time") %in% names(df))) {
    stop("Data frame shock_panel must contain 'dest' and 'time' columns plus 'h' or 'd'.")
  }
  t_ids <- as.character(sort(unique(df$time)))
  Nt <- length(t_ids)
  d_idx <- match(as.character(df$dest), dest_ids)
  t_idx <- match(as.character(df$time), t_ids)
  if (anyNA(d_idx)) stop("Unmatched destination IDs in shock_panel.")

  h_mat <- matrix(0.0, nrow = Nd, ncol = Nt)
  d_mat <- matrix(0.0, nrow = Nd, ncol = Nt)
  if ("h" %in% names(df)) {
    h_mat[cbind(d_idx, t_idx)] <- as.numeric(df$h)
    d_mat <- if ("d" %in% names(df)) {
      dm <- matrix(0.0, nrow = Nd, ncol = Nt)
      dm[cbind(d_idx, t_idx)] <- as.numeric(df$d)
      dm
    } else {
      sweep(h_mat, 1, tau0, "*")
    }
  } else if ("d" %in% names(df)) {
    d_mat[cbind(d_idx, t_idx)] <- as.numeric(df$d)
    h_mat <- sweep(d_mat, 1, ifelse(tau0 > 0, tau0, 1.0), "/")
    h_mat[tau0 <= 0, ] <- 0.0
  } else {
    stop("shock_panel must contain column 'h' or 'd'.")
  }
  list(h_mat = h_mat, d_mat = d_mat, time_ids = t_ids)
}

.map_dest_to_orig <- function(x_dest_mat, account) {
  out <- matrix(0.0, nrow = account$N_orig, ncol = ncol(x_dest_mat),
                dimnames = list(account$orig_ids, colnames(x_dest_mat)))
  m <- match(account$orig_ids, account$dest_ids)
  ok <- !is.na(m)
  if (any(ok)) out[ok, ] <- x_dest_mat[m[ok], , drop = FALSE]
  out
}

#' @method print urban_exposures
#' @export
print.urban_exposures <- function(x, ...) {
  cat(sprintf("<urban_exposures: design = '%s', %d origins, %d destinations, %d periods>\n",
              x$design, length(x$orig_ids), length(x$dest_ids), x$N_time))
  if (!is.null(x$expected_controls)) {
    cat("  Borusyak-Hull expected network controls: INCLUDED\n")
  }
  if (!is.null(x$iv_instruments)) {
    cat("  Shift-share 2SLS instruments: INCLUDED\n")
  }
  invisible(x)
}

#' @method summary urban_exposures
#' @export
summary.urban_exposures <- function(object, ...) {
  print.urban_exposures(object)
  cat(sprintf("  SD(d_it) = %.4f | SD(F^y_it) = %.4f | SD(F^v_nt) = %.4f\n",
              sd(as.vector(object$d)), sd(as.vector(object$Fy)), sd(as.vector(object$Fv))))
  if (length(object$horizons) > 0L) {
    cat("  Finite-round approximation max relative error vs. (I - M)^(-1):\n")
    for (nm in names(object$horizons)) {
      cat(sprintf("    Horizon h = %-2s : %.4e\n", nm, object$horizons[[nm]]$rel_err_Fy))
    }
  }
  invisible(object)
}
