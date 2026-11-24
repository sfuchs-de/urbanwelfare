#' Construct and Validate a Two-Sector Open-City Spatial Account
#'
#' Builds a validated \code{"urban_account"} object from initial bilateral
#' resident spending, commuting/ownership links, origin accounts, and
#' destination accounts. Verifies non-housing budget shares, destination sales
#' clearing, open-city factor leakage, and the spectral radius condition
#' \eqn{\rho(M) < 1} ensuring convergence of the spatial network multiplier
#' \eqn{(I - M)^{-1}}.
#'
#' @param spending_base Numeric matrix (\eqn{N_R \times N_I}) or 3-column data
#'   frame (\code{origin}, \code{dest}, \code{spending}) of initial bilateral
#'   resident non-housing expenditure \eqn{E_{ni}^{R,0}}.
#' @param commuting_base Numeric matrix (\eqn{N_R \times N_I}) or 3-column data
#'   frame (\code{origin}, \code{dest}, \code{flow}) of workplace payroll flows
#'   \eqn{W_{ni}^0} or direct factor-allocation shares \eqn{H_{ni}^m}.
#' @param origins Data frame with \eqn{N_R} rows containing at least \code{id}
#'   (optional, defaults to \code{1:nrow(origins)}), \code{v0} (initial
#'   disposable income), and optionally \code{pop} (population, default 1),
#'   \code{s_H} (housing/saving share \eqn{1 - \nu_n}, default 0.32),
#'   \code{s_out_NH} (conditional outside-good share \eqn{s_{n0}^{NH}}), and
#'   \code{is_city_resident} (logical flag separating core city residents from
#'   surrounding metropolitan origins, default \code{TRUE}).
#' @param destinations Data frame with \eqn{N_I} rows containing at least
#'   \code{id} (optional), \code{y0} (initial destination sales), \code{E_ext0}
#'   (initial external/visitor demand), and optionally \code{sector_group}
#'   (activity group for \eqn{\kappa_g^S}, default 1), \code{theta_G} (initial
#'   local factor share \eqn{G_i^0 / y_i^0}, default 0.45), \code{theta_L}
#'   (labor share of local factor income, default 0.65), and
#'   \code{resident_factor_share} (share of destination local factor income
#'   accruing to modeled origins \eqn{\sum_n H_{ni}^m \le 1}).
#' @param balance_method Character string: \code{"reconcile"} (default, computes
#'   implied outside shares and fixed demand residuals cleanly), \code{"ipf"}
#'   (iterative proportional fitting of \code{spending_base} when target origin
#'   and destination resident totals are supplied), or \code{"strict"} (errors
#'   if any accounting residual is negative beyond \code{tol}).
#' @param commuting_is_share Logical or \code{"auto"}: whether
#'   \code{commuting_base} already gives factor shares \eqn{H_{ni}^m} with
#'   column sums \eqn{\le 1} rather than nominal payroll levels.
#' @param tol Numeric tolerance for accounting inequalities (default \code{1e-7}).
#'
#' @return An S3 object of class \code{"urban_account"}.
#' @export
#' @examples
#' data(sim_city_160)
#' acct <- urban_account(
#'   spending_base  = sim_city_160$E_R0,
#'   commuting_base = sim_city_160$H_0,
#'   origins        = sim_city_160$origins,
#'   destinations   = sim_city_160$destinations
#' )
#' summary(acct)
urban_account <- function(spending_base,
                          commuting_base,
                          origins,
                          destinations,
                          balance_method = c("reconcile", "ipf", "strict"),
                          commuting_is_share = "auto",
                          tol = 1e-7) {
  balance_method <- match.arg(balance_method)
  origins <- as.data.frame(origins)
  destinations <- as.data.frame(destinations)

  No <- nrow(origins)
  Nd <- nrow(destinations)
  if (No < 1L || Nd < 1L) stop("origins and destinations must be non-empty.")

  orig_ids <- if ("id" %in% names(origins)) as.character(origins$id) else as.character(seq_len(No))
  dest_ids <- if ("id" %in% names(destinations)) as.character(destinations$id) else as.character(seq_len(Nd))
  origins$id <- orig_ids
  destinations$id <- dest_ids

  E_R0 <- .coerce_bip_matrix(spending_base, orig_ids, dest_ids, "spending_base")
  W_in <- .coerce_bip_matrix(commuting_base, orig_ids, dest_ids, "commuting_base")

  if (!("v0" %in% names(origins)) || any(!is.finite(origins$v0)) || any(origins$v0 <= 0)) {
    stop("origins$v0 must be positive and finite.")
  }
  if (!("y0" %in% names(destinations)) || any(!is.finite(destinations$y0)) || any(destinations$y0 <= 0)) {
    stop("destinations$y0 must be positive and finite.")
  }
  if (!("E_ext0" %in% names(destinations))) destinations$E_ext0 <- 0
  if (any(!is.finite(destinations$E_ext0)) || any(destinations$E_ext0 < -tol)) {
    stop("destinations$E_ext0 must be nonnegative and finite.")
  }

  pop <- if ("pop" %in% names(origins)) as.numeric(origins$pop) else rep(1.0, No)
  s_H <- if ("s_H" %in% names(origins)) as.numeric(origins$s_H) else rep(0.32, No)
  if (any(s_H < 0 | s_H >= 1)) stop("origins$s_H must lie in [0, 1).")
  nu <- 1.0 - s_H
  is_res <- if ("is_city_resident" %in% names(origins)) as.logical(origins$is_city_resident) else rep(TRUE, No)
  if (!any(is_res)) stop("At least one origin must have is_city_resident = TRUE.")

  v0 <- as.numeric(origins$v0)
  y0 <- as.numeric(destinations$y0)
  E_ext0 <- pmax(0.0, as.numeric(destinations$E_ext0))
  sector_group <- if ("sector_group" %in% names(destinations)) destinations$sector_group else rep(1L, Nd)
  theta_G <- if ("theta_G" %in% names(destinations)) as.numeric(destinations$theta_G) else rep(0.45, Nd)
  theta_L <- if ("theta_L" %in% names(destinations)) as.numeric(destinations$theta_L) else rep(0.65, Nd)
  if (any(theta_G <= 0 | theta_G >= 1)) stop("destinations$theta_G must lie in (0, 1).")
  G0 <- if ("G0" %in% names(destinations)) as.numeric(destinations$G0) else theta_G * y0
  theta_G <- G0 / y0

  # Optional IPF balancing of noisy resident spending matrix
  if (balance_method == "ipf") {
    target_row <- if ("s_out_NH" %in% names(origins)) {
      nu * v0 * (1.0 - as.numeric(origins$s_out_NH))
    } else {
      pmin(rowSums(E_R0), nu * v0 * 0.98)
    }
    A_fix_in <- if ("A_fixed0" %in% names(destinations)) as.numeric(destinations$A_fixed0) else pmax(0.0, y0 - colSums(E_R0) - E_ext0)
    target_col <- pmax(1e-10, y0 - E_ext0 - A_fix_in)
    # Scale target_row and target_col to common grand total
    grand <- 0.5 * (sum(target_row) + sum(target_col))
    target_row <- target_row * (grand / sum(target_row))
    target_col <- target_col * (grand / sum(target_col))
    E_R0 <- .ipf_balance(E_R0, target_row, target_col)
  }

  # Validate non-housing spending shares
  nh_budget <- nu * v0
  local_nh_spend <- rowSums(E_R0)
  if (any(local_nh_spend > nh_budget * (1 + tol))) {
    if (balance_method == "strict") {
      stop("Row sums of spending_base exceed non-housing budget nu * v0 under strict mode.")
    }
    scale_down <- pmin(1.0, (nh_budget * 0.995) / pmax(1e-12, local_nh_spend))
    E_R0 <- E_R0 * scale_down
    local_nh_spend <- rowSums(E_R0)
  }
  S_NH <- E_R0 / nh_budget
  s_out_NH <- pmax(0.0, 1.0 - rowSums(S_NH))
  S_tot <- E_R0 / v0

  # Destination sales clearing: y0 = sum_n E_R0[n, i] + E_ext0[i] + A_fixed0[i]
  res_sales <- colSums(E_R0)
  A_fixed0 <- y0 - res_sales - E_ext0
  if (any(A_fixed0 < -tol * y0)) {
    if (balance_method == "strict") {
      stop("Resident + external demand exceeds destination y0 under strict mode.")
    }
    A_fixed0 <- pmax(0.0, A_fixed0)
    y0 <- res_sales + E_ext0 + A_fixed0
    G0 <- theta_G * y0
  } else {
    A_fixed0 <- pmax(0.0, A_fixed0)
  }

  # Factor allocation matrix H_m (No x Nd), column sums <= 1
  col_w <- colSums(W_in)
  is_share <- if (identical(commuting_is_share, "auto")) {
    max(col_w) <= 1.0 + 1e-5 && !("resident_factor_share" %in% names(destinations))
  } else {
    isTRUE(commuting_is_share)
  }
  if (is_share) {
    H_m <- W_in
    if (any(colSums(H_m) > 1.0 + tol)) {
      stop("commuting_base treated as factor shares H_m, but column sums exceed 1.")
    }
  } else {
    res_f_share <- if ("resident_factor_share" %in% names(destinations)) {
      as.numeric(destinations$resident_factor_share)
    } else {
      pmin(1.0, col_w / pmax(1e-12, G0))
    }
    H_m <- sweep(W_in, 2, pmax(1e-12, col_w), "/")
    H_m <- sweep(H_m, 2, res_f_share, "*")
  }

  v_local0 <- as.numeric(H_m %*% G0)
  v_out <- v0 - v_local0
  if (any(v_out < -tol * v0) && balance_method == "strict") {
    stop("Implied local factor income H_m %*% G0 exceeds origins$v0 under strict mode.")
  }
  v_out <- pmax(0.0, v_out)

  # Build spatial network operators: Gamma_v (No x Nd), Omega_v (Nd x No), M = Omega_v %*% Gamma_v (Nd x Nd)
  Gamma_v <- sweep(H_m, 2, G0, "*")
  Gamma_v <- sweep(Gamma_v, 1, v0, "/")
  Omega_v <- sweep(t(E_R0), 1, y0, "/")
  M_op <- Omega_v %*% Gamma_v

  ev <- eigen(M_op, only.values = TRUE)$values
  rho_M <- max(Mod(ev))
  if (!is.finite(rho_M) || rho_M >= 1.0) {
    stop(sprintf("Spectral radius rho(M) = %.4f >= 1; network series (I - M)^(-1) does not converge.", rho_M))
  }

  # Leakage diagnostics (weighted over core city residents / all destinations)
  V_res0 <- sum(v0[is_res])
  Y_tot0 <- sum(y0)
  G_tot0 <- sum(G0)
  H_res_col <- colSums(H_m[is_res, , drop = FALSE])
  leakage <- list(
    flexible_inputs_share_sales = sum(y0 - G0) / Y_tot0,
    housing_saving_share_income = sum(s_H[is_res] * v0[is_res]) / V_res0,
    outside_goods_share_income  = sum((nu * s_out_NH * v0)[is_res]) / V_res0,
    external_factors_share_G    = sum((1.0 - H_res_col) * G0) / G_tot0,
    external_demand_share_sales = sum(E_ext0 + A_fixed0 + colSums(E_R0[!is_res, , drop = FALSE])) / Y_tot0,
    visitor_demand_share_sales  = sum(E_ext0) / Y_tot0
  )

  origins_clean <- origins
  origins_clean$pop <- pop
  origins_clean$v0 <- v0
  origins_clean$v_out <- v_out
  origins_clean$s_H <- s_H
  origins_clean$nu <- nu
  origins_clean$s_out_NH <- s_out_NH
  origins_clean$is_city_resident <- is_res

  destinations_clean <- destinations
  destinations_clean$y0 <- y0
  destinations_clean$G0 <- G0
  destinations_clean$E_ext0 <- E_ext0
  destinations_clean$A_fixed0 <- A_fixed0
  destinations_clean$tau0 <- E_ext0 / y0
  destinations_clean$theta_G <- theta_G
  destinations_clean$theta_L <- theta_L
  destinations_clean$sector_group <- sector_group
  destinations_clean$H_res <- H_res_col

  structure(
    list(
      N_orig = No,
      N_dest = Nd,
      orig_ids = orig_ids,
      dest_ids = dest_ids,
      origins = origins_clean,
      destinations = destinations_clean,
      E_R0 = E_R0,
      H_m = H_m,
      S_NH = S_NH,
      S_tot = S_tot,
      s_out_NH = s_out_NH,
      Gamma_v = Gamma_v,
      Omega_v = Omega_v,
      M = M_op,
      I_minus_M_inv = solve(diag(Nd) - M_op),
      rho_M = rho_M,
      multiplier = 1.0 / (1.0 - rho_M),
      leakage = leakage,
      group_levels = sort(unique(sector_group)),
      balance_method = balance_method
    ),
    class = "urban_account"
  )
}

.coerce_bip_matrix <- function(x, row_ids, col_ids, name) {
  if (is.matrix(x) || inherits(x, "Matrix")) {
    m <- as.matrix(x)
    if (nrow(m) != length(row_ids) || ncol(m) != length(col_ids)) {
      stop(sprintf("%s dimensions (%d x %d) do not match origins (%d) x destinations (%d).",
                   name, nrow(m), ncol(m), length(row_ids), length(col_ids)))
    }
    rownames(m) <- row_ids
    colnames(m) <- col_ids
    return(m)
  }
  if (is.data.frame(x) && ncol(x) >= 3L) {
    m <- matrix(0.0, nrow = length(row_ids), ncol = length(col_ids),
                dimnames = list(row_ids, col_ids))
    r_idx <- match(as.character(x[[1]]), row_ids)
    c_idx <- match(as.character(x[[2]]), col_ids)
    if (anyNA(r_idx) || anyNA(c_idx)) stop(sprintf("Unmatched origin/dest IDs in %s.", name))
    vals <- as.numeric(x[[3]])
    for (k in seq_along(vals)) {
      m[r_idx[k], c_idx[k]] <- m[r_idx[k], c_idx[k]] + vals[k]
    }
    return(m)
  }
  stop(sprintf("%s must be a matrix or a 3-column data.frame (origin, dest, value).", name))
}

.ipf_balance <- function(mat, row_targets, col_targets, max_iter = 500L, tol = 1e-10) {
  m <- as.matrix(mat)
  m[m < 0] <- 0.0
  for (it in seq_len(max_iter)) {
    rs <- rowSums(m)
    m <- sweep(m, 1, ifelse(rs > 0, row_targets / rs, 1.0), "*")
    cs <- colSums(m)
    m <- sweep(m, 2, ifelse(cs > 0, col_targets / cs, 1.0), "*")
    err <- max(abs(rowSums(m) - row_targets), abs(colSums(m) - col_targets))
    if (err < tol) break
  }
  m
}

#' @method print urban_account
#' @export
print.urban_account <- function(x, ...) {
  cat(sprintf("<urban_account: %d origins (%d city residents), %d destinations, %d sector groups>\n",
              x$N_orig, sum(x$origins$is_city_resident), x$N_dest, length(x$group_levels)))
  cat(sprintf("  Spectral radius rho(M) = %.4f  |  Network multiplier (1 - rho)^(-1) = %.3f\n",
              x$rho_M, x$multiplier))
  invisible(x)
}

#' @method summary urban_account
#' @export
summary.urban_account <- function(object, ...) {
  cat("=== Two-Sector Open-City Spatial Account Summary ===\n")
  print.urban_account(object)
  cat("\nLeakage & Openness Shares:\n")
  cat(sprintf("  Flexible traded inputs (share of sales)     : %6.2f%%\n", 100 * object$leakage$flexible_inputs_share_sales))
  cat(sprintf("  Housing & saving (share of resident income) : %6.2f%%\n", 100 * object$leakage$housing_saving_share_income))
  cat(sprintf("  Outside traded goods (share of res. income) : %6.2f%%\n", 100 * object$leakage$outside_goods_share_income))
  cat(sprintf("  External factor payments (share of G^0)     : %6.2f%%\n", 100 * object$leakage$external_factors_share_G))
  cat(sprintf("  Visitor/shocked demand (share of sales)     : %6.2f%%\n", 100 * object$leakage$visitor_demand_share_sales))
  cat(sprintf("  Total non-resident sales (share of sales)   : %6.2f%%\n", 100 * object$leakage$external_demand_share_sales))
  invisible(object$leakage)
}
