#' Build an \code{urban_account} Automatically from Tidy Panel Tables
#'
#' Convenience constructor for applied researchers starting from raw tidy CSVs
#' or data frames rather than pre-assembled matrices. Extracts initial-period
#' (\code{base_time}) resident bilateral spending and external demand, aligns
#' origin and destination IDs, infers consistent baseline budgets and sales
#' totals when \code{origins} or \code{destinations} are omitted, and returns a
#' ready-to-use \code{\link{urban_account}} plus normalized shock panel.
#'
#' @param bilateral_panel Data frame with columns \code{origin}, \code{dest},
#'   \code{time}, and positive \code{spending} (or \code{log_spending}).
#' @param commuting Data frame (\code{origin}, \code{dest}, \code{flow}) or
#'   matrix of initial commuting/factor-income flows.
#' @param external_demand Data frame with columns \code{dest}, \code{time}, and
#'   either \code{external_spending} (levels, converted automatically to
#'   \eqn{h_{it} = \log(E^{\mathrm{ext}}_{it} / E^{\mathrm{ext}}_{i,0})}) or
#'   pre-computed log shift \code{h}.
#' @param base_time Optional value of \code{time} used as the initial equilibrium
#'   period (defaults to the earliest period in \code{bilateral_panel}).
#' @param origins Optional origin metadata table (\code{id}, \code{pop},
#'   \code{v0}, \code{s_H}, \code{is_city_resident}). When omitted, \code{v0}
#'   is calibrated consistently from initial local spending and default
#'   \code{s_H = 0.32} and \code{s_out_NH = 0.12}.
#' @param destinations Optional destination metadata table (\code{id},
#'   \code{y0}, \code{sector_group}, \code{theta_G}). When omitted, \code{y0}
#'   is set to initial resident + external spending.
#' @param default_s_H Default housing/saving share when \code{origins} is omitted (0.32).
#' @param default_s_out_NH Default conditional outside-good share when \code{origins$v0} is omitted (0.12).
#' @param default_theta_G Default local factor share of sales when \code{destinations$theta_G} is omitted (0.45).
#' @param default_resident_factor_share Default share of local factor payments
#'   accruing to modeled origins when \code{commuting} is in raw counts/payroll (0.85).
#'
#' @return A named list with \code{account} (an \code{\link{urban_account}}),
#'   \code{shocks} (normalized shock data frame), and \code{base_time}.
#' @export
#' @examples
#' data(sim_city_160)
#' ext_df <- data.frame(
#'   dest = sim_city_160$shocks$dest,
#'   time = sim_city_160$shocks$time,
#'   external_spending = sim_city_160$shocks$E_ext
#' )
#' built <- urban_from_tables(
#'   bilateral_panel = sim_city_160$bilateral_panel,
#'   commuting       = sim_city_160$H_0,
#'   external_demand = ext_df,
#'   origins         = sim_city_160$origins,
#'   destinations    = sim_city_160$destinations
#' )
#' summary(built$account)
urban_from_tables <- function(bilateral_panel,
                              commuting,
                              external_demand,
                              base_time = NULL,
                              origins = NULL,
                              destinations = NULL,
                              default_s_H = 0.32,
                              default_s_out_NH = 0.12,
                              default_theta_G = 0.45,
                              default_resident_factor_share = 0.85) {
  bp <- as.data.frame(bilateral_panel)
  if (!all(c("origin", "dest", "time") %in% names(bp))) {
    stop("bilateral_panel must have columns 'origin', 'dest', and 'time'.")
  }
  if (!("spending" %in% names(bp))) {
    if ("log_spending" %in% names(bp)) bp$spending <- exp(bp$log_spending)
    else stop("bilateral_panel must have column 'spending' or 'log_spending'.")
  }

  times <- sort(unique(as.character(bp$time)))
  if (is.null(base_time)) base_time <- times[1L]
  base_time <- as.character(base_time)

  orig_ids <- if (!is.null(origins) && "id" %in% names(origins)) {
    as.character(origins$id)
  } else {
    sort(unique(as.character(bp$origin)))
  }
  dest_ids <- if (!is.null(destinations) && "id" %in% names(destinations)) {
    as.character(destinations$id)
  } else {
    sort(unique(c(as.character(bp$dest), as.character(external_demand$dest))))
  }

  bp0 <- bp[as.character(bp$time) == base_time, c("origin", "dest", "spending"), drop = FALSE]
  E_R0 <- .coerce_bip_matrix(bp0, orig_ids, dest_ids, "bilateral_panel[base_time]")

  ext <- as.data.frame(external_demand)
  if (!all(c("dest", "time") %in% names(ext))) {
    stop("external_demand must contain columns 'dest' and 'time'.")
  }
  ext0_sub <- ext[as.character(ext$time) == base_time, , drop = FALSE]
  d_ix0 <- match(as.character(ext0_sub$dest), dest_ids)

  if ("external_spending" %in% names(ext)) {
    E_ext0 <- numeric(length(dest_ids))
    ok0 <- !is.na(d_ix0)
    E_ext0[d_ix0[ok0]] <- as.numeric(ext0_sub$external_spending[ok0])
    # Compute h = log(E_ext_it / E_ext_i0) safely
    e_base_map <- setNames(E_ext0, dest_ids)
    eb <- e_base_map[as.character(ext$dest)]
    ext$h <- ifelse(!is.na(eb) & eb > 0 & ext$external_spending > 0,
                    log(as.numeric(ext$external_spending) / eb), 0.0)
  } else if ("h" %in% names(ext)) {
    E_ext0 <- if (!is.null(destinations) && "E_ext0" %in% names(destinations)) {
      as.numeric(destinations$E_ext0)
    } else {
      0.15 * colSums(E_R0)
    }
  } else {
    stop("external_demand must contain 'external_spending' or 'h'.")
  }

  if (is.null(origins)) {
    local_r <- pmax(1e-8, rowSums(E_R0))
    nu <- 1.0 - default_s_H
    v0_inf <- local_r / (nu * (1.0 - default_s_out_NH))
    origins <- data.frame(
      id = orig_ids,
      pop = local_r / mean(local_r),
      v0 = v0_inf,
      s_H = rep(default_s_H, length(orig_ids)),
      is_city_resident = rep(TRUE, length(orig_ids)),
      stringsAsFactors = FALSE
    )
  } else {
    origins <- as.data.frame(origins)
    if (!("v0" %in% names(origins))) {
      sH <- if ("s_H" %in% names(origins)) origins$s_H else default_s_H
      origins$v0 <- pmax(1e-8, rowSums(E_R0)) / ((1.0 - sH) * (1.0 - default_s_out_NH))
    }
  }

  if (is.null(destinations)) {
    destinations <- data.frame(
      id = dest_ids,
      y0 = pmax(1e-8, colSums(E_R0) + E_ext0),
      E_ext0 = E_ext0,
      sector_group = 1L,
      theta_G = default_theta_G,
      resident_factor_share = default_resident_factor_share,
      stringsAsFactors = FALSE
    )
  } else {
    destinations <- as.data.frame(destinations)
    if (!("E_ext0" %in% names(destinations))) destinations$E_ext0 <- E_ext0
    if (!("y0" %in% names(destinations))) destinations$y0 <- pmax(1e-8, colSums(E_R0) + destinations$E_ext0)
  }

  acct <- urban_account(
    spending_base = E_R0,
    commuting_base = commuting,
    origins = origins,
    destinations = destinations,
    balance_method = "reconcile"
  )

  shocks_out <- ext[as.character(ext$dest) %in% dest_ids, intersect(c("dest", "time", "h", "external_spending"), names(ext)), drop = FALSE]
  list(account = acct, shocks = shocks_out, base_time = base_time)
}

#' One-Call End-to-End Urban Welfare Incidence Pipeline
#'
#' Runs the complete five-step sufficient-statistics pipeline—spatial account
#' validation, network exposure construction, bilateral spending estimation, and
#' open-city anchored welfare recovery—from tidy panel tables in a single call,
#' returning all intermediate objects alongside a provenance manifest.
#'
#' @param bilateral_panel Tidy bilateral spending panel (\code{origin},
#'   \code{dest}, \code{time}, \code{spending} or \code{log_spending}).
#' @param commuting Commuting matrix or 3-column table (\code{origin}, \code{dest}, \code{flow}).
#' @param external_demand Destination demand table (\code{dest}, \code{time},
#'   \code{external_spending} or \code{h}).
#' @param origins Optional origin metadata table.
#' @param destinations Optional destination metadata table.
#' @param base_time Optional initial base period.
#' @param estimator Bilateral estimator (\code{"ols"}, \code{"ppml"}, or \code{"iv"}).
#' @param price_shares Price-share rule passed to \code{\link{recover_welfare}}.
#' @param period Focal shock period for welfare summary.
#' @param se_method Standard error method (\code{"bootstrap"}, \code{"delta"}, or \code{"none"}).
#' @param ... Additional arguments passed to \code{\link{recover_welfare}}.
#' @return A list of class \code{"urban_pipeline"} with \code{account},
#'   \code{exposures}, \code{fit}, \code{welfare}, and \code{manifest}.
#' @export
run_urban_welfare <- function(bilateral_panel,
                              commuting,
                              external_demand,
                              origins = NULL,
                              destinations = NULL,
                              base_time = NULL,
                              estimator = "ols",
                              price_shares = "estimate_by_group",
                              period = NULL,
                              se_method = "delta",
                              ...) {
  built <- urban_from_tables(
    bilateral_panel = bilateral_panel,
    commuting = commuting,
    external_demand = external_demand,
    base_time = base_time,
    origins = origins,
    destinations = destinations
  )
  exps <- build_exposures(built$account, built$shocks, design = "direct")
  fit <- estimate_bilateral(bilateral_panel, exps, estimator = estimator)
  welf <- recover_welfare(fit, built$account, price_shares = price_shares,
                          period = period, se_method = se_method, ...)

  manifest <- list(
    package_version = as.character(utils::packageVersion("urbanwelfare")),
    timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    n_origins = built$account$N_orig,
    n_city_residents = sum(built$account$origins$is_city_resident),
    n_destinations = built$account$N_dest,
    n_periods = exps$N_time,
    n_bilateral_obs = fit$n_obs,
    rho_M = built$account$rho_M,
    vartheta = welf$headline$vartheta,
    terms_of_trade_gap = welf$headline$terms_of_trade_gap,
    mean_welfare_lp = welf$headline$mean_welfare_lp
  )

  structure(
    list(
      account = built$account,
      exposures = exps,
      fit = fit,
      welfare = welf,
      manifest = manifest
    ),
    class = "urban_pipeline"
  )
}
