#' Restraint-First Tufte ggplot2 Theme for Urban Incidence Exhibits
#'
#' Clean, minimalist \code{ggplot2} theme removing non-data ink, heavy gridlines,
#' and background boxes.
#'
#' @param base_size Base font size in points (default \code{11}).
#' @return A \code{ggplot2} theme object.
#' @export
theme_urban_tufte <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.35),
      axis.line.x = ggplot2::element_line(color = "#374151", linewidth = 0.4),
      axis.ticks.x = ggplot2::element_line(color = "#374151", linewidth = 0.35),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(color = "#4B5563", size = base_size - 1),
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    )
}

#' Plot Resident Welfare Decompositions, Spatial Maps, and QSM Parity
#'
#' @param x An \code{\link{urban_welfare}} object.
#' @param type Character: \code{"decomposition"} (histograms/densities of
#'   Income, Living Costs, and Net Welfare), \code{"deciles"} (bar/dot profile
#'   by welfare decile), \code{"map"} (spatial tract map using coordinates or
#'   \code{sf_tracts}), or \code{"qsm_parity"} (sufficient-statistics vs.
#'   benchmark vector scatter with 45-degree line).
#' @param benchmark Optional numeric vector of ground-truth or nonlinear QSM
#'   tract welfare (in log points) used when \code{type = "qsm_parity"}.
#' @param sf_tracts Optional \code{sf} polygon data frame with column \code{id}
#'   for choropleth rendering when \code{type = "map"}.
#' @param ... Unused.
#' @return A \code{ggplot2} object.
#' @export
plot.urban_welfare <- function(x,
                               type = c("decomposition", "deciles", "map", "qsm_parity"),
                               benchmark = NULL,
                               sf_tracts = NULL,
                               ...) {
  type <- match.arg(type)
  df <- x$tract_summary

  if (type == "decomposition") {
    long_df <- rbind(
      data.frame(component = "1. Nominal Income (+)", value = df$income_lp),
      data.frame(component = "2. Living Costs (-)", value = df$living_cost_lp),
      data.frame(component = "3. Net Resident Welfare", value = df$welfare_lp)
    )
    long_df$component <- factor(long_df$component, levels = unique(long_df$component))
    return(
      ggplot2::ggplot(long_df, ggplot2::aes(x = .data$value, fill = .data$component)) +
        ggplot2::geom_histogram(bins = 28, color = "white", alpha = 0.88) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "#374151", linewidth = 0.4) +
        ggplot2::facet_wrap(~ component, scales = "free_y", nrow = 1) +
        ggplot2::scale_fill_manual(values = c("#1D4ED8", "#B91C1C", "#047857")) +
        ggplot2::labs(
          title = sprintf("Resident Welfare Incidence (%s)", x$headline$period),
          subtitle = sprintf("Population-weighted mean: %.2f log pts (%.2f income - %.2f living costs)",
                             x$headline$mean_welfare_lp, x$headline$mean_income_lp, x$headline$mean_living_cost_lp),
          x = "Log points (100 x log change)",
          y = "Number of tracts"
        ) +
        theme_urban_tufte() +
        ggplot2::theme(legend.position = "none")
    )
  }

  if (type == "deciles") {
    dec <- x$deciles
    return(
      ggplot2::ggplot(dec, ggplot2::aes(x = factor(.data$decile), y = .data$welfare_lp)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#6B7280") +
        ggplot2::geom_segment(ggplot2::aes(xend = factor(.data$decile), y = 0, yend = .data$welfare_lp),
                              color = "#047857", linewidth = 0.9) +
        ggplot2::geom_point(color = "#047857", size = 2.8) +
        ggplot2::labs(
          title = "Distributional Resident Welfare by Tract Decile",
          x = "Resident Welfare Decile",
          y = "Mean Welfare Change (log points)"
        ) +
        theme_urban_tufte()
    )
  }

  if (type == "map") {
    if (!is.null(sf_tracts) && requireNamespace("sf", quietly = TRUE)) {
      merged <- merge(sf_tracts, df, by = "id")
      return(
        ggplot2::ggplot(merged) +
          ggplot2::geom_sf(ggplot2::aes(fill = .data$welfare_lp), color = "white", linewidth = 0.1) +
          ggplot2::scale_fill_gradient2(low = "#B91C1C", mid = "#F3F4F6", high = "#1D4ED8", midpoint = 0) +
          ggplot2::labs(title = "Spatial Resident Welfare Incidence", fill = "Log pts") +
          theme_urban_tufte()
      )
    }
    if (!all(c("x", "y") %in% names(df))) {
      stop("type = 'map' requires either x/y coordinates in origins or an sf_tracts object.")
    }
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, color = .data$welfare_lp, size = .data$pop)) +
        ggplot2::geom_point(alpha = 0.88) +
        ggplot2::scale_color_gradient2(low = "#B91C1C", mid = "#E5E7EB", high = "#047857",
                                       midpoint = stats::median(df$welfare_lp)) +
        ggplot2::scale_size_continuous(range = c(1.8, 5.2), guide = "none") +
        ggplot2::coord_equal() +
        ggplot2::labs(
          title = "Neighborhood Resident Welfare Incidence",
          x = "East-West coordinate",
          y = "North-South coordinate",
          color = "Welfare (lp)"
        ) +
        theme_urban_tufte()
    )
  }

  if (is.null(benchmark)) {
    stop("type = 'qsm_parity' requires supplying a numeric 'benchmark' vector of tract welfare.")
  }
  cmp <- data.frame(empirical = df$welfare_lp, benchmark = as.numeric(benchmark), pop = df$pop)
  r_val <- stats::cor(cmp$empirical, cmp$benchmark)
  ggplot2::ggplot(cmp, ggplot2::aes(x = .data$benchmark, y = .data$empirical, size = .data$pop)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#374151") +
    ggplot2::geom_point(shape = 21, color = "#1D4ED8", fill = "#93C5FD", alpha = 0.75) +
    ggplot2::scale_size_continuous(range = c(1.5, 4.8), guide = "none") +
    ggplot2::labs(
      title = "Sufficient-Statistics Recovery vs. Benchmark Equilibrium",
      subtitle = sprintf("Cross-tract correlation r = %.3f (45-degree line dashed)", r_val),
      x = "Benchmark / QSM Welfare (log points)",
      y = "Estimated Sufficient-Statistics Welfare (log points)"
    ) +
    theme_urban_tufte()
}

#' Plot Optimal Spatial Tax Schedules and Ramsey Decompositions
#'
#' @param x An \code{\link{urban_tax}} object.
#' @param type Character: \code{"distribution"} (histogram of optimal
#'   place-specific taxes with uniform rate benchmark) or \code{"map"} (spatial
#'   schedule).
#' @param ... Unused.
#' @return A \code{ggplot2} object.
#' @export
plot.urban_tax <- function(x, type = c("distribution", "map"), ...) {
  type <- match.arg(type)
  if (!is.null(x$qsm)) {
    df <- x$qsm$decomposition
    df$rate_pct <- 100 * df$tau_opt
    unif_pct <- 100 * x$qsm$tau_uniform
  } else {
    df <- data.frame(
      dest = x$account$dest_ids,
      rate_pct = 100 * x$empirical$tau_place,
      x = if ("x" %in% names(x$account$destinations)) x$account$destinations$x else NA_real_,
      y = if ("y" %in% names(x$account$destinations)) x$account$destinations$y else NA_real_
    )
    unif_pct <- 100 * x$empirical$tau_uniform
  }

  if (type == "distribution") {
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = .data$rate_pct)) +
        ggplot2::geom_histogram(bins = 26, fill = "#1E3A8A", color = "white", alpha = 0.88) +
        ggplot2::geom_vline(xintercept = unif_pct, color = "#B91C1C", linetype = "dashed", linewidth = 0.7) +
        ggplot2::labs(
          title = "Optimal Place-Specific Visitor Tax Schedule",
          subtitle = sprintf("Dashed red line: optimal city-wide uniform rate (%.1f%%)", unif_pct),
          x = "Optimal Ad Valorem Tax Rate (%)",
          y = "Number of destinations"
        ) +
        theme_urban_tufte()
    )
  }

  if (!all(c("x", "y") %in% names(df)) || anyNA(df$x)) {
    stop("type = 'map' requires x/y coordinates in destinations.")
  }
  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, color = .data$rate_pct)) +
    ggplot2::geom_point(size = 3.0, alpha = 0.9) +
    ggplot2::scale_color_viridis_c(option = "C") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Geography of Optimal Place-Specific Visitor Taxes",
      x = "East-West coordinate",
      y = "North-South coordinate",
      color = "Tax (%)"
    ) +
    theme_urban_tufte()
}
