# urbanwelfare

**Sufficient-Statistics Welfare Incidence and Optimal Policy in Urban Networks**

`urbanwelfare` is a standalone R + C++ (`Rcpp`) package implementing the methodology of **Allen, Fuchs, Ganapati, Graziano, Madera, and Montoriol-Garriga (2026, *"Urban Welfare: Tourism in Barcelona"*)**.

Any researcher with an intra-city bilateral spending (or trip) panel, a workplace commuting/ownership matrix, and a localized demand shock can recover neighborhood nominal income growth, living-cost changes, and money-metric resident welfare **without unit price indices or parametric assumptions on preferences or production**—and optionally compare against a calibrated nonlinear Quantitative Spatial Model (QSM) and solve for optimal place-based taxes.

---

## Installation

```r
# Install from local clone or GitHub:
# remotes::install_github("sfuchs-de/urbanwelfare")
install.packages(".", repos = NULL, type = "source")
```

---

## End-to-End Workflow (~20 Lines of R)

```r
library(urbanwelfare)

# 0. Load built-in synthetic 160-tract city (Section 2.8 Monte Carlo DGP)
data(sim_city_160)

# 1. Validate two-sector open-city account & build (Gamma_v, Omega^v, M)
acct <- urban_account(
  spending_base  = sim_city_160$E_R0,
  commuting_base = sim_city_160$H_0,
  origins        = sim_city_160$origins,      # v0, pop, s_H, is_city_resident
  destinations   = sim_city_160$destinations, # y0, E_ext0, sector_group
  balance_method = "reconcile"                # "reconcile", "ipf", or "strict"
)
summary(acct) # prints leakage breakdown, spectral radius rho(M), multiplier

# 2. Construct direct & network exposures (supports "direct" shocks OR "shift_share")
exposures <- build_exposures(
  account     = acct,
  shock_panel = sim_city_160$shocks,
  design      = "direct",                     # or "shift_share" w/ Borusyak-Hull controls
  center_dest = TRUE,                         # C++ basket-centers \tilde{h}_{nit}, \tilde{F}^y_{nit}
  horizons    = c(0, 1, 2, Inf)
)

# 3. Step 1: Bilateral spending regression (wraps fixest::feols / fepois / 2SLS)
fit <- estimate_bilateral(
  panel     = sim_city_160$bilateral_panel,
  exposures = exposures,
  estimator = "ols",                          # "ols", "ppml", or "iv"
  cluster   = c("origin^time", "dest^time")
)

# 4. Steps 2 & 3: Producer price shares (\kappa^S_g), common level (\hat{c}_t), & EV
welf <- recover_welfare(
  fit          = fit,
  account      = acct,
  price_shares = "estimate_by_group",
  level_anchor = "open_city_factor_balance",  # Proposition 2.5
  se_method    = "bootstrap"
)
summary(welf)

# 5. Optional extensions & optimal spatial tax module
welf_ext <- welf |>
  add_shopping_wedge(sigma = c(2.0, 3.456, 8.235)) |>
  add_workplace_mobility(eps = 3.3) |>
  add_finite_bounds(sigma_ces = 3.456)

tax <- optimal_spatial_tax(
  welfare = welf,
  account = acct,
  demand_elasticities = list(sigma_T = 7.0, eta_out = 2.5),
  mode = "both" # empirical linear Ramsey + C++ nonlinear QSM adjoint optimum
)
summary(tax)

plot(welf, type = "decomposition")
```

---

## Core Modules & Mathematical Mapping

| Function | Role | Paper Reference |
| --- | --- | --- |
| `urban_account()` | Validates open-city accounting closure (or IPF-balances noisy matrices), computes leakage shares, builds $\Gamma_v$, $\Omega^v$, $M = \Omega^v\Gamma_v$, and verifies $\rho(M) < 1$ | Section 2.1–2.4 |
| `build_exposures()` | Computes $\mathbf F_t^y = (I - M)^{-1}\mathbf d_t$, $\mathbf F_t^v = \Gamma_v\mathbf F_t^y$, finite-round bounds, and Borusyak–Hull propagated expected shift-share controls | Proposition 2.3(a), Section 3.3 |
| `estimate_bilateral()` | C++ basket-centering (`cpp_basket_center`) + `fixest` HDFE estimation (`OLS`, `PPML`, `2SLS`), plus `rotemberg_weights()` and `loo_influence()` | Proposition 2.4 |
| `recover_welfare()` | Spatial factor-income matching for $\widehat{\boldsymbol\kappa}^S_g$, open-city common shift $\widehat c_t = -\psi_t / (1 - \vartheta)$, money-metric $\mathrm dEV_{nt}/v_n^0$, and open-city terms-of-trade check | Propositions 2.2 & 2.5 |
| `add_shopping_wedge()`, `add_workplace_mobility()`, `add_finite_bounds()` | Redding–Weinstein CES amenity residuals, Fréchet Le Chatelier wage/price rules, and Diewert + C++ BVLS producer-concave finite-change bounds | Section 5.2 |
| `solve_qsm()`, `optimal_spatial_tax()` | C++ damped-Newton nonlinear QSM solver + analytical adjoint Ramsey optimizer for uniform and place-specific taxes with 4-channel decomposition | Section 6 & Appendix B |
| `sim_urban_city()` | Synthetic multi-sector gravity cities (`CES`, non-CES `HSA` + capacity limits, heterogeneous elasticities) with known ground-truth welfare | Section 2.8 |

---

## Empirical Applications Beyond Tourism

By treating the localized demand shift $d_{it}$ generically, `urbanwelfare` applies directly to modern urban transaction and mobility panels (payment cards, Pix, scanner data, SafeGraph/Advan + LEHD/LODES):

- **Tourism, conventions, & mega-events**: Incidence of external visitor spending across neighborhoods.
- **Return-to-office / WFH shocks**: Reallocation of weekday commuter spending between CBDs and residential suburbs.
- **Transit openings, pedestrianization (*Superblocks*), & congestion pricing**: Neighborhood incidence of commercial corridor demand shifts.
- **Place-based taxation & short-term rental regulation**: Combining sufficient-statistics incidence with `optimal_spatial_tax()`.

---

## Citation

> Allen, Treb, Simon Fuchs, Sharat Ganapati, Alberto Graziano, Rocio Madera, and Judit Montoriol-Garriga (2026). *"Urban Welfare: Tourism in Barcelona."* Working Paper.
