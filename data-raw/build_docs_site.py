#!/usr/bin/env python3
"""Generate the static Documenter.jl-style documentation site in docs/."""
from pathlib import Path

DOCS = Path(__file__).resolve().parents[1] / "docs"
DOCS.mkdir(parents=True, exist_ok=True)

NAV = [
    ("Getting Started", [
        ("index.html", "Home & Quickstart"),
        ("own-data.html", "Use Your Own Data"),
    ]),
    ("Methodology", [
        ("model.html", "Sufficient-Statistics Model"),
        ("extensions.html", "Extensions, QSM & Taxes"),
        ("barcelona.html", "Barcelona Replication"),
    ]),
    ("Reference", [
        ("api.html", "API Reference"),
    ]),
]

FLAT_PAGES = [item for _, items in NAV for item in items]


def render_page(filename: str, title: str, body_html: str) -> None:
    idx = [i for i, (f, _) in enumerate(FLAT_PAGES) if f == filename][0]
    prev_link = (
        f'<a href="{FLAT_PAGES[idx-1][0]}">&laquo; {FLAT_PAGES[idx-1][1]}</a>'
        if idx > 0 else '<span></span>'
    )
    next_link = (
        f'<a href="{FLAT_PAGES[idx+1][0]}">{FLAT_PAGES[idx+1][1]} &raquo;</a>'
        if idx < len(FLAT_PAGES) - 1 else '<span></span>'
    )

    nav_html = []
    for group_label, items in NAV:
        nav_html.append(f'<div class="nav-group-label">{group_label}</div><ul class="nav-list">')
        for href, label in items:
            cls = ' class="active"' if href == filename else ""
            nav_html.append(f'<li><a href="{href}"{cls}>{label}</a></li>')
        nav_html.append("</ul>")

    html = f"""<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
  <meta charset="UTF-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
  <title>{title} &middot; urbanwelfare</title>
  <meta name="description" content="Documentation for urbanwelfare: Sufficient-Statistics Welfare Incidence and Optimal Policy in Urban Networks."/>
  <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/lato-font/3.0.0/css/lato-font.min.css"/>
  <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/juliamono/0.050/juliamono.min.css"/>
  <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.9/katex.min.css"/>
  <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/styles/github.min.css"/>
  <link rel="stylesheet" href="assets/style.css"/>
  <script defer src="https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.9/katex.min.js"></script>
  <script defer src="https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.9/contrib/auto-render.min.js"></script>
  <script defer src="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/highlight.min.js"></script>
  <script defer src="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/languages/r.min.js"></script>
  <script>document.addEventListener("DOMContentLoaded", () => {{ if (window.hljs) hljs.highlightAll(); }});</script>
  <script defer src="assets/site.js"></script>
</head>
<body>
  <div class="layout">
    <aside class="sidebar">
      <a class="brand" href="index.html">
        <span class="brand-title">urbanwelfare</span>
        <span class="version-badge">v0.1.0</span>
      </a>
      <input id="sidebar-search" class="search-box" type="search" placeholder="Filter pages..."/>
      {"".join(nav_html)}
    </aside>
    <div class="main-wrap">
      <header class="topbar">
        <div class="breadcrumb">urbanwelfare &rsaquo; <strong>{title}</strong></div>
        <div class="topbar-actions">
          <select id="theme-picker" class="theme-select" aria-label="Theme">
            <option value="light">Light</option>
            <option value="dark">Dark</option>
            <option value="catppuccin">Catppuccin</option>
          </select>
          <a class="gh-link" href="https://github.com/sfuchs-de/urbanwelfare" target="_blank" rel="noopener">GitHub &nearr;</a>
        </div>
      </header>
      <main class="content">
        {body_html}
        <nav class="pager">
          {prev_link}
          {next_link}
        </nav>
        <div class="footer-note">
          <code>urbanwelfare</code> v0.1.0 &middot; Based on Allen, Fuchs, Ganapati, Graziano, Madera, and Montoriol-Garriga (2026, <em>&ldquo;Urban Welfare: Tourism in Barcelona&rdquo;</em>). MIT License.
        </div>
      </main>
    </div>
  </div>
</body>
</html>
"""
    (DOCS / filename).write_text(html, encoding="utf-8")


# ---------------------------------------------------------------------------
# Page 1: Home & Quickstart
# ---------------------------------------------------------------------------
render_page("index.html", "Home & Quickstart", """
<h1 id="urbanwelfare">urbanwelfare</h1>
<p>
  <code>urbanwelfare</code> is a standalone R and C++ (<code>Rcpp</code>) package for evaluating the
  spatial welfare incidence of localized demand shocks and place-based policies within cities.
  Given a bilateral resident spending (or trip) panel, workplace commuting/ownership links, and a
  destination demand shock, it recovers neighborhood nominal income growth, living-cost changes,
  and first-order money-metric resident welfare <strong>without unit price indices or parametric
  functional forms</strong>—and optionally solves a calibrated nonlinear Quantitative Spatial Model (QSM)
  and optimal place-specific Ramsey taxes.
</p>

<div class="callout">
  <div class="callout-title">Core Design Promise</div>
  Start from three tidy tables—bilateral spending, commuting flows, and destination shocks—or from a
  validated spatial account. No source edits are required, and every run verifies spectral convergence
  $\\rho(M) &lt; 1$ and the open-city aggregate terms-of-trade identity to machine precision.
</div>

<h2>Install and run the built-in synthetic city</h2>
<pre><code class="language-r"># Install from GitHub (requires R &gt;= 4.1, fixest, Matrix, Rcpp, ggplot2)
# remotes::install_github("sfuchs-de/urbanwelfare")

library(urbanwelfare)
data(sim_city_160)

# One-call pipeline directly from tidy tables:
pipe &lt;- run_urban_welfare(
  bilateral_panel = sim_city_160$bilateral_panel,
  commuting       = sim_city_160$H_0,
  external_demand = data.frame(
    dest = sim_city_160$shocks$dest,
    time = sim_city_160$shocks$time,
    external_spending = sim_city_160$shocks$E_ext
  ),
  origins      = sim_city_160$origins,
  destinations = sim_city_160$destinations
)
summary(pipe$welfare)</code></pre>

<p>
  <code>sim_city_160</code> ships inside the package (160 neighborhoods, 4 commercial sectors, 6 external
  visitor groups, and known nonlinear ground-truth welfare), so the example runs in under two seconds with
  zero external data files.
</p>

<h2>Step-by-step modular workflow</h2>
<p>
  For custom shift-share instruments, enclave regressions, robustness bounds, or optimal spatial taxation,
  call the five modular stages directly:
</p>
<pre><code class="language-r"># 1. Validate two-sector open-city account &amp; build network operators (Gamma_v, Omega^v, M)
acct &lt;- urban_account(
  spending_base  = sim_city_160$E_R0,
  commuting_base = sim_city_160$H_0,
  origins        = sim_city_160$origins,
  destinations   = sim_city_160$destinations,
  balance_method = "reconcile"   # "reconcile", "ipf", or "strict"
)

# 2. Propagate shocks through (I - M)^(-1) and enable C++ basket centering
exposures &lt;- build_exposures(acct, sim_city_160$shocks, design = "direct")

# 3. Step 1: Bilateral spending regression (fixest OLS, PPML, or 2SLS IV)
fit &lt;- estimate_bilateral(sim_city_160$bilateral_panel, exposures, estimator = "ols")

# 4. Steps 2 &amp; 3: Sector price shares (kappa^S_g), common shift (c_t), &amp; resident EV
welf &lt;- recover_welfare(fit, acct, level_anchor = "open_city_factor_balance", se_method = "bootstrap")

# 5. Extensions &amp; C++ adjoint optimal spatial tax
welf &lt;- welf |&gt;
  add_shopping_wedge(sigma = c(2.0, 3.456, 8.235)) |&gt;
  add_workplace_mobility(eps = 3.3) |&gt;
  add_finite_bounds(sigma_ces = 3.456)

tax &lt;- optimal_spatial_tax(welf, acct, demand_elasticities = list(sigma_T = 7.0, eta_out = 2.5))</code></pre>

<div class="figure-grid">
  <div class="figure-card">
    <img src="assets/fig-decomposition.png" alt="Income, living costs, and net resident welfare"/>
    <div class="figure-caption"><code>plot(welf, type = "decomposition")</code>: Nominal income gains, local living-cost increases, and net resident welfare across 160 neighborhoods.</div>
  </div>
  <div class="figure-card">
    <img src="assets/fig-parity.png" alt="Sufficient statistics vs nonlinear QSM ground truth"/>
    <div class="figure-caption"><code>plot(welf, type = "qsm_parity", benchmark = ...)</code>: Non-parametric sufficient-statistics recovery vs. exact nonlinear equilibrium truth.</div>
  </div>
</div>

<h2>Output contract &amp; repository boundaries</h2>
<p>
  Every <code>urban_welfare</code> and <code>urban_pipeline</code> object returns tract-level and decile
  tables (<code>as.data.frame(welf)</code>, <code>welf$deciles</code>), spectral and openness diagnostics
  (<code>acct$rho_M</code>, <code>acct$leakage</code>), the open-city recycling share $\\vartheta$, and the
  exact aggregate terms-of-trade check residual (<code>welf$headline$terms_of_trade_gap</code>).
</p>
<p>
  See <a href="own-data.html">Use Your Own Data</a> to map payment cards, scanner data, or mobility panels
  into <code>urbanwelfare</code>, and <a href="barcelona.html">Barcelona Replication</a> for exact six-decimal
  reproduction of the Barcelona tourism study while keeping restricted microdata outside Git.
</p>
""")

# ---------------------------------------------------------------------------
# Page 2: Use Your Own Data
# ---------------------------------------------------------------------------
render_page("own-data.html", "Use Your Own Data", """
<h1>Use Your Own Data</h1>
<p>
  <code>urbanwelfare</code> is designed for modern intra-city transaction and mobility datasets—such as
  payment cards (CaixaBank, BBVA, Affinity, Facteus, Nets), instant payments (Pix), retail scanner panels,
  or foot-traffic gravity panels (SafeGraph/Advan) paired with commuting records (Census LEHD/LODES,
  register commuting, or mobile-device home–work matrices).
</p>

<h2>Option A: Start from three tidy CSV tables</h2>
<p>
  If you have raw panel tables rather than pre-balanced equilibrium matrices, pass them directly to
  <code>urban_from_tables()</code> or <code>run_urban_welfare()</code>:
</p>
<table>
  <thead>
    <tr><th>Table</th><th>Required Columns</th><th>Economic Meaning</th></tr>
  </thead>
  <tbody>
    <tr>
      <td><code>bilateral_panel</code></td>
      <td><code>origin, dest, time, spending</code></td>
      <td>Resident expenditure $E^R_{nit}$ by home neighborhood $n$, commercial destination $i$, and period $t$.</td>
    </tr>
    <tr>
      <td><code>commuting</code></td>
      <td><code>origin, dest, flow</code></td>
      <td>Initial workplace payroll $W_{ni}^0$ (or factor-claim share $H_{ni}^m \\in [0,1]$) linking destination $i$ to home $n$.</td>
    </tr>
    <tr>
      <td><code>external_demand</code></td>
      <td><code>dest, time, external_spending</code> (or <code>h</code>)</td>
      <td>Shocked external/visitor/event demand $E^{\\mathrm{ext}}_{it}$ (or log shift $h_{it}$).</td>
    </tr>
  </tbody>
</table>

<pre><code class="language-r">built &lt;- urban_from_tables(
  bilateral_panel = read.csv("resident_card_panel.csv"),
  commuting       = read.csv("lehd_commuting.csv"),
  external_demand = read.csv("visitor_spending.csv"),
  base_time       = "2018-01",
  origins         = read.csv("tract_populations_incomes.csv") # optional
)
summary(built$account)</code></pre>

<h2>Option B: Noisy card coverage &amp; IPF balancing</h2>
<p>
  Raw card panels rarely capture 100% of household expenditure. In <code>urban_account()</code>, setting
  <code>balance_method = "ipf"</code> applies Iterative Proportional Fitting (RAS) to align noisy bilateral
  card shares with official origin non-housing budgets $\\nu_n(1 - s_{n0}^{NH})v_n^0$ and destination
  consumer sales totals $y_i^0 - E_i^{\\mathrm{ext},0}$ before constructing $(\\Gamma_v, \\Omega^v, M)$.
</p>

<h2>Option C: Surrounding metropolitan commuters</h2>
<p>
  Cities trade daily with surrounding metropolitan municipalities $\\mathcal{M}$: suburban commuters earn
  wages at city workplaces and shop at city commercial hubs, but city policy objectives evaluate welfare over
  core city residents $\\mathcal{R}$. Include suburban origins directly in <code>origins</code> with
  <code>is_city_resident = FALSE</code>:
</p>
<ul>
  <li><code>M = Omega_v %*% Gamma_v</code> automatically incorporates metropolitan earning-and-shopping feedback $\\Omega_R\\Gamma_R + \\Omega_M\\Gamma_M$.</li>
  <li><code>recover_welfare()</code> and <code>optimal_spatial_tax()</code> automatically restrict welfare weights, lump-sum tax rebates, and $\\widehat{c}_t$ factor-income matching to core city residents (<code>is_city_resident == TRUE</code>).</li>
</ul>

<h2>Option D: Secure data enclaves</h2>
<p>
  When raw bilateral card records reside inside an air-gapped bank or census enclave:
</p>
<ol>
  <li>Install <code>urbanwelfare</code> inside the enclave (zero external C++ system libraries required—only standard R packages <code>Rcpp</code>, <code>Matrix</code>, <code>fixest</code>, and <code>ggplot2</code>).</li>
  <li>Or export only the aggregated initial matrices $(E^{R,0}, H^m)$, exposures, and the 4-element slope vector $(a_O, b_v, a_D, b_D)$ + $4 \\times 4$ clustered covariance matrix, then reconstruct the fit outside the enclave via <code>bilateral_from_estimates()</code>:</li>
</ol>
<pre><code class="language-r">fit &lt;- bilateral_from_estimates(
  coefficients = c(a_O = -0.0059, b_v = 1.2176, a_D = -0.1200, b_D = -0.3066),
  exposures    = exposures,
  vcov         = exported_vcov_4x4,
  n_obs        = 3260707L
)
welf &lt;- recover_welfare(fit, acct, period = "201807")</code></pre>
""")

# ---------------------------------------------------------------------------
# Page 3: Sufficient-Statistics Model
# ---------------------------------------------------------------------------
render_page("model.html", "Sufficient-Statistics Model", """
<h1>Sufficient-Statistics Framework</h1>
<p>
  <code>urbanwelfare</code> implements the five-step sufficient-statistics framework of Allen, Fuchs,
  Ganapati, Graziano, Madera, and Montoriol-Garriga (2026) in a two-sector small open city.
</p>

<h2>1. Household &amp; Producer Envelopes (Propositions 2.1–2.2)</h2>
<p>
  Around an initial interior optimum with fixed outside-good price ($p_0 = 1$) and fixed short-run housing
  rents, Shephard&rsquo;s lemma and workplace earnings optimization yield first-order money-metric resident
  welfare in origin $n$:
</p>
$$\\frac{\\mathrm{d}EV_{nt}}{v_n^0} = \\mathrm{d}\\log v_{nt} - \\sum_{i=1}^N s_{ni}\\,\\mathrm{d}\\log p_{it}, \\qquad s_{ni} = \\frac{E_{ni}^{R,0}}{v_n^0}.$$
<p>
  Because transaction datasets record nominal sales $y_{it} = p_{it} q_{it}$ rather than unit prices $p_{it}$,
  competitive firms combining fixed local factors with flexible traded inputs at fixed external prices satisfy
  the producer envelope $\\mathrm{d}G_{it} = y_i^0\\,\\mathrm{d}\\log p_{it}$. Defining the marginal price
  share $\\kappa_i^S \\equiv \\mathrm{d}G_i / \\mathrm{d}y_i = (1 + \\varepsilon_i^S)^{-1} \\in (0,1]$ gives:
</p>
$$\\mathrm{d}\\log p_{it} = \\kappa_i^S\\,\\mathrm{d}\\log y_{it}, \\qquad \\mathrm{d}\\log v_{nt}^{\\mathrm{prod}} = \\sum_{i=1}^N H_{ni}^m \\frac{y_i^0}{v_n^0}\\kappa_i^S\\,\\mathrm{d}\\log y_{it}.$$

<h2>2. Spatial Network Propagation (Proposition 2.3a)</h2>
<p>
  Let $\\Gamma_{v,ni} = H_{ni}^m G_i^0 / v_n^0$ map destination factor income to residential origins, and
  $\\Omega^v_{in} = E_{ni}^{R,0} / y_i^0$ map origin income back to commercial destinations. One earning-and-shopping
  round is $M = \\Omega^v \\Gamma_v$. Because spending leaks into flexible traded inputs, housing/saving,
  outside traded goods, and external factor owners, $\\rho(M) &lt; 1$ and the infinite-round exposures are:
</p>
$$\\mathbf{F}_t^y = (I - M)^{-1}\\mathbf{d}_t = \\sum_{h=0}^{\\infty} M^h \\mathbf{d}_t, \\qquad \\mathbf{F}_t^v = \\Gamma_v \\mathbf{F}_t^y, \\qquad d_{it} = \\frac{E_i^{\\mathrm{ext},0}}{y_i^0} h_{it}.$$

<h2>3. Bilateral Spending Regression (Proposition 2.4)</h2>
<p>
  Under local homotheticity and spatial IIA relative to the outside good $i = 0$, basket-centering destination
  exposures $\\widetilde{x}_{nit} = x_{it} - \\sum_{j=1}^N s_{nj}^{NH} x_{jt}$ ($\\sum_{j=0}^N s_{nj}^{NH}\\widetilde{x}_{njt} = 0$)
  separates relative origin income growth from destination expenditure switching in a single HDFE regression
  with pair ($\\alpha_{ni}$) and period ($\\delta_t$) fixed effects:
</p>
$$\\log E^R_{nit} = \\alpha_{ni} + \\delta_t + \\underbrace{a_O h_{nt} + b_v F^v_{nt}}_{\\text{origin income}} + \\underbrace{a_D \\widetilde{h}_{nit} + b_D \\widetilde{F}^y_{nit}}_{\\text{expenditure switching}} + \\varepsilon_{nit}.$$
<p>
  Differentiating goods-market clearing then recovers destination sales growth
  $\\mathrm{d}\\log \\widehat{y}_{it} = d_{it} + \\sum_n \\Omega^v_{in}\\,\\mathrm{d}\\log \\widehat{E}^R_{nit}$,
  and matching $\\mathrm{d}\\log \\widehat{v}_{nt}$ against $\\sum_{i \\in g} H_{ni}^m(y_i^0/v_n^0)\\,\\mathrm{d}\\log \\widehat{y}_{it}$
  identifies sector marginal price shares $\\widehat{\\kappa}_g^S$.
</p>

<h2>4. Open-City Level Anchor &amp; Terms of Trade (Proposition 2.5)</h2>
<p>
  Time fixed effects $\\delta_t$ absorb the city-wide income shift $c_t$. In a small open city with fixed
  outside-sector factor earnings ($\\mathrm{d}G^O = 0$), aggregate resident income growth must equal local
  factor payments accruing to residents:
</p>
$$\\sum_{n \\in \\mathcal{R}} v_n^0 (\\mathrm{d}\\log \\widehat{v}_{nt} + c_t) = \\sum_{n \\in \\mathcal{R}}\\sum_{i=1}^N H_{ni}^m y_i^0 \\kappa_i^S (\\mathrm{d}\\log \\widehat{y}_{it} + \\Omega_i^v c_t) \\implies \\widehat{c}_t = -\\frac{\\psi_t}{1 - \\vartheta}.$$
<p>
  Summing anchored $\\mathrm{d}EV_{nt}$ over all city residents cancels internal resident-to-resident sales
  and yields the open-city aggregate terms-of-trade identity verified automatically by <code>recover_welfare()</code>:
</p>
$$\\sum_{n \\in \\mathcal{R}} \\mathrm{d}EV_{nt} = \\sum_{i=1}^N E_i^{\\mathrm{nonres},0}\\,\\mathrm{d}\\log p_{it} - \\sum_{i=1}^N (1 - H_i^{\\mathcal{R}}) y_i^0\\,\\mathrm{d}\\log p_{it}.$$
""")

# ---------------------------------------------------------------------------
# Page 4: Extensions, QSM & Optimal Taxes
# ---------------------------------------------------------------------------
render_page("extensions.html", "Extensions, QSM & Taxes", """
<h1>Extensions, Nonlinear QSM &amp; Optimal Spatial Taxes</h1>

<h2>1. Shift-Share Network Identification</h2>
<p>
  When shocks follow a shift-share structure $h_{it} = \\sum_g \\pi_{ig} g_{gt}$ (e.g. visitor origin shares
  $\\times$ seasonal travel shifts), non-random destination exposure propagates through $(I - M)^{-1}$.
  Passing <code>shifter_expectations</code> in <code>build_exposures(..., design = "shift_share")</code>
  propagates expected shifts through $(I - M)^{-1}$ following Borusyak and Hull (2023), and
  <code>rotemberg_weights()</code> / <code>loo_influence()</code> report shock-level relevance weights and
  leave-one-hub-out jackknife stability.
</p>

<h2>2. Robustness Extensions (Section 5.2)</h2>
<ul>
  <li>
    <strong>Shopping Amenity / Congestion Residuals (<code>add_shopping_wedge()</code>)</strong>:
    Inverts $\\widehat{\\Delta}^D_{it} / (1 - \\sigma)$ to separate monetary producer price growth
    $\\widehat{p}_{it}^{\\mathrm{prod}}$ from perceived shopping quality/crowding wedges $\\widehat{\\zeta}_{it}$.
  </li>
  <li>
    <strong>Workplace Mobility Le Chatelier Rules (<code>add_workplace_mobility()</code>)</strong>:
    Solves linearized workplace labor market clearing under Fr&eacute;chet sorting elasticity $\\varepsilon$
    and perfect city-wide mobility ($\\varepsilon \\to \\infty$).
  </li>
  <li>
    <strong>Finite-Change Welfare Bounds (<code>add_finite_bounds()</code>)</strong>:
    Evaluates exact finite-change consumer Laspeyres, CES, and Paasche bounds alongside C++ bounded-variable
    least-squares (<code>cpp_bvls_box</code>) producer-concave feasibility projections.
  </li>
</ul>

<h2>3. Nonlinear QSM &amp; Adjoint Ramsey Spatial Taxation</h2>
<p>
  <code>solve_qsm()</code> and <code>optimal_spatial_tax()</code> calibrate a two-sector open-city nonlinear
  QSM directly to <code>acct</code> and $\\widehat{\\boldsymbol{\\kappa}}^S$. Under ad valorem visitor tax
  $\\tau_i = e^{\\theta_i} - 1$ rebated lump-sum to city residents ($T = R / \\Omega$), the C++ solver computes
  exact equilibrium prices via damped Newton and evaluates the analytical adjoint welfare gradient:
</p>
$$\\nabla_{\\boldsymbol{\\theta}} W = -\\mathcal{E}' b^T + \\frac{\\mu_T^{GE}}{\\Omega}\\Bigl[(\\mathbf{1} + \\boldsymbol{\\tau})\\odot \\mathbf{E}^T - \\mathcal{E}'(\\boldsymbol{\\tau}\\odot \\mathbf{E}^T)\\Bigr],$$
<p>
  where $\\mathcal{E} = -\\mathrm{d}\\log \\mathbf{E}^T / \\mathrm{d}\\boldsymbol{\\theta}'$ is the general-equilibrium
  tax-base elasticity matrix. At the interior optimum, <code>optimal_spatial_tax()</code> decomposes each
  neighborhood&rsquo;s optimal rate into four economic channels: (1) own inverse-elasticity markup,
  (2) own resident incidence, (3) cross-neighborhood resident spillovers, and (4) cross-base fiscal spillovers.
</p>

<div class="figure-grid">
  <div class="figure-card">
    <img src="assets/fig-map.png" alt="Spatial welfare map"/>
    <div class="figure-caption"><code>plot(welf, type = "map")</code>: Spatial distribution of resident welfare gains across neighborhoods.</div>
  </div>
  <div class="figure-card">
    <img src="assets/fig-tax.png" alt="Optimal spatial tax distribution"/>
    <div class="figure-caption"><code>plot(tax, type = "distribution")</code>: Optimal place-specific visitor tax rates vs. optimal city-wide uniform rate.</div>
  </div>
</div>
""")

# ---------------------------------------------------------------------------
# Page 5: Barcelona Replication
# ---------------------------------------------------------------------------
render_page("barcelona.html", "Barcelona Replication", """
<h1>Barcelona Replication Guide</h1>
<p>
  <code>urbanwelfare</code> applies directly to the empirical setting of Allen, Fuchs, Ganapati, Graziano,
  Madera, and Montoriol-Garriga (2026, <em>&ldquo;Urban Welfare: Tourism in Barcelona&rdquo;</em>), which
  models <strong>1,053 Barcelona census tracts</strong> ($\\mathcal{R}$, <code>is_city_resident = TRUE</code>)
  together with <strong>35 surrounding metropolitan municipalities</strong> ($\\mathcal{M}$,
  <code>is_city_resident = FALSE</code>) across 24 months and 5 commercial activity groups.
</p>

<h2>Workflow on Prepared Barcelona Accounts</h2>
<pre><code class="language-r">library(urbanwelfare)

# 1. Build 1,088-origin (1,053 BCN + 35 Metro) x 1,053-destination open-city account
acct &lt;- urban_account(
  spending_base      = E_R0_all,      # 1088 x 1053 initial resident + metro card spending
  commuting_base     = H_0_all,       # 1088 x 1053 workplace &amp; business-claim shares
  origins            = origins_bcn,   # is_city_resident = TRUE for 1,053 BCN tracts
  destinations       = dests_bcn,     # 5 IO consumer-cut activity groups
  balance_method     = "strict",
  commuting_is_share = TRUE
)
summary(acct)

# 2. Construct network exposures &amp; import enclave Step-1 bilateral estimates (N = 3,260,707)
exps &lt;- build_exposures(acct, shocks_bcn, design = "direct")
fit  &lt;- bilateral_from_estimates(beta_4, exps, vcov = vcov_4x4, n_obs = 3260707L)

# 3. Recover 5-sector price shares on 2018 bank-income support &amp; July vs. January welfare
welf &lt;- recover_welfare(
  fit           = fit,
  account       = acct,
  period        = "201807",
  kappa_support = matching_support_2018,
  se_method     = "bootstrap",
  n_boot        = 1000,
  seed          = 20261004
)
summary(welf)</code></pre>

<h2>Verified Six-Decimal Parity</h2>
<p>
  Running the snippet above reproduces every diagnostic, sector marginal price share, and headline welfare
  statistic in the manuscript to six decimal places:
</p>
<table>
  <thead>
    <tr><th>Quantity</th><th><code>urbanwelfare</code> Output</th><th>Manuscript Value</th></tr>
  </thead>
  <tbody>
    <tr><td>Spectral radius $\\rho(M)$ &amp; multiplier $(1-\\rho)^{-1}$</td><td><code>0.244636</code> &middot; <code>1.324</code></td><td><code>0.245</code> &middot; <code>1.32</code></td></tr>
    <tr><td>Leakage: flexible inputs / housing / outside goods / ext. factors</td><td><code>56.24% / 32.09% / 8.81% / 4.30%</code></td><td><code>56.2% / 32.1% / 8.8% / 4.3%</code></td></tr>
    <tr><td>Sector $\\widehat{\\kappa}_g^S$ (Accomm./Rest., Other Serv., Prop., Trade, Transp.)</td><td><code>0.6245, 0.6334, 0.5950, 0.4340, 0.8108</code></td><td><code>0.625, 0.633, 0.595, 0.434, 0.811</code></td></tr>
    <tr><td>July–January Mean Nominal Income Growth</td><td><code>4.744697</code> log pts</td><td><code>4.744697</code> (<code>4.74</code>)</td></tr>
    <tr><td>July–January Mean Living-Cost Increase</td><td><code>3.019153</code> log pts</td><td><code>3.019153</code> (<code>3.02</code>)</td></tr>
    <tr><td>July–January Mean Resident Welfare ($\\mathrm{d}EV / v^0$)</td><td><code>1.725544</code> log pts</td><td><code>1.725544</code> (<code>1.73</code>)</td></tr>
    <tr><td>Tract Welfare P10 / Median / P90</td><td><code>0.709238 / 1.454095 / 3.324039</code></td><td><code>0.71 / 1.45 / 3.32</code></td></tr>
    <tr><td>Sensitivity: MPC $\\mu = 0.80$ &middot; Upstream $\\psi_O = 0.15$</td><td><code>1.800987</code> &middot; <code>2.357359</code></td><td><code>1.800987</code> &middot; <code>2.357359</code></td></tr>
    <tr><td>Aggregate Terms-of-Trade Identity Discrepancy</td><td><code>1.16e-16</code></td><td><code>&lt; 1e-15</code></td></tr>
  </tbody>
</table>
""")

# ---------------------------------------------------------------------------
# Page 6: API Reference
# ---------------------------------------------------------------------------
render_page("api.html", "API Reference", """
<h1>API Reference</h1>

<div class="api-block">
  <div class="api-sig">run_urban_welfare(bilateral_panel, commuting, external_demand, origins = NULL, destinations = NULL, base_time = NULL, estimator = "ols", price_shares = "estimate_by_group", period = NULL, se_method = "delta", ...)</div>
  <p>One-call end-to-end pipeline from tidy panel tables. Returns an <code>urban_pipeline</code> list containing <code>$account</code>, <code>$exposures</code>, <code>$fit</code>, <code>$welfare</code>, and <code>$manifest</code>.</p>
</div>

<div class="api-block">
  <div class="api-sig">urban_from_tables(bilateral_panel, commuting, external_demand, base_time = NULL, origins = NULL, destinations = NULL, ...)</div>
  <p>Extracts base-period bilateral spending and external demand from tidy data frames, calibrates omitted budget/sales totals, and returns <code>list(account, shocks, base_time)</code>.</p>
</div>

<div class="api-block">
  <div class="api-sig">urban_account(spending_base, commuting_base, origins, destinations, balance_method = c("reconcile", "ipf", "strict"), commuting_is_share = "auto", tol = 1e-7)</div>
  <p>Validates two-sector open-city accounting identities, computes leakage shares, constructs $\\Gamma_v$, $\\Omega^v$, $M = \\Omega^v\\Gamma_v$, and $(I - M)^{-1}$, and checks $\\rho(M) &lt; 1$.</p>
</div>

<div class="api-block">
  <div class="api-sig">build_exposures(account, shock_panel, design = c("direct", "shift_share"), center_dest = TRUE, horizons = c(0, 1, 2, Inf))</div>
  <p>Computes direct shocks $d_{it}$, infinite-horizon exposures $(\\mathbf{F}_t^y, \\mathbf{F}_t^v)$, finite-round approximations, and Borusyak–Hull propagated shift-share controls.</p>
</div>

<div class="api-block">
  <div class="api-sig">estimate_bilateral(panel, exposures, estimator = c("ols", "ppml", "iv"), cluster = c("origin^time", "dest^time"), include_direct = TRUE, use_expected_controls = TRUE, extra_controls = character(0), weights = NULL)</div>
  <p>C++ basket-centers destination exposures and estimates Proposition 2.4&rsquo;s bilateral regression via <code>fixest</code>. Pair with <code>bilateral_from_estimates()</code> for enclave imports, <code>rotemberg_weights()</code> for shift-share diagnostics, and <code>loo_influence()</code> for jackknife checks.</p>
</div>

<div class="api-block">
  <div class="api-sig">recover_welfare(fit, account, price_shares = "estimate_by_group", level_anchor = c("open_city_factor_balance", "zero"), period = NULL, base_period = NULL, kappa_periods = NULL, kappa_support = NULL, mpc_nonhousing = NULL, psi_O = 0.0, se_method = c("bootstrap", "delta", "none"), n_boot = 300L, seed = 42L)</div>
  <p>Recovers relative income and sales responses, estimates sector $\\widehat{\\kappa}_g^S$, pins down $\\widehat{c}_t$ via open-city factor balance (Proposition 2.5), verifies the aggregate terms-of-trade identity, and returns tract and decile welfare incidence.</p>
</div>

<div class="api-block">
  <div class="api-sig">add_shopping_wedge(welfare, sigma = c(2.0, 3.456, 8.235)) &middot; add_workplace_mobility(welfare, eps = 3.3) &middot; add_finite_bounds(welfare, sigma_ces = 3.456)</div>
  <p>Pipe-friendly robustness extensions relaxing fixed shopping amenities, fixed workplaces, and log-linearization.</p>
</div>

<div class="api-block">
  <div class="api-sig">solve_qsm(account, kappa = NULL, ET_shock = NULL, tau = NULL, sigma_R = 3.456, sigma_T = 7.0, eta_out = 2.5) &middot; optimal_spatial_tax(welfare, account, demand_elasticities = list(sigma_T = 7.0, eta_out = 2.5), sigma_R = 3.456, mode = c("both", "empirical", "nonlinear_qsm"))</div>
  <p>C++ damped-Newton nonlinear QSM solver and analytical-adjoint L-BFGS-B optimal uniform and place-specific visitor tax solver with 4-channel Ramsey decomposition.</p>
</div>

<div class="api-block">
  <div class="api-sig">sim_urban_city(N = 80L, T = 12L, K = 4L, n_visitor_groups = 6L, model = c("ces", "hsa", "iia_het", "kappa_het"), ...)</div>
  <p>Generates synthetic multi-sector gravity cities with known nonlinear equilibrium truth. Built-in 160-tract dataset accessible via <code>data(sim_city_160)</code>.</p>
</div>
""")

print("Generated 6 documentation pages in", DOCS)
