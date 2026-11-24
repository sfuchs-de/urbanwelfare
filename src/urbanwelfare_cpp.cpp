#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Dense LU solver with partial pivoting: solves A * X = B in place (B -> X)
// A is n x n (column-major), B is n x nrhs (column-major). Returns false if singular.
// ---------------------------------------------------------------------------
static bool dense_lu_solve(std::vector<double> A, std::vector<double>& B, int n, int nrhs) {
  std::vector<int> piv(n);
  for (int i = 0; i < n; ++i) piv[i] = i;

  for (int k = 0; k < n; ++k) {
    int p = k;
    double max_val = std::fabs(A[k + k * n]);
    for (int i = k + 1; i < n; ++i) {
      double v = std::fabs(A[i + k * n]);
      if (v > max_val) {
        max_val = v;
        p = i;
      }
    }
    if (max_val < 1e-14) return false;
    if (p != k) {
      std::swap(piv[k], piv[p]);
      for (int j = 0; j < n; ++j) {
        std::swap(A[k + j * n], A[p + j * n]);
      }
    }
    double akk = A[k + k * n];
    for (int i = k + 1; i < n; ++i) {
      A[i + k * n] /= akk;
      double lik = A[i + k * n];
      if (lik != 0.0) {
        for (int j = k + 1; j < n; ++j) {
          A[i + j * n] -= lik * A[k + j * n];
        }
      }
    }
  }

  // Apply row permutations and forward/back substitution for each RHS
  std::vector<double> tmp(n);
  for (int r = 0; r < nrhs; ++r) {
    double* bcol = &B[r * n];
    for (int i = 0; i < n; ++i) tmp[i] = bcol[piv[i]];
    for (int i = 0; i < n; ++i) {
      double s = tmp[i];
      for (int j = 0; j < i; ++j) s -= A[i + j * n] * tmp[j];
      tmp[i] = s;
    }
    for (int i = n - 1; i >= 0; --i) {
      double s = tmp[i];
      for (int j = i + 1; j < n; ++j) s -= A[i + j * n] * tmp[j];
      tmp[i] = s / A[i + i * n];
    }
    for (int i = 0; i < n; ++i) bcol[i] = tmp[i];
  }
  return true;
}

// [[Rcpp::export]]
NumericVector cpp_basket_center(NumericMatrix x_it,
                                NumericMatrix S_NH,
                                IntegerVector origin_idx,
                                IntegerVector dest_idx,
                                IntegerVector time_idx) {
  int N_dest = x_it.nrow();
  int N_time = x_it.ncol();
  int N_orig = S_NH.nrow();
  if (S_NH.ncol() != N_dest) {
    stop("S_NH ncol must equal x_it nrow (number of destinations).");
  }
  int M = origin_idx.size();
  if (dest_idx.size() != M || time_idx.size() != M) {
    stop("origin_idx, dest_idx, and time_idx must have equal length.");
  }

  // Precompute bar_x (N_orig x N_time): bar_x(n, t) = sum_j S_NH(n, j) * x_it(j, t)
  std::vector<double> bar_x(N_orig * N_time, 0.0);
  for (int t = 0; t < N_time; ++t) {
    for (int j = 0; j < N_dest; ++j) {
      double xjt = x_it(j, t);
      if (xjt == 0.0) continue;
      for (int n = 0; n < N_orig; ++n) {
        bar_x[n + t * N_orig] += S_NH(n, j) * xjt;
      }
    }
  }

  NumericVector out(M);
  for (int m = 0; m < M; ++m) {
    int o = origin_idx[m] - 1;
    int d = dest_idx[m] - 1;
    int t = time_idx[m] - 1;
    if (o < 0 || o >= N_orig || d < 0 || d >= N_dest || t < 0 || t >= N_time) {
      stop("Panel index out of bounds in cpp_basket_center.");
    }
    out[m] = x_it(d, t) - bar_x[o + t * N_orig];
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_bvls_box(NumericMatrix A,
                           NumericVector b,
                           NumericVector lo,
                           NumericVector hi,
                           int max_iter = 2000,
                           double tol = 1e-11) {
  // Fast coordinate-descent / Barzilai-Borwein projected gradient for
  // min_{lo <= x <= hi} 0.5 * ||A x - b||_2^2
  int m = A.nrow();
  int n = A.ncol();
  NumericVector x(n);
  for (int j = 0; j < n; ++j) {
    x[j] = std::min(hi[j], std::max(lo[j], 0.0));
  }

  // Precompute column norms squared and residual r = A x - b
  std::vector<double> col_norm2(n, 0.0);
  std::vector<double> r(m, 0.0);
  for (int i = 0; i < m; ++i) r[i] = -b[i];
  for (int j = 0; j < n; ++j) {
    double s2 = 0.0;
    double xj = x[j];
    for (int i = 0; i < m; ++i) {
      double aij = A(i, j);
      s2 += aij * aij;
      if (xj != 0.0) r[i] += aij * xj;
    }
    col_norm2[j] = s2;
  }

  for (int iter = 0; iter < max_iter; ++iter) {
    double max_diff = 0.0;
    for (int j = 0; j < n; ++j) {
      if (col_norm2[j] <= 1e-20 || hi[j] <= lo[j]) continue;
      double grad_j = 0.0;
      for (int i = 0; i < m; ++i) {
        grad_j += A(i, j) * r[i];
      }
      double x_old = x[j];
      double x_new = x_old - grad_j / col_norm2[j];
      if (x_new < lo[j]) x_new = lo[j];
      if (x_new > hi[j]) x_new = hi[j];
      double delta = x_new - x_old;
      if (delta != 0.0) {
        x[j] = x_new;
        for (int i = 0; i < m; ++i) {
          r[i] += A(i, j) * delta;
        }
        double ad = std::fabs(delta);
        if (ad > max_diff) max_diff = ad;
      }
    }
    if (max_diff < tol) break;
  }
  return x;
}

// ---------------------------------------------------------------------------
// Exact nonlinear open-city QSM evaluation + Newton solver + adjoint Ramsey gradient
// ---------------------------------------------------------------------------
struct QSMData {
  int N_orig;
  int N_dest;
  const double* y0;
  const double* mu;       // flexible input share 1 - theta_G
  const double* kappa;    // marginal price share kappa^S
  const double* A_fixed;  // fixed other demand A_i^0
  const double* ET_base;  // external/tourist spending at zero tax (or shocked ET_it)
  const double* v_out;    // fixed outside income by origin
  const double* nu;       // nonhousing share by origin
  const double* s_out0;   // initial conditional outside share s_{n0}^{NH}
  const double* S_NH0;    // N_orig x N_dest initial conditional shares s_{ni}^{NH}
  const double* H_alloc;  // N_orig x N_dest factor allocation matrix H_{ni}^m
  const double* pop;      // origin population
  const double* is_res;   // 1 if city resident receiving tax rebate & entering welfare, 0 otherwise
  const double* w_welf;   // welfare weights omega_n (sum to 1 over city residents)
  double sigma_R;         // resident CES elasticity across destinations
  double sigma_T;         // tourist CES elasticity across destinations
  double eta_out;         // outer tourist demand elasticity
};

static void qsm_supply(const QSMData& d, const std::vector<double>& logp,
                       std::vector<double>& y, std::vector<double>& G) {
  for (int i = 0; i < d.N_dest; ++i) {
    double lp = logp[i];
    double p = std::exp(lp);
    double mu_i = d.mu[i];
    double kap_i = d.kappa[i];
    double s_prod = (1.0 / kap_i - 1.0) * (1.0 - mu_i) / std::max(1e-8, mu_i);
    if (std::fabs(s_prod - 1.0) < 1e-6) {
      // Smooth Cobb-Douglas limit (s_prod -> 1, mu_i = 1 - kap_i)
      y[i] = d.y0[i] * std::exp(lp / kap_i);
      G[i] = (1.0 - mu_i) * y[i];
    } else {
      double sh_m = mu_i * std::exp((s_prod - 1.0) * lp);
      // Guard domain so 1 - sh_m stays strictly positive during trial Newton steps
      if (sh_m > 0.999999) sh_m = 0.999999;
      double q_ratio = std::pow((1.0 - mu_i) / (1.0 - sh_m), s_prod / (s_prod - 1.0));
      y[i] = d.y0[i] * p * q_ratio;
      G[i] = (1.0 - sh_m) * y[i];
    }
  }
}

static void qsm_eval(const QSMData& d,
                     const std::vector<double>& logp,
                     const std::vector<double>& theta,
                     bool use_tax_demand,
                     std::vector<double>& f_excess,
                     std::vector<double>& y,
                     std::vector<double>& G,
                     std::vector<double>& v,
                     std::vector<double>& logP_NH,
                     std::vector<double>& ET,
                     double& R_rev,
                     double& W_obj) {
  int No = d.N_orig;
  int Nd = d.N_dest;

  qsm_supply(d, logp, y, G);

  // Tourist spending & tax revenue
  R_rev = 0.0;
  if (!use_tax_demand) {
    for (int i = 0; i < Nd; ++i) ET[i] = d.ET_base[i];
  } else {
    double ET0_tot = 0.0;
    for (int i = 0; i < Nd; ++i) ET0_tot += d.ET_base[i];
    if (ET0_tot <= 0.0) {
      for (int i = 0; i < Nd; ++i) ET[i] = 0.0;
    } else {
      double denom_T = 0.0;
      std::vector<double> wT(Nd);
      for (int i = 0; i < Nd; ++i) {
        double alpha0 = d.ET_base[i] / ET0_tot;
        double lq = theta[i] + logp[i];
        wT[i] = alpha0 * std::exp((1.0 - d.sigma_T) * lq);
        denom_T += wT[i];
      }
      double logPT = std::log(denom_T) / (1.0 - d.sigma_T);
      double XT_tot = ET0_tot * std::exp(-d.eta_out * logPT);
      for (int i = 0; i < Nd; ++i) {
        double alpha_i = wT[i] / denom_T;
        double XT_i = alpha_i * XT_tot;
        ET[i] = std::exp(-theta[i]) * XT_i;
        R_rev += (1.0 - std::exp(-theta[i])) * XT_i;
      }
    }
  }

  double Omega_pop = 0.0;
  for (int n = 0; n < No; ++n) {
    if (d.is_res[n] > 0.5) Omega_pop += d.pop[n];
  }
  double T_reb = (Omega_pop > 0.0) ? (R_rev / Omega_pop) : 0.0;

  // Origin incomes v_n
  for (int n = 0; n < No; ++n) {
    double inc = d.v_out[n];
    for (int i = 0; i < Nd; ++i) {
      inc += d.H_alloc[n + i * No] * G[i];
    }
    if (d.is_res[n] > 0.5) {
      inc += d.pop[n] * T_reb;
    }
    v[n] = inc;
  }

  // Resident nonhousing CES shares & goods market clearing
  std::vector<double> exp_p(Nd);
  for (int i = 0; i < Nd; ++i) exp_p[i] = std::exp((1.0 - d.sigma_R) * logp[i]);

  std::vector<double> res_spend_dest(Nd, 0.0);
  W_obj = 0.0;
  for (int n = 0; n < No; ++n) {
    double denom = d.s_out0[n];
    for (int i = 0; i < Nd; ++i) {
      denom += d.S_NH0[n + i * No] * exp_p[i];
    }
    logP_NH[n] = std::log(denom) / (1.0 - d.sigma_R);
    double x_n = d.nu[n] * v[n];
    for (int i = 0; i < Nd; ++i) {
      double s_ni = (d.S_NH0[n + i * No] * exp_p[i]) / denom;
      res_spend_dest[i] += x_n * s_ni;
    }
    if (d.is_res[n] > 0.5 && d.w_welf[n] > 0.0) {
      // Baseline v0 is evaluated at logp=0, theta=0 outside; log(v_n) - nu_n * logP_NH[n] up to constant
      W_obj += d.w_welf[n] * (std::log(std::max(1e-12, v[n])) - d.nu[n] * logP_NH[n]);
    }
  }

  for (int i = 0; i < Nd; ++i) {
    f_excess[i] = (y[i] - res_spend_dest[i] - ET[i] - d.A_fixed[i]) / d.y0[i];
  }
}

static bool qsm_newton_solve(const QSMData& d,
                             const std::vector<double>& theta,
                             bool use_tax_demand,
                             std::vector<double>& logp,
                             int max_iter,
                             double tol,
                             int& iter_used,
                             double& max_resid) {
  int Nd = d.N_dest;
  int No = d.N_orig;
  std::vector<double> f(Nd), f_up(Nd), f_dn(Nd), y(Nd), G(Nd), v(No), logP(No), ET(Nd);
  double R_rev = 0.0, W_obj = 0.0;
  double eps_fd = 1e-6;

  qsm_eval(d, logp, theta, use_tax_demand, f, y, G, v, logP, ET, R_rev, W_obj);
  max_resid = 0.0;
  for (int i = 0; i < Nd; ++i) max_resid = std::max(max_resid, std::fabs(f[i]));
  iter_used = 0;
  if (max_resid < tol) return true;

  std::vector<double> J(Nd * Nd);
  std::vector<double> lp_pert = logp;

  for (int it = 1; it <= max_iter; ++it) {
    iter_used = it;
    for (int j = 0; j < Nd; ++j) {
      double orig = logp[j];
      lp_pert[j] = orig + eps_fd;
      qsm_eval(d, lp_pert, theta, use_tax_demand, f_up, y, G, v, logP, ET, R_rev, W_obj);
      lp_pert[j] = orig - eps_fd;
      qsm_eval(d, lp_pert, theta, use_tax_demand, f_dn, y, G, v, logP, ET, R_rev, W_obj);
      lp_pert[j] = orig;
      for (int i = 0; i < Nd; ++i) {
        J[i + j * Nd] = (f_up[i] - f_dn[i]) / (2.0 * eps_fd);
      }
    }

    std::vector<double> step(Nd);
    for (int i = 0; i < Nd; ++i) step[i] = -f[i];
    if (!dense_lu_solve(J, step, Nd, 1)) return false;

    // Backtracking line search on merit 0.5 * ||f||_2^2
    double norm0 = 0.0;
    for (int i = 0; i < Nd; ++i) norm0 += f[i] * f[i];
    double alpha = 1.0;
    std::vector<double> lp_cand(Nd);
    bool accepted = false;
    for (int ls = 0; ls < 12; ++ls) {
      for (int i = 0; i < Nd; ++i) {
        lp_cand[i] = std::min(3.5, std::max(-3.5, logp[i] + alpha * step[i]));
      }
      qsm_eval(d, lp_cand, theta, use_tax_demand, f_up, y, G, v, logP, ET, R_rev, W_obj);
      double norm1 = 0.0;
      for (int i = 0; i < Nd; ++i) norm1 += f_up[i] * f_up[i];
      if (std::isfinite(norm1) && norm1 < (1.0 - 1e-4 * alpha) * norm0) {
        logp = lp_cand;
        f = f_up;
        accepted = true;
        break;
      }
      alpha *= 0.5;
    }
    if (!accepted) {
      for (int i = 0; i < Nd; ++i) logp[i] += 0.1 * step[i];
      qsm_eval(d, logp, theta, use_tax_demand, f, y, G, v, logP, ET, R_rev, W_obj);
    }

    max_resid = 0.0;
    for (int i = 0; i < Nd; ++i) max_resid = std::max(max_resid, std::fabs(f[i]));
    if (max_resid < tol) return true;
  }
  return max_resid < 1e-8;
}

// [[Rcpp::export]]
List cpp_solve_qsm(NumericVector y0,
                   NumericVector mu,
                   NumericVector kappa,
                   NumericVector A_fixed,
                   NumericVector ET_shock,
                   NumericVector v_out,
                   NumericVector nu,
                   NumericVector s_out0,
                   NumericMatrix S_NH0,
                   NumericMatrix H_alloc,
                   NumericVector pop,
                   NumericVector is_res,
                   NumericVector w_welf,
                   NumericVector theta,
                   bool use_tax_demand,
                   double sigma_R,
                   double sigma_T,
                   double eta_out,
                   NumericVector logp_init,
                   bool compute_jacobian = false,
                   int max_iter = 40,
                   double tol = 1e-11) {
  int No = S_NH0.nrow();
  int Nd = S_NH0.ncol();
  QSMData d{No, Nd, y0.begin(), mu.begin(), kappa.begin(), A_fixed.begin(),
            ET_shock.begin(), v_out.begin(), nu.begin(), s_out0.begin(),
            S_NH0.begin(), H_alloc.begin(), pop.begin(), is_res.begin(),
            w_welf.begin(), sigma_R, sigma_T, eta_out};

  std::vector<double> th(theta.begin(), theta.end());
  std::vector<double> lp(logp_init.begin(), logp_init.end());

  // Baseline welfare offset at logp=0, theta=0
  std::vector<double> zero_Nd(Nd, 0.0), f0(Nd), y_b(Nd), G_b(Nd), v_b(No), logP_b(No), ET_b(Nd);
  double R0 = 0.0, W0 = 0.0;
  qsm_eval(d, zero_Nd, zero_Nd, use_tax_demand, f0, y_b, G_b, v_b, logP_b, ET_b, R0, W0);

  int iter_used = 0;
  double max_resid = 0.0;
  bool ok = qsm_newton_solve(d, th, use_tax_demand, lp, max_iter, tol, iter_used, max_resid);
  if (!ok) {
    // Homotopy continuation from zero
    lp.assign(Nd, 0.0);
    std::vector<double> th_step(Nd, 0.0);
    std::vector<double> et_orig(ET_shock.begin(), ET_shock.end());
    for (int s = 1; s <= 4; ++s) {
      double lam = 0.25 * s;
      for (int i = 0; i < Nd; ++i) th_step[i] = lam * th[i];
      ok = qsm_newton_solve(d, th_step, use_tax_demand, lp, max_iter, tol, iter_used, max_resid);
    }
  }

  std::vector<double> f(Nd), y(Nd), G(Nd), v(No), logP(No), ET(Nd);
  double R_rev = 0.0, W_obj = 0.0;
  qsm_eval(d, lp, th, use_tax_demand, f, y, G, v, logP, ET, R_rev, W_obj);

  NumericVector out_logp(lp.begin(), lp.end());
  NumericVector out_y(y.begin(), y.end());
  NumericVector out_G(G.begin(), G.end());
  NumericVector out_v(v.begin(), v.end());
  NumericVector out_v0(v_b.begin(), v_b.end());
  NumericVector out_logP(logP.begin(), logP.end());
  NumericVector out_ET(ET.begin(), ET.end());
  NumericVector out_welf(No);
  for (int n = 0; n < No; ++n) {
    out_welf[n] = std::log(v[n] / v_b[n]) - nu[n] * logP[n];
  }

  // Optional adjoint gradient dW/dtheta and GE tax-base elasticity matrix E = -d log ET / d theta'
  NumericVector grad_theta(Nd, 0.0);
  NumericMatrix E_ge(Nd, Nd);
  NumericVector b_T(Nd, 0.0);
  double mu_T_ge = 0.0;

  if (use_tax_demand) {
    double eps_fd = 1e-6;
    std::vector<double> J_p(Nd * Nd), f_u(Nd), f_d(Nd), y_t(Nd), G_t(Nd), v_t(No), lp_t(No), ET_u(Nd), ET_d(Nd);
    std::vector<double> g_p(Nd, 0.0);
    std::vector<double> dET_dp(Nd * Nd, 0.0); // d log ET_i / d log p_j
    double R_u = 0.0, R_d = 0.0, W_u = 0.0, W_d = 0.0;
    std::vector<double> lp_pert = lp;

    for (int j = 0; j < Nd; ++j) {
      double orig = lp[j];
      lp_pert[j] = orig + eps_fd;
      qsm_eval(d, lp_pert, th, true, f_u, y_t, G_t, v_t, lp_t, ET_u, R_u, W_u);
      lp_pert[j] = orig - eps_fd;
      qsm_eval(d, lp_pert, th, true, f_d, y_t, G_t, v_t, lp_t, ET_d, R_d, W_d);
      lp_pert[j] = orig;
      g_p[j] = (W_u - W_d) / (2.0 * eps_fd);
      for (int i = 0; i < Nd; ++i) {
        J_p[i + j * Nd] = (f_u[i] - f_d[i]) / (2.0 * eps_fd);
        if (ET[i] > 0.0) {
          dET_dp[i + j * Nd] = ( std::log(std::max(1e-18, ET_u[i])) - std::log(std::max(1e-18, ET_d[i])) ) / (2.0 * eps_fd);
        }
      }
    }

    // Adjoint solve: J_p' * lambda = g_p
    std::vector<double> J_pT(Nd * Nd);
    for (int i = 0; i < Nd; ++i) {
      for (int j = 0; j < Nd; ++j) {
        J_pT[i + j * Nd] = J_p[j + i * Nd];
      }
    }
    std::vector<double> lam = g_p;
    dense_lu_solve(J_pT, lam, Nd, 1);

    // Partial derivatives wrt theta_j (holding logp fixed)
    std::vector<double> th_pert = th;
    std::vector<double> J_th(Nd * Nd, 0.0);
    std::vector<double> dET_dth_partial(Nd * Nd, 0.0);
    for (int j = 0; j < Nd; ++j) {
      double orig = th[j];
      th_pert[j] = orig + eps_fd;
      qsm_eval(d, lp, th_pert, true, f_u, y_t, G_t, v_t, lp_t, ET_u, R_u, W_u);
      th_pert[j] = orig - eps_fd;
      qsm_eval(d, lp, th_pert, true, f_d, y_t, G_t, v_t, lp_t, ET_d, R_d, W_d);
      th_pert[j] = orig;

      double g_th_j = (W_u - W_d) / (2.0 * eps_fd);
      double dot_lam = 0.0;
      for (int i = 0; i < Nd; ++i) {
        double jf = (f_u[i] - f_d[i]) / (2.0 * eps_fd);
        J_th[i + j * Nd] = jf;
        dot_lam += jf * lam[i];
        if (ET[i] > 0.0) {
          dET_dth_partial[i + j * Nd] = (std::log(std::max(1e-18, ET_u[i])) - std::log(std::max(1e-18, ET_d[i]))) / (2.0 * eps_fd);
        }
      }
      grad_theta[j] = g_th_j - dot_lam;
    }

    if (compute_jacobian) {
      // Solve J_p * (d logp / d theta') = -J_th
      std::vector<double> dp_dth(Nd * Nd);
      for (int k = 0; k < Nd * Nd; ++k) dp_dth[k] = -J_th[k];
      dense_lu_solve(J_p, dp_dth, Nd, Nd);
      for (int j = 0; j < Nd; ++j) {
        for (int i = 0; i < Nd; ++i) {
          double tot = dET_dth_partial[i + j * Nd];
          for (int k = 0; k < Nd; ++k) {
            tot += dET_dp[i + k * Nd] * dp_dth[k + j * Nd];
          }
          E_ge(i, j) = -tot;
        }
      }

      // Also compute conditional welfare derivatives b_T = dW/d log E^T |_T and mu_T^GE = dW/dT |_{E^T}
      // In conditional system F^c(logp, e, T):
      // df_i / d e_j = -ET_j / y0_i * 1_{i=j}; dW/de_j|_p = 0 => b_T_j = - (-ET_j / y0_j) * lam_c_j
      // Let's evaluate conditional Jacobian J_p^c (holding ET and T fixed!)
      std::vector<double> J_pc(Nd * Nd), g_pc(Nd, 0.0);
      QSMData dc = d;
      dc.ET_base = ET.data(); // fixed producer-price tourist spending ET at solution
      double Omega_pop = 0.0;
      for (int n = 0; n < No; ++n) if (d.is_res[n] > 0.5) Omega_pop += d.pop[n];
      double T_reb = (Omega_pop > 0.0) ? (R_rev / Omega_pop) : 0.0;
      std::vector<double> v_out_fixedT(No);
      for (int n = 0; n < No; ++n) {
        v_out_fixedT[n] = d.v_out[n] + ((d.is_res[n] > 0.5) ? d.pop[n] * T_reb : 0.0);
      }
      dc.v_out = v_out_fixedT.data();

      for (int j = 0; j < Nd; ++j) {
        double orig = lp[j];
        lp_pert[j] = orig + eps_fd;
        qsm_eval(dc, lp_pert, zero_Nd, false, f_u, y_t, G_t, v_t, lp_t, ET_u, R_u, W_u);
        lp_pert[j] = orig - eps_fd;
        qsm_eval(dc, lp_pert, zero_Nd, false, f_d, y_t, G_t, v_t, lp_t, ET_d, R_d, W_d);
        lp_pert[j] = orig;
        g_pc[j] = (W_u - W_d) / (2.0 * eps_fd);
        for (int i = 0; i < Nd; ++i) {
          J_pc[i + j * Nd] = (f_u[i] - f_d[i]) / (2.0 * eps_fd);
        }
      }
      std::vector<double> J_pcT(Nd * Nd);
      for (int i = 0; i < Nd; ++i) for (int j = 0; j < Nd; ++j) J_pcT[i + j * Nd] = J_pc[j + i * Nd];
      std::vector<double> lam_c = g_pc;
      dense_lu_solve(J_pcT, lam_c, Nd, 1);

      for (int j = 0; j < Nd; ++j) {
        // df^c_j / d log E^T_j = -ET_j / y0_j, partial dW/d log E^T_j = 0
        b_T[j] = (ET[j] / y0[j]) * lam_c[j];
      }

      // Perturb T_reb holding ET fixed for mu_T^GE
      double eps_T = 1e-6 * std::max(1.0, std::fabs(T_reb));
      std::vector<double> v_out_Tp = v_out_fixedT, v_out_Tm = v_out_fixedT;
      for (int n = 0; n < No; ++n) {
        if (d.is_res[n] > 0.5) {
          v_out_Tp[n] += d.pop[n] * eps_T;
          v_out_Tm[n] -= d.pop[n] * eps_T;
        }
      }
      dc.v_out = v_out_Tp.data();
      qsm_eval(dc, lp, zero_Nd, false, f_u, y_t, G_t, v_t, lp_t, ET_u, R_u, W_u);
      dc.v_out = v_out_Tm.data();
      qsm_eval(dc, lp, zero_Nd, false, f_d, y_t, G_t, v_t, lp_t, ET_d, R_d, W_d);
      double dW_dT_partial = (W_u - W_d) / (2.0 * eps_T);
      double dot_T = 0.0;
      for (int i = 0; i < Nd; ++i) {
        double df_dT = (f_u[i] - f_d[i]) / (2.0 * eps_T);
        dot_T += df_dT * lam_c[i];
      }
      mu_T_ge = dW_dT_partial - dot_T;
    }
  }

  return List::create(
    _["converged"] = ok,
    _["iterations"] = iter_used,
    _["max_residual"] = max_resid,
    _["logp"] = out_logp,
    _["y"] = out_y,
    _["G"] = out_G,
    _["v"] = out_v,
    _["v0"] = out_v0,
    _["logP_NH"] = out_logP,
    _["ET"] = out_ET,
    _["revenue"] = R_rev,
    _["objective"] = W_obj - W0,
    _["grad_theta"] = grad_theta,
    _["E_ge"] = E_ge,
    _["b_T"] = b_T,
    _["mu_T_ge"] = mu_T_ge
  );
}
