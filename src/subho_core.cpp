// Compiled core of GGMsubho.
//
// Implements the two computational steps of SUBHO (Step-Up Bayesian HOrseshoe):
// the horseshoe Gibbs sampler for one node-wise regression and the evaluation of
// the working-model posterior along the ranked, nested model path.
//
// Method: Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional
// Gaussian Graphical Model Estimation. arXiv:2308.02713.
// https://arxiv.org/abs/2308.02713
//
// Sampler: Makalic and Schmidt (2016), IEEE Signal Processing Letters 23(1),
// with the beta update of Bhattacharya, Chakraborty and Mallick (2016), Biometrika,
// when the number of columns is not small relative to the number of rows.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

namespace {

// 1 / Gamma(shape, rate = rate)
inline double r_inv_gamma(double shape, double rate) {
  return 1.0 / R::rgamma(shape, 1.0 / rate);
}

inline arma::vec r_std_normal(arma::uword n) {
  Rcpp::NumericVector z = Rcpp::rnorm(n);
  return arma::vec(z.begin(), n, true, false);
}

// Keeps the scale parameters inside the range where 1/x and x are both finite.
inline double clamp_scale(double x) {
  const double lo = 1e-150, hi = 1e150;
  return x < lo ? lo : (x > hi ? hi : x);
}

} // namespace

// Horseshoe Gibbs sampler
//
//   y | beta, sigma2 ~ N(X beta, sigma2 I)
//   beta_j | lambda_j, tau, sigma2 ~ N(0, lambda_j^2 tau^2 sigma2)
//   lambda_j, tau ~ C+(0, 1), through the inverse-gamma augmentation of
//   Makalic and Schmidt; sigma2 ~ IG(a_sigma, b_sigma).
//
// Returns draws saved at iterations burn + thin, burn + 2 thin, ...
//
// [[Rcpp::export]]
Rcpp::List hs_gibbs_cpp(const arma::mat& X, const arma::vec& y, int n_iter,
                        int burn, int thin, double a_sigma, double b_sigma,
                        double sigma2_init, double tau2_init, bool fix_tau) {
  const arma::uword n = X.n_rows, p = X.n_cols;
  const int n_save = (n_iter - burn) / thin;
  const bool use_primal = static_cast<double>(p) < 0.8 * static_cast<double>(n);

  arma::mat beta_out(p, n_save, arma::fill::zeros);
  arma::vec sigma2_out(n_save, arma::fill::zeros), tau2_out(n_save, arma::fill::zeros);

  arma::vec beta(p, arma::fill::ones), lambda2(p, arma::fill::ones), nu(p, arma::fill::ones);
  arma::vec lambda_star(p);
  double sigma2 = sigma2_init, tau2 = tau2_init, xi = 1.0;
  lambda_star = tau2 * lambda2;

  arma::mat XtX;
  arma::vec Xty;
  if (use_primal) {
    XtX = X.t() * X;
    Xty = X.t() * y;
  }
  arma::mat A, U, XD;
  arma::vec resid;

  for (int iter = 1; iter <= n_iter; ++iter) {
    if (iter % 100 == 0) Rcpp::checkUserInterrupt();

    // beta | rest
    if (use_primal) {
      // A = X'X + diag(1 / lambda_star);  beta ~ N(A^{-1} X'y, sigma2 A^{-1})
      A = XtX;
      A.diag() += 1.0 / lambda_star;
      if (!arma::chol(U, A, "upper"))
        Rcpp::stop("Cholesky factorisation failed in the horseshoe beta update.");
      arma::vec z = arma::solve(arma::trimatl(U.t()), Xty, arma::solve_opts::fast);
      z += std::sqrt(sigma2) * r_std_normal(p);
      beta = arma::solve(arma::trimatu(U), z, arma::solve_opts::fast);
    } else {
      // Bhattacharya et al.: Phi = X / sigma, D = sigma2 diag(lambda_star)
      const double sigma = std::sqrt(sigma2);
      arma::vec d = sigma2 * lambda_star;
      arma::vec u = arma::sqrt(d) % r_std_normal(p);
      arma::vec v = (X * u) / sigma + r_std_normal(n);
      // Phi D Phi' + I = X diag(lambda_star) X' + I
      XD = X;
      XD.each_row() %= arma::sqrt(lambda_star).t();
      A = XD * XD.t();
      A.diag() += 1.0;
      if (!arma::chol(U, A, "upper"))
        Rcpp::stop("Cholesky factorisation failed in the horseshoe beta update.");
      arma::vec w = arma::solve(arma::trimatl(U.t()), y / sigma - v, arma::solve_opts::fast);
      w = arma::solve(arma::trimatu(U), w, arma::solve_opts::fast);
      beta = u + d % (X.t() * w) / sigma;
    }

    // xi | tau2
    if (!fix_tau) xi = r_inv_gamma(1.0, 1.0 + 1.0 / tau2);

    // nu_j | lambda_j, then lambda_j | nu_j, beta_j, tau2, sigma2
    for (arma::uword j = 0; j < p; ++j) {
      nu[j] = r_inv_gamma(1.0, 1.0 + 1.0 / lambda2[j]);
      lambda2[j] = clamp_scale(
        r_inv_gamma(1.0, 1.0 / nu[j] + beta[j] * beta[j] / (2.0 * sigma2 * tau2)));
    }

    // tau2 | rest
    if (!fix_tau) {
      double s = 1.0 / xi + arma::accu(arma::square(beta) / lambda2) / (2.0 * sigma2);
      tau2 = clamp_scale(r_inv_gamma((static_cast<double>(p) + 1.0) / 2.0, s));
    }
    lambda_star = lambda2 * tau2;
    lambda_star.transform([](double x) { return clamp_scale(x); });

    // sigma2 | rest
    resid = y - X * beta;
    double rate = 0.5 * (arma::dot(resid, resid) +
                         arma::accu(arma::square(beta) / lambda_star)) + b_sigma;
    sigma2 = clamp_scale(
      r_inv_gamma(0.5 * (static_cast<double>(n) + static_cast<double>(p)) + a_sigma, rate));

    if (iter > burn && (iter - burn) % thin == 0) {
      const int k = (iter - burn) / thin - 1;
      if (k < n_save) {
        beta_out.col(k) = beta;
        sigma2_out[k] = sigma2;
        tau2_out[k] = tau2;
      }
    }
  }

  return Rcpp::List::create(Rcpp::Named("beta") = beta_out,
                            Rcpp::Named("sigma2") = sigma2_out,
                            Rcpp::Named("tau2") = tau2_out);
}

// Working-model log posterior along a nested path.
//
// `gram` is the cross-product matrix of a design whose columns include those in
// `path` (0-based); path[0] is the intercept, path[1..K] the ranked predictors.
// `xty` is the cross-product of the same design with the response.
//
// For each k the working prior is beta | sigma2 ~ N(0, sigma2 v I) with
// sigma2 ~ IG(c0, d0). The Cholesky factor is extended one column at a time, so the
// whole path costs O(K^3), and everything is computed on the log scale.
//
// [[Rcpp::export]]
Rcpp::List step_up_path_cpp(const arma::mat& gram, const arma::vec& xty, double yty,
                            int n, const arma::uvec& path, double c0, double d0,
                            double v, const arma::vec& log_prior) {
  const arma::uword m = path.n_elem;
  arma::mat L(m, m, arma::fill::zeros);
  arma::vec w(m, arma::fill::zeros);
  arma::vec log_marg(m);
  log_marg.fill(R_NegInf);

  const double half_n = 0.5 * static_cast<double>(n);
  const double cst = c0 * std::log(d0) + R::lgammafn(c0 + half_n) - R::lgammafn(c0) -
    half_n * std::log(2.0 * M_PI);
  const double log_v = std::log(v);
  double log_det = 0.0, quad = 0.0;
  arma::uword used = 0;

  for (arma::uword j = 0; j < m; ++j) {
    const arma::uword cj = path[j];
    double pivot = gram(cj, cj) + 1.0 / v, wj = xty[cj];
    for (arma::uword i = 0; i < j; ++i) {
      double s = gram(path[i], cj);
      for (arma::uword l = 0; l < i; ++l) s -= L(i, l) * L(j, l);
      s /= L(i, i);
      L(j, i) = s;
      pivot -= s * s;
      wj -= s * w[i];
    }
    if (pivot <= 0.0) break;  // numerically singular: longer models are not evaluated
    const double ljj = std::sqrt(pivot);
    L(j, j) = ljj;
    w[j] = wj / ljj;
    log_det += 2.0 * std::log(ljj);
    quad += w[j] * w[j];
    const double dstar = d0 + 0.5 * std::max(yty - quad, 0.0);
    log_marg[j] = cst - 0.5 * log_det - 0.5 * static_cast<double>(j + 1) * log_v -
      (c0 + half_n) * std::log(dstar);
    used = j + 1;
  }
  if (used == 0) Rcpp::stop("The intercept-only working model could not be evaluated.");

  arma::vec log_post = log_marg + log_prior.head(m);
  const arma::uword best = log_post.index_max();

  // Posterior mean of the coefficients under the selected working model
  arma::vec coef = arma::solve(arma::trimatu(L.submat(0, 0, best, best).t()),
                               w.head(best + 1), arma::solve_opts::fast);

  return Rcpp::List::create(Rcpp::Named("log_marginal") = log_marg,
                            Rcpp::Named("log_posterior") = log_post,
                            Rcpp::Named("k") = static_cast<int>(best),
                            Rcpp::Named("coef") = coef);
}
