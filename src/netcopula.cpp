#include "netcopula.h"

// TO DO:
//   1. complete help pages for all commands
//   2. check if speed can be improved

// [[Rcpp::depends("RcppArmadillo")]]

//' Log-likelihood of copula based model for a multivariate NMA.
//'
//' Evaluation of the log-likelihood.
//'
//' @param y pippo
//' @param n pippo
//' @param x pippo
//' @param trt pippo
//' @param mu pippo
//' @param delta pippo
//' @param Gamma pippo
//'
//' @return A length-one numeric vector.
//' @export
//'
//' @examples
//' # nothing for now!
// [[Rcpp::export]]
double nc_loglik(const Rcpp::NumericMatrix& y, const Rcpp::NumericMatrix& n, const Rcpp::NumericMatrix& x, const Rcpp::IntegerVector& trt, const Rcpp::NumericMatrix& mu, const Rcpp::NumericMatrix& delta, const Rcpp::List& Gamma) {
  int n_datapoints = y.nrow(), M = y.ncol(), nGamma = Gamma.size();

  arma::mat Gamma_k(M, M);
  Rcpp::NumericVector x_ik(M);
  double llik = 0;
  int ind_a_b = 0;
  for (unsigned int ik = 0; ik < n_datapoints; ik++) {
    if (nGamma > 1) {
      Gamma_k = Rcpp::as<arma::mat>(Gamma(trt(ik) - 1));
    } else {
      Gamma_k = Rcpp::as<arma::mat>(Gamma(0));
    }
    for (unsigned int m = 0; m < M; m++) {
      x_ik(m) = x(ik, m);
      ind_a_b = indic_a_b(y(ik, m), n(ik, m), x(ik, m), mu(ik, m), delta(ik, m));
      if (!ind_a_b) {
      //   return NA_REAL;
      }
    }
    llik += gausscopdens(x_ik, Gamma_k, false, true);
  }

  return llik;
}

//' @export
// [[Rcpp::export]]
int indic_a_b(const double& y_ikm, const int& n_ikm, const double& x_ikm, const double& mu_ikm, const double& delta_ikm) {
  int out = 0;
  double theta = mu_ikm + delta_ikm, pi = expit_double(theta);

  double phi_inv_a = 0.0, phi_inv_b = 0.0;
  phi_inv_a = R::qnorm(R::pbinom(y_ikm - 1, n_ikm, pi, 1, 1), 0.0, 1.0, 1, 1);
  phi_inv_b = R::qnorm(R::pbinom(y_ikm, n_ikm, pi, 1, 1), 0.0, 1.0, 1, 1);
  if ((phi_inv_a < x_ikm) && (x_ikm <= phi_inv_b)) {
    out = 1;
  }

  return out;
}

//' @export
// [[Rcpp::export]]
double nc_logprior(const Rcpp::NumericMatrix& mu, const double& mu_sigma, const arma::mat& d, const double& d_sigma, const arma::mat& Sigma_M, const double& sigma_r, const int& ref_trt) {
  int M = mu.ncol(), n_trt = d.n_rows, count = 0;

  arma::uvec row_id(n_trt - 1);
  for (unsigned int k = 0; k < n_trt; k++) {
    if (k != (ref_trt - 1)) {
      row_id(count) = k;
      count++;
    }
  }

  double logprior = arma::sum(dmvn_arma(Rcpp::as<arma::mat>(mu), arma::zeros<arma::vec>(M), pow(mu_sigma, 2)*arma::eye(M, M), true));
  arma::mat d_tmp(n_trt - 1, M);
  d_tmp = d.rows(row_id);
  logprior += arma::sum(dmvn_arma(d_tmp, arma::zeros<arma::vec>(M), pow(d_sigma, 2)*arma::eye(M, M), true));
  logprior += dlogchol_arma(Sigma_M, sigma_r, true);

  return logprior;
}

//' Log-likelihood of copula based model for a multivariate NMA.
//'
//' Evaluation of the log-likelihood.
//'
//' @param data pippo
//' @param init pippo
//' @param totiter pippo
//' @param prior pippo
//' @param prop pippo
//' @param verbose pippo
//'
//' @return A list containing the output of the MCMC simulation.
//' @export
//'
//' @examples
//' # nothing for now!
// [[Rcpp::export]]
Rcpp::NumericMatrix x_imputed(const Rcpp::NumericMatrix& x, const Rcpp::List& Gamma, const Rcpp::IntegerVector& trt) {
  int n = x.nrow(), M = x.ncol(), nGamma = Gamma.size(), n_q, n_q_vec, n_q_obs, n_q_miss;
  Rcpp::List s2t_split = split_iv(Rcpp::seq_len(n), trt);
  Rcpp::List x2t_split = split_nm(x, trt);
  Rcpp::IntegerVector split_s;
  Rcpp::NumericMatrix x_imp = Rcpp::clone(x), x_q, x_q_imp;
  Rcpp::NumericVector x_q_vec, x_q_vec_obs, x_q_vec_imp;
  arma::uvec x_q_vec_miss_i, x_q_vec_obs_i;
  arma::vec x_q_vec_miss, mu;
  arma::mat Gamma_q_tilde, Gamma_q_miss, Gamma_q_obs, Gamma_q_12, Gamma_q_21, Gamma_q_obs_inv, Omega;
  if (nGamma > 1) {
    for (int q = 0; q < nGamma; q++) {
      split_s = s2t_split[q];
      x_q = Rcpp::as<Rcpp::NumericMatrix>(x2t_split[q]);
      n_q = x_q.nrow();
      x_q_vec = nm_stack(x_q);
      x_q_vec_obs = nv_omit(x_q_vec);
      n_q_vec = x_q_vec.size();
      n_q_obs = x_q_vec_obs.size();
      n_q_miss = n_q_vec - n_q_obs;
      x_q_vec_miss_i = nv_na_index(x_q_vec, n_q_miss, true);
      x_q_vec_obs_i = nv_na_index(x_q_vec, n_q_obs, false);
      Gamma_q_tilde = mat_block_diag(Rcpp::as<arma::mat>(Gamma[q]), n_q);
      Gamma_q_miss = Gamma_q_tilde.submat(x_q_vec_miss_i, x_q_vec_miss_i);
      Gamma_q_obs = Gamma_q_tilde.submat(x_q_vec_obs_i, x_q_vec_obs_i);
      Gamma_q_21 = Gamma_q_tilde.submat(x_q_vec_obs_i, x_q_vec_miss_i);
      Gamma_q_12 = trans(Gamma_q_21);
      Gamma_q_obs_inv = arma::inv(Gamma_q_obs);
      mu = Gamma_q_12*Gamma_q_obs_inv*Rcpp::as<arma::vec>(x_q_vec_obs);
      Omega = Gamma_q_miss - Gamma_q_12*Gamma_q_obs_inv*Gamma_q_21;
      if (!is_positive_definite(Omega, 116)) {
        Rprintf("x_imputed\n");
        Omega = make_positive_definite(Omega);
      }
      x_q_vec_miss = arma::vectorise(rmvn_arma(1, mu, Omega));
      x_q_vec_imp = nv_miss_replace(x_q_vec, x_q_vec_miss, x_q_vec_miss_i);
      x_q_imp = nv_unstack(x_q_vec_imp, M);
      for (unsigned int k = 0; k < n_q; k++) {
        x_imp(split_s[k] - 1, Rcpp::_) = x_q_imp(k, Rcpp::_);
      }
    }
  } else {
    x_q_vec = nm_stack(x);
    x_q_vec_obs = nv_omit(x_q_vec);
    n_q_vec = x_q_vec.size();
    n_q_obs = x_q_vec_obs.size();
    n_q_miss = n_q_vec - n_q_obs;
    x_q_vec_miss_i = nv_na_index(x_q_vec, n_q_miss, true);
    x_q_vec_obs_i = nv_na_index(x_q_vec, n_q_obs, false);
    Gamma_q_tilde = mat_block_diag(Rcpp::as<arma::mat>(Gamma[0]), n);
    Gamma_q_miss = Gamma_q_tilde.submat(x_q_vec_miss_i, x_q_vec_miss_i);
    Gamma_q_obs = Gamma_q_tilde.submat(x_q_vec_obs_i, x_q_vec_obs_i);
    Gamma_q_21 = Gamma_q_tilde.submat(x_q_vec_obs_i, x_q_vec_miss_i);
    Gamma_q_12 = trans(Gamma_q_21);
    Gamma_q_obs_inv = arma::inv(Gamma_q_obs);
    mu = Gamma_q_12*Gamma_q_obs_inv*Rcpp::as<arma::vec>(x_q_vec_obs);
    Omega = Gamma_q_miss - Gamma_q_12*Gamma_q_obs_inv*Gamma_q_21;
    if (!is_positive_definite(Omega, 143)) {
      Rprintf("x_imputed\n");
      Omega = make_positive_definite(Omega);
    }
    x_q_vec_miss = arma::vectorise(rmvn_arma(1, mu, Omega));
    x_q_vec_imp = nv_miss_replace(x_q_vec, x_q_vec_miss, x_q_vec_miss_i);
    x_q_imp = nv_unstack(x_q_vec_imp, M);
    for (unsigned int k = 0; k < n; k++) {
      x_imp(k, Rcpp::_) = x_q_imp(k, Rcpp::_);
    }
  }

  return x_imp;
}

//' Log-likelihood of copula based model for a multivariate NMA.
//'
//' Evaluation of the log-likelihood.
//'
//' @param data pippo
//' @param init pippo
//' @param totiter pippo
//' @param prior pippo
//' @param prop pippo
//' @param verbose pippo
//'
//' @return A list containing the output of the MCMC simulation.
//' @export
//'
//' @examples
//' # nothing for now!
// [[Rcpp::export]]
Rcpp::NumericMatrix n_imputed(const Rcpp::NumericMatrix& n_data) {
  int n_datapoints = n_data.nrow(), M = n_data.ncol(), n_miss;
  Rcpp::NumericMatrix n_imp = Rcpp::clone(n_data);
  for (unsigned int i = 0; i < n_datapoints; i++) {
    n_miss = arma::as_scalar(arma::min(nm_omit(n_data, i))); // use either min, max or mean
    for (unsigned int j = 0; j < M; j++) {
      if (ISNAN(n_imp(i, j))) {
        n_imp(i, j) = n_miss;
      }
    }
  }

  return n_imp;
}

//' Log-likelihood of copula based model for a multivariate NMA.
//'
//' Evaluation of the log-likelihood.
//'
//' @param data pippo
//' @param init pippo
//' @param totiter pippo
//' @param prior pippo
//' @param prop pippo
//' @param verbose pippo
//'
//' @return A list containing the output of the MCMC simulation.
//' @export
//'
//' @examples
//' # nothing for now!
// [[Rcpp::export]]
Rcpp::NumericMatrix y_imputed(const Rcpp::NumericMatrix& y, const Rcpp::NumericMatrix& x_imp, const Rcpp::IntegerVector& narms, const Rcpp::NumericMatrix& mu, const Rcpp::NumericMatrix& delta, const Rcpp::NumericMatrix& n_imp) {
  int n_datapoints = y.nrow(), M = y.ncol();
  Rcpp::NumericMatrix y_imp = Rcpp::clone(y), mu_long = param_long(mu, narms, false);
  double theta, pi;
  for (unsigned int i = 0; i < n_datapoints; i++) {
    for (unsigned int j = 0; j < M; j++) {
      if (ISNAN(y_imp(i, j))) {
        theta = mu_long(i, j) + delta(i, j);
        pi = expit_double(theta);
        y_imp(i, j) = R::qbinom(R::pnorm(x_imp(i, j), 0.0, 1.0, 1, 0), n_imp(i, j), pi, 1, 0);
      }
    }
  }

  return y_imp;
}

//' Log-likelihood of copula based model for a multivariate NMA.
//'
//' Evaluation of the log-likelihood.
//'
//' @param nc_data pippo
//' @param x pippo
//' @param mu pippo
//' @param delta pippo
//' @param Gamma pippo
//'
//' @return A length-one numeric vector.
//' @export
//'
//' @examples
//' # nothing for now!
// [[Rcpp::export]]
double Gamma_logpost(const arma::mat& Gamma, const arma::mat& x, const double& eta) {
  int M = Gamma.n_rows, n = x.n_rows;
  double tmp = 0;
  Rcpp::NumericVector x_row_i(M);

  for (unsigned i = 0; i < n; i++) {
    for (unsigned m = 0; m < M; m++) {
      x_row_i(m) = x(i, m);
    }
    tmp += gausscopdens(x_row_i, Gamma, false, true);
  }

  double lpost = tmp + dlkj_arma(Gamma, eta, true);

  return lpost;
}

//' @export
// [[Rcpp::export]]
double logpost(const double& mu, const double& delta, const double& y, const double& n, const double& w, const double& gamma, const double& eps, const double& eps_ab) {
  double theta = mu + delta;
  double p = expit_double(theta);

  double a = R::pbinom(y - 1, n, p, 1, 1);
  double b = R::pbinom(y, n, p, 1, 1);
  // if (a > (1 - eps_ab)) {
  //   a = 1 - eps_ab;
  // } else if (a < eps_ab) {
  //   a = eps_ab;
  // }
  // if (b > (1 - eps_ab)) {
  //   b = 1 - eps_ab;
  // } else if (b < eps_ab) {
  //   b = eps_ab;
  // }
  // double eps_a = a/2.0, eps_b = (1.0 - b)/2.0;
  // if (a == b) {
  //   if (a == 0) {
  //     b += eps;
  //   } else if (b == 1) {
  //     a -= eps;
  //   } else {
  //     a -= fmin(eps/2.0, eps_a);
  //     b += fmin(eps/2.0, eps_b);
  //   }
  // }
  double phi_inv_a = R::qnorm(a, 0.0, 1.0, 1, 1);
  double phi_inv_b = R::qnorm(b, 0.0, 1.0, 1, 1);
  double tmp = log(R::pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 1, 0) - R::pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 1, 0));
  if (!R_FINITE(tmp)) {
    tmp = log(R::pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 0, 0) - R::pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 0, 0));
  }
  if (!R_FINITE(tmp)) {
    double tmp_a = R::pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 1, 1);
    double tmp_b = R::pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 1, 1);
    tmp = log(exp(tmp_b - tmp_a) - 1) + tmp_a;
  }
  if (!R_FINITE(tmp)) {
    double tmp_a = R::pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 0, 1);
    double tmp_b = R::pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 0, 1);
    tmp = log(exp(tmp_a - tmp_b) - 1) + tmp_b;
  }
  if (!R_FINITE(tmp)) {
    // Rprintf("mu = %4.10f - delta = %4.10f - y = %4.10f - n = %4.10f - w = %4.10f - gamma = %4.10f - theta = %4.10f - p = %4.10f - a = %4.10f - b) = %4.10f - phi_inv_a = %4.10f - phi_inv_b = %4.10f\n", mu, delta, y, n, w, gamma, theta, p, a, b, phi_inv_a, phi_inv_b);
    return NA_REAL;
  }

  return tmp;
}

//' Log-likelihood of copula based model for a multivariate NMA.
//'
//' Evaluation of the log-likelihood.
//'
//' @param nc_data pippo
//' @param x pippo
//' @param mu pippo
//' @param delta pippo
//' @param Gamma pippo
//'
//' @return A length-one numeric vector.
//' @export
//'
//' @examples
//' # nothing for now!
// [[Rcpp::export]]
double mu_logpost(const double& mu, const double& delta, const double& y, const double& n, const double& w, const double& gamma, const double& mu_sigma, const double& eps, const double& eps_ab) {
  double tmp = logpost(mu, delta, y, n, w, gamma, eps, eps_ab);
  if (ISNAN(tmp)) {
    // Rprintf("mu_logpost ==> NA\n");
    return NA_REAL;
  }

  double lpost = tmp + R::dnorm(mu, 0, mu_sigma, 1);

  return lpost;
}

//' Log-likelihood of copula based model for a multivariate NMA.
//'
//' Evaluation of the log-likelihood.
//'
//' @param nc_data pippo
//' @param x pippo
//' @param mu pippo
//' @param delta pippo
//' @param Gamma pippo
//'
//' @return A length-one numeric vector.
//' @export
//'
//' @examples
//' # nothing for now!
// [[Rcpp::export]]
double delta_logpost(const double& delta, const double& mu, const double& tau, const double& eta, const double& y, const double& n, const double& w, const double& gamma, const double& eps, const double& eps_ab) {
  double tmp = logpost(mu, delta, y, n, w, gamma, eps, eps_ab);
  if (ISNAN(tmp)) {
    // Rprintf("delta_logpost ==> NA\n");
    return NA_REAL;
  }

  double lpost = tmp + R::dnorm(delta, tau, eta, 1);

  return lpost;
}

//' Log-likelihood of copula based model for a multivariate NMA.
//'
//' Evaluation of the log-likelihood.
//'
//' @param nc_data pippo
//' @param x pippo
//' @param mu pippo
//' @param delta pippo
//' @param Gamma pippo
//'
//' @return A length-one numeric vector.
//' @export
//'
//' @examples
//' # nothing for now!
// [[Rcpp::export]]
double delta_logprior(const arma::mat& delta, const arma::mat& d, const arma::mat& Sigma_M, const Rcpp::IntegerVector& trt_arms, const Rcpp::IntegerVector& baseline, const Rcpp::IntegerVector& narms) {
  int M = d.n_cols, n_study = narms.size(), ns_count = 0;
  double dens = 0;
  arma::mat delta_tmp, d_tmp, Sigma, d_base(1, M), delta_mat_tmp;
  arma::vec delta_vec_tmp, d_vec_tmp, index_tmp, trt_arms_arma = Rcpp::as<arma::vec>(trt_arms);
  arma::uvec index, which_base;

  for (unsigned int s = 0; s < n_study; s++) {
    delta_tmp.set_size(narms(s), M);
    delta_tmp = delta.submat(ns_count, 0, ns_count + narms(s) - 1, M - 1);
    index_tmp.set_size(narms(s));
    index_tmp = trt_arms_arma.subvec(ns_count, ns_count + narms(s) - 1);
    which_base.set_size(narms(s));
    which_base = arma::find(index_tmp == baseline(s));
    delta_tmp.shed_row(which_base(0));
    delta_vec_tmp.set_size(M*(narms(s) - 1));
    delta_vec_tmp = mat_to_vec(delta_tmp, false, 0);

    d_base = d.row(baseline(s) - 1);
    index.set_size(narms(s));
    index = arma::conv_to<arma::uvec>::from(index_tmp - 1);
    d_tmp = d.rows(index);
    d_tmp.shed_row(which_base(0));
    d_tmp.each_row() -= d_base;
    d_vec_tmp.set_size(M*(narms(s) - 1));
    d_vec_tmp = mat_to_vec(d_tmp, false, 0);

    Sigma.set_size(M*(narms(s) - 1), M*(narms(s) - 1));
    Sigma = Sigma_block(Sigma_M, narms(s) - 1);

    delta_mat_tmp.set_size(1, M*(narms(s) - 1));
    delta_mat_tmp.row(0) = arma::trans(delta_vec_tmp);
    dens += dmvn_arma(delta_mat_tmp, d_vec_tmp, Sigma, true)(0);

    ns_count += narms(s);
  }

  return dens;
}

//' Log-likelihood of copula based model for a multivariate NMA.
//'
//' Evaluation of the log-likelihood.
//'
//' @param nc_data pippo
//' @param x pippo
//' @param mu pippo
//' @param delta pippo
//' @param Gamma pippo
//'
//' @return A length-one numeric vector.
//' @export
//'
//' @examples
//' # nothing for now!
// [[Rcpp::export]]
double d_logprior(const arma::mat& d, const double& d_sigma, const int& ref_trt) {
  int n_trt = d.n_rows, M = d.n_cols;
  double d_prior = 0;

  for (unsigned int k = 0; k < n_trt; k++) {
    if (k != (ref_trt - 1)) {
      for (unsigned int m = 0; m < M; m++) {
        d_prior += R::dnorm(d(k, m), 0, d_sigma, 1);
      }
    }
  }

  return d_prior;
}

//' Log-likelihood of copula based model for a multivariate NMA.
//'
//' Evaluation of the log-likelihood.
//'
//' @param nc_data pippo
//' @param x pippo
//' @param mu pippo
//' @param delta pippo
//' @param Gamma pippo
//'
//' @return A length-one numeric vector.
//' @export
//'
//' @examples
//' # nothing for now!
// [[Rcpp::export]]
double d_logpost_multi(const arma::mat& d, const arma::mat& delta, const arma::mat& Sigma_M, const Rcpp::IntegerVector& trt_arms, const Rcpp::IntegerVector& baseline, const Rcpp::IntegerVector& narms, const double& d_sigma, const int& ref_trt) {

  double delta_prior = delta_logprior(delta, d, Sigma_M, trt_arms, baseline, narms);

  double d_prior = d_logprior(d, d_sigma, ref_trt);

  return (delta_prior + d_prior);
}

//' Log-likelihood of copula based model for a multivariate NMA.
//'
//' Evaluation of the log-likelihood.
//'
//' @param nc_data pippo
//' @param x pippo
//' @param mu pippo
//' @param delta pippo
//' @param Gamma pippo
//'
//' @return A length-one numeric vector.
//' @export
//'
//' @examples
//' # nothing for now!
// [[Rcpp::export]]
double d_logpost(const arma::mat& d, const arma::mat& delta, const arma::mat& Sigma_M, const Rcpp::IntegerVector& trt, const Rcpp::IntegerVector& baseline, const Rcpp::IntegerVector& narms_study, const double& d_sigma, const int& ref_trt) {
  int M = d.n_cols, n_datapoints = delta.n_rows;
  double delta_prior = 0;
  arma::mat d_base(1, M), d_trt(1, M), d_k(1, M), delta_k(1, M), d_k_tmp(1, M);
  arma::mat Sigma_M_k(Sigma_M);
  int count_trt = 1;

  for (unsigned int ik = 0; ik < n_datapoints; ik++) {
    if (trt(ik) != baseline(ik)) {
      if (narms_study(ik) == 2) {
        d_base = d.row(baseline(ik) - 1);
        d_trt = d.row(trt(ik) - 1);
        d_k = d_trt - d_base;
        delta_k = delta.row(ik);
        delta_prior += dmvn_arma(delta_k, arma::trans(d_k), Sigma_M, true)(0);
      } else {
        if (count_trt == 1) {
          d_base = d.row(baseline(ik) - 1);
          d_trt = d.row(trt(ik) - 1);
          d_k = d_trt - d_base;
          delta_k = delta.row(ik);
          delta_prior += dmvn_arma(delta_k, arma::trans(d_k), Sigma_M, true)(0);
        } else {
          d_base = d.row(baseline(ik) - 1);
          d_trt = d.row(trt(ik) - 1);
          d_k = d_trt - d_base;
          for (unsigned int k = 1; k < count_trt; k++) {
            d_base = d.row(baseline(ik) - 1);
            d_trt = d.row(trt(ik - k) - 1);
            d_k_tmp = d_trt - d_base;
            delta_k = delta.row(ik - k);
            d_k += (delta_k - d_k_tmp)/((double)count_trt);
          }
          Sigma_M_k = Sigma_M*(count_trt + 1)/((double)(2*count_trt));
          delta_k = delta.row(ik);
          if (!is_positive_definite(Sigma_M_k, 528)) {
            Rprintf("d_logpost\n");
            Sigma_M_k = make_positive_definite(Sigma_M_k);
          }
          delta_prior += dmvn_arma(delta_k, arma::trans(d_k), Sigma_M_k, true)(0);
        }
        count_trt++;
        if (count_trt == narms_study(ik)) {
          count_trt = 1;
        }
      }
    }
  }

  double d_prior = d_logprior(d, d_sigma, ref_trt);

  return (delta_prior + d_prior);
}

//' Log-likelihood of copula based model for a multivariate NMA.
//'
//' Evaluation of the log-likelihood.
//'
//' @param nc_data pippo
//' @param x pippo
//' @param mu pippo
//' @param delta pippo
//' @param Gamma pippo
//'
//' @return A length-one numeric vector.
//' @export
//'
//' @examples
//' # nothing for now!
// [[Rcpp::export]]
double Sigma_M_logpost(const arma::mat& d, const arma::mat& delta, const arma::mat& Sigma_M, const Rcpp::IntegerVector& trt, const Rcpp::IntegerVector& baseline, const Rcpp::IntegerVector& narms_study, const double& sigma_r) {
  int M = d.n_cols, n_datapoints = delta.n_rows;
  double delta_prior = 0;
  arma::mat d_base(1, M), d_trt(1, M), d_k(1, M), delta_k(1, M), d_k_tmp(1, M);
  arma::mat Sigma_M_k(Sigma_M);
  int count_trt = 1;

  for (unsigned int ik = 0; ik < n_datapoints; ik++) {
    if (trt(ik) != baseline(ik)) {
      if (narms_study(ik) == 2) {
        d_base = d.row(baseline(ik) - 1);
        d_trt = d.row(trt(ik) - 1);
        d_k = d_trt - d_base;
        delta_k = delta.row(ik);
        delta_prior += dmvn_arma(delta_k, arma::trans(d_k), Sigma_M, true)(0);
      } else {
        if (count_trt == 1) {
          d_base = d.row(baseline(ik) - 1);
          d_trt = d.row(trt(ik) - 1);
          d_k = d_trt - d_base;
          delta_k = delta.row(ik);
          delta_prior += dmvn_arma(delta_k, arma::trans(d_k), Sigma_M, true)(0);
        } else {
          d_base = d.row(baseline(ik) - 1);
          d_trt = d.row(trt(ik) - 1);
          d_k = d_trt - d_base;
          for (unsigned int k = 1; k < count_trt; k++) {
            d_base = d.row(baseline(ik) - 1);
            d_trt = d.row(trt(ik - k) - 1);
            d_k_tmp = d_trt - d_base;
            delta_k = delta.row(ik - k);
            d_k += (delta_k - d_k_tmp)/((double)count_trt);
          }
          Sigma_M_k = Sigma_M*(count_trt + 1)/((double)(2*count_trt));
          delta_k = delta.row(ik);
          if (!is_positive_definite(Sigma_M_k, 598)) {
            Rprintf("Sigma_M_logpost\n");
            Sigma_M_k = make_positive_definite(Sigma_M_k);
          }
          delta_prior += dmvn_arma(delta_k, arma::trans(d_k), Sigma_M_k, true)(0);
        }
        count_trt++;
        if (count_trt == narms_study(ik)) {
          count_trt = 1;
        }
      }
    }
  }

  double Sigma_M_prior = dlogchol_arma(Sigma_M, sigma_r, true);

  return (delta_prior + Sigma_M_prior);
}

//' Log-likelihood of copula based model for a multivariate NMA.
//'
//' Evaluation of the log-likelihood.
//'
//' @param data pippo
//' @param init pippo
//' @param totiter pippo
//' @param prior pippo
//' @param prop pippo
//' @param verbose pippo
//'
//' @return A list containing the output of the MCMC simulation.
//' @export
//'
//' @examples
//' # nothing for now!
// [[Rcpp::export]]
Rcpp::List nc_mcmc_opt(const Rcpp::RObject& data, const Rcpp::List& init, const int& totiter, const Rcpp::List& prior, const Rcpp::List& prop, const Rcpp::List& tuning, const Rcpp::List& adapt, const bool& verbose) {

  int niter = 0, adapt_k = 0, n_datapoints = data.slot("n_datapoints"), M = data.slot("n_outcomes"), n_trt = data.slot("n_treatments"), n_study = data.slot("n_study"), nn = M*(M + 1)/2;

  double eps = Rcpp::as<double>(tuning["eps"]), eps_ab = Rcpp::as<double>(tuning["eps_ab"]);//, eps_a = 0.0, eps_b = 0.0;

  int every = Rcpp::as<int>(adapt["every"]), maxiter = Rcpp::as<int>(adapt["maxiter"]); //, miniter = Rcpp::as<int>(adapt["miniter"]);

  double alpha = Rcpp::as<double>(adapt["alpha"]), beta = Rcpp::as<double>(adapt["beta"]), gamma = Rcpp::as<double>(adapt["gamma"]), tar = Rcpp::as<double>(adapt["tar"]), tol = Rcpp::as<double>(adapt["tol"]);

  Rcpp::DataFrame study_data(data.slot("study_data"));
  Rcpp::DataFrame study_id(data.slot("study_id"));
  Rcpp::IntegerVector cols_y(M), cols_n(M);
  for (unsigned int i = 0; i < M; i++) {
    cols_y(i) = 6 + M + i;
    cols_n(i) = 6 + i;
  }
  Rcpp::List Gamma_init = Rcpp::as<Rcpp::List>(init["Gamma"]);
  Rcpp::List Gamma = Rcpp::clone(Gamma_init);
  Rcpp::IntegerVector trt = study_data["trt"];
  int ref_trt = data.slot("ref_trt");
  Rcpp::IntegerVector narms = study_id["narms"];
  Rcpp::IntegerVector narms_study = study_data["narms"];

  Rcpp::NumericMatrix mu_init = init["mu"], delta_init = init["delta"], x_init = init["x"];
  Rcpp::NumericMatrix mu = Rcpp::clone(mu_init), delta = Rcpp::clone(delta_init), x = Rcpp::clone(x_init);
  arma::mat Sigma_M = Rcpp::as<arma::mat>(init["Sigma_M"]), d = Rcpp::as<arma::mat>(init["d"]);
  arma::mat delta_arma = Rcpp::as<arma::mat>(delta), delta_arma_prop = Rcpp::as<arma::mat>(delta);
  Rcpp::NumericMatrix y = df_nm(study_data, cols_y);
  Rcpp::NumericMatrix n_data = df_nm(study_data, cols_n);
  Rcpp::NumericMatrix n_imp = n_imputed(n_data);
  Rcpp::NumericMatrix x_imp(n_datapoints, M), y_imp(n_datapoints, M);
  arma::mat x_imp_arma(n_datapoints, M);

  Rcpp::List D_list_init = Rcpp::as<Rcpp::List>(init["D"]);
  Rcpp::List D_list = Rcpp::clone(D_list_init);
  arma::mat D(M, M), D_inv(M, M), S_q = arma::zeros<arma::mat>(M, M), D_prop(M, M);
  Rcpp::List x2t_split(n_trt);
  arma::mat x_q;
  int n_q;

  int nGamma = prior["nGamma"];
  arma::mat Gamma_q_prop(M, M), D_prop_inv(M, M), Gamma_q_curr(M, M), Sigma_q_prop(M, M);
  double ran_unif, Gamma_rate, mu_rate, delta_rate, d_rate, Sigma_M_rate;
  arma::vec accept_Gamma_q = arma::zeros<arma::vec>(nGamma);
  arma::cube Gamma_chain(nGamma, M*(M - 1)/2, totiter);

  arma::mat Gamma_k(M, M), Gamma_k_m_m(1, M), Gamma_k_m(M, M), x_ik_m(1, M);
  double a_ikm = 0.0, b_ikm = 0.0, w_ikm = 0.0, gamma_ikm = 0.0, theta_ikm = 0.0, p_ikm = 0.0;
  arma::cube x_chain(n_datapoints, M, totiter);

  Rcpp::Environment nc_env("package:netcopula");
  Rcpp::Function mu_logpost_R = nc_env["mu_logpost_func"];
  Rcpp::Function mu_logpost_R_quad = nc_env["mu_logpost_quad"];
  Rcpp::Function delta_logpost_R = nc_env["delta_logpost_func"];
  Rcpp::Function delta_logpost_R_quad = nc_env["delta_logpost_quad"];
  Rcpp::NumericVector parscale(1), ndeps(1);
  double mu_sigma = sqrt(Rcpp::as<double>(prior["mu_sigma2"]));
  parscale(0) = 1;
  ndeps(0) = 0.001;
  arma::vec coef(3);
  Rcpp::List args_mu = Rcpp::List::create(
      Rcpp::Named("delta") = 0.0,   // temporary value
      Rcpp::Named("y") = 0,         // temporary value
      Rcpp::Named("n") = 1,         // temporary value
      Rcpp::Named("w") = 0.0,       // temporary value
      Rcpp::Named("gamma") = 1.0,   // temporary value
      Rcpp::Named("mu_sigma") = mu_sigma,
      Rcpp::Named("eps") = eps,
      Rcpp::Named("eps_ab") = eps_ab);
  Rcpp::List args_delta = Rcpp::List::create(
      Rcpp::Named("mu") = 0,        // temporary value
      Rcpp::Named("tau") = 0.0,     // temporary value
      Rcpp::Named("eta") = 1.0,     // temporary value
      Rcpp::Named("y") = 0,         // temporary value
      Rcpp::Named("n") = 1,         // temporary value
      Rcpp::Named("w") = 0.0,       // temporary value
      Rcpp::Named("gamma") = 1.0,   // temporary value
      Rcpp::Named("eps") = eps,
      Rcpp::Named("eps_ab") = eps_ab);
  Rcpp::List args_coef = Rcpp::List::create(Rcpp::Named("coef") = coef);
  Rcpp::List options = Rcpp::List::create(
      Rcpp::Named("trace") = 0,
      Rcpp::Named("fnscale") = -1,
      Rcpp::Named("parscale") = Rcpp::as<Rcpp::NumericVector>(parscale),
      Rcpp::Named("ndeps") = Rcpp::as<Rcpp::NumericVector>(ndeps),
      Rcpp::Named("maxit") = 500,
      Rcpp::Named("abstol") = R_NegInf,
      Rcpp::Named("reltol") = sqrt(machine_eps),
      Rcpp::Named("alpha") = 1.0,
      Rcpp::Named("beta") = 0.5,
      Rcpp::Named("gamma") = 2.0);
  double xmin = Rcpp::as<double>(tuning["xmin"]), xmax = Rcpp::as<double>(tuning["xmax"]);
  double tol_laplace = pow(machine_eps, .25);
  Rcpp::List fit_mu(2), fit_delta(2);
  arma::vec mode_mu(1), mode_delta(1);
  arma::mat var_mu(1, 1), var_delta(1, 1);
  double tau_ikm = 0.0, eta_ikm = 0.0;
  arma::mat Sigma_M_m_m(1, M), Sigma_M_m(M, M), delta_ik_m(1, M);
  arma::mat d_1(1, M), d_k(1, M), d_1k(1, M), d_1k_m(1, M);
  Rcpp::IntegerVector baseline = study_data["baseline"];
  Rcpp::IntegerVector baseline_id = study_id["baseline"];
  Rcpp::NumericMatrix mu_long = param_long(mu, narms, false);

  double mu_curr = 0, mu_prop = 0;
  arma::mat mu_curr_arma(1, 1), mu_prop_arma(1, 1);
  double mu_A = 0.0, mu_B = 0.0, target_mu_prop = 0.0, target_mu_curr = 0.0;
  int accept_mu = 0;
  arma::cube mu_chain(n_study, M, totiter);

  double delta_curr = 0, delta_prop = 0;
  arma::mat delta_curr_arma(1, 1), delta_prop_arma(1, 1);
  double delta_A = 0.0, delta_B = 0.0, target_delta_prop = 0.0, target_delta_curr = 0.0;
  int accept_delta = 0;
  arma::cube delta_chain(n_datapoints, M, totiter);

  int accept_d = 0;
  arma::cube d_chain(n_trt, M, totiter);

  int accept_Sigma_M = 0;
  arma::mat Sigma_M_chain(totiter, nn);
  arma::mat beta_chain(totiter, nn);

  double rho_d = 2.38/sqrt(M*(n_trt - 1));
  arma::mat cov_d = arma::eye(M*(n_trt - 1), M*(n_trt - 1));
  arma::vec d_curr(M*(n_trt - 1));
  arma::mat d_prop(1, M*(n_trt - 1)), d_curr_ref(n_trt, M), d_prop_ref(n_trt, M);
  double d_sigma = sqrt(Rcpp::as<double>(prior["d_sigma2"]));
  double d_A = 0.0, d_B = 0.0, target_d_prop = 0.0, target_d_curr = 0.0;

  double rho_beta = 2.38/sqrt(nn);
  arma::mat cov_beta = arma::eye(nn, nn);
  arma::mat Sigma_M_curr(Sigma_M), Sigma_M_prop(Sigma_M);
  arma::vec beta_curr(nn), beta_prop(nn);
  double sigma_r = Rcpp::as<double>(prior["beta_sigma"]);
  double Sigma_M_A = 0.0, Sigma_M_B = 0.0, target_Sigma_M_prop = 0.0, target_Sigma_M_curr = 0.0;

  arma::vec mu_d(M*(n_trt - 1));
  arma::cube theta_d(n_trt, M, every);
  arma::mat theta_d_reshaped(every, M*(n_trt - 1));
  arma::vec ar_d_vec(totiter), ar_d(2);
  ar_d(0) = 0.0;
  Rcpp::List prop_d(4);

  arma::vec mu_beta(nn);
  arma::mat theta_beta(every, nn);
  arma::vec ar_Sigma_M_vec(totiter), ar_Sigma_M(2);
  ar_Sigma_M(0) = 0.0;
  Rcpp::List prop_beta(4);

  arma::vec loglik(totiter, arma::fill::zeros), logprior(totiter, arma::fill::zeros), logpost(totiter, arma::fill::zeros);

  GetRNGstate();

  while (niter < totiter) {
    // imputing x and y variables
    x_imp = x_imputed(x, Gamma, trt);
    x_imp_arma = Rcpp::as<arma::mat>(x_imp);
    y_imp = y_imputed(y, x_imp, narms, mu, delta, n_imp);

    // updating Gamma (correlation of latent variables)
    if (nGamma > 1) {
      x2t_split = split_nm(x_imp, trt);
      for (unsigned q = 0; q < n_trt; q++) {
        x_q = Rcpp::as<arma::mat>(x2t_split[q]);
        n_q = x_q.n_rows;
        D = Rcpp::as<arma::mat>(D_list(q));
        // D_inv = arma::inv_sympd(D);
        for (unsigned int k = 0; k < n_q; k++) {
          S_q += D*arma::trans(x_q.row(k))*x_q.row(k)*D;
        }
        if (!is_positive_definite(S_q, 812)) {
          Rprintf("S_q - nGamma > 1\n");
          S_q = make_positive_definite(S_q);
        }
        // in the following, degrees of freedom are set to (n_q > M) ? n_q : M)
        // otherwise, when n_q < M the IW routine returns an error
        // (indeterminate system)
        Sigma_q_prop = rinvwish_arma(((n_q > M) ? n_q : M), S_q);
        D_prop = arma::sqrt(arma::diagmat(Sigma_q_prop));
        D_prop_inv = arma::inv_sympd(D_prop);
        Gamma_q_prop = cov2cor_rcpp(Sigma_q_prop);

        Gamma_q_curr = Rcpp::as<arma::mat>(Gamma(q));
        Gamma_rate = 0.5*(M + 1)*(log(arma::det(Gamma_q_prop)) - log(arma::det(Gamma_q_curr)));
        ran_unif = R::runif(0.0, 1.0);
        if (ran_unif < exp(Gamma_rate)) {
          Gamma_q_curr = Gamma_q_prop;
          D = D_prop;
          accept_Gamma_q(q)++;
        }
        Gamma(q) = Gamma_q_curr;
        D_list(q) = D;

        S_q = arma::zeros<arma::mat>(M, M);
      }
    } else {
      D = Rcpp::as<arma::mat>(D_list(0));
      // D_inv = arma::inv_sympd(D);
      for (unsigned int k = 0; k < n_datapoints; k++) {
        S_q += D*arma::trans(x_imp_arma.row(k))*x_imp_arma.row(k)*D;
      }
      if (!is_positive_definite(S_q, 843)) {
        Rprintf("S_q - nGamma = 1\n");
        S_q = make_positive_definite(S_q);
      }
      // in the following, degrees of freedom are set to (n_datapoints > M) ?
      // n_datapoints : M) otherwise, when n_datapoints < M the IW routine
      // returns an error (indeterminate system)
      Sigma_q_prop = rinvwish_arma(((n_datapoints > M) ? n_datapoints : M), S_q);
      D_prop = arma::sqrt(arma::diagmat(Sigma_q_prop));
      D_prop_inv = arma::inv_sympd(D_prop);
      Gamma_q_prop = cov2cor_rcpp(Sigma_q_prop);

      Gamma_q_curr = Rcpp::as<arma::mat>(Gamma(0));
      Gamma_rate = 0.5*(M + 1)*(log(arma::det(Gamma_q_prop)) - log(arma::det(Gamma_q_curr)));
      ran_unif = R::runif(0.0, 1.0);
      if (ran_unif < exp(Gamma_rate)) {
        Gamma_q_curr = Gamma_q_prop;
        D = D_prop;
        accept_Gamma_q(0)++;
      }
      Gamma(0) = Gamma_q_curr;
      D_list(0) = D;

      S_q = arma::zeros<arma::mat>(M, M);
    }

    Gamma_chain.slice(niter) = list_mat(Gamma);

    // updating mu (study-specific baseline effects) and delta (study-specific
    // [random] treatment effects)
    for (unsigned int ik = 0; ik < n_datapoints; ik++) {
      if (nGamma > 1) {
        Gamma_k = Rcpp::as<arma::mat>(Gamma(trt(ik) - 1));
      } else {
        Gamma_k = Rcpp::as<arma::mat>(Gamma(0));
      }
      d_1 = d.row(baseline(ik) - 1);
      d_k = d.row(trt(ik) - 1);
      d_1k = d_k - d_1;
      if (trt(ik) == baseline(ik)) {
        // update mu (study-specific baseline effect)
        for (unsigned int m = 0; m < M; m++) {
          Gamma_k_m_m = Gamma_k.row(m);
          Gamma_k_m_m.shed_col(m);
          Gamma_k_m = Gamma_k;
          Gamma_k_m.shed_col(m);
          Gamma_k_m.shed_row(m);
          x_ik_m = x_imp_arma.row(ik);
          x_ik_m.shed_col(m);
          w_ikm = arma::as_scalar(Gamma_k_m_m*arma::inv(Gamma_k_m)*arma::trans(x_ik_m));
          gamma_ikm = sqrt(arma::as_scalar(Gamma_k(m, m) - Gamma_k_m_m*arma::inv(Gamma_k_m)*arma::trans(Gamma_k_m_m)));

          args_mu["delta"] = delta(ik, m);
          args_mu["y"] = y_imp(ik, m);
          args_mu["n"] = n_imp(ik, m);
          args_mu["w"] = w_ikm;
          args_mu["gamma"] = gamma_ikm;
          fit_mu = laplace_u_rcpp(mu_logpost_R, xmin, xmax, tol_laplace, args_mu, options);
          mode_mu(0) = Rcpp::as<double>(fit_mu["mode"]);
          var_mu(0, 0) = Rcpp::as<double>(fit_mu["var"]);
          if (ISNAN(var_mu(0, 0))) {
            coef = ols_coef(xmin, xmax, args_mu, false);
            args_coef["coef"] = coef;
            fit_mu = laplace_u_rcpp(mu_logpost_R_quad, xmin, xmax, tol_laplace, args_coef, options);
            mode_mu(0) = Rcpp::as<double>(fit_mu["mode"]);
            var_mu(0, 0) = Rcpp::as<double>(fit_mu["var"]);
          }
          if (!is_positive_definite(var_mu, 910)) {
            Rprintf("var_mu\n");
            var_mu = make_positive_definite(var_mu);
          }
          if (var_mu(0, 0) <= 0 || ISNAN(var_mu(0, 0))) {
            Rprintf("mu\n");
            mode_mu.print();
            var_mu.print();
            Rprintf("mu = %4.10f\n", mu_long(ik, m));
            Rprintf("delta = %4.10f\n", delta(ik, m));
            Rprintf("tau = %4.10f\n", tau_ikm);
            Rprintf("eta = %4.10f\n", eta_ikm);
            Rprintf("y = %4.0f\n", y_imp(ik, m));
            Rprintf("n = %4.0f\n", n_imp(ik, m));
            Rprintf("w = %4.10f\n", w_ikm);
            Rprintf("gamma = %4.10f\n", gamma_ikm);
          }

          mu_prop = rmvt_arma(1, mode_mu, var_mu, 7)(0, 0);
          mu_prop_arma(0, 0) = mu_prop;
          mu_curr = mu_long(ik, m);
          mu_curr_arma(0, 0) = mu_curr;

          target_mu_curr = mu_logpost(mu_curr, delta(ik, m), static_cast<double>(y_imp(ik, m)), static_cast<double>(n_imp(ik, m)), w_ikm, gamma_ikm, mu_sigma, eps, eps_ab);
          target_mu_prop = mu_logpost(mu_prop, delta(ik, m), static_cast<double>(y_imp(ik, m)), static_cast<double>(n_imp(ik, m)), w_ikm, gamma_ikm, mu_sigma, eps, eps_ab);

          mu_A = target_mu_prop - target_mu_curr;
          mu_B = dmvt_arma(mu_curr_arma, mode_mu, var_mu, 7, true)(0) - dmvt_arma(mu_prop_arma, mode_mu, var_mu, 7, true)(0);
          mu_rate = exp(mu_A + mu_B);
          ran_unif = R::runif(0.0, 1.0);
          if (!ISNAN(target_mu_prop) && !ISNAN(target_mu_curr)) {
            if (ran_unif < mu_rate) {
              mu_long(ik, m) = mu_prop;
              mu = param_wide(mu_long, narms, trt, baseline);
              mu_long = param_long(mu, narms, false); // is this necessary?
              accept_mu++;
            }
          } else {
            mu_rate = 0;
          }
        }
      } else {
        // update delta (study-specific [random] treatment effects)
        for (unsigned int m = 0; m < M; m++) {
          Gamma_k_m_m = Gamma_k.row(m);
          Gamma_k_m_m.shed_col(m);
          Gamma_k_m = Gamma_k;
          Gamma_k_m.shed_col(m);
          Gamma_k_m.shed_row(m);
          x_ik_m = x_imp_arma.row(ik);
          x_ik_m.shed_col(m);
          w_ikm = arma::as_scalar(Gamma_k_m_m*arma::inv(Gamma_k_m)*arma::trans(x_ik_m));
          gamma_ikm = sqrt(arma::as_scalar(Gamma_k(m, m) - Gamma_k_m_m*arma::inv(Gamma_k_m)*arma::trans(Gamma_k_m_m)));

          Sigma_M_m_m = Sigma_M.row(m);
          Sigma_M_m_m.shed_col(m);
          Sigma_M_m = Sigma_M;
          Sigma_M_m.shed_col(m);
          Sigma_M_m.shed_row(m);
          delta_ik_m = delta_arma.row(ik);
          delta_ik_m.shed_col(m);
          d_1k_m = d_1k;
          d_1k_m.shed_col(m);
          tau_ikm = d_1k(m) + arma::as_scalar(Sigma_M_m_m*arma::inv(Sigma_M_m)*arma::trans(delta_ik_m - d_1k_m));
          eta_ikm = sqrt(arma::as_scalar(Sigma_M(m, m) - Sigma_M_m_m*arma::inv(Sigma_M_m)*arma::trans(Sigma_M_m_m)));

          args_delta["mu"] = mu_long(ik, m);
          args_delta["tau"] = tau_ikm;
          args_delta["eta"] = eta_ikm;
          args_delta["y"] = y_imp(ik, m);
          args_delta["n"] = n_imp(ik, m);
          args_delta["w"] = w_ikm;
          args_delta["gamma"] = gamma_ikm;
          fit_delta = laplace_u_rcpp(delta_logpost_R, xmin, xmax, tol_laplace, args_delta, options);
          mode_delta(0) = Rcpp::as<double>(fit_delta["mode"]);
          var_delta(0, 0) = Rcpp::as<double>(fit_delta["var"]);
          if (ISNAN(var_delta(0, 0))) {
            coef = ols_coef(xmin, xmax, args_delta, true);
            args_coef["coef"] = coef;
            fit_delta = laplace_u_rcpp(delta_logpost_R_quad, xmin, xmax, tol_laplace, args_coef, options);
            mode_delta(0) = Rcpp::as<double>(fit_delta["mode"]);
            var_delta(0, 0) = Rcpp::as<double>(fit_delta["var"]);
          }
          if (!is_positive_definite(var_delta, 991)) {
            Rprintf("var_delta\n");
            var_delta = make_positive_definite(var_delta);
          }
          if (var_delta(0, 0) <= 0 || ISNAN(var_delta(0, 0))) {
            Rprintf("delta\n");
            mode_delta.print();
            var_delta.print();
            Rprintf("mu = %4.10f\n", mu_long(ik, m));
            Rprintf("delta = %4.10f\n", delta(ik, m));
            Rprintf("tau = %4.10f\n", tau_ikm);
            Rprintf("eta = %4.10f\n", eta_ikm);
            Rprintf("y = %4.0f\n", y_imp(ik, m));
            Rprintf("n = %4.0f\n", n_imp(ik, m));
            Rprintf("w = %4.10f\n", w_ikm);
            Rprintf("gamma = %4.10f\n", gamma_ikm);
          }

          delta_prop = rmvt_arma(1, mode_delta, var_delta, 7)(0, 0);
          delta_prop_arma(0, 0) = delta_prop;
          delta_curr = delta(ik, m);
          delta_curr_arma(0, 0) = delta_curr;

          // use here the quadratic approximation to delta_logpost()?

          target_delta_curr = delta_logpost(delta_curr, mu_long(ik, m), tau_ikm, eta_ikm, static_cast<double>(y_imp(ik, m)), static_cast<double>(n_imp(ik, m)), w_ikm, gamma_ikm, eps, eps_ab);
          target_delta_prop = delta_logpost(delta_prop, mu_long(ik, m), tau_ikm, eta_ikm, static_cast<double>(y_imp(ik, m)), static_cast<double>(n_imp(ik, m)), w_ikm, gamma_ikm, eps, eps_ab);

          delta_A = target_delta_prop - target_delta_curr;
          delta_B = dmvt_arma(delta_curr_arma, mode_delta, var_delta, 7, true)(0) - dmvt_arma(delta_prop_arma, mode_delta, var_delta, 7, true)(0);
          delta_rate = exp(delta_A + delta_B);
          ran_unif = R::runif(0.0, 1.0);
          if (!ISNAN(target_delta_prop) && !ISNAN(target_delta_curr)) {
            if (ran_unif < delta_rate) {
              delta(ik, m) = delta_prop;
              accept_delta++;
            }
          } else {
            delta_rate = 0;
          }
        }
      }
    }
    mu_chain.slice(niter) = Rcpp::as<arma::mat>(mu);
    delta_chain.slice(niter) = Rcpp::as<arma::mat>(delta);
    delta_arma = Rcpp::as<arma::mat>(delta);

    // updating x (latent variables)
    for (unsigned int ik = 0; ik < n_datapoints; ik++) {
      if (nGamma > 1) {
        Gamma_k = Rcpp::as<arma::mat>(Gamma(trt(ik) - 1));
      } else {
        Gamma_k = Rcpp::as<arma::mat>(Gamma(0));
      }
      for (unsigned int m = 0; m < M; m++) {
        if (!ISNAN(y(ik, m))) {
          theta_ikm = mu_long(ik, m) + delta(ik, m);
          p_ikm = expit_double(theta_ikm);
          a_ikm = R::pbinom(y_imp(ik, m) - 1, n_imp(ik, m), p_ikm, 1, 1); // log
          b_ikm = R::pbinom(y_imp(ik, m), n_imp(ik, m), p_ikm, 1, 1); // log
          // eps_a = a_ikm/2.0;
          // eps_b = (1.0 - b_ikm)/2.0;
          if (a_ikm == b_ikm) {
          //   if (a_ikm == 0) {
          //     b_ikm += eps;
          //   } else if (b_ikm == 1) {
          //     a_ikm -= eps;
          //   } else {
          //     a_ikm -= fmin(eps/2.0, eps_a);
          //     b_ikm += fmin(eps/2.0, eps_b);
          //   }
            Rprintf("niter = %d a = %4.5f - b = %4.5f - y = %4.0f - n = %4.0f - p = %4.10f\n", niter, a_ikm, b_ikm, y_imp(ik, m), n_imp(ik, m), p_ikm);
          } else {
            Gamma_k_m_m = Gamma_k.row(m);
            Gamma_k_m_m.shed_col(m);
            Gamma_k_m = Gamma_k;
            Gamma_k_m.shed_col(m);
            Gamma_k_m.shed_row(m);
            x_ik_m = x_imp_arma.row(ik);
            x_ik_m.shed_col(m);
            w_ikm = arma::as_scalar(Gamma_k_m_m*arma::inv(Gamma_k_m)*arma::trans(x_ik_m));
            gamma_ikm = sqrt(arma::as_scalar(Gamma_k(m, m) - Gamma_k_m_m*arma::inv(Gamma_k_m)*arma::trans(Gamma_k_m_m)));
            x(ik, m) = rtruncnorm_rcpp(1, R::qnorm(a_ikm, 0.0, 1.0, 1, 1), R::qnorm(b_ikm, 0.0, 1.0, 1, 1), w_ikm, gamma_ikm)[0];
            x_imp(ik, m) = x(ik, m);
            x_imp_arma(ik, m) = x(ik, m);
          }
        }
      }
    }
    x_chain.slice(niter) = Rcpp::as<arma::mat>(x);

    // updating d (pooled treatment effects across trials)
    d_curr = mat_to_vec(d, true, ref_trt);
    d_curr_ref = d;
    d_prop = rmvt_arma(1, d_curr, pow(rho_d, 2)*cov_d, 7);
    d_prop_ref = vec_to_mat(arma::vectorise(d_prop), M, true, ref_trt);
    target_d_curr = d_logpost(d_curr_ref, delta_arma, Sigma_M, trt, baseline, narms_study, d_sigma, ref_trt);
    target_d_prop = d_logpost(d_prop_ref, delta_arma, Sigma_M, trt, baseline, narms_study, d_sigma, ref_trt);
    d_A = target_d_prop - target_d_curr;
    // d_B dovrebbe essere la differenza della proposal distribution, ma e'
    // simmetrica, quindi sparisce; la prior e' gia' dentro la target
    d_B = 0.0; //d_logprior(d_prop_ref, d_sigma, ref_trt) - d_logprior(d_curr_ref, d_sigma, ref_trt);
    d_rate = exp(d_A + d_B);
    // Rprintf("d_A = %4.8f - d_rate = %4.8f\n", d_A, d_rate);
    ran_unif = R::runif(0.0, 1.0);
    if (ran_unif < d_rate) {
      d = d_prop_ref;
      accept_d++;
    }
    ar_d_vec(niter) = fmin(d_rate, 1);
    d_chain.slice(niter) = d;

    // updating Sigma_M (common between-study covariance structure)
    Sigma_M_curr = Sigma_M;
    beta_curr = Sigma_M_to_beta(Sigma_M_curr);
    beta_prop = arma::vectorise(rmvt_arma(1, beta_curr, pow(rho_beta, 2)*cov_beta, 7));
    Sigma_M_prop = beta_to_Sigma_M(beta_prop, M);
    if (!is_positive_definite(Sigma_M_prop, 1104)) {
      Rprintf("Sigma_M_prop\n");
      Sigma_M_prop = make_positive_definite(Sigma_M_prop);
    }
    target_Sigma_M_prop = Sigma_M_logpost(d, delta_arma, Sigma_M_prop, trt, baseline, narms_study, sigma_r);
    target_Sigma_M_curr = Sigma_M_logpost(d, delta_arma, Sigma_M_curr, trt, baseline, narms_study, sigma_r);
    Sigma_M_A = target_Sigma_M_prop - target_Sigma_M_curr;
    // Sigma_M_B dovrebbe essere la differenza della proposal distribution, ma
    // e' simmetrica, quindi sparisce; la prior invece e' gia' dentro la target
    Sigma_M_B = 0.0; //dlogchol_arma(Sigma_M_prop, sigma_r, true) - dlogchol_arma(Sigma_M_curr, sigma_r, true);
    Sigma_M_rate = exp(Sigma_M_A + Sigma_M_B);
    ran_unif = R::runif(0.0, 1.0);
    if (ran_unif < Sigma_M_rate) {
      Sigma_M = Sigma_M_prop;
      accept_Sigma_M++;
    }
    ar_Sigma_M_vec(niter) = fmin(Sigma_M_rate, 1);
    Sigma_M_chain.row(niter) = diag_tri(Sigma_M);
    beta_chain.row(niter) = arma::trans(Sigma_M_to_beta(Sigma_M));

    // adaptation step
    if (((adapt_k + 1) <= maxiter) && (((niter + 1) % every) == 0)) {
      if (verbose) {
        Rprintf("      --> performing adaptation [step %d/%d]\n", adapt_k + 1, maxiter);
      }

      // d proposal parameters
      if (fabs(ar_d(0) - tar) > tol) {
        theta_d = d_chain.slices(niter - every + 1, niter);
        theta_d_reshaped = cube_to_mat(theta_d, true, ref_trt);
        if (adapt_k == 0) {
          mu_d = arma::vectorise(arma::mean(theta_d_reshaped, 0)); // means by columns
        }
        ar_d(1) = arma::as_scalar(arma::mean(ar_d_vec.subvec(niter - every + 1, niter)));
        prop_d = rwmh_adapt(theta_d_reshaped, mu_d, rho_d, cov_d, ar_d, alpha, beta, gamma, tar, adapt_k, false, 5);
        rho_d = Rcpp::as<double>(prop_d["rho"]);
        cov_d = Rcpp::as<arma::mat>(prop_d["covariance"]);
        if (!is_positive_definite(cov_d, 1141)) {
          Rprintf("d\n");
          cov_d = make_positive_definite(cov_d);
        }
        mu_d = Rcpp::as<arma::vec>(prop_d["mu"]);
        ar_d(0) = Rcpp::as<double>(prop_d["ar"]);
      }

      // log Cholesky betas (Sigma_M) proposal parameters
      if (fabs(ar_Sigma_M(0) - tar) > tol) {
        theta_beta = beta_chain.rows(niter - every + 1, niter);
        if (adapt_k == 0) {
          mu_beta = arma::vectorise(arma::mean(theta_beta, 0)); // means by columns
        }
        ar_Sigma_M(1) = arma::as_scalar(arma::mean(ar_Sigma_M_vec.subvec(niter - every + 1, niter)));
        prop_beta = rwmh_adapt(theta_beta, mu_beta, rho_beta, cov_beta, ar_Sigma_M, alpha, beta, gamma, tar, adapt_k, false, 6);
        rho_beta = Rcpp::as<double>(prop_beta["rho"]);
        cov_beta = Rcpp::as<arma::mat>(prop_beta["covariance"]);
        if (!is_positive_definite(cov_beta, 1159)) {
          Rprintf("cov_beta\n");
          cov_beta = make_positive_definite(cov_beta);
        }
        mu_beta = Rcpp::as<arma::vec>(prop_beta["mu"]);
        ar_Sigma_M(0) = Rcpp::as<double>(prop_beta["ar"]);
      }

      adapt_k++;
    }

    // calculate the loglikelihood, logprior and logposterior
    loglik(niter) = nc_loglik(y_imp, n_imp, x_imp, trt, mu, delta, Gamma); // c'e' il problema degli indicatori nella definizione della augmented likelihood
    logprior(niter) = nc_logprior(mu, mu_sigma, d, d_sigma, Sigma_M, sigma_r, ref_trt);

    // print the information
    if ((((niter + 1) % 100) == 0) && verbose) {
      Rprintf("   iter. %d/%d ==> d: %1.3f - Gamma: %1.3f - mu: %1.3f - delta: %1.3f - Sigma_M: %1.3f\n", (niter + 1), totiter, accept_d/(static_cast<double>(totiter)), arma::mean(accept_Gamma_q)/(static_cast<double>(totiter)), accept_mu/(static_cast<double>(totiter*n_study*M)), accept_delta/(static_cast<double>(totiter*(n_datapoints - n_study)*M)), accept_Sigma_M/(static_cast<double>(totiter)));
    }

    niter++;

    R_CheckUserInterrupt();
  }

  PutRNGstate();

  Rcpp::List accept = Rcpp::List::create(
    Rcpp::Named("d") = accept_d/static_cast<double>(totiter),
    Rcpp::Named("mu") = accept_mu/static_cast<double>(totiter*n_study*M),
    Rcpp::Named("delta") = accept_delta/static_cast<double>(totiter*(n_datapoints - n_study)*M),
    Rcpp::Named("Gamma") = accept_Gamma_q/static_cast<double>(totiter),
    Rcpp::Named("Sigma") = accept_Sigma_M/static_cast<double>(totiter));

  return Rcpp::List::create(Rcpp::Named("mu") = mu_chain,
                            Rcpp::Named("delta") = delta_chain,
                            Rcpp::Named("d") = d_chain,
                            Rcpp::Named("Sigma") = Sigma_M_chain,
                            Rcpp::Named("Gamma") = Gamma_chain,
                            Rcpp::Named("x") = x_chain,
                            Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("logprior") = logprior,
                            Rcpp::Named("logpost") = (loglik + logprior),
                            Rcpp::Named("accept") = accept);
}

//' Log-likelihood of copula based model for a multivariate NMA.
//'
//' Evaluation of the log-likelihood.
//'
//' @param y pippo
//' @param n pippo
//' @param x pippo
//' @param trt pippo
//' @param mu pippo
//' @param delta pippo
//' @param Gamma pippo
//'
//' @return A length-one numeric vector.
//' @export
//'
//' @examples
//' # nothing for now!
// [[Rcpp::export]]
Rcpp::List rwmh_adapt(const arma::mat& theta, const arma::vec& mu, const double& rho, const arma::mat& cov, const arma::vec& ar, const double& alpha, const double& beta, const double& gamma, const double& tar, const int& k, const bool& iter_cols, const int& what) {

  // [See 'Stata Bayesian Analysis Reference Manual', page 158]

  double ar_new = (1 - alpha)*ar(0) + alpha*ar(1);
  double beta_k = beta/pow(k + 1, gamma);
  double rho_new = rho*exp(beta_k*(R::qnorm(ar_new/2.0, 0.0, 1.0, 1, 0) - R::qnorm(tar/2.0, 0.0, 1.0, 1, 0)));
  // if (ISNAN(rho_new)) {
  //   Rprintf("what = %d - ar(0) = %4.5f - ar(1) = %4.5f - ar_new = %4.5f - beta_k = %4.5f - rho_new = %4.5f\n", what, ar(0), ar(1), ar_new, beta_k, rho_new);
  // }
  arma::mat theta_k;
  if (iter_cols) {
    theta_k.set_size(arma::size(theta));
    theta_k = theta;
  } else {
    theta_k.set_size(arma::size(arma::trans(theta)));
    theta_k = arma::trans(theta);
  }
  int niter = theta_k.n_cols;
  arma::vec theta_hat = arma::vectorise(arma::mean(theta_k, 1)); // means by rows
  arma::vec mu_new(arma::size(mu));
  if (k > 0) {
    mu_new = mu + beta_k*(theta_hat - mu);
  } else {
    mu_new = mu;
  }
  arma::mat theta_k_c = theta_k.each_col() - theta_hat;
  arma::mat cov_hat = theta_k_c*arma::trans(theta_k_c)/niter;
  arma::mat cov_new = (1 - beta_k)*cov + beta_k*cov_hat;
  if ((what == 1 || what == 2)  && (cov_new(0, 0) <= 0)) {
    cov_new = 1e-1;
    // Rprintf("cov = %4.5f - cov_hat = %4.5f\n", cov(0, 0), cov_hat(0, 0));
  }

  return Rcpp::List::create(Rcpp::Named("rho") = rho_new,
                            Rcpp::Named("covariance") = cov_new,
                            Rcpp::Named("mu") = mu_new,
                            Rcpp::Named("ar") = ar_new);
}

//' Log-likelihood of copula based model for a multivariate NMA.
//'
//' Evaluation of the log-likelihood.
//'
//' @param data pippo
//' @param init pippo
//' @param totiter pippo
//' @param prior pippo
//' @param prop pippo
//' @param verbose pippo
//'
//' @return A list containing the output of the MCMC simulation.
//' @export
//'
//' @examples
//' # nothing for now!
// [[Rcpp::export]]
Rcpp::List nc_mcmc_mh(const Rcpp::RObject& data, const Rcpp::List& init, const int& totiter, const Rcpp::List& prior, const Rcpp::List& prop, const Rcpp::List& tuning, const Rcpp::List& adapt, const bool& verbose) {

  int niter = 0, adapt_k = 0, n_datapoints = data.slot("n_datapoints"), M = data.slot("n_outcomes"), n_trt = data.slot("n_treatments"), n_study = data.slot("n_study"), nn = M*(M + 1)/2;

  double eps = Rcpp::as<double>(tuning["eps"]), eps_ab = Rcpp::as<double>(tuning["eps_ab"]);//, eps_a = 0.0, eps_b = 0.0;

  int print_every = 500, every = Rcpp::as<int>(adapt["every"]), maxiter = Rcpp::as<int>(adapt["maxiter"]); //, miniter = Rcpp::as<int>(adapt["miniter"]);

  double alpha = Rcpp::as<double>(adapt["alpha"]), beta = Rcpp::as<double>(adapt["beta"]), gamma = Rcpp::as<double>(adapt["gamma"]), tar = Rcpp::as<double>(adapt["tar"]), tol = Rcpp::as<double>(adapt["tol"]);

  Rcpp::DataFrame study_data(data.slot("study_data"));
  Rcpp::DataFrame study_id(data.slot("study_id"));
  Rcpp::IntegerVector cols_y(M), cols_n(M);
  for (unsigned int i = 0; i < M; i++) {
    cols_y(i) = 6 + M + i;
    cols_n(i) = 6 + i;
  }
  Rcpp::List Gamma_init = Rcpp::as<Rcpp::List>(init["Gamma"]);
  Rcpp::List Gamma = Rcpp::clone(Gamma_init);
  Rcpp::IntegerVector trt = study_data["trt"];
  int ref_trt = data.slot("ref_trt");
  Rcpp::IntegerVector narms = study_id["narms"];
  Rcpp::IntegerVector narms_study = study_data["narms"];
  Rcpp::IntegerVector study_index = study_data["study"];

  Rcpp::NumericMatrix mu_init = init["mu"], delta_init = init["delta"], x_init = init["x"];
  Rcpp::NumericMatrix mu = Rcpp::clone(mu_init), delta = Rcpp::clone(delta_init), x = Rcpp::clone(x_init);
  arma::mat Sigma_M = Rcpp::as<arma::mat>(init["Sigma_M"]), d = Rcpp::as<arma::mat>(init["d"]);
  arma::mat delta_arma = Rcpp::as<arma::mat>(delta), delta_arma_prop = Rcpp::as<arma::mat>(delta);
  Rcpp::NumericMatrix y = df_nm(study_data, cols_y);
  Rcpp::NumericMatrix n_data = df_nm(study_data, cols_n);
  Rcpp::NumericMatrix n_imp = n_imputed(n_data);
  Rcpp::NumericMatrix x_imp(n_datapoints, M), y_imp(n_datapoints, M);
  arma::mat x_imp_arma(n_datapoints, M), x_star_arma(n_datapoints, M);

  Rcpp::List D_list_init = Rcpp::as<Rcpp::List>(init["D"]);
  Rcpp::List D_list = Rcpp::clone(D_list_init);
  arma::mat D(M, M), D_inv(M, M), S_q = arma::zeros<arma::mat>(M, M), D_prop(M, M);
  Rcpp::List x2t_split(n_trt);
  arma::mat x_q;
  int n_q;

  int nGamma = prior["nGamma"];
  std::string Gamma_update = prop["Gamma_update"];
  double eta_prior = prior["eta_prior"];
  double eta_prop = prop["eta_prop"];
  arma::mat Gamma_q_prop(M, M), D_prop_inv(M, M), Gamma_q_curr(M, M), Sigma_q_prop(M, M);
  double ran_unif, Gamma_rate, mu_rate, delta_rate, d_rate, Sigma_M_rate;
  double Gamma_A = 0.0, Gamma_B = 0.0, target_Gamma_prop = 0.0, target_Gamma_curr = 0.0;
  arma::vec accept_Gamma_q = arma::zeros<arma::vec>(nGamma);
  arma::cube D_chain(nGamma, M, totiter), S_q_chain(nGamma, nn, totiter), Sigma_q_prop_chain(nGamma, nn, totiter);
  arma::cube Gamma_chain(nGamma, M*(M - 1)/2, totiter);

  arma::mat Gamma_k(M, M), Gamma_k_m_m(1, M), Gamma_k_m(M, M), x_ik_m(1, M);
  double a_ikm = 0.0, b_ikm = 0.0, w_ikm = 0.0, gamma_ikm = 0.0, theta_ikm = 0.0, p_ikm = 0.0;
  arma::cube a_chain(n_datapoints, M, totiter), b_chain(n_datapoints, M, totiter), x_chain(n_datapoints, M, totiter);//, x_adj_chain(n_datapoints, M, totiter), x_unadj_chain(n_datapoints, M, totiter);
  arma::vec D_tmp(M);
  arma::mat Gamma_tmp(M, M), Gamma_inv(M, M), Gamma_new(M, M), S(M, M), W(n_datapoints, M);

  double mu_sigma = sqrt(Rcpp::as<double>(prior["mu_sigma2"]));
  arma::vec mu_coef(3), delta_coef(3);
  Rcpp::List args_mu = Rcpp::List::create(
      Rcpp::Named("delta") = 0.0,   // temporary value
      Rcpp::Named("y") = 0,         // temporary value
      Rcpp::Named("n") = 1,         // temporary value
      Rcpp::Named("w") = 0.0,       // temporary value
      Rcpp::Named("gamma") = 1.0,   // temporary value
      Rcpp::Named("mu_sigma") = mu_sigma,
      Rcpp::Named("eps") = eps,
      Rcpp::Named("eps_ab") = eps_ab);
  Rcpp::List args_delta = Rcpp::List::create(
      Rcpp::Named("mu") = 0,        // temporary value
      Rcpp::Named("tau") = 0.0,     // temporary value
      Rcpp::Named("eta") = 1.0,     // temporary value
      Rcpp::Named("y") = 0,         // temporary value
      Rcpp::Named("n") = 1,         // temporary value
      Rcpp::Named("w") = 0.0,       // temporary value
      Rcpp::Named("gamma") = 1.0,   // temporary value
      Rcpp::Named("eps") = eps,
      Rcpp::Named("eps_ab") = eps_ab);
  // double xmin = Rcpp::as<double>(tuning["xmin"]), xmax = Rcpp::as<double>(tuning["xmax"]);

  arma::vec mode_mu(1), mode_delta(1);
  arma::mat var_mu(1, 1), var_delta(1, 1);
  double tau_ikm = 0.0, eta_ikm = 0.0;
  arma::mat Sigma_M_m_m(1, M), Sigma_M_m(M, M), delta_ik_m(1, M);
  arma::mat d_1(1, M), d_k(1, M), d_1k(1, M), d_1k_m(1, M);
  Rcpp::IntegerVector baseline = study_data["baseline"];
  Rcpp::IntegerVector baseline_id = study_id["baseline"];
  Rcpp::NumericMatrix mu_long = param_long(mu, narms, false);

  arma::mat rho_mu(n_study, M);
  rho_mu.fill(2.38/sqrt(M*n_study));
  arma::mat cov_mu(n_study, M, arma::fill::ones), cov_mu_mat(1, 1);
  double mu_curr = 0.0, mu_prop = 0.0;
  // arma::mat mu_curr_arma(1, 1), mu_prop_arma(1, 1);
  double mu_A = 0.0, mu_B = 0.0, target_mu_prop = 0.0, target_mu_curr = 0.0;
  int accept_mu = 0;
  arma::vec mu_mu_vec(1);
  mu_mu_vec(0) = 0.0;
  arma::mat mu_mu(n_study, M, arma::fill::zeros);
  arma::vec theta_mu(every);
  arma::cube ar_mu_vec(n_study, M, totiter), ar_mu(2, n_study, M);
  arma::vec ar_mu_sub(2);
  ar_mu_sub(0) = 0.0;
  for (unsigned int s = 0; s < n_study; s++) {
    for (unsigned int m = 0; m < M; m++) {
      ar_mu(0, s, m) = 0.0;
    }
  }
  Rcpp::List prop_mu(4);
  arma::cube cov_mu_chain(n_study, M, totiter), rho_mu_chain(n_study, M, totiter), mu_prop_chain(n_study, M, totiter), mu_rate_chain(n_study, M, totiter), tp_mu_chain(n_study, M, totiter), mu_chain(n_study, M, totiter);

  arma::mat rho_delta(n_datapoints, M);
  rho_delta.fill(2.38/sqrt(M*n_datapoints));
  arma::mat cov_delta(n_datapoints, M, arma::fill::ones), cov_delta_mat(1, 1);
  double delta_curr = 0.0, delta_prop = 0.0;
  // arma::mat delta_curr_arma(1, 1), delta_prop_arma(1, 1);
  double delta_A = 0.0, delta_B = 0.0, target_delta_prop = 0.0, target_delta_curr = 0.0;
  int accept_delta = 0;
  arma::vec mu_delta_vec(1);
  mu_delta_vec(0) = 0.0;
  arma::mat mu_delta(n_datapoints, M, arma::fill::zeros);
  arma::vec theta_delta(every);
  arma::cube ar_delta_vec(n_datapoints, M, totiter), ar_delta(2, n_datapoints, M);
  arma::vec ar_delta_sub(2);
  ar_delta_sub(0) = 0.0;
  for (unsigned int ik = 0; ik < n_datapoints; ik++) {
    for (unsigned int m = 0; m < M; m++) {
      ar_delta(0, ik, m) = 0.0;
    }
  }
  Rcpp::List prop_delta(4);
  arma::cube cov_delta_chain(n_datapoints, M, totiter), rho_delta_chain(n_datapoints, M, totiter), delta_prop_chain(n_datapoints, M, totiter), delta_chain(n_datapoints, M, totiter);

  int accept_d = 0;
  arma::cube d_chain(n_trt, M, totiter);

  int accept_Sigma_M = 0;
  arma::mat Sigma_M_chain(totiter, nn);
  arma::mat beta_chain(totiter, nn);

  double rho_d = 2.38/sqrt(M*(n_trt - 1)), sd_multplier = 1;
  arma::mat cov_d = arma::eye(M*(n_trt - 1), M*(n_trt - 1));
  arma::vec d_curr(M*(n_trt - 1));
  arma::mat d_prop(1, M*(n_trt - 1)), d_curr_ref(n_trt, M), d_prop_ref(n_trt, M);
  double d_sigma = sqrt(Rcpp::as<double>(prior["d_sigma2"]));
  double d_A = 0.0, d_B = 0.0, target_d_prop = 0.0, target_d_curr = 0.0;
  arma::vec mu_d(M*(n_trt - 1));
  arma::cube theta_d(n_trt, M, every);
  arma::mat theta_d_reshaped(every, M*(n_trt - 1));
  arma::vec ar_d_vec(totiter), ar_d(2);
  ar_d(0) = 0.0;
  Rcpp::List prop_d(4);

  double rho_beta = 2.38/sqrt(nn);
  arma::mat cov_beta = arma::eye(nn, nn);
  arma::mat Sigma_M_curr(Sigma_M), Sigma_M_prop(Sigma_M);
  arma::vec beta_curr(nn), beta_prop(nn);
  double sigma_r = Rcpp::as<double>(prior["beta_sigma"]);
  double Sigma_M_A = 0.0, Sigma_M_B = 0.0, target_Sigma_M_prop = 0.0, target_Sigma_M_curr = 0.0;
  arma::vec mu_beta(nn);
  arma::mat theta_beta(every, nn);
  arma::vec ar_Sigma_M_vec(totiter), ar_Sigma_M(2);
  ar_Sigma_M(0) = 0.0;
  Rcpp::List prop_beta(4);

  arma::vec loglik(totiter, arma::fill::zeros), logprior(totiter, arma::fill::zeros), logpost(totiter, arma::fill::zeros);

  arma::vec tmp(every);

  GetRNGstate();

  // Set the Armadillo seed from R's
  // int seed = (int)Rf_runif(0.0, 10000.0);
  // std::srand(seed);

  while (niter < totiter) {
    // if (niter > 0) {
    //   for (unsigned int ik = 0; ik < n_datapoints; ik++) {
    //     for (unsigned int m = 0; m < M; m++) {
    //       x(ik, m) = x_unadj_chain(ik, m, niter - 1);
    //     }
    //   }
    // }

    // imputing x and y variables
    // x_imp = x_imputed(x, Gamma, trt);
    x_imp = Rcpp::clone(x); // CHECK THAT THIS PART IS OK WHEN IMPUTING!!!
    x_imp_arma = Rcpp::as<arma::mat>(x_imp);
    y_imp = y_imputed(y, x_imp, narms, mu, delta, n_imp);

    // updating mu (study-specific baseline effects) and delta (study-specific
    // [random] treatment effects)
    // for (unsigned int ik = 0; ik < n_datapoints; ik++) {
    //   if (nGamma > 1) {
    //     Gamma_k = Rcpp::as<arma::mat>(Gamma(trt(ik) - 1));
    //   } else {
    //     Gamma_k = Rcpp::as<arma::mat>(Gamma(0));
    //   }
    //   d_1 = d.row(baseline(ik) - 1);
    //   d_k = d.row(trt(ik) - 1);
    //   d_1k = d_k - d_1;
    //   if (trt(ik) == baseline(ik)) {
    //     // update mu (study-specific baseline effect)
    //     for (unsigned int m = 0; m < M; m++) {
    //       Gamma_k_m_m = Gamma_k.row(m);
    //       Gamma_k_m_m.shed_col(m);
    //       Gamma_k_m = Gamma_k;
    //       Gamma_k_m.shed_col(m);
    //       Gamma_k_m.shed_row(m);
    //       x_ik_m = x_imp_arma.row(ik);
    //       x_ik_m.shed_col(m);
    //       w_ikm = arma::as_scalar(Gamma_k_m_m*arma::inv(Gamma_k_m)*arma::trans(x_ik_m));
    //       gamma_ikm = sqrt(arma::as_scalar(Gamma_k(m, m) - Gamma_k_m_m*arma::inv(Gamma_k_m)*arma::trans(Gamma_k_m_m)));

    //       mu_curr = mu_long(ik, m);
    //       // mu_curr_arma(0, 0) = mu_curr;
    //       mode_mu(0) = mu_curr;
    //       var_mu(0, 0) = pow(rho_mu(study_index(ik) - 1, m), 2)*cov_mu(study_index(ik) - 1, m);
    //       mu_prop = rmvt_arma(1, mode_mu, var_mu, 7)(0, 0);
    //       mu_prop_chain(study_index(ik) - 1, m, niter) = mu_prop;
    //       rho_mu_chain(study_index(ik) - 1, m, niter) = rho_mu(study_index(ik) - 1, m);
    //       cov_mu_chain(study_index(ik) - 1, m, niter) = cov_mu(study_index(ik) - 1, m);
    //       // mu_prop_arma(0, 0) = mu_prop;

    //       args_mu["delta"] = delta(ik, m);
    //       args_mu["y"] = y_imp(ik, m);
    //       args_mu["n"] = n_imp(ik, m);
    //       args_mu["w"] = w_ikm;
    //       args_mu["gamma"] = gamma_ikm;
    //       target_mu_curr = mu_logpost(mu_curr, delta(ik, m), static_cast<double>(y_imp(ik, m)), static_cast<double>(n_imp(ik, m)), w_ikm, gamma_ikm, mu_sigma, eps, eps_ab);
    //       // if (ISNA(target_mu_curr)) {
    //       //   mu_coef = ols_coef(xmin, xmax, args_mu, false);
    //       //   if (!ISNA(mu_coef(0))) {
    //       //     target_mu_curr = ols_pred(mu_coef, mu_curr);
    //       //   } else {
    //       //     target_mu_curr = NA_REAL;
    //       //   }
    //       // }
    //       target_mu_prop = mu_logpost(mu_prop, delta(ik, m), static_cast<double>(y_imp(ik, m)), static_cast<double>(n_imp(ik, m)), w_ikm, gamma_ikm, mu_sigma, eps, eps_ab);
    //       // if (ISNA(target_mu_prop)) {
    //       //   mu_coef = ols_coef(xmin, xmax, args_mu, false);
    //       //   if (!ISNA(mu_coef(0))) {
    //       //     target_mu_prop = ols_pred(mu_coef, mu_prop);
    //       //   } else {
    //       //     target_mu_prop = NA_REAL;
    //       //   }
    //       // }
    //       tp_mu_chain(study_index(ik) - 1, m, niter) = target_mu_prop;

    //       mu_A = target_mu_prop - target_mu_curr;
    //       // mu_B dovrebbe essere la differenza della proposal distribution, ma
    //       // e' simmetrica, quindi sparisce; la prior e' gia' dentro la target
    //       mu_B = 0.0; //dmvt_arma(mu_curr_arma, mode_mu, var_mu, 7, true)(0) - dmvt_arma(mu_prop_arma, mode_mu, var_mu, 7, true)(0);
    //       mu_rate = exp(mu_A + mu_B);
    //       mu_rate_chain(study_index(ik) - 1, m, niter) = mu_rate;
    //       ran_unif = R::runif(0.0, 1.0);
    //       if (!ISNAN(target_mu_prop) && !ISNAN(target_mu_curr)) {
    //         if (ran_unif < mu_rate) {
    //           mu_long(ik, m) = mu_prop;
    //           mu = param_wide(mu_long, narms, trt, baseline);
    //           mu_long = param_long(mu, narms, false); // is this necessary?
    //           accept_mu++;
    //         }
    //       } else {
    //         mu_rate = 0;
    //       }
    //       ar_mu_vec(study_index(ik) - 1, m, niter) = fmin(mu_rate, 1);
    //     }
    //   } else {
    //     // update delta (study-specific [random] treatment effects)
    //     for (unsigned int m = 0; m < M; m++) {
    //       Gamma_k_m_m = Gamma_k.row(m);
    //       Gamma_k_m_m.shed_col(m);
    //       Gamma_k_m = Gamma_k;
    //       Gamma_k_m.shed_col(m);
    //       Gamma_k_m.shed_row(m);
    //       x_ik_m = x_imp_arma.row(ik);
    //       x_ik_m.shed_col(m);
    //       w_ikm = arma::as_scalar(Gamma_k_m_m*arma::inv(Gamma_k_m)*arma::trans(x_ik_m));
    //       gamma_ikm = sqrt(arma::as_scalar(Gamma_k(m, m) - Gamma_k_m_m*arma::inv(Gamma_k_m)*arma::trans(Gamma_k_m_m)));

    //       Sigma_M_m_m = Sigma_M.row(m);
    //       Sigma_M_m_m.shed_col(m);
    //       Sigma_M_m = Sigma_M;
    //       Sigma_M_m.shed_col(m);
    //       Sigma_M_m.shed_row(m);
    //       delta_ik_m = delta_arma.row(ik);
    //       delta_ik_m.shed_col(m);
    //       d_1k_m = d_1k;
    //       d_1k_m.shed_col(m);
    //       tau_ikm = d_1k(m) + arma::as_scalar(Sigma_M_m_m*arma::inv(Sigma_M_m)*arma::trans(delta_ik_m - d_1k_m));
    //       eta_ikm = sqrt(arma::as_scalar(Sigma_M(m, m) - Sigma_M_m_m*arma::inv(Sigma_M_m)*arma::trans(Sigma_M_m_m)));

    //       delta_curr = delta(ik, m);
    //       // delta_curr_arma(0, 0) = delta_curr;
    //       mode_delta(0) = delta_curr;
    //       var_delta(0, 0) = pow(rho_delta(ik, m), 2)*cov_delta(ik, m);
    //       delta_prop = rmvt_arma(1, mode_delta, var_delta, 7)(0, 0);
    //       delta_prop_chain(ik, m, niter) = delta_prop;
    //       rho_delta_chain(ik, m, niter) = rho_delta(ik, m);
    //       cov_delta_chain(ik, m, niter) = cov_delta(ik, m);
    //       // delta_prop_arma(0, 0) = delta_prop;

    //       args_delta["mu"] = mu_long(ik, m);
    //       args_delta["tau"] = tau_ikm;
    //       args_delta["eta"] = eta_ikm;
    //       args_delta["y"] = y_imp(ik, m);
    //       args_delta["n"] = n_imp(ik, m);
    //       args_delta["w"] = w_ikm;
    //       args_delta["gamma"] = gamma_ikm;
    //       target_delta_curr = delta_logpost(delta_curr, mu_long(ik, m), tau_ikm, eta_ikm, static_cast<double>(y_imp(ik, m)), static_cast<double>(n_imp(ik, m)), w_ikm, gamma_ikm, eps, eps_ab);
    //       // if (ISNA(target_delta_curr)) {
    //       //   delta_coef = ols_coef(xmin, xmax, args_delta, true);
    //       //   if (!ISNA(delta_coef(0))) {
    //       //     target_delta_curr = ols_pred(delta_coef, delta_curr);
    //       //   } else {
    //       //     target_delta_curr = NA_REAL;
    //       //   }
    //       // }
    //       target_delta_prop = delta_logpost(delta_prop, mu_long(ik, m), tau_ikm, eta_ikm, static_cast<double>(y_imp(ik, m)), static_cast<double>(n_imp(ik, m)), w_ikm, gamma_ikm, eps, eps_ab);
    //       // if (ISNA(target_delta_prop)) {
    //       //   delta_coef = ols_coef(xmin, xmax, args_delta, true);
    //       //   if (!ISNA(delta_coef(0))) {
    //       //     target_delta_prop = ols_pred(delta_coef, delta_prop);
    //       //   } else {
    //       //     target_delta_prop = NA_REAL;
    //       //   }
    //       // }

    //       delta_A = target_delta_prop - target_delta_curr;
    //       // delta_B dovrebbe essere la differenza della proposal distribution,
    //       // ma e' simmetrica, quindi sparisce; la prior e' gia' dentro la
    //       // target
    //       delta_B = 0.0; //dmvt_arma(delta_curr_arma, mode_delta, var_delta, 7, true)(0) - dmvt_arma(delta_prop_arma, mode_delta, var_delta, 7, true)(0);
    //       delta_rate = exp(delta_A + delta_B);
    //       ran_unif = R::runif(0.0, 1.0);
    //       if (!ISNAN(target_delta_prop) && !ISNAN(target_delta_curr)) {
    //         if (ran_unif < delta_rate) {
    //           delta(ik, m) = delta_prop;
    //           accept_delta++;
    //         }
    //       } else {
    //         delta_rate = 0;
    //       }
    //       ar_delta_vec(ik, m, niter) = fmin(delta_rate, 1);
    //       // if (ISNAN(delta_rate)) {
    //       //   Rprintf("%4.10f\n", fmin(delta_rate, 1));
    //       // }
    //     }
    //   }
    // }
    mu_chain.slice(niter) = Rcpp::as<arma::mat>(mu);
    delta_chain.slice(niter) = Rcpp::as<arma::mat>(delta);
    delta_arma = Rcpp::as<arma::mat>(delta);

    // updating x (latent variables)
    for (unsigned int ik = 0; ik < n_datapoints; ik++) {
      if (nGamma > 1) {
        Gamma_k = Rcpp::as<arma::mat>(Gamma(trt(ik) - 1));
      } else {
        Gamma_k = Rcpp::as<arma::mat>(Gamma(0));
      }
      for (unsigned int m = 0; m < M; m++) {
        // if (!ISNAN(y(ik, m))) {
          theta_ikm = mu_long(ik, m) + delta(ik, m);
          p_ikm = expit_double(theta_ikm);
          a_ikm = R::pbinom(y_imp(ik, m) - 1, n_imp(ik, m), p_ikm, 1, 1); // log
          b_ikm = R::pbinom(y_imp(ik, m), n_imp(ik, m), p_ikm, 1, 1); // log
          a_chain(ik, m, niter) = a_ikm;
          b_chain(ik, m, niter) = b_ikm;
          // eps_a = a_ikm/2.0;
          // eps_b = (1.0 - b_ikm)/2.0;
          if (a_ikm == b_ikm) {
            // if (a_ikm == 0) {
            //   b_ikm += eps;
            // } else if (b_ikm == 1) {
            //   a_ikm -= eps;
            // } else {
            //   a_ikm -= fmin(eps/2.0, eps_a);
            //   b_ikm += fmin(eps/2.0, eps_b);
            // }
            Rprintf("niter = %d a = %4.5f - b = %4.5f - y = %4.0f - n = %4.0f - p = %4.10f\n", niter, a_ikm, b_ikm, y_imp(ik, m), n_imp(ik, m), p_ikm);
          } else {
            Gamma_k_m_m = Gamma_k.row(m);
            Gamma_k_m_m.shed_col(m);
            Gamma_k_m = Gamma_k;
            Gamma_k_m.shed_col(m);
            Gamma_k_m.shed_row(m);
            x_ik_m = x_imp_arma.row(ik);
            x_ik_m.shed_col(m);
            w_ikm = arma::as_scalar(Gamma_k_m_m*arma::inv(Gamma_k_m)*arma::trans(x_ik_m));
            gamma_ikm = sqrt(arma::as_scalar(Gamma_k(m, m) - Gamma_k_m_m*arma::inv(Gamma_k_m)*arma::trans(Gamma_k_m_m)));
            x(ik, m) = rtruncnorm_rcpp(1, R::qnorm(a_ikm, 0.0, 1.0, 1, 1), R::qnorm(b_ikm, 0.0, 1.0, 1, 1), w_ikm, gamma_ikm)[0];
            x_imp(ik, m) = x(ik, m);
            x_imp_arma(ik, m) = x(ik, m);
          }
        // }
      }
    }
    x_chain.slice(niter) = Rcpp::as<arma::mat>(x);
    // x_unadj_chain.slice(niter) = Rcpp::as<arma::mat>(x);

    // updating Gamma (correlation of latent variables)
    if (Gamma_update == "PX-RPMH") {
      // if (niter > 0) {
      //   for (unsigned int ik = 0; ik < n_datapoints; ik++) {
      //     for (unsigned int m = 0; m < M; m++) {
      //       x(ik, m) = x_adj_chain(ik, m, niter - 1);
      //     }
      //   }
      //   x_imp = Rcpp::clone(x); // CHECK THAT THIS PART IS OK WHEN IMPUTING!!!
      //   x_imp_arma = Rcpp::as<arma::mat>(x_imp);
      // }

      if (nGamma > 1) {
        x2t_split = split_nm(x_imp, trt);
        for (unsigned q = 0; q < n_trt; q++) {
          x_q = Rcpp::as<arma::mat>(x2t_split[q]);
          n_q = x_q.n_rows;
          D = Rcpp::as<arma::mat>(D_list(q));
          // D_inv = arma::inv_sympd(D);
          for (unsigned int k = 0; k < n_q; k++) {
            S_q += D*arma::trans(x_q.row(k))*x_q.row(k)*D;
          }
          if (!is_positive_definite(S_q, 1423)) {
            Rprintf("S_q - nGamma > 1\n");
            S_q = make_positive_definite(S_q);
          }
          S_q_chain(arma::span(q), arma::span::all, arma::span(niter)) = diag_tri(S_q);
          // in the following, degrees of freedom are set to (n_q > M) ? n_q : M)
          // otherwise, when n_q < M the IW routine returns an error
          // (indeterminate system)
          Sigma_q_prop = rinvwish_arma(((n_q > M) ? n_q : M), S_q);
          Sigma_q_prop_chain(arma::span(q), arma::span::all, arma::span(niter)) = diag_tri(Sigma_q_prop);
          D_prop = arma::sqrt(arma::diagmat(Sigma_q_prop));
          D_prop_inv = arma::inv_sympd(D_prop);
          Gamma_q_prop = cov2cor_rcpp(Sigma_q_prop);

          Gamma_q_curr = Rcpp::as<arma::mat>(Gamma(q));
          Gamma_rate = 0.5*(M + 1)*(log(arma::det(Gamma_q_prop)) - log(arma::det(Gamma_q_curr)));
          ran_unif = R::runif(0.0, 1.0);
          if (ran_unif < exp(Gamma_rate)) {
            Gamma_q_curr = Gamma_q_prop;
            D = D_prop;
            accept_Gamma_q(q)++;
          }
          Gamma(q) = Gamma_q_curr;
          D_list(q) = D;

          S_q = arma::zeros<arma::mat>(M, M);
          D_chain(arma::span(q), arma::span::all, arma::span(niter)) = arma::diagvec(D);
        }
      } else {
        D = arma::diagmat(arma::sqrt(1/arma::sum(arma::pow(x_imp_arma, 2), 0)));
        x_star_arma = x_imp_arma * D;
        S_q = arma::trans(x_star_arma)*x_star_arma;
        // for (unsigned int k = 0; k < n_datapoints; k++) {
        //   S_q += D*arma::trans(x_imp_arma.row(k))*x_imp_arma.row(k)*D;
        // }
        if (!is_positive_definite(S_q, 1456)) {
          Rprintf("S_q - nGamma = 1\n");
          S_q = make_positive_definite(S_q);
        }
        S_q_chain(arma::span(0), arma::span::all, arma::span(niter)) = diag_tri(S_q);
        // in the following, degrees of freedom are set to (n_datapoints > M) ?
        // n_datapoints : M) otherwise, when n_datapoints < M the IW routine
        // returns an error (indeterminate system)
        Sigma_q_prop = rinvwish_arma(((n_datapoints > M) ? n_datapoints : M), S_q);
        Sigma_q_prop_chain(arma::span(0), arma::span::all, arma::span(niter)) = diag_tri(Sigma_q_prop);
        D_prop = arma::sqrt(arma::diagmat(Sigma_q_prop));
        Gamma_q_prop = cov2cor_rcpp(Sigma_q_prop);

        Gamma_q_curr = Rcpp::as<arma::mat>(Gamma(0));
        Gamma_rate = 0.5*(M + 1)*(log(arma::det(Gamma_q_prop)) - log(arma::det(Gamma_q_curr)));
        ran_unif = R::runif(0.0, 1.0);
        if (ran_unif < exp(Gamma_rate)) {
          // adjust latent variables x according to the new scales in 'D_prop'
          D_prop_inv = arma::inv_sympd(D_prop);
          x_imp_arma = x_star_arma * D_prop_inv;
          // x_imp_arma = x_imp_arma * D_prop_inv;
          for (unsigned int ik = 0; ik < n_datapoints; ik++) {
            for (unsigned int m = 0; m < M; m++) {
              x_imp(ik, m) = x_imp_arma(ik, m);
              x(ik, m) = x_imp_arma(ik, m);
            }
          }

          D = D_prop;
          Gamma_q_curr = Gamma_q_prop;
          accept_Gamma_q(0)++;
        }
        Gamma(0) = Gamma_q_curr;

        S_q = arma::zeros<arma::mat>(M, M);
      }
    } else if (Gamma_update == "IMH") {
      if (nGamma > 1) {
        x2t_split = split_nm(x_imp, trt);
        for (unsigned q = 0; q < n_trt; q++) {
          x_q = Rcpp::as<arma::mat>(x2t_split[q]);
          Gamma_q_curr = Rcpp::as<arma::mat>(Gamma(q));
          Gamma_q_prop = rlkj_arma(M, eta_prop);

          target_Gamma_curr = Gamma_logpost(Gamma_q_curr, x_q, eta_prior);
          target_Gamma_prop = Gamma_logpost(Gamma_q_prop, x_q, eta_prior);
          Gamma_A = target_Gamma_prop - target_Gamma_curr;
          Gamma_B = dlkj_arma(Gamma_q_curr, eta_prop, true) - dlkj_arma(Gamma_q_prop, eta_prop, true);
          Gamma_rate = exp(Gamma_A + Gamma_B);
          ran_unif = R::runif(0.0, 1.0);
          if (ran_unif < exp(Gamma_rate)) {
            Gamma_q_curr = Gamma_q_prop;
            accept_Gamma_q(q)++;
          }
          Gamma(q) = Gamma_q_curr;
        }
      } else {
        Gamma_q_curr = Rcpp::as<arma::mat>(Gamma(0));
        Gamma_q_prop = rlkj_arma(M, eta_prop);

        target_Gamma_curr = Gamma_logpost(Gamma_q_curr, x_imp_arma, eta_prior);
        target_Gamma_prop = Gamma_logpost(Gamma_q_prop, x_imp_arma, eta_prior);
        Gamma_A = target_Gamma_prop - target_Gamma_curr;
        Gamma_B = dlkj_arma(Gamma_q_curr, eta_prop, true) - dlkj_arma(Gamma_q_prop, eta_prop, true);
        Gamma_rate = exp(Gamma_A + Gamma_B);
        ran_unif = R::runif(0.0, 1.0);
        if (ran_unif < exp(Gamma_rate)) {
          Gamma_q_curr = Gamma_q_prop;
          accept_Gamma_q(0)++;
        }
        Gamma(0) = Gamma_q_curr;
      }
    } else if (Gamma_update == "Talhouketal") {
      if (nGamma > 1) {
        Rprintf("STILL TO DO!\n");
      } else {
        Gamma_tmp = Rcpp::as<arma::mat>(Gamma(0));
        // Gamma_tmp.print();
        Gamma_inv = arma::inv_sympd(Gamma_tmp);
        // Gamma_inv.print();
        // Rprintf("\n");
        for (unsigned int m = 0; m < M; m++) {
          D_tmp(m) = rinvgamma_rcpp(1, (M + 1)/2.0, Gamma_inv(m, m)/2.0)(0);
        }
        D = arma::diagmat(arma::sqrt(D_tmp));
        W = x_imp_arma * D;
        S = arma::trans(W) * W + arma::eye(M, M);
        Gamma_new = rinvwish_arma(2 + n_datapoints, S); // the MATLAB code reports (n + M + 1) + M + 1 as the degrees of freedom
        D = arma::sqrt(arma::diagmat(Gamma_new));
        Gamma(0) = cov2cor_rcpp(Gamma_new);
      }
    } else {
      error("the specified update method for the Gamma parameter is not available. use either 'IMH' or 'PX-RPMH'.");
    }
    D_chain(arma::span(0), arma::span::all, arma::span(niter)) = arma::diagvec(D);
    Gamma_chain.slice(niter) = list_mat(Gamma);
    x_chain.slice(niter) = Rcpp::as<arma::mat>(x);
    // x_adj_chain.slice(niter) = Rcpp::as<arma::mat>(x);

    // updating d (pooled treatment effects across trials)
    // d_curr = mat_to_vec(d, true, ref_trt);
    // d_curr_ref = d;
    // d_prop = rmvt_arma(1, d_curr, pow(sd_multplier*rho_d, 2)*cov_d, 7);
    // d_prop_ref = vec_to_mat(arma::vectorise(d_prop), M, true, ref_trt);
    // target_d_curr = d_logpost(d_curr_ref, delta_arma, Sigma_M, trt, baseline, narms_study, d_sigma, ref_trt);
    // target_d_prop = d_logpost(d_prop_ref, delta_arma, Sigma_M, trt, baseline, narms_study, d_sigma, ref_trt);
    // d_A = target_d_prop - target_d_curr;
    // // d_B dovrebbe essere la differenza della proposal distribution, ma e'
    // // simmetrica, quindi sparisce; la prior e' gia' dentro la target
    // d_B = 0.0; //d_logprior(d_prop_ref, d_sigma, ref_trt) - d_logprior(d_curr_ref, d_sigma, ref_trt);
    // d_rate = exp(d_A + d_B);
    // ran_unif = R::runif(0.0, 1.0);
    // if (ran_unif < d_rate) {
    //   d = d_prop_ref;
    //   accept_d++;
    // }
    ar_d_vec(niter) = fmin(d_rate, 1);
    d_chain.slice(niter) = d;

    // updating Sigma_M (common between-study covariance structure)
    // Sigma_M_curr = Sigma_M;
    // beta_curr = Sigma_M_to_beta(Sigma_M_curr);
    // beta_prop = arma::vectorise(rmvt_arma(1, beta_curr, pow(sd_multplier*rho_beta, 2)*cov_beta, 7));
    // // if (niter == 0) {
    // //   beta_prop.print();
    // // }
    // Sigma_M_prop = beta_to_Sigma_M(beta_prop, M);
    // if (!is_positive_definite(Sigma_M_prop, 1749)) {
    //   Rprintf("Sigma_M_prop\n");
    //   Sigma_M_prop = make_positive_definite(Sigma_M_prop);
    // }
    // target_Sigma_M_prop = Sigma_M_logpost(d, delta_arma, Sigma_M_prop, trt, baseline, narms_study, sigma_r);
    // target_Sigma_M_curr = Sigma_M_logpost(d, delta_arma, Sigma_M_curr, trt, baseline, narms_study, sigma_r);
    // Sigma_M_A = target_Sigma_M_prop - target_Sigma_M_curr;
    // // Sigma_M_B dovrebbe essere la differenza della proposal distribution, ma
    // // e' simmetrica, quindi sparisce; la prior invece e' gia' dentro la target
    // Sigma_M_B = 0.0; //dlogchol_arma(Sigma_M_prop, sigma_r, true) - dlogchol_arma(Sigma_M_curr, sigma_r, true);
    // Sigma_M_rate = exp(Sigma_M_A + Sigma_M_B);
    // ran_unif = R::runif(0.0, 1.0);
    // if (ran_unif < Sigma_M_rate) {
    //   Sigma_M = Sigma_M_prop;
    //   accept_Sigma_M++;
    // }
    ar_Sigma_M_vec(niter) = fmin(Sigma_M_rate, 1);
    Sigma_M_chain.row(niter) = diag_tri(Sigma_M);
    beta_chain.row(niter) = arma::trans(Sigma_M_to_beta(Sigma_M));

    // adaptation step
    // if (((adapt_k + 1) <= maxiter) && (((niter + 1) % every) == 0)) {
    //   if (verbose) {
    //     Rprintf("      --> performing adaptation [step %d/%d]\n", adapt_k + 1, maxiter);
    //   }

    //   // mu proposal parameters
    //   for (unsigned int s = 0; s < n_study; s++) {
    //     for (unsigned int m = 0; m < M; m++) {
    //       if (fabs(ar_mu(0, s, m) - tar) > tol) {
    //         theta_mu = mu_chain.subcube(s, m, niter - every + 1, s, m, niter);
    //         if (adapt_k == 0) {
    //           mu_mu(s, m) = arma::mean(theta_mu);
    //         }
    //         mu_mu_vec(0) = mu_mu(s, m);
    //         tmp = ar_mu_vec.subcube(s, m, niter - every + 1, s, m, niter);
    //         ar_mu(1, s, m) = arma::as_scalar(arma::mean(tmp));
    //         ar_mu_sub = ar_mu.subcube(0, s, m, 1, s, m);
    //         cov_mu_mat(0, 0) = cov_mu(s, m);
    //         prop_mu = rwmh_adapt(theta_mu, mu_mu_vec, rho_mu(s, m), cov_mu_mat, ar_mu_sub, alpha, beta, gamma, tar, adapt_k, false, 1);
    //         rho_mu(s, m) = Rcpp::as<double>(prop_mu["rho"]);
    //         cov_mu_mat = Rcpp::as<arma::mat>(prop_mu["covariance"]);
    //         if (!is_positive_definite(cov_mu_mat, 1791)) {
    //           Rprintf("cov_mu_mat\n");
    //           cov_mu_mat = make_positive_definite(cov_mu_mat);
    //         }
    //         cov_mu(s, m) = cov_mu_mat(0, 0);
    //         mu_mu_vec = Rcpp::as<arma::vec>(prop_mu["mu"]);
    //         mu_mu(s, m) = mu_mu_vec(0);
    //         ar_mu(0, s, m) = Rcpp::as<double>(prop_mu["ar"]);
    //       }
    //     }
    //   }

    //   // delta proposal parameters
    //   for (unsigned int ik = 0; ik < n_datapoints; ik++) {
    //     if (trt(ik) != baseline(ik)) {
    //       for (unsigned int m = 0; m < M; m++) {
    //         if (fabs(ar_delta(0, ik, m) - tar) > tol) {
    //           theta_delta = delta_chain.subcube(ik, m, niter - every + 1, ik, m, niter);
    //           if (adapt_k == 0) {
    //             mu_delta(ik, m) = arma::mean(theta_delta);
    //           }
    //           mu_delta_vec(0) = mu_delta(ik, m);
    //           tmp = ar_delta_vec.subcube(ik, m, niter - every + 1, ik, m, niter);
    //           ar_delta(1, ik, m) = arma::as_scalar(arma::mean(tmp));
    //           ar_delta_sub(0) = ar_delta(0, ik, m);
    //           ar_delta_sub(1) = ar_delta(1, ik, m);
    //           cov_delta_mat(0, 0) = cov_delta(ik, m);
    //           prop_delta = rwmh_adapt(theta_delta, mu_delta_vec, rho_delta(ik, m), cov_delta_mat, ar_delta_sub, alpha, beta, gamma, tar, adapt_k, false, 2);
    //           rho_delta(ik, m) = Rcpp::as<double>(prop_delta["rho"]);
    //           cov_delta_mat = Rcpp::as<arma::mat>(prop_delta["covariance"]);
    //           if (!is_positive_definite(cov_delta_mat, 1820)) {
    //             Rprintf("cov_delta_mat\n");
    //             cov_delta_mat = make_positive_definite(cov_delta_mat);
    //           }
    //           cov_delta(ik, m) = cov_delta_mat(0, 0);
    //           mu_delta_vec = Rcpp::as<arma::vec>(prop_delta["mu"]);
    //           mu_delta(ik, m) = mu_delta_vec(0);
    //           ar_delta(0, ik, m) = Rcpp::as<double>(prop_delta["ar"]);
    //         }
    //       }
    //     }
    //   }

    //   // d proposal parameters
    //   if (fabs(ar_d(0) - tar) > tol) {
    //     theta_d = d_chain.slices(niter - every + 1, niter);
    //     theta_d_reshaped = cube_to_mat(theta_d, true, ref_trt);
    //     if (adapt_k == 0) {
    //       mu_d = arma::vectorise(arma::mean(theta_d_reshaped, 0)); // means by columns
    //     }
    //     ar_d(1) = arma::as_scalar(arma::mean(ar_d_vec.subvec(niter - every + 1, niter)));
    //     prop_d = rwmh_adapt(theta_d_reshaped, mu_d, rho_d, cov_d, ar_d, alpha, beta, gamma, tar, adapt_k, false, 3);
    //     rho_d = Rcpp::as<double>(prop_d["rho"]);
    //     cov_d = Rcpp::as<arma::mat>(prop_d["covariance"]);
    //     if (!is_positive_definite(cov_d, 1843)) {
    //       Rprintf("cov_d\n");
    //       cov_d = make_positive_definite(cov_d);
    //     }
    //     mu_d = Rcpp::as<arma::vec>(prop_d["mu"]);
    //     ar_d(0) = Rcpp::as<double>(prop_d["ar"]);
    //   }

    //   // log Cholesky betas (Sigma_M) proposal parameters
    //   if (fabs(ar_Sigma_M(0) - tar) > tol) {
    //     theta_beta = beta_chain.rows(niter - every + 1, niter);
    //     if (adapt_k == 0) {
    //       mu_beta = arma::vectorise(arma::mean(theta_beta, 0)); // means by columns
    //     }
    //     ar_Sigma_M(1) = arma::as_scalar(arma::mean(ar_Sigma_M_vec.subvec(niter - every + 1, niter)));
    //     prop_beta = rwmh_adapt(theta_beta, mu_beta, rho_beta, cov_beta, ar_Sigma_M, alpha, beta, gamma, tar, adapt_k, false, 4);
    //     rho_beta = Rcpp::as<double>(prop_beta["rho"]);
    //     cov_beta = Rcpp::as<arma::mat>(prop_beta["covariance"]);
    //     if (!is_positive_definite(cov_beta, 1861)) {
    //       Rprintf("cov_beta\n");
    //       cov_beta = make_positive_definite(cov_beta);
    //     }
    //     mu_beta = Rcpp::as<arma::vec>(prop_beta["mu"]);
    //     ar_Sigma_M(0) = Rcpp::as<double>(prop_beta["ar"]);
    //   }

    //   adapt_k++;
    // }

    // calculate the loglikelihood, logprior and logposterior
    // loglik(niter) = nc_loglik(y_imp, n_imp, x_imp, trt, mu, delta, Gamma); // c'e' il problema degli indicatori nella augmented likelihood
    // logprior(niter) = nc_logprior(mu, mu_sigma, d, d_sigma, Sigma_M, sigma_r, ref_trt);

    // print the information
    if ((((niter + 1) % print_every) == 0) && verbose) {
      Rprintf("   iter. %d/%d ==> d: %1.3f - Gamma: %1.3f - mu: %1.3f - delta: %1.3f - Sigma_M: %1.3f\n", (niter + 1), totiter, accept_d/(static_cast<double>(totiter)), arma::mean(accept_Gamma_q)/(static_cast<double>(totiter)), accept_mu/(static_cast<double>(totiter*n_study*M)), accept_delta/(static_cast<double>(totiter*(n_datapoints - n_study)*M)), accept_Sigma_M/(static_cast<double>(totiter)));
    }

    niter++;

    R_CheckUserInterrupt();
  }

  PutRNGstate();

  Rcpp::List accept = Rcpp::List::create(
    Rcpp::Named("d") = accept_d/static_cast<double>(totiter),
    Rcpp::Named("mu") = accept_mu/static_cast<double>(totiter*n_study*M),
    Rcpp::Named("delta") = accept_delta/static_cast<double>(totiter*(n_datapoints - n_study)*M),
    Rcpp::Named("Gamma") = accept_Gamma_q/static_cast<double>(totiter),
    Rcpp::Named("Sigma") = accept_Sigma_M/static_cast<double>(totiter));

  return Rcpp::List::create(Rcpp::Named("mu") = mu_chain,
                            Rcpp::Named("delta") = delta_chain,
                            Rcpp::Named("d") = d_chain,
                            Rcpp::Named("Sigma") = Sigma_M_chain,
                            Rcpp::Named("Gamma") = Gamma_chain,
                            // Rcpp::Named("x") = x_adj_chain,
                            Rcpp::Named("x") = x_chain,
                            // Rcpp::Named("x_unadj") = x_unadj_chain,
                            Rcpp::Named("x_unadj") = x_chain,
                            // Rcpp::Named("a") = a_chain,
                            // Rcpp::Named("b") = b_chain,
                            Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("logprior") = logprior,
                            Rcpp::Named("logpost") = (loglik + logprior),
                            Rcpp::Named("accept") = accept);
                            // Rcpp::Named("D") = D_chain,
                            // Rcpp::Named("S_q") = S_q_chain,
                            // Rcpp::Named("Sigma_q_prop") = Sigma_q_prop_chain);
                            // Rcpp::Named("rho_mu") = rho_mu_chain,
                            // Rcpp::Named("cov_mu") = cov_mu_chain,
                            // Rcpp::Named("ar_mu_vec") = ar_mu_vec,
                            // Rcpp::Named("mu_rate") = mu_rate_chain,
                            // Rcpp::Named("mu_prop") = mu_prop_chain,
                            // Rcpp::Named("tp_mu") = tp_mu_chain,
                            // Rcpp::Named("delta_prop") = delta_prop_chain);
}
