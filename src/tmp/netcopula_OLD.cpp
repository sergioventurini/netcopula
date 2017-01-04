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
  // int ind_a_b = 0;
  for (unsigned int ik = 0; ik < n_datapoints; ik++) {
    if (nGamma > 1) {
      Gamma_k = Rcpp::as<arma::mat>(Gamma(trt(ik) - 1));
    } else {
      Gamma_k = Rcpp::as<arma::mat>(Gamma(0));
    }
    for (unsigned int m = 0; m < M; m++) {
      x_ik(m) = x(ik, m);
      // ind_a_b = indic_a_b(y(ik, m), n(ik, m), x(ik, m), mu(ik, m), delta(ik, m));
      // if (!ind_a_b) {
      //   return NA_REAL;
      // }
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
  phi_inv_a = R::qnorm(R::pbinom(y_ikm - 1, n_ikm, pi, 1, 0), 0.0, 1.0, 1, 0);
  phi_inv_b = R::qnorm(R::pbinom(y_ikm, n_ikm, pi, 1, 0), 0.0, 1.0, 1, 0);
  if ((phi_inv_a <= x_ikm) && (x_ikm < phi_inv_b)) {
    out = 1;
  }

  return out;
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
      if (!is_positive_definite(Omega)) {
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
    if (!is_positive_definite(Omega)) {
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
double mu_delta_logpost(const double& mu, const double& delta, const double& tau, const double& eta, const double& y, const double& n, const double& w, const double& gamma, const double& mu_sigma, const double& eps, const double& eps_ab) {
  double theta = mu + delta;
  double p = expit_double(theta);

  double a = R::pbinom(y - 1, n, p, 1, 0);
  double b = R::pbinom(y, n, p, 1, 0);
  if (a > (1 - eps_ab)) {
    a = 1 - eps_ab;
  } else if (a < eps_ab) {
    a = eps_ab;
  }
  if (b > (1 - eps_ab)) {
    b = 1 - eps_ab;
  } else if (b < eps_ab) {
    b = eps_ab;
  }
  double eps_a = a/2.0, eps_b = (1.0 - b)/2.0;
  if (a == b) {
    if (a == 0) {
      b += eps;
    } else if (b == 1) {
      a -= eps;
    } else {
      a -= fmin(eps/2.0, eps_a);
      b += fmin(eps/2.0, eps_b);
    }
  }
  double phi_inv_a = R::qnorm(a, 0.0, 1.0, 1, 0);
  double phi_inv_b = R::qnorm(b, 0.0, 1.0, 1, 0);
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
  // if (!R_FINITE(tmp)) {
  //   double tmp_div = 700;
  //   double tmp_a = R::pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 0, 1);
  //   double tmp_b = R::pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 0, 1);
  //   double tmp_diff = tmp_a - tmp_b;
  //   double k = ((tmp_diff > tmp_div) ? (tmp_diff - tmp_div) : 0);
  //   tmp = k + log(exp(tmp_div) - exp(-k)) + tmp_b;
  // }

  if (!R_FINITE(tmp)) {
    return NA_REAL;
  }

  double lpost = tmp + R::dnorm(delta, tau, eta, 1) + R::dnorm(mu, 0, mu_sigma, 1);

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
          if (!is_positive_definite(Sigma_M_k)) {
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
          if (!is_positive_definite(Sigma_M_k)) {
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
Rcpp::List nc_mcmc(const Rcpp::RObject& data, const Rcpp::List& init, const int& totiter, const Rcpp::List& prior, const Rcpp::List& prop, const Rcpp::List& tuning, const Rcpp::List& adapt, const bool& verbose) {

  int niter = 0, adapt_k = 0, n_datapoints = data.slot("n_datapoints"), M = data.slot("n_outcomes"), n_trt = data.slot("n_treatments"), n_study = data.slot("n_study"), nn = M*(M + 1)/2;

  double eps = Rcpp::as<double>(tuning["eps"]), eps_ab = Rcpp::as<double>(tuning["eps_ab"]), eps_a = 0.0, eps_b = 0.0;

  int every = Rcpp::as<int>(adapt["every"]), maxiter = Rcpp::as<int>(adapt["maxiter"]); //, miniter = Rcpp::as<int>(adapt["miniter"]);

  double alpha = Rcpp::as<double>(adapt["alpha"]), beta = Rcpp::as<double>(adapt["beta"]), gamma = Rcpp::as<double>(adapt["gamma"]), tar = Rcpp::as<double>(adapt["tar"]), tol = Rcpp::as<double>(adapt["tol"]);

  Rcpp::DataFrame study_data(data.slot("study_data"));
  Rcpp::DataFrame study_id(data.slot("study_id"));
  Rcpp::IntegerVector cols_y(M), cols_n(M);
  for (unsigned int i = 0; i < M; i++) {
    cols_y(i) = 6 + M + i;
    cols_n(i) = 6 + i;
  }
  Rcpp::List Gamma = Rcpp::as<Rcpp::List>(init["Gamma"]);
  Rcpp::IntegerVector trt = study_data["trt"];
  int ref_trt = data.slot("ref_trt");
  Rcpp::IntegerVector narms = study_id["narms"];
  Rcpp::IntegerVector narms_study = study_data["narms"];

  Rcpp::NumericMatrix mu = init["mu"], delta = init["delta"], x = init["x"];
  arma::mat delta_arma = Rcpp::as<arma::mat>(delta), delta_arma_prop = Rcpp::as<arma::mat>(delta);
  Rcpp::NumericMatrix y = df_nm(study_data, cols_y);
  Rcpp::NumericMatrix n_data = df_nm(study_data, cols_n);
  Rcpp::NumericMatrix n_imp = n_imputed(n_data);
  Rcpp::NumericMatrix x_imp(n_datapoints, M), y_imp(n_datapoints, M);;
  arma::mat x_imp_arma(n_datapoints, M);

  Rcpp::List D_list = Rcpp::as<Rcpp::List>(init["D"]);
  arma::mat D(M, M), D_inv(M, M), S_q = arma::zeros<arma::mat>(M, M), D_prop(M, M);
  Rcpp::List x2t_split(n_trt);
  arma::mat x_q;
  int n_q;

  int nGamma = prior["nGamma"];
  arma::mat Gamma_q_prop(M, M), Gamma_q_prop_inv(M, M), Gamma_q_curr(M, M);
  double ran_unif, Gamma_rate, mu_delta_rate, d_rate, Sigma_M_rate;
  arma::vec accept_Gamma_q = arma::zeros<arma::vec>(nGamma);
  arma::cube Gamma_chain(nGamma, M*(M - 1)/2, totiter);

  Rcpp::NumericMatrix mu_long = param_long(mu, narms, false);
  arma::mat Gamma_k(M, M), Gamma_k_m_m(1, M), Gamma_k_m(M, M), x_ik_m(1, M);
  double a_ikm = 0.0, b_ikm = 0.0, w_ikm = 0.0, gamma_ikm = 0.0, theta_ikm = 0.0, p_ikm = 0.0;
  arma::cube x_chain(n_datapoints, M, totiter);

  Rcpp::Environment nc_env("package:netcopula");
  Rcpp::Function mu_delta_logpost_R = nc_env["mu_delta_logpost_func"];
  Rcpp::NumericVector par(2), parscale(2), ndeps(2);
  double mu_sigma = sqrt(Rcpp::as<double>(prior["mu_sigma2"]));
  for (unsigned int h = 0; h < 2; h++) {
    parscale(h) = 1;
    ndeps(h) = 0.001;
  }
  Rcpp::List args = Rcpp::List::create(
      Rcpp::Named("tau") = 0.0,     // temporary value
      Rcpp::Named("eta") = 1.0,     // temporary value
      Rcpp::Named("y") = 0,         // temporary value
      Rcpp::Named("n") = 1,         // temporary value
      Rcpp::Named("w") = 0.0,       // temporary value
      Rcpp::Named("gamma") = 1.0,   // temporary value
      Rcpp::Named("mu_sigma") = mu_sigma,
      Rcpp::Named("eps") = eps,
      Rcpp::Named("eps_ab") = eps_ab);
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
  Rcpp::List fit(4);
  arma::vec mode(2);
  mode.fill(0);
  arma::mat var(2, 2, arma::fill::eye), var_old(2, 2);
  double tau_ikm = 0.0, eta_ikm = 0.0;
  arma::mat Sigma_M = Rcpp::as<arma::mat>(init["Sigma_M"]);
  arma::mat d = Rcpp::as<arma::mat>(init["d"]);
  arma::mat Sigma_M_m_m(1, M), Sigma_M_m(M, M), delta_ik_m(1, M);
  arma::mat d_1(1, M), d_k(1, M), d_1k(1, M), d_1k_m(1, M);
  Rcpp::IntegerVector baseline = study_data["baseline"];
  Rcpp::IntegerVector baseline_id = study_id["baseline"];
  arma::mat mu_delta_prop(1, 2), mu_delta_curr(1, 2);
  Rcpp::NumericMatrix mu_prop = Rcpp::clone(mu);
  double mu_delta_A = 0.0, mu_delta_B = 0.0, target_mu_delta_prop = 0.0, target_mu_delta_curr = 0.0;
  Rcpp::NumericMatrix mu_long_prop(n_datapoints, M), delta_prop(n_datapoints, M);
  int accept_mu_delta = 0;
  arma::cube mu_chain(n_study, M, totiter);
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
  arma::mat d_tmp(n_trt - 1, M);

  int count = 0;
  arma::uvec row_id(n_trt - 1);
  for (unsigned int k = 0; k < n_trt; k++) {
    if (k != (ref_trt - 1)) {
      row_id(count) = k;
      count++;
    }
  }

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
        D_inv = arma::inv_sympd(D);
        for (unsigned int k = 0; k < n_q; k++) {
          S_q += D_inv*arma::trans(x_q.row(k))*x_q.row(k)*D_inv;
        }
        if (!is_positive_definite(S_q)) {
          S_q = make_positive_definite(S_q);
        }
        // in the following, degrees of freedom are set to (n_q > M) ? n_q : M)
        // otherwise, when n_q < M the IW routine returns an error
        // (indeterminate system)
        Gamma_q_prop = rinvwish_arma(((n_q > M) ? n_q : M), S_q);
        D_prop = arma::sqrt(arma::diagmat(Gamma_q_prop));
        Gamma_q_prop_inv = arma::inv_sympd(D_prop);
        Gamma_q_prop = Gamma_q_prop_inv*Gamma_q_prop*Gamma_q_prop_inv;

        Gamma_q_curr = Rcpp::as<arma::mat>(Gamma(q));
        Gamma_rate = (M + 1)*(log(arma::det(Gamma_q_prop)) - log(arma::det(Gamma_q_curr)))/2.0;
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
      D_inv = arma::inv_sympd(D);
      for (unsigned int k = 0; k < n_datapoints; k++) {
        S_q += D_inv*arma::trans(x_imp_arma.row(k))*x_imp_arma.row(k)*D_inv;
      }
      if (!is_positive_definite(S_q)) {
        S_q = make_positive_definite(S_q);
      }
      // in the following, degrees of freedom are set to (n_datapoints > M) ?
      // n_datapoints : M) otherwise, when n_datapoints < M the IW routine
      // returns an error (indeterminate system)
      Gamma_q_prop = rinvwish_arma(((n_datapoints > M) ? n_datapoints : M), S_q);
      D_prop = arma::sqrt(arma::diagmat(Gamma_q_prop));
      Gamma_q_prop_inv = arma::inv_sympd(D_prop);
      Gamma_q_prop = Gamma_q_prop_inv*Gamma_q_prop*Gamma_q_prop_inv;

      Gamma_q_curr = Rcpp::as<arma::mat>(Gamma(0));
      Gamma_rate = (M + 1)*(log(arma::det(Gamma_q_prop)) - log(arma::det(Gamma_q_curr)))/2.0;
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
          a_ikm = R::pbinom(y_imp(ik, m) - 1, n_imp(ik, m), p_ikm, 1, 0);
          b_ikm = R::pbinom(y_imp(ik, m), n_imp(ik, m), p_ikm, 1, 0);
          eps_a = a_ikm/2.0;
          eps_b = (1.0 - b_ikm)/2.0;
          if (a_ikm == b_ikm) {
            if (a_ikm == 0) {
              b_ikm += eps;
            } else if (b_ikm == 1) {
              a_ikm -= eps;
            } else {
              a_ikm -= fmin(eps/2.0, eps_a);
              b_ikm += fmin(eps/2.0, eps_b);
            }
          }
          Gamma_k_m_m = Gamma_k.row(m);
          Gamma_k_m_m.shed_col(m);
          Gamma_k_m = Gamma_k;
          Gamma_k_m.shed_col(m);
          Gamma_k_m.shed_row(m);
          x_ik_m = x_imp_arma.row(ik);
          x_ik_m.shed_col(m);
          w_ikm = arma::as_scalar(Gamma_k_m_m*arma::inv(Gamma_k_m)*arma::trans(x_ik_m));
          gamma_ikm = sqrt(arma::as_scalar(Gamma_k(m, m) - Gamma_k_m_m*arma::inv(Gamma_k_m)*arma::trans(Gamma_k_m_m)));
          x(ik, m) = rtruncnorm_rcpp(1, R::qnorm(a_ikm, 0.0, 1.0, 1, 0), R::qnorm(b_ikm, 0.0, 1.0, 1, 0), w_ikm, gamma_ikm)[0];
          x_imp(ik, m) = x(ik, m);
          x_imp_arma(ik, m) = x(ik, m);
        }
      }
    }
    x_chain.slice(niter) = Rcpp::as<arma::mat>(x);

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

        par(0) = mu_long(ik, m);
        par(1) = delta(ik, m); // par() should be the marginal model estimates
        args["y"] = y_imp(ik, m);
        args["n"] = n_imp(ik, m);
        args["tau"] = tau_ikm;
        args["eta"] = eta_ikm;
        args["w"] = w_ikm;
        args["gamma"] = gamma_ikm;
        fit = laplace_rcpp(mu_delta_logpost_R, par, args, options);
        if (Rcpp::as<bool>(fit["converge"])) {
          mode = Rcpp::as<arma::vec>(fit["mode"]);
          var = Rcpp::as<arma::mat>(fit["var"]);
          if (ISNAN(var(0, 0)) || ISNAN(var(0, 1)) || ISNAN(var(1, 0)) || ISNAN(var(1, 1))) {
            var = var_old;
          }
          if (!is_positive_definite(var)) {
            var = make_positive_definite(var);
          }
        }
        var_old = var;

        mu_delta_prop = rmvt_arma(1, mode, var, 7);
        mu_delta_curr(0, 0) = mu_long(ik, m);
        mu_delta_curr(0, 1) = delta(ik, m);
        mu_long_prop = Rcpp::clone(mu_long);
        delta_prop = Rcpp::clone(delta);
        // if (trt(ik) == baseline(ik)) {
          mu_long_prop(ik, m) = mu_delta_prop(0, 0);
        //   delta_prop(ik, m) = 0;
        // } else {
          delta_prop(ik, m) = mu_delta_prop(0, 1);
        // }
        // mu_prop = param_wide(mu_long_prop, narms, trt, baseline);
        // mu_long_prop = param_long(mu_prop, narms, false);
        delta_arma_prop = Rcpp::as<arma::mat>(delta_prop);

        target_mu_delta_curr = mu_delta_logpost(mu_long(ik, m), delta(ik, m), tau_ikm, eta_ikm, static_cast<double>(y_imp(ik, m)), static_cast<double>(n_imp(ik, m)), w_ikm, gamma_ikm, mu_sigma, eps, eps_ab);

        delta_ik_m = delta_arma_prop.row(ik);
        delta_ik_m.shed_col(m);
        d_1k_m = d_1k;
        d_1k_m.shed_col(m);
        tau_ikm = d_1k(m) + arma::as_scalar(Sigma_M_m_m*arma::inv(Sigma_M_m)*arma::trans(delta_ik_m - d_1k_m));

        target_mu_delta_prop = mu_delta_logpost(mu_long_prop(ik, m), delta_prop(ik, m), tau_ikm, eta_ikm, static_cast<double>(y_imp(ik, m)), static_cast<double>(n_imp(ik, m)), w_ikm, gamma_ikm, mu_sigma, eps, eps_ab);

        mu_delta_A = target_mu_delta_prop - target_mu_delta_curr;
        mu_delta_B = dmvt_arma(mu_delta_curr, mode, var, 7, true)(0) - dmvt_arma(mu_delta_prop, mode, var, 7, true)(0);
        mu_delta_rate = exp(mu_delta_A + mu_delta_B);
        ran_unif = R::runif(0.0, 1.0);
        if (!ISNAN(target_mu_delta_prop) && !ISNAN(target_mu_delta_curr)) {
          if (ran_unif < mu_delta_rate) {
            mu_long(ik, m) = mu_long_prop(ik, m);
            // mu = param_wide(mu_long, narms, trt, baseline);
            // mu_long = param_long(mu, narms, false);
            delta(ik, m) = delta_prop(ik, m);
            accept_mu_delta++;
          }
        }
      }
    }
    mu = param_wide(mu_long, narms, trt, baseline);
    // mu_long = param_long(mu, narms, false);
    // for (unsigned int ik = 0; ik < n_datapoints; ik++) {
    //   for (unsigned int m = 0; m < M; m++) {
    //     if (trt(ik) == baseline(ik)) {
    //       delta(ik, m) = 0;
    //     }
    //   }
    // }
    delta_arma = Rcpp::as<arma::mat>(delta);
    mu_chain.slice(niter) = Rcpp::as<arma::mat>(mu);
    delta_chain.slice(niter) = Rcpp::as<arma::mat>(delta);

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
    if (!is_positive_definite(Sigma_M_prop)) {
      Sigma_M_prop = make_positive_definite(Sigma_M_prop);
    }
    target_Sigma_M_prop = Sigma_M_logpost(d, delta_arma, Sigma_M_prop, trt, baseline, narms_study, sigma_r);
    target_Sigma_M_curr = Sigma_M_logpost(d, delta_arma, Sigma_M_curr, trt, baseline, narms_study, sigma_r);
    Sigma_M_A = target_Sigma_M_prop - target_Sigma_M_curr;
    // Sigma_M_B dovrebbe essere la differenza della proposal distribution, ma
    // e' simmetrica, quindi sparisce; la prior invece e' gia' dentro la target
    Sigma_M_B = 0.0; //dlogchol_arma(Sigma_M_prop, sigma_r, true) - dlogchol_arma(Sigma_M_curr, sigma_r, true);
    Sigma_M_rate = exp(Sigma_M_A + Sigma_M_B);
    // Rprintf("Sigma_M_A = %4.8f - Sigma_M_rate = %4.8f\n", Sigma_M_A, Sigma_M_rate);
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
        prop_d = rwmh_adapt(theta_d_reshaped, mu_d, rho_d, cov_d, ar_d, alpha, beta, gamma, tar, adapt_k, false);
        rho_d = Rcpp::as<double>(prop_d["rho"]);
        cov_d = Rcpp::as<arma::mat>(prop_d["covariance"]);
        if (!is_positive_definite(cov_d)) {
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
        prop_beta = rwmh_adapt(theta_beta, mu_beta, rho_beta, cov_beta, ar_Sigma_M, alpha, beta, gamma, tar, adapt_k, false);
        rho_beta = Rcpp::as<double>(prop_beta["rho"]);
        cov_beta = Rcpp::as<arma::mat>(prop_beta["covariance"]);
        if (!is_positive_definite(cov_beta)) {
          cov_beta = make_positive_definite(cov_beta);
        }
        mu_beta = Rcpp::as<arma::vec>(prop_beta["mu"]);
        ar_Sigma_M(0) = Rcpp::as<double>(prop_beta["ar"]);
      }

      adapt_k++;
    }

    for (unsigned int it = 0; it < totiter; it++) {
      for (unsigned int ik = 0; ik < n_datapoints; ik++) {
        for (unsigned int m = 0; m < M; m++) {
          if (trt(ik) == baseline(ik)) {
            delta_chain(ik, m, it) = 0;
          }
        }
      }
    }

    // calculate the loglikelihood, logprior and logposterior
    loglik(niter) = nc_loglik(y_imp, n_imp, x_imp, trt, mu, delta, Gamma); // c'e' il problema degli indicatori nella definizione della augmented likelihood
    logprior(niter) += arma::sum(dmvn_arma(Rcpp::as<arma::mat>(mu), arma::zeros<arma::vec>(M), pow(mu_sigma, 2)*arma::eye(M, M), true));
    d_tmp = d.rows(row_id);
    logprior(niter) += arma::sum(dmvn_arma(d_tmp, arma::zeros<arma::vec>(M), pow(d_sigma, 2)*arma::eye(M, M), true));
    logprior(niter) += dlogchol_arma(Sigma_M, sigma_r, true);

    // print the information
    if ((((niter + 1) % 100) == 0) && verbose) {
      Rprintf("   iter. %d/%d ==> acc. d: %1.3f - acc. Gamma: %1.3f - acc. mu/delta: %1.3f - acc. Sigma_M: %1.3f\n", (niter + 1), totiter, accept_d/(static_cast<double>(totiter)), arma::mean(accept_Gamma_q)/(static_cast<double>(totiter)), accept_mu_delta/(static_cast<double>(totiter*n_datapoints*M)), accept_Sigma_M/(static_cast<double>(totiter)));
    }

    niter++;

    R_CheckUserInterrupt();
  }

  PutRNGstate();

  Rcpp::List accept = Rcpp::List::create(
    Rcpp::Named("d") = accept_d/static_cast<double>(totiter),
    Rcpp::Named("mu_delta") = accept_mu_delta/static_cast<double>(totiter*n_datapoints*M),
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
Rcpp::List rwmh_adapt(const arma::mat& theta, const arma::vec& mu, const double& rho, const arma::mat& cov, const arma::vec& ar, const double& alpha, const double& beta, const double& gamma, const double& tar, const int& k, const bool& iter_cols) {
  double ar_new = (1 - alpha)*ar(0) + alpha*ar(1);
  double beta_k = beta/pow(k + 1, gamma);
  double rho_new = rho*exp(beta_k*(R::qnorm(ar_new/2.0, 0.0, 1.0, 1, 0) - R::qnorm(tar/2.0, 0.0, 1.0, 1, 0)));
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

  return Rcpp::List::create(Rcpp::Named("rho") = rho_new,
                            Rcpp::Named("covariance") = cov_new,
                            Rcpp::Named("mu") = mu_new,
                            Rcpp::Named("ar") = ar_new);
}
