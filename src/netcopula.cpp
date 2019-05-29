#include "netcopula.h"

// [[Rcpp::depends("RcppArmadillo")]]

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
  if (Rcpp::any(Rcpp::is_na(x))) {
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
        Omega = make_positive_definite(Omega);
      }
      x_q_vec_miss = arma::vectorise(rmvn_arma(1, mu, Omega));
      x_q_vec_imp = nv_miss_replace(x_q_vec, x_q_vec_miss, x_q_vec_miss_i);
      x_q_imp = nv_unstack(x_q_vec_imp, M);
      for (unsigned int k = 0; k < n; k++) {
        x_imp(k, Rcpp::_) = x_q_imp(k, Rcpp::_);
      }
    }
  }

  return x_imp;
}

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

Rcpp::NumericMatrix y_imputed(const Rcpp::NumericMatrix& y, const Rcpp::NumericMatrix& x_imp, const Rcpp::IntegerVector& narms, const Rcpp::NumericMatrix& mu, const Rcpp::NumericMatrix& delta, const Rcpp::NumericMatrix& n_imp) {
  int n_datapoints = y.nrow(), M = y.ncol();
  Rcpp::NumericMatrix y_imp = Rcpp::clone(y);
  Rcpp::NumericMatrix mu_long = param_long(mu, narms, false);
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

double logpost(const double& mu, const double& delta, const double& y, const double& n, const double& w, const double& gamma, const double& eps, const double& eps_ab) {
  double theta = mu + delta;
  double p = expit_double(theta);

  double a = R::pbinom(y - 1, n, p, 1, 1);
  double b = R::pbinom(y, n, p, 1, 1);
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
    return NA_REAL;
  }
  
  return tmp;
}

double mu_logpost(const double& mu, const arma::vec& delta, const arma::vec& y, const arma::vec& n, const arma::vec& w, const double& gamma, const double& mu_sigma, const double& eps, const double& eps_ab) {
  int n_arms = y.n_elem;
  double llik = 0;

  for (unsigned k = 0; k < n_arms; k++) {
    llik += logpost(mu, delta(k), y(k), n(k), w(k), gamma, eps, eps_ab);
  }
  if (ISNAN(llik)) {
    return NA_REAL;
  }

  double lpost = llik + R::dnorm(mu, 0, mu_sigma, 1);

  return lpost;
}

double delta_logpost(const double& delta, const double& mu, const double& tau, const double& eta, const double& y, const double& n, const double& w, const double& gamma, const double& eps, const double& eps_ab) {
  double llik = logpost(mu, delta, y, n, w, gamma, eps, eps_ab);
  if (ISNAN(llik)) {
    return NA_REAL;
  }

  double lpost = llik + R::dnorm(delta, tau, eta, 1);

  return lpost;
}

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

Rcpp::List rwmh_adapt(const arma::mat& theta, const arma::vec& mu, const double& rho, const arma::mat& cov, const arma::vec& ar, const double& alpha, const double& beta, const double& gamma, const double& tar, const int& k, const bool& iter_cols, const int& what, const bool& diagonal) {

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
  int every = theta_k.n_cols;
  arma::vec theta_hat = arma::vectorise(arma::mean(theta_k, 1)); // means by rows
  arma::vec mu_new(arma::size(mu));
  if (k > 0) {
    mu_new = mu + beta_k*(theta_hat - mu);
  } else {
    mu_new = mu;
  }
  arma::mat theta_k_c = theta_k.each_col() - theta_hat;
  arma::mat cov_hat = theta_k_c*arma::trans(theta_k_c)/every;
  if (diagonal) {
    cov_hat = arma::diagmat(cov_hat);
  }
  arma::mat cov_new = (1 - beta_k)*cov + beta_k*cov_hat;
  if ((what == 1 || what == 2)  && (cov_new(0, 0) <= 0)) {
    cov_new = 1e-1;
  }

  return Rcpp::List::create(Rcpp::Named("rho") = rho_new,
                            Rcpp::Named("covariance") = cov_new,
                            Rcpp::Named("mu") = mu_new,
                            Rcpp::Named("ar") = ar_new);
}

double mudelta_logpost(const arma::vec& mudelta, const arma::mat& delta_arma, const arma::mat& y, const arma::mat& n, const arma::mat& x, const Rcpp::IntegerVector& baseline, const Rcpp::IntegerVector& trt, const Rcpp::IntegerVector& study, const Rcpp::IntegerVector& narms, const arma::mat& d, const arma::mat& Sigma_M, const Rcpp::List& Gamma, const int& i, const int& m, const double& mu_sigma, const double& eps, const double& eps_ab) {
  int n_datapoints = y.n_rows, M = y.n_cols, nGamma = Gamma.size();
  double mu = 0.0, delta = 0.0, llik = 0.0, mu_lprior = 0.0, delta_lprior = 0.0;
  arma::mat d_1(1, M), d_k(1, M), d_1k(1, M), d_1k_m(1, M);
  arma::mat Sigma_M_m_m(1, M), Sigma_M_m(M, M), delta_ik_m(1, M);
  arma::mat Gamma_k(M, M), Gamma_k_m_m(1, M), Gamma_k_m(M, M), x_ik_m(1, M);
  double w_ikm = 0.0, gamma_ikm = 0.0, tau_ikm = 0.0, eta_ikm = 0.0;
  int count_arms = 0;

  mu = mudelta(0);
  mu_lprior = R::dnorm(mu, 0, mu_sigma, 1);

  for (unsigned int ik = 0; ik < n_datapoints; ik++) {
    if (study(ik) == (i + 1)) {
      if (nGamma > 1) {
        Gamma_k = Rcpp::as<arma::mat>(Gamma(trt(ik) - 1));
      } else {
        Gamma_k = Rcpp::as<arma::mat>(Gamma(0));
      }
      Gamma_k_m_m = Gamma_k.row(m);
      Gamma_k_m_m.shed_col(m);
      Gamma_k_m = Gamma_k;
      Gamma_k_m.shed_col(m);
      Gamma_k_m.shed_row(m);
      x_ik_m = x.row(ik);
      x_ik_m.shed_col(m);
      w_ikm = arma::as_scalar(Gamma_k_m_m*arma::inv(Gamma_k_m)*arma::trans(x_ik_m));
      gamma_ikm = sqrt(arma::as_scalar(Gamma_k(m, m) - Gamma_k_m_m*arma::inv(Gamma_k_m)*arma::trans(Gamma_k_m_m)));

      if (baseline(ik) == trt(ik)) {
        delta = 0.0;
        count_arms++;
      } else {
        delta = mudelta(count_arms);
        count_arms++;

        Sigma_M_m_m = Sigma_M.row(m);
        Sigma_M_m_m.shed_col(m);
        Sigma_M_m = Sigma_M;
        Sigma_M_m.shed_col(m);
        Sigma_M_m.shed_row(m);
        delta_ik_m = delta_arma.row(ik);
        delta_ik_m.shed_col(m);
        d_1 = d.row(baseline(ik) - 1);
        d_k = d.row(trt(ik) - 1);
        d_1k = d_k - d_1;
        d_1k_m = d_1k;
        d_1k_m.shed_col(m);
        tau_ikm = d_1k(m) + arma::as_scalar(Sigma_M_m_m*arma::inv(Sigma_M_m)*arma::trans(delta_ik_m - d_1k_m));
        eta_ikm = sqrt(arma::as_scalar(Sigma_M(m, m) - Sigma_M_m_m*arma::inv(Sigma_M_m)*arma::trans(Sigma_M_m_m)));

        delta_lprior += R::dnorm(delta, tau_ikm, eta_ikm, 1);
      }

      llik += logpost(mu, delta, y(ik, m), n(ik, m), w_ikm, gamma_ikm, eps, eps_ab);
    }
  }
  if (ISNAN(llik)) {
    return NA_REAL;
  }

  double lpost = llik + mu_lprior + delta_lprior;

  return lpost;
}

double r_logpost(const arma::mat& Gamma, const arma::mat& x) {
  int M = Gamma.n_rows, n = x.n_rows;

  arma::mat S(M, M);
  S = arma::trans(x) * x;

  double logdet, logdet_sgn;
  arma::log_det(logdet, logdet_sgn, Gamma);

  double out = -0.5*(n*logdet + arma::trace(S * arma::inv_sympd(Gamma)));

  return out;
}

//' Internal functions for MCMC simulation.
//'
//' For internal use only.
//'
//' @param data internal data structure
//' @param init internal data structure
//' @param totiter internal data structure
//' @param prior internal data structure
//' @param prop internal data structure
//' @param tuning internal data structure
//' @param adapt internal data structure
//' @param verbose internal data structure
//'
//' @aliases netcopula-internal
//' @aliases netcopula_internal
//'
// [[Rcpp::export]]
Rcpp::List nc_mcmc_mh(const Rcpp::RObject& data, const Rcpp::List& init, const int& totiter, const Rcpp::List& prior, const Rcpp::List& prop, const Rcpp::List& tuning, const Rcpp::List& adapt, const bool& verbose) {

  int niter = 0, adapt_k = 0, n_datapoints = data.slot("n_datapoints"), M = data.slot("n_outcomes"), n_trt = data.slot("n_treatments"), n_study = data.slot("n_study"), nn = M*(M + 1)/2;

  double eps = Rcpp::as<double>(tuning["eps"]), eps_ab = Rcpp::as<double>(tuning["eps_ab"]); //, eps_a = 0.0, eps_b = 0.0;

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
  arma::mat delta_arma = Rcpp::as<arma::mat>(delta);
  arma::mat Sigma_M = Rcpp::as<arma::mat>(init["Sigma_M"]), d = Rcpp::as<arma::mat>(init["d"]);
  Rcpp::NumericMatrix y = df_nm(study_data, cols_y);
  Rcpp::NumericMatrix n_data = df_nm(study_data, cols_n);
  Rcpp::NumericMatrix n_imp = n_imputed(n_data);
  arma::mat n_imp_arma = Rcpp::as<arma::mat>(n_imp);
  Rcpp::NumericMatrix x_imp(n_datapoints, M), y_imp(n_datapoints, M);
  arma::mat y_imp_arma(n_datapoints, M), x_imp_arma(n_datapoints, M), x_star_arma(n_datapoints, M);

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
  double ran_unif, mu_delta_rate, Gamma_rate, d_rate, Sigma_M_rate, r_rate;
  double Gamma_A = 0.0, Gamma_B = 0.0, target_Gamma_prop = 0.0, target_Gamma_curr = 0.0;
  arma::vec accept_Gamma_q = arma::zeros<arma::vec>(nGamma);
  arma::cube D_chain(nGamma, M, totiter), S_q_chain(nGamma, nn, totiter);
  int nn_minus = 0;
  if (M > 1) {
    nn_minus = M*(M - 1)/2;
  } else {
    nn_minus = 1;
  }
  arma::cube Gamma_chain(nGamma, nn_minus, totiter);

  arma::mat r(nGamma, nn_minus);
  arma::mat pippo;
  for (unsigned q = 0; q < nGamma; q++) {
    if (M > 1) {
      r.row(q) = get_upper_tri(arma::chol(arma::inv_sympd(Rcpp::as<arma::mat>(Gamma(q)))), false);
    } else {
      r(q, 0) = 1;
    }
  }
  arma::cube r_chain(nGamma, nn_minus, totiter);
  double r_curr = 0.0, r_prop = 0.0;
  double r_A = 0.0, r_B = 0.0, target_r_prop = 0.0, target_r_curr = 0.0;
  double sigma_r_prop = prop["sigma_r_prop"];
  arma::mat R_q_curr(M, M), R_q_prop(M, M), Sigma_q_curr(M, M);
  arma::rowvec r_curr_vec(nn_minus), r_prop_vec(nn_minus);
  arma::vec orig_idx(nn_minus), new_idx(nn_minus);
  orig_idx = arma::linspace(0, nn_minus - 1, nn_minus);

  arma::mat Gamma_k(M, M), Gamma_k_m_m(1, M), Gamma_k_m(M, M), x_ik_m(1, M);
  double a_ikm = 0.0, b_ikm = 0.0, w_ikm = 0.0, gamma_ikm = 0.0, theta_ikm = 0.0, p_ikm = 0.0;
  arma::cube x_chain(n_datapoints, M, totiter);
  arma::vec D_tmp(M);
  arma::mat Gamma_tmp(M, M), Gamma_inv(M, M), Gamma_new(M, M), S(M, M), W(n_datapoints, M);

  double mu_sigma = sqrt(Rcpp::as<double>(prior["mu_sigma2"]));
  arma::vec mu_coef(3), delta_coef(3);

  arma::mat Sigma_M_m_m(1, M), Sigma_M_m(M, M), delta_ik_m(1, M);
  arma::mat d_1(1, M), d_k(1, M), d_1k(1, M), d_1k_m(1, M);
  Rcpp::IntegerVector baseline = study_data["baseline"];
  Rcpp::IntegerVector baseline_id = study_id["baseline"];
  Rcpp::NumericMatrix mu_long = param_long(mu, narms, false);

  arma::cube mu_chain(n_study, M, totiter);
  arma::cube delta_chain(n_datapoints, M, totiter);

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
  ar_d(0) = 1.0;
  Rcpp::List prop_d(4);

  arma::vec rho_r = 2.38/sqrt(nn_minus)*arma::ones(nn_minus);
  arma::vec mu_r_arma(1);
  arma::vec cov_r = pow(sigma_r_prop, 2)*arma::ones(nn_minus);
  arma::mat cov_r_arma(1, 1);
  arma::vec mu_r(nn_minus);
  arma::cube theta_r(nGamma, nn_minus, every);
  arma::mat theta_r_reshaped(every, nn_minus*nGamma);
  arma::mat ar_r_single_vec(nn_minus, totiter), ar_r(nn_minus, 2);
  ar_r.col(0) = arma::ones<arma::vec>(nn_minus);
  arma::mat accept_r_q = arma::zeros<arma::mat>(nGamma, nn_minus);

  double rho_r_vec = 2.38/sqrt(nn_minus);
  arma::mat cov_r_vec = arma::eye(nn_minus, nn_minus);
  arma::vec mu_r_vec(nn_minus);
  arma::cube theta_r_vec(nGamma, nn_minus, every);
  arma::mat theta_r_vec_reshaped(every, nn_minus*nGamma);
  arma::vec ar_r_vec_vec(totiter), ar_r_vec(2);
  ar_r_vec(0) = 1.0;
  Rcpp::List prop_r(4);

  double rho_beta = 2.38/sqrt(nn);
  arma::mat cov_beta = arma::eye(nn, nn);
  arma::mat Sigma_M_curr(Sigma_M), Sigma_M_prop(Sigma_M);
  arma::vec beta_curr(nn), beta_prop(nn);
  double sigma_b = Rcpp::as<double>(prior["beta_sigma"]);
  double Sigma_M_A = 0.0, Sigma_M_B = 0.0, target_Sigma_M_prop = 0.0, target_Sigma_M_curr = 0.0;
  arma::vec mu_beta(nn);
  arma::mat theta_beta(every, nn);
  arma::vec ar_Sigma_M_vec(totiter), ar_Sigma_M(2);
  ar_Sigma_M(0) = 1.0;
  Rcpp::List prop_beta(4);

  arma::vec loglik(totiter, arma::fill::zeros), logprior(totiter, arma::fill::zeros), logpost(totiter, arma::fill::zeros);

  arma::field<arma::vec> mu_delta_chain(n_study, M, totiter);
  int count_arms = 0;
  arma::mat rho_mu_delta(n_study, M);
  for (unsigned int i = 0; i < n_study; i++) {
    for (unsigned int m = 0; m < M; m++) {
      rho_mu_delta(i, m) = 2.38/sqrt(narms(i));
    }
  }
  arma::vec mu_delta_curr, mu_delta_prop;
  double mu_delta_A = 0.0, mu_delta_B = 0.0, target_mu_delta_curr = 0.0, target_mu_delta_prop = 0.0;
  int accept_mu_delta = 0;
  arma::field<arma::vec> mu_mu_delta(n_study, M);
  for (unsigned int i = 0; i < n_study; i++) {
    for (unsigned int m = 0; m < M; m++) {
      mu_mu_delta(i, m) = arma::zeros(narms(i));
    }
  }
  arma::vec mu_mu_delta_vec;
  arma::mat theta_mu_delta;
  arma::cube ar_mu_delta_vec(n_study, M, totiter);
  arma::cube ar_mu_delta(2, n_study, M, arma::fill::ones);
  arma::vec ar_mu_delta_sub(2, arma::fill::zeros);
  Rcpp::List prop_mu_delta(4);
  arma::field<arma::mat> cov_mu_delta(n_study, M);
  for (unsigned int i = 0; i < n_study; i++) {
    for (unsigned int m = 0; m < M; m++) {
      cov_mu_delta(i, m) = arma::eye(narms(i), narms(i));
    }
  }
  arma::mat cov_mu_delta_mat;

  arma::vec tmp(every);

  GetRNGstate();

  while (niter < totiter) {
    // imputing x and y variables
    x_imp = x_imputed(x, Gamma, trt);
    x_imp_arma = Rcpp::as<arma::mat>(x_imp);
    y_imp = y_imputed(y, x_imp, narms, mu, delta, n_imp);
    y_imp_arma = Rcpp::as<arma::mat>(y_imp);

    for (unsigned int i = 0; i < n_study; i++) {
    // updating mu (study-specific baseline effects) and delta (study-specific
    // [random] treatment effects); the update occurs simultaneously for mu and delta
    // in each study and separately for each outcome
      mu_delta_curr.set_size(narms(i));
      mu_delta_prop.set_size(narms(i));
      for (unsigned int m = 0; m < M; m++) {
        mu_delta_curr(0) = mu(i, m);
        for (unsigned int k = 1; k < narms(i); k++) {
          mu_delta_curr(k) = delta(k + count_arms, m);
        }
        mu_delta_prop = arma::vectorise(rmvn_arma(1, mu_delta_curr, pow(sd_multplier*rho_mu_delta(i, m), 2)*arma::diagmat(cov_mu_delta(i, m))));

        target_mu_delta_curr = mudelta_logpost(mu_delta_curr, delta_arma, y_imp_arma, n_imp_arma, x_imp_arma, baseline, trt, study_index, narms_study, d, Sigma_M, Gamma, i, m, mu_sigma, eps, eps_ab);
        target_mu_delta_prop = mudelta_logpost(mu_delta_prop, delta_arma, y_imp_arma, n_imp_arma, x_imp_arma, baseline, trt, study_index, narms_study, d, Sigma_M, Gamma, i, m, mu_sigma, eps, eps_ab);

        mu_delta_A = target_mu_delta_prop - target_mu_delta_curr;
        mu_delta_B = 0.0;
        mu_delta_rate = exp(mu_delta_A + mu_delta_B);
        ran_unif = R::runif(0.0, 1.0);
        if (!ISNAN(target_mu_delta_prop) && !ISNAN(target_mu_delta_curr)) {
          if (ran_unif < mu_delta_rate) {
            mu(i, m) = mu_delta_prop(0);
            mu_delta_curr(0) = mu_delta_prop(0);
            for (unsigned int k = 1; k < narms(i); k++) {
              delta(k + count_arms, m) = mu_delta_prop(k);
              mu_delta_curr(k) = mu_delta_prop(k);
            }
            mu_long = param_long(mu, narms, false);
            delta_arma = Rcpp::as<arma::mat>(delta);
            accept_mu_delta++;
          }
        } else {
          mu_delta_rate = 0.0;
        }
        ar_mu_delta_vec(i, m, niter) = fmin(mu_delta_rate, 1);

        // updating x (latent variables)
        for (unsigned int ik = 0; ik < n_datapoints; ik++) {
          if (study_index(ik) == (i + 1)) {
            if (nGamma > 1) {
              Gamma_k = Rcpp::as<arma::mat>(Gamma(trt(ik) - 1));
            } else {
              Gamma_k = Rcpp::as<arma::mat>(Gamma(0));
            }
            if (!ISNAN(y(ik, m))) {
              theta_ikm = mu_long(ik, m) + delta(ik, m);
              p_ikm = expit_double(theta_ikm);
              a_ikm = R::pbinom(y_imp(ik, m) - 1, n_imp(ik, m), p_ikm, 1, 1); // log
              b_ikm = R::pbinom(y_imp(ik, m), n_imp(ik, m), p_ikm, 1, 1); // log
              if (a_ikm == b_ikm) {
                // if (a_ikm == 0) {
                //   b_ikm += eps;
                // } else if (b_ikm == 1) {
                //   a_ikm -= eps;
                // } else {
                //   a_ikm -= fmin(eps/2.0, eps_a);
                //   b_ikm += fmin(eps/2.0, eps_b);
                // }
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
        mu_delta_chain(i, m, niter) = mu_delta_curr;
      }
      count_arms += narms(i);
    }
    count_arms = 0;
    mu_chain.slice(niter) = Rcpp::as<arma::mat>(mu);
    delta_chain.slice(niter) = delta_arma;
    x_chain.slice(niter) = Rcpp::as<arma::mat>(x);

    // updating Gamma (correlation of latent variables)
    if (M > 1) {
      if (Gamma_update == "PX-RPMH") {
        if (nGamma > 1) {
          x2t_split = split_nm(x_imp, trt);
          for (unsigned q = 0; q < n_trt; q++) {
            x_q = Rcpp::as<arma::mat>(x2t_split[q]);
            n_q = x_q.n_rows;
            D = Rcpp::as<arma::mat>(D_list(q));
            for (unsigned int k = 0; k < n_q; k++) {
              S_q += D*arma::trans(x_q.row(k))*x_q.row(k)*D;
            }
            if (!is_positive_definite(S_q, 360)) {
              S_q = make_positive_definite(S_q);
            }
            S_q_chain(arma::span(q), arma::span::all, arma::span(niter)) = diag_tri(S_q);
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
            D_chain(arma::span(q), arma::span::all, arma::span(niter)) = arma::diagvec(D);
          }
        } else {
          D = arma::diagmat(arma::sqrt(1/arma::sum(arma::pow(x_imp_arma, 2), 0)));
          x_star_arma = x_imp_arma * D;
          S_q = arma::trans(x_star_arma)*x_star_arma;
          if (!is_positive_definite(S_q, 394)) {
            S_q = make_positive_definite(S_q);
          }
          S_q_chain(arma::span(0), arma::span::all, arma::span(niter)) = diag_tri(S_q);
          Sigma_q_prop = rinvwish_arma(((n_datapoints > M) ? n_datapoints : M), S_q);
          D_prop = arma::sqrt(arma::diagmat(Sigma_q_prop));
          Gamma_q_prop = cov2cor_rcpp(Sigma_q_prop);

          Gamma_q_curr = Rcpp::as<arma::mat>(Gamma(0));
          Gamma_rate = 0.5*(M + 1)*(log(arma::det(Gamma_q_prop)) - log(arma::det(Gamma_q_curr)));
          ran_unif = R::runif(0.0, 1.0);
          if (ran_unif < exp(Gamma_rate)) {
            // adjust latent variables x according to the new scales in 'D_prop'
            D_prop_inv = arma::inv_sympd(D_prop);
            x_imp_arma = x_star_arma * D_prop_inv;
            for (unsigned int ik = 0; ik < n_datapoints; ik++) {
              for (unsigned int m = 0; m < M; m++) {
                x_imp(ik, m) = x_imp_arma(ik, m);
                if (!ISNAN(y(ik, m))) {
                  x(ik, m) = x_imp_arma(ik, m);
                }
              }
            }

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
            if (ran_unif < Gamma_rate) {
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
          if (ran_unif < Gamma_rate) {
            Gamma_q_curr = Gamma_q_prop;
            accept_Gamma_q(0)++;
          }
          Gamma(0) = Gamma_q_curr;
        }
      } else if (Gamma_update == "DanaherSmith") {
        new_idx = arma::shuffle(orig_idx);
        for (unsigned int v = 0; v < nn_minus; v++) {
          // r_curr = r(0, new_idx(v));
          r_curr = r(0, v);
          // R_q_curr = put_upper_tri(r.row(0), false);
          // R_q_curr.diag().ones();
          // Sigma_q_curr = arma::inv_sympd(R_q_curr.t() * R_q_curr);
          // Gamma_q_curr = cov2cor_rcpp(Sigma_q_curr);
          Gamma_q_curr = Rcpp::as<arma::mat>(Gamma(0));

          r_prop = Rcpp::rnorm(1, r_curr, sigma_r_prop)(0);
          // r_prop = Rcpp::rnorm(1, r_curr, sd_multplier*rho_r(v)*sqrt(cov_r(v)))(0);
          r_prop_vec = r.row(0);
          // r_prop_vec(0, new_idx(v)) = r_prop;
          r_prop_vec(0, v) = r_prop;
          R_q_prop = put_upper_tri(r_prop_vec, false);
          R_q_prop.diag().ones();
          Sigma_q_prop = arma::inv_sympd(R_q_prop.t() * R_q_prop);
          Gamma_q_prop = cov2cor_rcpp(Sigma_q_prop);

          target_r_curr = r_logpost(Gamma_q_curr, x_imp_arma);
          target_r_prop = r_logpost(Gamma_q_prop, x_imp_arma);

          r_A = target_r_prop - target_r_curr;
          r_B = 0;
          r_rate = exp(r_A + r_B);
          ran_unif = R::runif(0.0, 1.0);
          if (ran_unif < r_rate) {
            // r(0, new_idx(v)) = r_prop;
            r(0, v) = r_prop;
            Gamma_q_curr = Gamma_q_prop;
            accept_r_q(0, v)++;
            accept_Gamma_q(0)++;
          }
          Gamma(0) = Gamma_q_curr;

          ar_r_single_vec(v, niter) = fmin(r_rate, 1);
        }
      } else if (Gamma_update == "JointR") {
        // the difference with "DanaherSmith" is that here we propose all the r_ij elements jointly
        r_curr_vec = r.row(0);
        Gamma_q_curr = Rcpp::as<arma::mat>(Gamma(0));

        r_prop_vec = arma::trans(arma::vectorise(rmvt_arma(1, arma::trans(r_curr_vec), pow(sd_multplier*rho_r_vec, 2)*cov_r_vec, 7)));
        R_q_prop = put_upper_tri(r_prop_vec, false);
        R_q_prop.diag().ones();
        Sigma_q_prop = arma::inv_sympd(R_q_prop.t() * R_q_prop);
        Gamma_q_prop = cov2cor_rcpp(Sigma_q_prop);

        target_r_curr = r_logpost(Gamma_q_curr, x_imp_arma);
        target_r_prop = r_logpost(Gamma_q_prop, x_imp_arma);

        r_A = target_r_prop - target_r_curr;
        r_B = 0;
        r_rate = exp(r_A + r_B);
        ran_unif = R::runif(0.0, 1.0);
        if (ran_unif < r_rate) {
          r.row(0) = r_prop_vec;
          Gamma(0) = Gamma_q_prop;
          accept_Gamma_q(0)++;
        }

        ar_r_vec_vec(niter) = fmin(r_rate, 1);
      } else {
        error("the specified method for updating the Gamma parameter is not available. only 'Danaher-Smith', 'IMH', 'JointR' or 'PX-RPMH' are allowed.");
      }
    }
    D_chain(arma::span(0), arma::span::all, arma::span(niter)) = arma::diagvec(D);
    Gamma_chain.slice(niter) = list_mat(Gamma);
    r_chain.slice(niter) = r;
    x_chain.slice(niter) = Rcpp::as<arma::mat>(x);

    // updating d (pooled treatment effects across trials)
    d_curr = mat_to_vec(d, true, ref_trt);
    d_curr_ref = d;
    d_prop = rmvt_arma(1, d_curr, pow(sd_multplier*rho_d, 2)*cov_d, 7);
    d_prop_ref = vec_to_mat(arma::vectorise(d_prop), M, true, ref_trt);
    target_d_curr = d_logpost(d_curr_ref, delta_arma, Sigma_M, trt, baseline, narms_study, d_sigma, ref_trt);
    target_d_prop = d_logpost(d_prop_ref, delta_arma, Sigma_M, trt, baseline, narms_study, d_sigma, ref_trt);
    d_A = target_d_prop - target_d_curr;
    d_B = 0.0; //d_logprior(d_prop_ref, d_sigma, ref_trt) - d_logprior(d_curr_ref, d_sigma, ref_trt);
    d_rate = exp(d_A + d_B);
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
    beta_prop = arma::vectorise(rmvt_arma(1, beta_curr, pow(sd_multplier*rho_beta, 2)*cov_beta, 7));
    Sigma_M_prop = beta_to_Sigma_M(beta_prop, M);
    if (!is_positive_definite(Sigma_M_prop, 515)) {
      Sigma_M_prop = make_positive_definite(Sigma_M_prop);
    }
    target_Sigma_M_prop = Sigma_M_logpost(d, delta_arma, Sigma_M_prop, trt, baseline, narms_study, sigma_b);
    target_Sigma_M_curr = Sigma_M_logpost(d, delta_arma, Sigma_M_curr, trt, baseline, narms_study, sigma_b);
    Sigma_M_A = target_Sigma_M_prop - target_Sigma_M_curr;
    Sigma_M_B = 0.0; //dlogchol_arma(Sigma_M_prop, sigma_b, true) - dlogchol_arma(Sigma_M_curr, sigma_b, true);
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

      // mu/delta proposal parameters
      for (unsigned int i = 0; i < n_study; i++) {
        for (unsigned int m = 0; m < M; m++) {
          if (fabs(ar_mu_delta(0, i, m) - tar) > tol) {
            theta_mu_delta = field_to_matrix(mu_delta_chain.subfield(i, m, niter - every + 1, i, m, niter));
            if (adapt_k == 0) {
              // means by rows
              mu_mu_delta(i, m) = arma::mean(theta_mu_delta, 1);
            }
            mu_mu_delta_vec = mu_mu_delta(i, m);
            tmp = arma::vectorise(ar_mu_delta_vec.subcube(i, m, niter - every + 1, i, m, niter));
            ar_mu_delta(1, i, m) = arma::as_scalar(arma::mean(tmp));
            ar_mu_delta_sub(0) = ar_mu_delta(0, i, m);
            ar_mu_delta_sub(1) = ar_mu_delta(1, i, m);
            cov_mu_delta_mat = cov_mu_delta(i, m);
            prop_mu_delta = rwmh_adapt(theta_mu_delta, mu_mu_delta_vec, rho_mu_delta(i, m), cov_mu_delta_mat, ar_mu_delta_sub, alpha, beta, gamma, tar, adapt_k, true, 2, true); // NOTE: the diagonal argument is set
                                                                                 //       to true
            rho_mu_delta(i, m) = Rcpp::as<double>(prop_mu_delta["rho"]);
            cov_mu_delta_mat = Rcpp::as<arma::mat>(prop_mu_delta["covariance"]);
            if (!is_positive_definite(cov_mu_delta_mat, 573)) {
              cov_mu_delta_mat = make_positive_definite(cov_mu_delta_mat);
            }
            cov_mu_delta(i, m) = cov_mu_delta_mat;
            mu_mu_delta_vec = Rcpp::as<arma::vec>(prop_mu_delta["mu"]);
            mu_mu_delta(i, m) = mu_mu_delta_vec;
            ar_mu_delta(0, i, m) = Rcpp::as<double>(prop_mu_delta["ar"]);
          }
        }
      }

      // R (elements of Cholesky decomposition of Gamma) proposal parameters
      if (Gamma_update == "JointR") {
        if (fabs(ar_r_vec(0) - tar) > tol) {
          theta_r_vec = r_chain.slices(niter - every + 1, niter);
          theta_r_vec_reshaped = cube_to_mat(theta_r_vec, false, ref_trt);
          if (adapt_k == 0) {
            // means by columns
            mu_r_vec = arma::vectorise(arma::mean(theta_r_vec_reshaped, 0));
          }
          ar_r_vec(1) = arma::as_scalar(arma::mean(ar_r_vec_vec.subvec(niter - every + 1, niter)));
          prop_r = rwmh_adapt(theta_r_vec_reshaped, mu_r_vec, rho_r_vec, cov_r_vec, ar_r_vec, alpha, beta, gamma, tar, adapt_k, false, 4, false);
          rho_r_vec = Rcpp::as<double>(prop_r["rho"]);
          cov_r_vec = Rcpp::as<arma::mat>(prop_r["covariance"]);
          if (!is_positive_definite(cov_r_vec, 602)) {
            cov_r_vec = make_positive_definite(cov_r_vec);
          }
          mu_r_vec = Rcpp::as<arma::vec>(prop_r["mu"]);
          ar_r_vec(0) = Rcpp::as<double>(prop_r["ar"]);
        }
      }

      // d proposal parameters
      if (fabs(ar_d(0) - tar) > tol) {
        theta_d = d_chain.slices(niter - every + 1, niter);
        theta_d_reshaped = cube_to_mat(theta_d, true, ref_trt);
        if (adapt_k == 0) {
          // means by columns
          mu_d = arma::vectorise(arma::mean(theta_d_reshaped, 0));
        }
        ar_d(1) = arma::as_scalar(arma::mean(ar_d_vec.subvec(niter - every + 1, niter)));
        prop_d = rwmh_adapt(theta_d_reshaped, mu_d, rho_d, cov_d, ar_d, alpha, beta, gamma, tar, adapt_k, false, 3, false);
        rho_d = Rcpp::as<double>(prop_d["rho"]);
        cov_d = Rcpp::as<arma::mat>(prop_d["covariance"]);
        if (!is_positive_definite(cov_d, 583)) {
          cov_d = make_positive_definite(cov_d);
        }
        mu_d = Rcpp::as<arma::vec>(prop_d["mu"]);
        ar_d(0) = Rcpp::as<double>(prop_d["ar"]);
      }

      // log Cholesky betas (Sigma_M) proposal parameters
      if (fabs(ar_Sigma_M(0) - tar) > tol) {
        theta_beta = beta_chain.rows(niter - every + 1, niter);
        if (adapt_k == 0) {
          // means by columns
          mu_beta = arma::vectorise(arma::mean(theta_beta, 0));
        }
        ar_Sigma_M(1) = arma::as_scalar(arma::mean(ar_Sigma_M_vec.subvec(niter - every + 1, niter)));
        prop_beta = rwmh_adapt(theta_beta, mu_beta, rho_beta, cov_beta, ar_Sigma_M, alpha, beta, gamma, tar, adapt_k, false, 4, false);
        rho_beta = Rcpp::as<double>(prop_beta["rho"]);
        cov_beta = Rcpp::as<arma::mat>(prop_beta["covariance"]);
        if (!is_positive_definite(cov_beta, 602)) {
          cov_beta = make_positive_definite(cov_beta);
        }
        mu_beta = Rcpp::as<arma::vec>(prop_beta["mu"]);
        ar_Sigma_M(0) = Rcpp::as<double>(prop_beta["ar"]);
      }

      adapt_k++;
    }

    // compute loglikelihood, logprior and logposterior
    // loglik(niter) = nc_loglik(y_imp, n_imp, x_imp, trt, mu, delta, Gamma);
    // logprior(niter) = nc_logprior(mu, mu_sigma, d, d_sigma, Sigma_M, sigma_b, ref_trt);

    // print the information
    if ((((niter + 1) % print_every) == 0) && verbose) {
      Rprintf("   iter. %d/%d ==> d: %1.3f - Gamma: %1.3f - mu/delta: %1.3f - Sigma_M: %1.3f\n", (niter + 1), totiter, accept_d/(static_cast<double>(totiter)), arma::mean(accept_Gamma_q)/(static_cast<double>(totiter*nn_minus)), accept_mu_delta/(static_cast<double>(totiter*n_study*M)), accept_Sigma_M/(static_cast<double>(totiter)));
    }

    niter++;

    R_CheckUserInterrupt();
  }

  PutRNGstate();

  Rcpp::List accept = Rcpp::List::create(
    Rcpp::Named("d") = accept_d/static_cast<double>(totiter),
    Rcpp::Named("mudelta") = accept_mu_delta/static_cast<double>(totiter*n_study*M),
    Rcpp::Named("Gamma") = accept_Gamma_q/static_cast<double>(totiter),
    Rcpp::Named("Sigma") = accept_Sigma_M/static_cast<double>(totiter));

  return Rcpp::List::create(Rcpp::Named("mu") = mu_chain,
                            Rcpp::Named("delta") = delta_chain,
                            Rcpp::Named("d") = d_chain,
                            Rcpp::Named("Sigma") = Sigma_M_chain,
                            Rcpp::Named("Gamma") = Gamma_chain,
                            Rcpp::Named("r") = r_chain,
                            Rcpp::Named("x") = x_chain,
                            Rcpp::Named("x_unadj") = x_chain,
                            Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("logprior") = logprior,
                            Rcpp::Named("logpost") = (loglik + logprior),
                            Rcpp::Named("accept") = accept);
}
