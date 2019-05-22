#ifndef NETCOPULA_H
#define NETCOPULA_H

// #define ARMA_NO_DEBUG // this macro disables bounds checks in the Armadillo library making the code faster (but also more frail!); the suggestion is to disable it only for the final release (see http://arma.sourceforge.net/docs.html)

#include <R.h>
#include <RcppArmadillo.h>
#include <string>

static const double machine_eps = 2.220446049250313080847e-16;
static const double log_pi = std::log(M_PI);
static const double log_2pi = std::log(2.0 * M_PI);
static const double log_two = std::log(2.0);

// UTILITY FUNCTIONS ----------------------------------------------------------
void tableC(int* counts, const int* x, int nelem, int ndistelem);
arma::mat mat_block_diag(const arma::mat& A, const int& n);
arma::mat Sigma_block(const arma::mat& Sigma_M, const int& n_arms); 
Rcpp::NumericMatrix df_nm(const Rcpp::DataFrame& x, const Rcpp::IntegerVector& cols);
Rcpp::NumericVector nm_stack(const Rcpp::NumericMatrix& x);
Rcpp::NumericVector nv_omit(const Rcpp::NumericVector& x);
arma::uvec nv_na_index(const Rcpp::NumericVector& x, const int& dim, const bool& type);
Rcpp::NumericMatrix nv_unstack(const Rcpp::NumericVector& x, const int& nc);
Rcpp::NumericVector nv_miss_replace(const Rcpp::NumericVector& x, const arma::vec& miss, const arma::uvec& miss_i);
Rcpp::List split_iv(const Rcpp::IntegerVector& x, const Rcpp::IntegerVector& f);
Rcpp::List split_nm(const Rcpp::NumericMatrix& x, const Rcpp::IntegerVector& f);
Rcpp::NumericVector logit_rcpp(Rcpp::NumericVector p);
Rcpp::NumericVector expit_rcpp(Rcpp::NumericVector x);
double logit_double(const double& p);
double expit_double(const double& x);
arma::vec nm_omit(const Rcpp::NumericMatrix& x, const int& rownum);
Rcpp::NumericMatrix param_long(const Rcpp::NumericMatrix& prm_wide, const Rcpp::IntegerVector& narms, const bool& rowindex);
Rcpp::NumericMatrix param_wide(const Rcpp::NumericMatrix& prm_long, const Rcpp::IntegerVector& narms, const Rcpp::IntegerVector& trt, const Rcpp::IntegerVector& baseline);
arma::mat list_mat(const Rcpp::List& X);
arma::mat diag_tri(const arma::mat& A);
arma::rowvec get_upper_tri(const arma::mat& A, const bool& main);
arma::mat put_upper_tri(const arma::rowvec& a, const bool& main);
arma::mat cube_to_mat(const arma::cube& X, const bool& is_d, const int& ref);
arma::vec mat_to_vec(const arma::mat& X, const bool& is_d, const int& ref);
arma::mat vec_to_mat(const arma::vec& x, const int& nc, const bool& is_d, const int& ref);
arma::vec Sigma_M_to_beta(const arma::mat& A);
arma::mat beta_to_Sigma_M(const arma::vec& beta, const int& M);
bool is_symmetric(const arma::mat& m);
bool is_correlation(const arma::mat& A);
bool is_positive_definite(const arma::mat& A, const int& line);
arma::mat make_positive_definite(const arma::mat& A);
bool is_singular(const arma::mat& A);
arma::vec ols_coef(const double& xmin, const double& xmax, const Rcpp::List& args, const bool& delta_par);
double ols_pred(const arma::vec& coef, const double& x);
arma::mat cov2cor_rcpp(const arma::mat& V);
arma::mat spearman_mcmc(const arma::cube& Gamma_chain, const double& n, const double& M);

// DISTRIBUTION FUNCTIONS -----------------------------------------------------
arma::vec dmvn_arma(const arma::mat& x, const arma::vec& mean, const arma::mat& sigma, const  bool& logd);
arma::mat rmvn_arma(const int& n, const arma::vec& mu, const arma::mat& sigma);
double dinvwish_arma(const arma::mat& IW, const int& nu, const arma::mat& S, const bool& logd);
arma::mat rinvwish_arma(const int& nu, const arma::mat& S);
double dlkj_arma(const arma::mat& R, const double& eta, const bool& logd);
arma::mat rlkj_arma(const int& K, const double& eta);
Rcpp::NumericVector rtruncnorm_rcpp(const int& n, const double& a, const double& b, const double& mean, const double& sd);
Rcpp::NumericVector rtruncnorm2_rcpp(const int& n, const double& a, const double& b, const double& mean, const double& sd);
Rcpp::NumericVector dtruncnorm_rcpp(const Rcpp::NumericVector& x, const double& a, const double& b, const double& mean, const double& sd);
double dlogchol_arma(const arma::mat& A, const double& sigma_r, const bool& logd);
arma::mat rlogchol_arma(const int& M, const double& sigma_r);
Rcpp::NumericVector dmvt_arma(const arma::mat& x, const arma::vec& mean, const arma::mat& sigma, const int& df, const bool& logd);
arma::mat rmvt_arma(const int& n, const arma::vec& mean, const arma::mat& sigma, const int& df);
double dmatvn_arma(const arma::mat& X, const arma::mat& M, const arma::mat& U, const arma::mat& V, const bool& logd);
arma::mat rmatvn_arma(const arma::mat& M, const arma::mat& U, const arma::mat& V);
double dmatvt_arma(const arma::mat& X, const arma::mat& M, const arma::mat& Omega, const arma::mat& Sigma, const int& df, const bool& logd);
arma::mat rmatvt_arma(const arma::mat& M, const arma::mat& Sigma, const arma::mat& Omega, const int& df);
arma::vec dinvgamma_rcpp(const arma::vec& x, const double& alpha, const double& beta, const bool& logd);
arma::vec rinvgamma_rcpp(const int& n, const double& alpha, const double& beta);

// GAUSSIAN COPULA FUNCTIONS --------------------------------------------------
double gausscopdens(const Rcpp::NumericVector& u, const arma::mat& Gamma, const bool& is_u, const bool& logd);

// MCMC FUNCTIONS -------------------------------------------------------------
double nc_loglik(const Rcpp::NumericMatrix& y, const Rcpp::NumericMatrix& n, const Rcpp::NumericMatrix& x, const Rcpp::IntegerVector& trt, const Rcpp::NumericMatrix& mu, const Rcpp::NumericMatrix& delta, const Rcpp::List& Gamma);
double nc_logprior(const Rcpp::NumericMatrix& mu, const double& mu_sigma, const arma::mat& d, const double& d_sigma, const arma::mat& Sigma_M, const double& sigma_r, const int& ref_trt);
int indic_a_b(const double& y_ikm, const int& n_ikm, const double& x_ikm, const double& mu_ikm, const double& delta_ikm);
Rcpp::NumericMatrix x_imputed(const Rcpp::NumericMatrix& x, const Rcpp::List& Gamma, const Rcpp::IntegerVector& trt);
Rcpp::NumericMatrix n_imputed(const Rcpp::NumericMatrix& n_data);Rcpp::NumericMatrix y_imputed(const Rcpp::NumericMatrix& y, const Rcpp::NumericMatrix& x_imp, const Rcpp::IntegerVector& narms, const Rcpp::NumericMatrix& mu, const Rcpp::NumericMatrix& delta, const Rcpp::NumericMatrix& n_imp);
double Gamma_logpost(const arma::mat& Gamma, const arma::mat& x, const double& eta);
double r_logpost(const arma::mat& Gamma, const arma::mat& x);
double logpost(const double& mu, const double& delta, const double& y, const double& n, const double& w, const double& gamma, const double& eps, const double& eps_ab);
double mudelta_logpost(const arma::vec& mudelta, const arma::mat& delta_arma, const arma::mat& y, const arma::mat& n, const arma::mat& x, const Rcpp::IntegerVector& baseline, const Rcpp::IntegerVector& trt, const Rcpp::IntegerVector& study, const Rcpp::IntegerVector& narms, const arma::mat& d, const arma::mat& Sigma_M, const Rcpp::List& Gamma, const int& m, const double& mu_sigma, const double& eps, const double& eps_ab);
double mudelta_logpost2(const arma::vec& mudelta, const arma::mat& delta_arma, const arma::mat& y, const arma::mat& n, const arma::mat& x, const Rcpp::IntegerVector& baseline, const Rcpp::IntegerVector& trt, const Rcpp::IntegerVector& study, const Rcpp::IntegerVector& narms, const arma::mat& d, const arma::mat& Sigma_M, const Rcpp::List& Gamma, const int& i, const int& m, const double& mu_sigma, const double& eps, const double& eps_ab);
double mu_logpost(const double& mu, const arma::vec& delta, const arma::vec& y, const arma::vec& n, const arma::vec& w, const double& gamma, const double& mu_sigma, const double& eps, const double& eps_ab);
double mu_logpost2(const double& mu, const arma::vec& delta, const arma::vec& y, const arma::vec& n, const arma::vec& w, const double& gamma, const double& mu_sigma, const double& eps, const double& eps_ab);
double delta_logpost(const double& delta, const double& mu, const double& tau, const double& eta, const double& y, const double& n, const double& w, const double& gamma, const double& eps, const double& eps_ab);
double delta_logprior(const arma::mat& delta, const arma::mat& d, const arma::mat& Sigma_M, const Rcpp::IntegerVector& trt_arms, const Rcpp::IntegerVector& baseline, const Rcpp::IntegerVector& narms);
double d_logprior(const arma::mat& d, const double& d_sigma, const int& ref_trt);
double d_logpost_multi(const arma::mat& d, const arma::mat& delta, const arma::mat& Sigma_M, const Rcpp::IntegerVector& trt_arms, const Rcpp::IntegerVector& baseline, const Rcpp::IntegerVector& narms, const double& d_sigma, const int& ref_trt);
double d_logpost(const arma::mat& d, const arma::mat& delta, const arma::mat& Sigma_M, const Rcpp::IntegerVector& trt, const Rcpp::IntegerVector& baseline, const Rcpp::IntegerVector& narms_study, const double& d_sigma, const int& ref_trt);
double Sigma_M_logpost(const arma::mat& d, const arma::mat& delta, const arma::mat& Sigma_M, const Rcpp::IntegerVector& trt, const Rcpp::IntegerVector& baseline, const Rcpp::IntegerVector& narms_study, const double& sigma_r);
Rcpp::List nc_mcmc_opt(const Rcpp::List& data, const Rcpp::List& init, const int& totiter, const Rcpp::List& prior, const Rcpp::List& prop, const Rcpp::List& tuning, const Rcpp::List& adapt, const bool& verbose);
Rcpp::List nc_mcmc_mh(const Rcpp::List& data, const Rcpp::List& init, const int& totiter, const Rcpp::List& prior, const Rcpp::List& prop, const Rcpp::List& tuning, const Rcpp::List& adapt, const bool& verbose);
Rcpp::List nc_mcmc_mh_new(const Rcpp::List& data, const Rcpp::List& init, const int& totiter, const Rcpp::List& prior, const Rcpp::List& prop, const Rcpp::List& tuning, const Rcpp::List& adapt, const bool& verbose);
Rcpp::List nc_mcmc_mh_new2(const Rcpp::List& data, const Rcpp::List& init, const int& totiter, const Rcpp::List& prior, const Rcpp::List& prop, const Rcpp::List& tuning, const Rcpp::List& adapt, const bool& verbose);
Rcpp::List rwmh_adapt(const arma::mat& theta, const arma::vec& mu, const double& rho, const arma::mat& cov, const arma::vec& ar, const double& alpha, const double& beta, const double& gamma, const double& tar, const int& k, const bool& iter_cols, const int& what, const bool& diagonal);

// OPTIM FUNCTIONS ------------------------------------------------------------
void nmmin_rcpp(int n, Rcpp::NumericVector Bvec, Rcpp::NumericVector X, double *Fmin, int *fail, double abstol, double intol, const Rcpp::List& ex, double alpha, double bet, double gamm, int trace, int *fncount, int maxit);
Rcpp::List optim_rcpp(const Rcpp::NumericVector& par, const Rcpp::Function& fn, const Rcpp::List& args, const Rcpp::List& options, const bool& hessian);
Rcpp::NumericMatrix optimhess_rcpp(const Rcpp::NumericVector& par, const Rcpp::Function& fn, const Rcpp::List& args, const Rcpp::List& options);
Rcpp::List laplace_rcpp(const Rcpp::Function& logpost, const Rcpp::NumericVector& guess, const Rcpp::List& args, const Rcpp::List& options);
Rcpp::List optim_rcpp_example();
double optimize_rcpp(const Rcpp::Function& fn, const double& xmin, const double& xmax, const double& tol, const Rcpp::List& args);
Rcpp::List laplace_u_rcpp(const Rcpp::Function& fn, const double& xmin, const double& xmax, const double& tol, const Rcpp::List& args, const Rcpp::List& options);

// BOOST FUNCTIONS ------------------------------------------------------------
Rcpp::NumericVector qnorm_boost(const Rcpp::NumericVector& p, const double& mean, const double& sd, const bool& lower_tail);
Rcpp::NumericVector pnorm_boost(const Rcpp::NumericVector& x, const double& mean, const double& sd, const bool& lower_tail);
Rcpp::NumericVector qbinom_boost(const Rcpp::NumericVector& p, const int& n, const double& prob, const bool& lower_tail);
Rcpp::NumericVector pbinom_boost(const Rcpp::NumericVector& x, const int& n, const double& prob, const bool& lower_tail);

// MANIPULATION FUNCTIONS -----------------------------------------------------
arma::vec get_elements(const arma::mat& x, const arma::uvec& row_ind, const arma::uvec& col_ind);
arma::mat rev_col_subset(arma::mat x, unsigned int start, unsigned int end);
arma::mat rev_row_subset(arma::mat x, unsigned int start, unsigned int end);
arma::vec reverse_vec(arma::vec x);
arma::mat field_to_matrix(arma::field<arma::vec> x);
double sum_field_vec(const arma::field<arma::vec>& x);

#endif
