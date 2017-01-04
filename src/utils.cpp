// utils.cpp

/* do this first to get the right options for math.h */
#include <R_ext/Arith.h>

#include <algorithm>
#include <Rmath.h>
#include "netcopula.h"

#define both_non_NA(a,b) (!ISNAN(a) && !ISNAN(b))

// [[Rcpp::depends("RcppArmadillo")]]

void tableC(int* counts, const int* x, int nelem, int ndistelem){
  int* xtmp = new int[nelem];
  for(int i = 0; i < nelem; i++){
    xtmp[i] = x[i];
  }
  // R_isort(xtmp, nelem);
  R_qsort_int(xtmp, 1, nelem);
  
  int xlead = 1;
  counts[xlead - 1] = 0;
  for(int i = 0; i < nelem; i++){
    if(xtmp[i] != xlead){
      xlead++;
      counts[xlead - 1] = 0;
    }
    counts[xlead - 1]++;
  }

  delete[] xtmp;
}

// THE FOLLOWING FUNCTIONS ARE EXPORTED ONLY FOR CHECKING.
// IN THE FINAL RELEASE THEY WILL NOT BE EXPORTED.

//' @export
// [[Rcpp::export]]
arma::mat mat_block_diag(const arma::mat& A, const int& n) {
  int d = A.n_cols;
  arma::mat A_block(n*d, n*d, arma::fill::zeros);
  for (unsigned int i = 0; i < n; i++) {
    A_block.submat(d*i, d*i, d*(i + 1) - 1, d*(i + 1) - 1) = A;
  }
  return A_block;
}

//' @export
// [[Rcpp::export]]
arma::mat Sigma_block(const arma::mat& Sigma_M, const int& n) {
  int M = Sigma_M.n_cols;
  arma::mat Sigma_block(n*M, n*M, arma::fill::zeros);
  for (unsigned int i = 0; i < n; i++) {
    Sigma_block.submat(M*i, M*i, M*(i + 1) - 1, M*(i + 1) - 1) = Sigma_M;
    for (unsigned int j = (i + 1); j < n; j++) {
      Sigma_block.submat(M*i, M*j, M*(i + 1) - 1, M*(j + 1) - 1) = 0.5*Sigma_M;
      Sigma_block.submat(M*j, M*i, M*(j + 1) - 1, M*(i + 1) - 1) = 0.5*Sigma_M;
    }
  }
  return Sigma_block;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix df_nm(const Rcpp::DataFrame& x, const Rcpp::IntegerVector& cols) {
  int nrows = x.nrows(), ncols = cols.size();
  Rcpp::NumericMatrix y(nrows, ncols);
  for (unsigned int j = 0; j < ncols; j++) {
    y(Rcpp::_, j) = Rcpp::NumericVector(x[cols(j)]);
  }
  return y;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector nm_stack(const Rcpp::NumericMatrix& x) {
  int nr = x.nrow(), nc = x.ncol();
  Rcpp::NumericVector y(nr*nc);
  for (unsigned int i = 0; i < nr; i++) {
    for (unsigned int j = 0; j < nc; j++) {
      y(j + nc*i) = x(i, j);
    }
  }
  return y;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix nv_unstack(const Rcpp::NumericVector& x, const int& nc) {
  int n = x.size(), nr = n/nc;
  Rcpp::NumericMatrix y(nr, nc);
  for (unsigned int i = 0; i < nr; i++) {
    for (unsigned int j = 0; j < nc; j++) {
      y(i, j) = x(j + nc*i);
    }
  }
  return y;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector nv_omit(const Rcpp::NumericVector& x) {
  return x[!is_na(x)];
}

//' @export
// [[Rcpp::export]]
arma::uvec nv_na_index(const Rcpp::NumericVector& x, const int& dim, const bool& type) {
  int n = x.size(), cum = 0;
  Rcpp::LogicalVector x_na = is_na(x);
  arma::uvec y(dim);
  for (unsigned int i = 0; i < n; i++) {
    if (x_na(i) == type) {
      y(cum) = i;
      cum++;
    }
  }
  return y;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector nv_miss_replace(const Rcpp::NumericVector& x, const arma::vec& miss, const arma::uvec& miss_i) {
  int n_miss = miss_i.n_elem;
  Rcpp::NumericVector y = x;
  for (unsigned int i = 0; i < n_miss; i++) {
    y(miss_i(i)) = miss(i);
  }
  return y;
}

//' @export
// [[Rcpp::export]]
Rcpp::List split_iv(const Rcpp::IntegerVector& x, const Rcpp::IntegerVector& f) {
  if (f.size() <= 0) {
    Rcpp::stop("f must be a vector with at least one element.");
  }
  if (x.size() <= 0) {
    Rcpp::stop("x must be a vector with at least one element.");
  }
  if (x.size() != f.size()) {
    Rcpp::stop("x and f must be of the same size.");
  }

  Rcpp::IntegerVector f_tbl = table(f);
  int nf = f_tbl.size(), cum = 0, nx = x.size();
  Rcpp::List res(nf);
  std::vector<int> tmp;
  for (unsigned int i = 0; i < nf; i++) {
    tmp.resize(f_tbl[i]);
    for (unsigned int j = 0; j < nx; j++) {
      if (f[j] == (i + 1)) {
        tmp[cum] = x[j];
        cum++;
      }
    }
    res[i] = tmp;
    cum = 0;
    tmp.clear();
  }

  return res;
}

//' @export
// [[Rcpp::export]]
Rcpp::List split_nm(const Rcpp::NumericMatrix& x, const Rcpp::IntegerVector& f) {
  if (f.size() <= 0) {
    Rcpp::stop("f must be a vector with at least one element.");
  }
  if (x.nrow() <= 0) {
    Rcpp::stop("x must be a data frame with at least one row.");
  }
  if (x.nrow() != f.size()) {
    Rcpp::stop("number of rows of x and length of f must be equal.");
  }

  Rcpp::IntegerVector f_tbl = table(f), idx;
  int nf = f_tbl.size(), xr = x.nrow(), xc = x.ncol();
  arma::mat tmp;
  Rcpp::List res(nf), split_indices = split_iv(Rcpp::seq_len(xr), f);
  for (unsigned int k = 0; k < nf; k++) {
    idx = Rcpp::as<Rcpp::IntegerVector>(split_indices[k]);
    tmp.resize(idx.size(), xc);
    for (unsigned int i = 0; i < idx.size(); i++) {
      for (unsigned int j = 0; j < xc; j++) {
        tmp(i, j) = x(idx[i] - 1, j);
      }
    }
    res[k] = tmp;
  }

  return res;
}

double logit_double(const double& p) {
  return log(p/(1 - p));
}

double expit_double(const double& x) {
  return 1/(1 + exp(-x));
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector logit_rcpp(Rcpp::NumericVector p) {
  return log(p/(1 - p));
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector expit_rcpp(Rcpp::NumericVector x) {
  return 1/(1 + exp(-x));
}

//' @export
// [[Rcpp::export]]
arma::vec nm_omit(const Rcpp::NumericMatrix& x, const int& rownum) {
  int n_col = x.ncol(), cum = 0;
  arma::vec y(n_col);
  for (unsigned int j = 0; j < n_col; j++) {
    if (!ISNAN(x(rownum, j))) {
      y(cum) = x(rownum, j);
      cum++;
    }
  }
  y.resize(cum);

  return y;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix param_long(const Rcpp::NumericMatrix& prm_wide, const Rcpp::IntegerVector& narms, const bool& rowindex) {
  int n = narms.size(), M = prm_wide.ncol(), cumul = 0, cumul_arms = 1;
  Rcpp::NumericMatrix prm_long(sum(narms), M);

  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int k = 0; k < narms(i); k++) {
      for (unsigned int j = 0; j < M; j++) {
        if (rowindex) {
          prm_long(cumul, j) = cumul_arms;
        } else {
          prm_long(cumul, j) = prm_wide(i, j);
        }
      }
      cumul++;
    }
    cumul_arms += narms(i);
  }

  return prm_long;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix param_wide(const Rcpp::NumericMatrix& prm_long, const Rcpp::IntegerVector& narms, const Rcpp::IntegerVector& trt, const Rcpp::IntegerVector& baseline) {
  int n = narms.size(), ndata = prm_long.nrow(), M = prm_long.ncol(), cum = 0;
  Rcpp::NumericMatrix prm_wide(n, M);

  for (unsigned int i = 0; i < ndata; i++) {
    if (trt(i) == baseline(i)) {
      for (unsigned int j = 0; j < M; j++) {
        prm_wide(cum, j) = prm_long(i, j);
      }
      cum++;
    }
  }

  return prm_wide;
}

//' @export
// [[Rcpp::export]]
arma::mat list_mat(const Rcpp::List& X) {
  arma::mat tmp = Rcpp::as<arma::mat>(X(0));
  int nr = X.size(), M = tmp.n_rows, nc = M*(M - 1)/2, count = 0;
  arma::mat Y(nr, nc);

  for (unsigned int k = 0; k < nr; k++) {
    tmp = Rcpp::as<arma::mat>(X(k));
    for (unsigned int j = 0; j < (M - 1); j++) {
      for (unsigned int i = (j + 1); i < M; i++) {
        Y(k, count) = tmp(i, j);
        count++;
      }
    }
    count = 0;
  }

  return Y;
}

//' @export
// [[Rcpp::export]]
arma::mat diag_tri(const arma::mat& A) {
  int M = A.n_rows, count = 0;
  arma::mat out(1, M*(M + 1)/2);

  for (unsigned int j = 0; j < M; j++) {
    for (unsigned int i = j; i < M; i++) {
      out(0, count) = A(i, j);
      count++;
    }
  }

  return out;
}

//' @export
// [[Rcpp::export]]
arma::mat cube_to_mat(const arma::cube& X, const bool& is_d, const int& ref) {
  int nr = X.n_rows, nc = X.n_cols, ns = X.n_slices, nc_out = (nr - 1)*nc, count = 0;
  if (!is_d) {
    nc_out += nc;
  }
  arma::mat out(ns, nc_out);
  arma::mat tmp(nr, nc);
  arma::uvec row_id(nr - 1);

  if (is_d) {
    for (unsigned int k = 0; k < nr; k++) {
      if (k != (ref - 1)) {
        row_id(count) = k;
        count++;
      }
    }
    for (unsigned int s = 0; s < ns; s++) {
      tmp = X.slice(s);
      out.row(s) = arma::vectorise(tmp.rows(row_id), 1);
    }
  } else {
    for (unsigned int s = 0; s < ns; s++) {
      tmp = X.slice(s);
      out.row(s) = arma::vectorise(tmp, 1);
    }
  }

  return out;
}

//' @export
// [[Rcpp::export]]
arma::vec mat_to_vec(const arma::mat& X, const bool& is_d, const int& ref) {
  int nr = X.n_rows, nc = X.n_cols, n_out = (nr - 1)*nc, count = 0;
  if (!is_d) {
    n_out += nc;
  }
  arma::vec out(n_out);
  arma::uvec row_id(nr - 1);

  if (is_d) {
    for (unsigned int k = 0; k < nr; k++) {
      if (k != (ref - 1)) {
        row_id(count) = k;
        count++;
      }
    }
    out = arma::trans(arma::vectorise(X.rows(row_id), 1));
  } else {
    out = arma::trans(arma::vectorise(X, 1));
  }

  return out;
}

//' @export
// [[Rcpp::export]]
arma::mat vec_to_mat(const arma::vec& x, const int& nc, const bool& is_d, const int& ref) {
  int n = x.n_elem, nr_out = (n/nc);
  arma::mat out(nr_out, nc);

  for (unsigned int k = 0; k < nr_out; k++) {
    out.row(k) = arma::trans(x.subvec(nc*k, nc*(k + 1) - 1));
  }

  if (is_d) {
    out.insert_rows(ref - 1, 1);
  }

  return out;
}

//' @export
// [[Rcpp::export]]
arma::vec Sigma_M_to_beta(const arma::mat& A) {
  int M = A.n_cols, nn = M*(M + 1)/2;
  
  arma::mat A_inv = arma::solve(A, arma::eye(M, M));
  // if (!is_positive_definite(A_inv, 391)) {
  //   Rprintf("Sigma_M_to_beta\n");
  //   A_inv = make_positive_definite(A_inv);
  // }
  arma::mat R = arma::chol(A_inv);
  R.diag() = arma::log(R.diag());
  arma::vec beta(nn);
  int position = 0;
  for (int m = 0; m < M; m++) {
    beta(arma::span(position, position + m)) = R(arma::span(0, m), m);
    position += (m + 1);
  }

  return beta;
}

//' @export
// [[Rcpp::export]]
arma::mat beta_to_Sigma_M(const arma::vec& beta, const int& M) {
  int position = 0;
  arma::mat ret = arma::zeros(M, M);  
  for (unsigned int m = 0; m < M; m++) {
    ret(arma::span(0, m), m) = beta(arma::span(position, position + m));
    position += (m + 1);
  }
  ret.diag() = arma::exp(ret.diag());
  arma::mat R = arma::trimatu(ret);
  arma::mat R_inv = arma::solve(R, arma::eye(M, M));

  return (R_inv * arma::trans(R_inv));
}

//' @export
// [[Rcpp::export]]
bool is_symmetric(const arma::mat& A) {
  arma::umat Z = (A == arma::trans(A));

  return arma::all(arma::vectorise(Z) == 1);
}

//' @export
// [[Rcpp::export]]
bool is_correlation(const arma::mat& A) {
  return (arma::all(arma::diagvec(A) == 1) && is_positive_definite(A, 434));
}

//' @export
// [[Rcpp::export]]
bool is_positive_definite(const arma::mat& A, const int& line) {
  int d = A.n_rows;

  // Rprintf("is_positive_definite_IN_%d\n", line);
  arma::vec eigval = arma::eig_sym(A);;
  // Rprintf("is_positive_definite_OUT_%d\n", line);
  double tol = d*arma::max(arma::abs(eigval))*machine_eps;

  return arma::all(eigval > tol);
}

//' @export
// [[Rcpp::export]]
arma::mat make_positive_definite(const arma::mat& A) {
  // (from corpcor::make.positive.definite)
  // (see also sfsmisc::posdefify and Matrix::nearPD)
  int d = A.n_rows;
  arma::mat dA(A), eigvec;
  arma::vec eigval, tau(d);

  arma::eig_sym(eigval, eigvec, A);

  double tol_eig = d*arma::max(arma::abs(eigval))*machine_eps;
  for (unsigned int i = 0; i < d; i++) {
    tau(i) = fmax(0, 2*tol_eig - eigval(i));
  }
  dA = eigvec*arma::diagmat(tau)*arma::trans(eigvec);

  return A + dA;
}

//' @export
// [[Rcpp::export]]
bool is_singular(const arma::mat& A) {
  return ((arma::det(A) == 0) ? true : false);
}

//' @export
// [[Rcpp::export]]
arma::vec ols_coef(const double& xmin, const double& xmax, const Rcpp::List& args, const bool& delta_par) {
  int n_val = 1e+02, count = 0;
  double step = (xmax - xmin)/static_cast<double>(n_val);
  arma::vec x_grid(n_val + 1), lp(n_val + 1);
  double mu = 0.0, delta = 0.0, tau = 0.0, eta = 0.0, y = 0.0, n = 0.0, w = 0.0, gamma = 0.0, mu_sigma = 0.0, eps = 0.0, eps_ab = 0.0;
  arma::vec y_val(1);
  arma::mat X(1, 3);
  if (delta_par) {
    mu = Rcpp::as<double>(args["mu"]);
    tau = Rcpp::as<double>(args["tau"]);
    eta = Rcpp::as<double>(args["eta"]);
    y = Rcpp::as<double>(args["y"]);
    n = Rcpp::as<double>(args["n"]);
    w = Rcpp::as<double>(args["w"]);
    gamma = Rcpp::as<double>(args["gamma"]);
    eps = Rcpp::as<double>(args["eps"]);
    eps_ab = Rcpp::as<double>(args["eps_ab"]);
  } else {
    delta = Rcpp::as<double>(args["delta"]);
    y = Rcpp::as<double>(args["y"]);
    n = Rcpp::as<double>(args["n"]);
    w = Rcpp::as<double>(args["w"]);
    gamma = Rcpp::as<double>(args["gamma"]);
    mu_sigma = Rcpp::as<double>(args["mu_sigma"]);
    eps = Rcpp::as<double>(args["eps"]);
    eps_ab = Rcpp::as<double>(args["eps_ab"]);
  }
  for (unsigned int i = 0; i < (n_val + 1); i++) {
    x_grid(i) = xmin + i*step;
    if (delta_par) {
      lp(i) = delta_logpost(x_grid(i), mu, tau, eta, y, n, w, gamma, eps, eps_ab);
    } else {
      lp(i) = mu_logpost(x_grid(i), delta, y, n, w, gamma, mu_sigma, eps, eps_ab);
    }
    if (!ISNA(lp(i))) {
      X.resize(count + 1, 3);
      y_val.resize(count + 1);
      X(count, 0) = 1;
      X(count, 1) = x_grid(i);
      X(count, 2) = pow(x_grid(i), 2);
      y_val(count) = lp(i);
      count++;
    }
  }

  arma::vec coef(3);
  if (count > 1) {
    coef = arma::solve(X, y_val);
  } else {
    coef.fill(NA_REAL);
  }

  return coef;
}

//' @export
// [[Rcpp::export]]
double ols_pred(const arma::vec& coef, const double& x) {
  arma::rowvec X(3);
  X(0) = 1;
  X(1) = x;
  X(2) = pow(x, 2);

  double pred = arma::as_scalar(X*coef);
  if (ISNAN(pred)) {
    coef.print();
    X.print();
  }

  return pred;
}

//' @export
// [[Rcpp::export]]
arma::mat cov2cor_rcpp(const arma::mat& V) {
  int p = V.n_rows, d = V.n_cols;
  if (p != d) {
      Rcpp::stop("'V' is not a square numeric matrix.");
  }
  arma::mat Is = arma::diagmat(arma::sqrt(1/arma::diagvec(V)));
  arma::mat R = Is * V * Is;
  R.diag().ones();

  return R;
}
