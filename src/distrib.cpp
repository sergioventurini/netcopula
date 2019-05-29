#include "netcopula.h"

#ifndef MAX
#define MAX(A, B) ((A > B) ? (A) : (B))
#endif

static const double t1 = 0.15; // this comes from Geweke (1991), but...
// static const double t1 = 0.375; // ...this is reported in Geweke (pag. 113)
static const double t2 = 2.18;
static const double t3 = 0.725;
static const double t4 = 0.45;

// [[Rcpp::depends("RcppArmadillo")]]

//' Multivariate normal density.
//'
//' Evaluate the multivariate normal density.
//'
//' @param x A numeric matrix whose rows represent the points at which to
//' evaluate the density.
//' @param mean A numeric vector containing the univariate means.
//' @param sigma A positive definite numeric matrix representing the covariance
//' matrix of the distribution.
//' @param logd Boolean length-one vector; if TRUE the log density is returned.
//'
//' @return A numeric vector.
//' @export
//'
//' @references
//' \url{http://gallery.rcpp.org/articles/dmvnorm_arma/}
//' 
//' @examples
//' set.seed(123)
//' ### Covariance matrix and mean vector
//' sigma <- bayesm::rwishart(10, diag(8))$IW
//' mean <- rnorm(8)
//' ### Benchmarking
//' n <- 1e+4
//' X <- mvtnorm::rmvnorm(n, mean, sigma)
//' require(rbenchmark)
//' benchmark(mvtnorm::dmvnorm(X, mean, sigma, log = FALSE), 
//'           dmvn_arma(X, mean, sigma, FALSE),
//'           order = "relative", replications = 100)[, 1:4]
// [[Rcpp::export]]
arma::vec dmvn_arma(const arma::mat& x, const arma::vec& mean, const arma::mat& sigma, const bool& logd = false) {
  int n = x.n_rows, d = x.n_cols;

  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(arma::trimatu(arma::chol(sigma))));
  double rootisum = sum(log(rooti.diag()));
  double constant = -(d/2.0) * log_2pi;

  for (unsigned int i = 0; i < n; i++) {
    arma::vec z = rooti * arma::trans(x.row(i) - arma::trans(mean));
    out(i) = constant - 0.5 * sum(z%z) + rootisum;
  }

  if (!logd) {
    out = exp(out);
  }

  return out;
}

//' Multivariate normal variates generation.
//'
//' Generation of multivariate normal variates.
//'
//' @param n A length-one numeric vector providing the number of draws to
//' generate.
//' @param mean A numeric vector containing the univariate means.
//' @param sigma A positive definite numeric matrix representing the covariance
//' matrix of the distribution.
//'
//' @return A numeric matrix.
//' @export
//'
//' @references
//' \url{http://gallery.rcpp.org/articles/simulate-multivariate-normal/}
//' 
//' @examples
//' set.seed(123)
//' ### Covariance matrix and mean vector
//' sigma <- matrix(c(1, 0.9, -0.3, 0.9, 1, -0.4, -0.3, -0.4, 1), ncol = 3)
//' mu <- c(10, 5, -3)
//' ### Benchmarking
//' n <- 1e+4
//' require(rbenchmark)
//' benchmark(mvtnorm::rmvnorm(n, mu, sigma),
//'           MASS::mvrnorm(n, mu, sigma),
//'           rmvn_arma(n, mu, sigma),
//'           columns = c("test", "replications", "relative", "elapsed"),
//'           order = "relative")
// [[Rcpp::export]]
arma::mat rmvn_arma(const int& n, const arma::vec& mean, const arma::mat& sigma) {
  int d = sigma.n_cols;
  // if (!is_positive_definite(sigma, 99)) {
  //   Rcpp::stop("non-positive definite sigma in rmvn_arma().\n");
  // }

  // GetRNGstate();  // http://gallery.rcpp.org/articles/random-number-generation/
  arma::mat Y = Rcpp::rnorm(n*d);
  Y = arma::reshape(Y, n, d);
  // PutRNGstate();

  // arma::mat Y = arma::randn(n, d);
  return arma::repmat(mean, 1, n).t() + Y * arma::chol(sigma);
}

//' Inverse Wishart density.
//'
//' Evaluation of the inverse Wishart density.
//'
//' @param IW A numeric matrix at which to evaluate the density.
//' @param nu Integer length-one vector providing the number of degrees of
//' freedom.
//' @param S Symmetric, positive definite numeric matrix representing the
//' distribution scale matrix.
//' @param logd Boolean length-one vector; if TRUE the log density is returned.
//'
//' @return A length-one numeric vector.
//' @export
//'
//' @examples
//' set.seed(123)
//' ### Parameters
//' k <- 8
//' S <- bayesm::rwishart(10, diag(k))$IW
//' nu <- round(runif(1, min = (k + 1), max = 30))
//' ### Benchmarking
//' n <- 1e+2
//' X <- array(NA, dim = c(n, k, k))
//' for (i in 1:n) X[i, , ] <- MCMCpack::riwish(nu, S)
//' require(rbenchmark)
//' benchmark(apply(X, 1, MCMCpack::diwish, nu, S),
//'           apply(X, 1, dinvwish_arma, nu, S, FALSE),
//'           order = "relative", replications = 100)[, 1:4]
// [[Rcpp::export]]
double dinvwish_arma(const arma::mat& IW, const int& nu, const arma::mat& S, const bool& logd = false) {
  int k = S.n_cols;

  double ldet_S = log(arma::det(S));
  double ldet_IW = log(arma::det(IW));

  long double sum_lgamma = 0;
  for (unsigned int i = 0; i < k; i++) {
    sum_lgamma += R::lgammafn((nu + 1 - (i + 1))/2.0);
  }
  double cnst = ((nu*k)/2.0)*log(2.0) + ((k*(k - 1))/4.0)*log(M_PI) + sum_lgamma;

  arma::mat IW_inv = arma::inv(IW);

  double out = -cnst + (nu*ldet_S - (nu + k + 1)*ldet_IW - arma::trace((S * IW_inv)))/2.0;

  if (!logd) {
    out = exp(out);
  }

  return out;
}

//' Inverse Wishart variates generation.
//'
//' Generate a single inverse Wishart variate.
//'
//' @param nu Integer length-one vector providing the number of degrees of
//' freedom.
//' @param S Symmetric, positive definite numeric matrix representing the
//' distribution scale matrix.
//'
//' @return A length-one numeric vector.
//' @export
//'
//' @references
//' \url{https://en.wikipedia.org/wiki/
//' Wishart_distribution#Bartlett_decomposition}
//'
//' @examples
//' set.seed(123)
//' ### Parameters
//' k <- 8
//' S <- bayesm::rwishart(10, diag(k))$IW
//' nu <- round(runif(1, min = (k + 1), max = 30))
//' ### Benchmarking
//' n <- 1e+2
//' require(rbenchmark)
//' benchmark(apply(rWishart(n, nu, S), 3, solve),
//'           replicate(n, MCMCpack::riwish(nu, S), simplify = "array"),
//'           replicate(n, rinvwish_arma(nu, S), simplify = "array"),
//'           columns = c("test", "replications", "relative", "elapsed"),
//'           order = "relative")
// [[Rcpp::export]]
arma::mat rinvwish_arma(const int& nu, const arma::mat& S) {
  int k = S.n_rows;

  if (S.n_rows != S.n_cols) {
    Rcpp::stop("S is not square in rinvwish_arma().\n");
  }
  if (nu < k) {
    Rcpp::stop("nu is less than the dimension of S in rinvwish_arma().\n");
  }

  // A has sqrt of chi-squares on diagonal and normals below diagonal
  arma::mat A = arma::zeros(k, k);
  for (unsigned int j = 0; j < k; j++) {
    A(j, j) = sqrt(R::rchisq(nu - j));
  }

  for (unsigned int j = 0; j < (k - 1); j++) {
    for (unsigned int i = (j + 1); i < k; i++) {
      A(i, j) = R::rnorm(0.0, 1.0);
    }
  }
  
  // L is the upper triangular root of a Wishart therefore, W = L'L
  // This is the UL (not LU!) decomposition W^-1 = (L^-1)(L^-1)'
  arma::mat L = arma::trans(A) * arma::chol(arma::inv_sympd(S));
  arma::mat L_inv = arma::solve(arma::trimatu(L), arma::eye(k, k)); // trimatu interprets L as a triangular matrix and allows a faster and more accurate inversion (http://arma.sourceforge.net/docs.html#solve)
  
  return (L_inv * arma::trans(L_inv));
}

//' Compute the density function of a Lewandowski-Kurowicka-Joe distribution.
//'
//' Compute the density function of a Lewandowski-Kurowicka-Joe distribution
//' for a correlation matrix.
//'
//' @param R A numeric matrix that representing the correlation matrix.
//' @param eta Length-one numeric vector representing the parameter of the
//' distribution.
//'
//' @return A length-one numeric vector.
//' @export
//'
//' @references
//' Lewandowski, D., Kurowicka, D., and Joe, H. (2009). Generating random
//' correlation matrices based on vines and extended onion method. Journal of
//' Multivariate Analysis, 100:1989-2001.
//' 
//' @details
//' This function is a port of the corresponding Stan (\url{http://mc-stan.
//' org}) function developed using the Eigen library (\url{http://eigen.
//' tuxfamily.org}).
//' 
//' @examples
//' set.seed(123)
//' K <- 5
//' R <- rlkj_arma(K, .5) # eta = 0.5
//' dlkj_arma(R, .5)
//' R <- rlkj_arma(K, 1) # eta = 1
//' dlkj_arma(R, 1)
//' R <- rlkj_arma(K, 2) # eta = 2
//' dlkj_arma(R, 2)
// [[Rcpp::export]]
double dlkj_arma(const arma::mat& R, const double& eta, const bool& logd = false) {
  double constant = 0, out = 0;
  int K = R.n_rows, Km1 = K - 1;
  if (eta == 1.0) {
    // C++ integer division is appropriate in this block
    arma::vec numerator(Km1/2);
    for (unsigned int k = 1; k <= numerator.n_rows; k++) {
      numerator(k - 1) = R::lgammafn(2 * k);
    }
    constant = arma::sum(numerator);
    if ((K % 2) == 1) {
      constant += 0.25 * (K * K - 1) * log_pi - 0.25 * (Km1 * Km1) * log_two - Km1 * R::lgammafn((K + 1) / 2);
    } else {
      constant += 0.25 * K * (K - 2) * log_pi + 0.25 * (3 * K * K - 4 * K) * log_two + K * R::lgammafn(K / 2) - Km1 * R::lgammafn(K);
    }
  } else {
    constant = -Km1 * R::lgammafn(eta + 0.5 * Km1);
    for (unsigned int k = 1; k <= Km1; k++) {
      constant += 0.5 * k * log_pi + R::lgammafn(eta + 0.5 * (Km1 - k));
    }
  }

  // double constant2 = 0;
  // arma::vec part1(Km1), part2(Km1);
  // for (unsigned k = 1; k <= Km1; k++) {
  //   part1(k - 1) = (2 * eta - 2 + K - k) * (K - k);
  //   part2(k - 1) = (K - k) * R::lbeta(eta + 0.5 * (K - k - 1), eta + 0.5 * (K - k - 1));
  // }
  // constant2 = arma::sum(part1) * log_two + arma::sum(part2);

  // double logdet, logdet_sgn;
  // arma::log_det(logdet, logdet_sgn, R);
  // out += constant + (eta - 1.0) * logdet;

  arma::vec values = arma::log(arma::eig_sym(R));
  out += constant + (eta - 1.0) * arma::sum(values);

  if (!logd) {
    out = exp(out);
  }

  return out;
}

//' Generate random variate from a Lewandowski-Kurowicka-Joe distribution.
//'
//' Generation of a single draw from a Lewandowski-Kurowicka-Joe distribution
//' for a correlation matrix.
//'
//' @param K Length-one numeric vector that represents the size of the matrix.
//' @param eta Length-one numeric vector representing the parameter of the
//' distribution.
//'
//' @return A numeric (correlation) matrix.
//' @export
//'
//' @references
//' Lewandowski, D., Kurowicka, D., and Joe, H. (2009). Generating random
//' correlation matrices based on vines and extended onion method. Journal of
//' Multivariate Analysis, 100:1989-2001.
//' 
//' @details
//' This function is a port of the corresponding Stan (\url{http://mc-stan.
//' org}) function developed using the Eigen library (\url{http://eigen.
//' tuxfamily.org}).
//' 
//' @examples
//' set.seed(123)
//' K <- 5
//' rlkj_arma(K, .5) # eta = 0.5
//' rlkj_arma(K, 1) # eta = 1
//' rlkj_arma(K, 2) # eta = 2
// [[Rcpp::export]]
arma::mat rlkj_arma(const int& K, const double& eta) {
  arma::vec CPCs(K*(K - 1)/2);
  double alpha = eta + 0.5 * (K - 1);
  unsigned int count = 0;
  for (size_t i = 0; i < (K - 1); i++) {
    alpha -= 0.5;
    for (size_t j = i + 1; j < K; j++) {
      CPCs(count) = 2.0*R::rbeta(alpha, alpha) - 1.0;
      count++;
    }
  }

  arma::vec temp(K - 1);
  arma::vec acc(K - 1);
  acc.ones();
  // Cholesky factor of correlation matrix
  arma::mat L(K, K);
  L.zeros();

  size_t position = 0;
  size_t pull = (K - 1);

  L(0, 0) = 1.0;
  L.col(0).tail(pull) = temp = CPCs.head(pull);
  acc.tail(pull) = 1.0 - arma::square(temp);
  for (size_t i = 1; i < (K - 1); i++) {
    position += pull;
    pull--;
    temp.resize(pull);
    temp = CPCs(arma::span(position, position + pull - 1));
    L(i, i) = sqrt(acc(i - 1));
    L.col(i).tail(pull) = temp % sqrt(acc.tail(pull));
    acc.tail(pull) = acc.tail(pull) % (1.0 - arma::square(temp));
  }
  L(K - 1, K - 1) = sqrt(acc(K - 2));

  arma::mat L_lower_tri = arma::trimatl(L);
  return L_lower_tri * L_lower_tri.t();
}

/* Exponential rejection sampling (a,inf) */
static R_INLINE double ers_a_inf(double a) {
  const double ainv = 1.0 / a;
  double x, rho;
  do {
    x = R::rexp(ainv) + a; /* rexp works with 1/lambda */
    rho = exp(-0.5 * pow((x - a), 2));
  } while (R::runif(0.0, 1.0) > rho);
  return x;
}

/* Exponential rejection sampling (a,b) */
static R_INLINE double ers_a_b(double a, double b) {
  const double ainv = 1.0 / a;
  double x, rho;
  do {
    x = R::rexp(ainv) + a; /* rexp works with 1/lambda */
    rho = exp(-0.5 * pow((x-a), 2));
  } while (R::runif(0.0, 1.0) > rho || x > b);
  return x;
}

/* Normal rejection sampling (a,b) */
static R_INLINE double nrs_a_b(double a, double b){
  double x = -DBL_MAX;
  while (x < a || x > b) {
    x = R::rnorm(0.0, 1.0);
  }
  return x;
}

/* Normal rejection sampling (a,inf) */
static R_INLINE double nrs_a_inf(double a){
  double x = -DBL_MAX;
  while (x < a) {
    x = R::rnorm(0.0, 1.0);
  }
  return x;
}

/* Half-normal rejection sampling */
double hnrs_a_b(double a, double b){
  double x = a - 1.0;
  while (x < a || x > b) {
    x = R::rnorm(0.0, 1.0);
    x = fabs(x);
  }
  return x;
}

/* Uniform rejection sampling */
static R_INLINE double urs_a_b(double a, double b){
  const double phi_a = R::dnorm(a, 0.0, 1.0, false);
  double x = 0.0;
  
  /* Upper bound of normal density on [a, b] */
  const double ub = a < 0 && b > 0 ? M_1_SQRT_2PI : phi_a;
  do {
    x = R::runif(a, b);
  } while (R::runif(0.0, 1.0) * ub > R::dnorm(x, 0.0, 1.0, false));
  return x;
}

/* Previously this was referred to as type 1 sampling: */
static inline double r_lefttruncnorm(double a, double mean, double sd) {
  const double alpha = (a - mean) / sd;
  if (alpha < t4) {
    return mean + sd * nrs_a_inf(alpha);
  } else {
    return mean + sd * ers_a_inf(alpha);
  }
}

static R_INLINE double r_righttruncnorm(double b, double mean, double sd) {
  const double beta = (b - mean) / sd;
  /* Exploit symmetry: */
  return mean - sd * r_lefttruncnorm(-beta, 0.0, 1.0);
}

static R_INLINE double r_truncnorm(double a, double b, double mean, double sd) {
  const double alpha = (a - mean)/sd;
  const double beta = (b - mean)/sd;
  const double phi_a = R::dnorm(alpha, 0.0, 1.0, false);
  const double phi_b = R::dnorm(beta, 0.0, 1.0, false);

  if (beta <= alpha) {
    return NA_REAL;
  } else if (alpha <= 0 && 0 <= beta) { /* 2 */
    if (phi_a <= t1 || phi_b <= t1) { /* 2 (a) */
      return mean + sd * nrs_a_b(alpha, beta);
    } else { /* 2 (b) */
      return mean + sd * urs_a_b(alpha, beta);
    }
  } else if (alpha > 0) { /* 3 */
    if (phi_a / phi_b <= t2) { /* 3 (a) */
      return mean + sd * urs_a_b(alpha, beta);
    } else {
      if (alpha < t3) { /* 3 (b) */                
        return mean + sd * hnrs_a_b(alpha, beta);
      } else { /* 3 (c) */
        return mean + sd * ers_a_b(alpha, beta);
      }
    }
  } else { /* 3s */
    if (phi_b / phi_a <= t2) { /* 3s (a) */
      return mean - sd * urs_a_b(-beta, -alpha);
    } else {
      if (beta > -t3) { /* 3s (b) */
        return mean - sd * hnrs_a_b(-beta, -alpha);
      } else { /* 3s (c) */
        return mean - sd * ers_a_b(-beta, -alpha);
      }
    }
  }
}

//' The truncated univariate normal distribution.
//'
//' Generation of draws from a truncated normal distribution with mean equal
//' to \code{mean} and standard deviation equal to \code{sd}.
//'
//' @param n Length-one numeric vector representing the number of draws.
//' @param a A numeric vector providing the lower bounds. These may be
//' \code{-Inf}.
//' @param b A numeric vector providing the upper bounds. These may be
//' \code{Inf}.
//' @param mean A numeric vector of means.
//' @param sd A numeric vector of standard deviations.
//'
//' @return A numeric vector.
//' @export
//'
//' @references
//' Geweke, J. (1991). Efficient simulation from the multivariate normal and
//' student-t distributions subject to linear constraints. In Computing Science
//' and Statistics: Proceedings of the 23rd Symposium on the Interface, ed. E.
//' Keramidas and S. Kaufman, pp. 571-8. Fairfax Station, VA: Interface
//' Foundation of North America.
//' 
//' @details
//' The values of \code{a}, \code{b}, \code{mean} and \code{sd} are recycled as
//' needed.
//'
//' This function is a Rcpp port of the corresponding 
//' \code{\link[truncnorm]{rtruncnorm}} function from the \code{truncnorm}
//' package.
//' 
//' @examples
//' n <- 20
//' a <- 0
//' b <- Inf
//' mean <- 0
//' sd <- 1
//' rng <- round(runif(1, 1, 10000))
//' set.seed(rng)
//' rtruncnorm_rcpp(n, a, b, mean, sd)
// [[Rcpp::export]]
Rcpp::NumericVector rtruncnorm_rcpp(const int& n, const double& a, const double& b, const double& mean, const double& sd) {
  Rcpp::NumericVector ret(n);
  
  for (unsigned int i = 0; i < n; i++) {
    if (R_FINITE(a) && R_FINITE(b)) {
      ret(i) = r_truncnorm(a, b, mean, sd);
    } else if (R_NegInf == a && R_FINITE(b)) {
      ret(i) = r_righttruncnorm(b, mean, sd);
    } else if (R_FINITE(a) && R_PosInf == b) {
      ret(i) = r_lefttruncnorm(a, mean, sd);
    } else if (R_NegInf == a && R_PosInf == b) {
      ret(i) = R::rnorm(mean, sd);
    } else {
      ret(i) = NA_REAL;
    }
    R_CheckUserInterrupt();
  }

  return ret;
}

//' The truncated univariate normal distribution.
//'
//' Density of a truncated normal distribution with mean equal to \code{mean}
//' and standard deviation equal to \code{sd}.
//'
//' @param x A numeric vector containing the values at which to evaluate the
//' density.
//' @param a A numeric vector providing the lower bounds. These may be
//' \code{-Inf}.
//' @param b A numeric vector providing the upper bounds. These may be
//' \code{Inf}.
//' @param mean A numeric vector of means.
//' @param sd A numeric vector of standard deviations.
//'
//' @return A numeric vector.
//' @export
//'
//' @details
//' The values of \code{a}, \code{b}, \code{mean} and \code{sd} are recycled as
//' needed.
//'
//' This function is a Rcpp port of the corresponding 
//' \code{\link[truncnorm]{dtruncnorm}} function from the \code{truncnorm}
//' package.
//' 
//' @examples
//' a <- 0
//' b <- Inf
//' mean <- 0
//' sd <- 1
//' rng <- round(runif(1, 1, 10000))
//' set.seed(rng)
//' n <- 20
//' x <- rtruncnorm_rcpp(n, a, b, mean, sd)
//' dtruncnorm_rcpp(x, a, b, mean, sd)
// [[Rcpp::export]]
Rcpp::NumericVector dtruncnorm_rcpp(const Rcpp::NumericVector& x, const double& a, const double& b, const double& mean, const double& sd) {
  int n = x.size();
  Rcpp::NumericVector ret(n);

  for (unsigned int i = 0; i < n; ++i) {
    const double cx = x(i);
    if (a <= cx && cx <= b) { /* In range: */
      const double c1 = R::pnorm(a, mean, sd, true, false);
      const double c2 = R::pnorm(b, mean, sd, true, false);
      const double c3 = sd * (c2 - c1);
      const double c4 = R::dnorm((cx - mean)/sd, 0.0, 1.0, false);
      ret(i) = c4/c3;
    } else { /* Truncated: */
      ret(i) = 0.0;
    }
    R_CheckUserInterrupt();
  }

  return ret;
}

//' Log Cholesky prior distribution.
//'
//' Evaluation of the prior distribution based on the log Cholesky
//' factorization for a covariance matrix \eqn{\Sigma_M} .
//'
//' @param A A numeric matrix at which to evaluate the prior density.
//' @param sigma_r Length-one numeric vector providing the standard deviation
//' of the underlying normal distributions used to generate the Cholesky
//' factors.
//'
//' @return A length-one numeric vector.
//' @export
//'
//' @examples
//' M <- 5
//' sigma_r <- 10^-(1/3)
//' rng <- round(runif(1, 1, 10000))
//' set.seed(rng)
//' Sigma <- rlogchol_arma(M, sigma_r)
//' dlogchol_arma(Sigma, sigma_r)
// [[Rcpp::export]]
double dlogchol_arma(const arma::mat& A, const double& sigma_r, const bool& logd = false) {
  int M = A.n_cols, nn = M*(M + 1)/2;
  
  arma::vec beta = Sigma_M_to_beta(A);

  double out = 0;
  for (unsigned int i = 0; i < nn; i++) {
    out += R::dnorm(beta(i), 0.0, sigma_r, true);
  }

  if (!logd) {
    out = exp(out);
  }

  return out;
}

//' Log Cholesky prior distribution.
//'
//' Generation of a single draw from the prior distribution of the
//' \eqn{\Sigma_M} matrix based on the log Cholesky factorization.
//'
//' @param M Length-one numeric vector representing the covariance matrix size.
//' @param sigma_r Length-one numeric vector providing the standard deviation
//' of the underlying normal distributions used to generate the Cholesky
//' factors.
//'
//' @return A numeric matrix.
//' @export
//'
//' @examples
//' M <- 5
//' sigma_r <- 10^-(1/3)
//' rng <- round(runif(1, 1, 10000))
//' set.seed(rng)
//' Sigma <- rlogchol_arma(M, sigma_r)
//' all(eigen(Sigma)$values > 0)
//' Sigma_sqrt <- solve(diag(diag(Sigma)))^0.5
//' Sigma_sqrt %*% Sigma %*% Sigma_sqrt
// [[Rcpp::export]]
arma::mat rlogchol_arma(const int& M, const double& sigma_r) {
  int nn = M*(M + 1)/2;
  arma::vec beta = Rcpp::as<arma::vec>(Rcpp::rnorm(nn, 0.0, sigma_r));

  return beta_to_Sigma_M(beta, M);
}

//' Multivariate t density.
//'
//' Evaluate the multivariate t density.
//'
//' @param x A numeric matrix whose rows represent the points at which to
//' evaluate the density.
//' @param mean A numeric vector containing the univariate means.
//' @param sigma A positive definite numeric matrix representing the covariance
//' matrix of the distribution.
//' @param df Integer length-one vector providing the number of degrees of
//' freedom.
//' @param logd Boolean length-one vector; if TRUE the log density is returned.
//'
//' @return A numeric vector.
//' @export
//'
//' @references
//' \url{http://gallery.rcpp.org/articles/dmvnorm_arma/}
//' 
//' @examples
//' set.seed(123)
//' ### Covariance matrix and mean vector
//' sigma <- bayesm::rwishart(10, diag(8))$IW
//' mean <- rnorm(8)
//' ### Benchmarking
//' n <- 1e+4
//' df <- 20
//' X <- mvtnorm::rmvt(n, sigma, df, mean)
//' require(rbenchmark)
//' benchmark(mvtnorm::dmvt(X, mean, sigma, df, log = FALSE), 
//'           LearnBayes::dmt(X, mean, sigma, df, FALSE),
//'           dmvt_arma(X, mean, sigma, df, FALSE),
//'           order = "relative", replications = 100)[, 1:4]
// [[Rcpp::export]]
Rcpp::NumericVector dmvt_arma(const arma::mat& x, const arma::vec& mean, const arma::mat& sigma, const int& df, const bool& logd = false) {
  int n = x.n_rows, d = sigma.n_cols;

  arma::mat X = arma::trans(x);
  X.each_col() -= mean;
  arma::vec Q = arma::vectorise(sum((arma::inv(sigma) * X) % X, 0), 0);

  double logdet, logdet_sgn;
  arma::log_det(logdet, logdet_sgn, sigma);
  double constant = R::lgammafn((df + d)/2.0) - (d*log(M_PI*df) + logdet)/2.0 - R::lgammafn(df/2.0);

  Rcpp::NumericVector out(n);
  for (unsigned int i = 0; i < n; i++) {
    out(i) = constant - 0.5*(df + d)*log(1 + Q(i)/static_cast<double>(df));
  }

  if (!logd) {
    out = exp(out);
  }

  return out;
}

//' Multivariate t variates generation.
//'
//' Generation of multivariate t variates.
//'
//' @param n A length-one numeric vector providing the number of draws to
//' generate.
//' @param mean A numeric vector containing the univariate means.
//' @param sigma A positive definite numeric matrix representing the covariance
//' matrix of the distribution.
//' @param df Integer length-one vector providing the number of degrees of
//' freedom.
//'
//' @return A numeric matrix.
//' @export
//'
//' @references
//' \url{http://gallery.rcpp.org/articles/simulate-multivariate-normal/}
//' 
//' @examples
//' set.seed(123)
//' ### Covariance matrix and mean vector
//' sigma <- matrix(c(1, 0.9, -0.3, 0.9, 1, -0.4, -0.3, -0.4, 1), ncol = 3)
//' mu <- c(10, 5, -3)
//' df <- 7
//' ### Benchmarking
//' n <- 1e+4
//' require(rbenchmark)
//' benchmark(mvtnorm::rmvt(n, sigma, df, mu),
//'           LearnBayes::rmt(n, mu, sigma, df),
//'           rmvt_arma(n, mu, sigma, df),
//'           columns = c("test", "replications", "relative", "elapsed"),
//'           order = "relative")
// [[Rcpp::export]]
arma::mat rmvt_arma(const int& n, const arma::vec& mean, const arma::mat& sigma, const int& df) {
  int d = sigma.n_cols;

  arma::vec x = Rcpp::as<arma::vec>(Rcpp::rchisq(n, df));
  x = arma::sqrt(x/df);
  arma::vec zero_mean = arma::zeros(d);
  arma::mat t = rmvn_arma(n, zero_mean, sigma);
  t.each_col() /= x;
  t.each_row() += mean.t();

  return t;
}

//' Inverse gamma density.
//'
//' Evaluate the inverse gamma density.
//'
//' @param x A numeric vector whose representing the points at which to evaluate
//' the density.
//' @param alpha Inverse gamma shape parameter. Must be strictly positive.
//' @param beta Inverse gamma scale parameter. Must be strictly positive.
//' @param logd Boolean length-one vector; if TRUE the log density is returned.
//'
//' @return A numeric vector.
//' @export
//'
//' @examples
//' x <- seq(.01, 50, length.out = 1000)
//' alpha <- 1
//' beta <- 5
//' res <- dinvgamma_rcpp(x, alpha, beta)
//' plot(x, res, type = "l", main = "Inverse Gamma density")
// [[Rcpp::export]]
arma::vec dinvgamma_rcpp(const arma::vec& x, const double& alpha, const double& beta, const bool& logd = false) {
  if ((alpha <= 0) || (beta <= 0)) {
    error("alpha (shape) and/or beta (scale) parameters in dinvgamma_rcpp() need to be both strictly positive.\n");
  }

  int n = x.n_elem;
  arma::vec out(n);
  double lbeta = log(beta);
  double lgalpha = R::lgammafn(alpha);

  for (unsigned int i = 0; i < n; i++) {
    out(i) = alpha * lbeta - lgalpha - (alpha + 1) * log(x(i)) - (beta/x(i));
    if (!logd) {
      out(i) = exp(out(i));
    }
  }

  return out;
}

//' Inverse gamma variates generation.
//'
//' Generation of inverse gamma variates.
//'
//' @param n A length-one numeric vector providing the number of draws to
//' generate.
//' @param alpha Inverse gamma shape parameter. Must be strictly positive.
//' @param beta Inverse gamma scale parameter. Must be strictly positive.
//'
//' @return A numeric vector.
//' @export
//'
//' @examples
//' set.seed(123)
//' n <- 1e4
//' alpha <- 3
//' beta <- 5
//' x <- sort(rinvgamma_rcpp(n, alpha, beta))
//' hist(x, breaks = 30, xlab = "x", freq = FALSE,
//'      main = "Inverse gamma variates")
//' lines(x, dinvgamma_rcpp(x, alpha, beta))
// [[Rcpp::export]]
arma::vec rinvgamma_rcpp(const int& n, const double& alpha, const double& beta) {
  if ((alpha <= 0) || (beta <= 0)) {
    error("alpha (shape) and/or beta (scale) parameters in rinvgamma() need to be both strictly positive.\n");
  }

  arma::vec out(n);

  for (unsigned int i = 0; i < n; i++) {
    out(i) = 1.0/R::rgamma(alpha, 1.0/beta);
  }

  return out;
}

//' Gaussian copula density evaluation.
//'
//' Evaluates a Gaussian copula density.
//'
//' @param u A numeric vector providing the point at which to evaluate the
//' copula density.
//' @param Gamma A numeric matrix providing the Gaussian copula correlation
//' matrix. These elements need to be provided row-wise (i.e., g_12, g_13,
//' \ldots, g_1M, g_23, g_24,\ldots, g_2M,\ldots,g_(M-1)(M-1)).
//' @param is_u A length-one logical vector indicating whether the \code{u}
//' argument provides the \eqn{u} or the \eqn{\Phi^{-1}(u)} values.
//' @param logd Boolean length-one vector; if TRUE the log density is returned.
//'
//' @return A length-one numeric vector.
//' @export
//'
//' @examples
//' u <- rep(0.1, 3)
//' Gamma <- c(.3, .5, .2)
//' gausscopdens(u, Gamma)
//' 
//' Gamma <- .3
//' len <- 100
//' u1 <- u2 <- seq(.01, .99, length.out = len)
//' dens <- matrix(NA, nrow = len, ncol = len)
//' for (i in 1:len) {
//'   for (j in 1:len) {
//'     dens[i, j] <- gausscopdens(c(u1[i], u2[j]), Gamma)
//'   }
//' }
//' persp(u1, u2, dens, theta = 120, phi = 25)
// [[Rcpp::export]]
double gausscopdens(const Rcpp::NumericVector& u, const arma::mat& Gamma, const bool& is_u, const bool& logd = false) {
  int M = u.size();

  arma::vec x(M);
  if (is_u) {
    x = Rcpp::qnorm(u, 0.0, 1.0);
  }
  else {
    x = Rcpp::as<arma::vec>(u);
  }
  double logdet, logdet_sgn;
  arma::log_det(logdet, logdet_sgn, Gamma);
  arma::mat H = arma::eye(M, M) - arma::inv(Gamma);
  double udens = arma::as_scalar(x.t() * H * x);

  double out = -0.5 * (logdet - udens);
  if (!logd) {
    out = exp(out);
  }

  return out;
}
