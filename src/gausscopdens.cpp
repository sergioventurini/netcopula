#include "netcopula.h"

// [[Rcpp::depends(RcppArmadillo)]]

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
