// boost.cpp

#include "netcopula.h"
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/normal.hpp>

// [[Rcpp::depends("BH")]]

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector qnorm_boost(const Rcpp::NumericVector& p, const double& mean = 0.0, const double& sd = 1.0, const bool& lower_tail = true) {
  int N = p.size();
  Rcpp::NumericVector out(N);
  boost::math::normal norm(mean, sd);

  for (unsigned int i = 0; i < N; i++) {
    if (lower_tail) {
      out(i) = boost::math::quantile(norm, p(i));
    } else {
      out(i) = boost::math::quantile(boost::math::complement(norm, p(i)));
    }
  }

  return out;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pnorm_boost(const Rcpp::NumericVector& x, const double& mean = 0.0, const double& sd = 1.0, const bool& lower_tail = true) {
  int N = x.size();
  Rcpp::NumericVector out(N);
  boost::math::normal norm(mean, sd);

  for (unsigned int i = 0; i < N; i++) {
    if (lower_tail) {
      out(i) = boost::math::cdf(norm, x(i));
    } else {
      out(i) = boost::math::cdf(boost::math::complement(norm, x(i)));
    }
  }

  return out;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector qbinom_boost(const Rcpp::NumericVector& p, const int& n, const double& prob, const bool& lower_tail = true) {
  int N = p.size();
  Rcpp::NumericVector out(N);
  boost::math::binomial_distribution<double> binom(n, prob);

  for (unsigned int i = 0; i < N; i++) {
    if (lower_tail) {
      out(i) = boost::math::quantile(binom, p(i));
    } else {
      out(i) = boost::math::quantile(boost::math::complement(binom, p(i)));
    }
  }

  return out;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pbinom_boost(const Rcpp::NumericVector& x, const int& n, const double& prob, const bool& lower_tail = true) {
  int N = x.size();
  Rcpp::NumericVector out(N);
  boost::math::binomial_distribution<double> binom(n, prob);

  for (unsigned int i = 0; i < N; i++) {
    if (lower_tail) {
      out(i) = boost::math::cdf(binom, x(i));
    } else {
      out(i) = boost::math::cdf(boost::math::complement(binom, x(i)));
    }
  }

  return out;
}
