#include <netcopula.h>

//' Transform an Armadillo field<vec> to a matrix
//'
//' Unlists vectors in a field and places them into a matrix
//' @param x A \code{field<vec>}.
//' @return A \code{mat} containing the field elements within a column.
//' @author https://github.com/coatless/r-to-armadillo
// [[Rcpp::export]]
arma::mat field_to_matrix(arma::field<arma::vec> x) {
  unsigned int nx = x.n_elem;
  unsigned int row = x(0).n_elem;
  arma::mat A(row, nx);

  for (unsigned int i = 0; i < nx; i++) {
    A.col(i) = x(i);
  }

  return A;
}
