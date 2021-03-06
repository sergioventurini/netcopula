% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rlogchol_arma}
\alias{rlogchol_arma}
\title{Log Cholesky prior distribution.}
\usage{
rlogchol_arma(M, sigma_r)
}
\arguments{
\item{M}{Length-one numeric vector representing the covariance matrix size.}

\item{sigma_r}{Length-one numeric vector providing the standard deviation
of the underlying normal distributions used to generate the Cholesky
factors.}
}
\value{
A numeric matrix.
}
\description{
Generation of a single draw from the prior distribution of the
\eqn{\Sigma_M} matrix based on the log Cholesky factorization.
}
\examples{
M <- 5
sigma_r <- 10^-(1/3)
rng <- round(runif(1, 1, 10000))
set.seed(rng)
Sigma <- rlogchol_arma(M, sigma_r)
all(eigen(Sigma)$values > 0)
Sigma_sqrt <- solve(diag(diag(Sigma)))^0.5
Sigma_sqrt \%*\% Sigma \%*\% Sigma_sqrt
}
