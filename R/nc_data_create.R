#' Create an instance of the \code{nc_data} class.
#'
#' @param type A length-one character vector representing the data type and
#' must be either \code{"binary"}, \code{"count"} or \code{"continuous"}.
#' @param n_study A length-one numeric vector representing the number of
#' studies.
#' @param n_outcomes A length-one numeric vector representing the number of
#' studies.
#' @param n_datapoints A length-one numeric vector representing the number of
#' data points in the sample (i.e., the overall number of arms for the trials
#' in the sample).
#' @param n_treatments A length-one numeric vector representing the total
#' number of treatments involved in the analysis.
#' @param study_id A four-columns data frame providing studies-specific data.
#' The four columns must be named \code{studyid}, \code{study}, \code{narms}
#' and \code{baseline}, with the second one that needs to be a progressive
#' number that match with those in the \code{study_data} slot. If NULL, a
#' default set of IDs will be created.
#' @param study_data A data frame containing the outcome data for each study.
#' Its first five columns must be named \code{study}, \code{arm}, \code{trt},
#' \code{baseline} and \code{narms}. Values in the \code{study} column must
#' match with those in the same column of the \code{study_id} slot.
#' @param ref_trt A length-one integer vector providing the reference
#' treatment (which should not be confused with the baseline treatment for each
#' study).
#'
#' @return Returns an instance of the \code{\link{nc_data-class}}.
#'
#' @name nc_data_create
#' @rdname nc_data_create
#'
#' @export
#'
#' @examples
#' data(homesafety)
#'
#' n_study <- max(homesafety$study)
#' n_outcomes <- 3
#' n_datapoints <- nrow(homesafety)
#' n_treatments <- max(homesafety$trt)
#' studyid <- data.frame(studyid = unique(homesafety$studyid),
#'                       study = 1:n_study)
#'
#' nc_data_create("binary", n_study, n_outcomes, n_datapoints, n_treatments,
#'                studyid, homesafety[, 2:12])

nc_data_create <- function(type = character(1), n_study = numeric(1), n_outcomes = numeric(1), n_datapoints = numeric(1), n_treatments = numeric(1), study_id = data.frame(), study_data = data.frame(), ref_trt = integer(1)) {
  out <- new("nc_data",
    type = type,
    n_study = n_study,
    n_outcomes = n_outcomes,
    n_datapoints = n_datapoints,
    n_treatments = n_treatments,
    study_id = study_id,
    study_data = study_data,
    ref_trt = ref_trt)
  return(out)
}

#' Estimation of copula based models for multivariate network meta-analysis.
#'
#' Estimation of copula based models for multivariate network meta-analysis.
#'
#' @param data Object of class \code{nc_data} containing the data to analyze.
#' @param prm_init Named list of parameters for the proposal distributions of
#' the \eqn{Z} latent positions and \eqn{\alpha} parameters; the elements'
#' names must be \code{"z"} and \code{"alpha"} and both need to be a single
#' number representing the standard deviations of the normal proposal
#' distributions for each \eqn{z_{ih}^{g}}{z_ih^g} and \eqn{\alpha_g}
#' respectively.
#' @param random_start Length-one logical vector; if \code{TRUE} (default), it
#' allows the algorithm to start from a random starting point.
#' @param init Named list providing user-defined starting values. If provided,
#' it overrides the value given for \code{random_start}.
#'
#' @return \code{netcopula} returns an object of class \code{\link{nc_mcmc}}.
#' @export
#'
#' @examples
#' \dontrun{
#' # netcopula(...)
#' }
#'
#' @references
#' Venturini, S. and Graziani, R. (2016), "A Bayesian Copula Based Model for
#' Multivariate Network Meta-Analysis". Technical report.
#'
#' @author Sergio Venturini \email{sergio.venturini@unibocconi.it}
#'
#' @seealso
#' \code{\link{nc_data-class}},
#' \code{\link{nc_mcmc}}
nc_data_simulate <- function(d, sigma_vec, rho_vec, gamma_vec, n_studies, n_trt, n_bin, n_outcomes, mu_sigma2, ref_trt) {
  # code for simulating data according to the 'netcopula' model for the case of
  # two arms studies the first half comparing 1 to 2 and second half 2 to 3,
  # treat 1 taken as baseline

  n_datapoints <- 2*n_studies
  study <- sort(rep(1:n_studies, 2))
  arm <- rep(c(1, 2), n_studies)
  trt <- c(rep(c(1, 2), n_studies/2), rep(c(2, 3), n_studies/2))
  baseline <- c(rep(1, n_studies), rep(2, n_studies))
  narms <- rep(2, n_datapoints)
  n <- kronecker(matrix(1, n_datapoints, 1), t(n_bin))
  colnames(n) <- paste0("n", 1:n_outcomes)
  data <- data.frame(study, arm, trt, baseline, narms, n)
  study_id <- data.frame(baseline = c(rep(1, n_studies/2), rep(2, n_studies/2)), narms = rep(2, n_studies))

  Corr_M <- diag(n_outcomes)
  Corr_M[lower.tri(Corr_M)] <- Corr_M[upper.tri(Corr_M)] <- rho_vec
  Sigma_M <- diag(sigma_vec)
  Sigma_M <- Sigma_M %*% Corr_M %*% Sigma_M

  Gamma <- diag(n_outcomes)
  Gamma[lower.tri(Gamma)] <- Gamma[upper.tri(Gamma)] <- gamma_vec

  d_trt <- matrix(NA, nrow = n_trt, ncol = n_outcomes)
  d_trt[ref_trt, ] <- 0
  trt_count <- 1
  for (q in seq_along(1:n_trt)[-ref_trt]) {
    d_trt[q, ] <- d[trt_count, ]
    trt_count <- trt_count + 1
  }

  mu_mean <- rep(0, n_outcomes)
  mu_Sigma <- mu_sigma2*diag(n_outcomes)
  mu_tmp <- rmvn_arma(n_studies, mu_mean, mu_Sigma)
  mu <- matrix(0, nrow = n_datapoints, ncol = n_outcomes)
  mu[seq(1, n_datapoints, by = 2), ] <- mu_tmp # this is not general

  delta <- matrix(NA, nrow = n_datapoints, ncol = n_outcomes)
  trt_cum <- 1
  for (i in 1:n_studies) {
    trt_study <- data[data[, "study"] == i, "trt"]
    base_study <- study_id[i, "baseline"]
    trt_no_base_study <- setdiff(trt_study, base_study)
    d_study <- d_trt[trt_no_base_study, ] - d_trt[base_study, ]
    Sigma <- Sigma_block(Sigma_M, study_id[i, "narms"] - 1)
    delta_tmp <- rmvn_arma(1, d_study, Sigma)
    delta_tmp <- matrix(delta_tmp, nrow = length(trt_no_base_study), ncol = n_outcomes, byrow = TRUE)
    delta[trt_cum:(trt_cum + length(trt_study) - 1), ] <- rbind(rep(0, n_outcomes), delta_tmp) # this is not general
    trt_cum <- trt_cum + length(trt_study)
  }

  x <- rmvn_arma(n = n_datapoints, rep(0, n_outcomes), Gamma)
  y <- matrix(NA, nrow = n_datapoints, ncol = n_outcomes)
  theta <- mu + delta
  x_inv <- pnorm(x)
  for (i in 1:n_datapoints) {
    for (j in 1:n_outcomes) {
      y[i, j] <- qbinom(x_inv[i, j], n[i, j], expit_rcpp(theta[i, j]))
    }
  }
  colnames(y) <- paste0("y", 1:n_outcomes)

  data <- cbind(data, y)
  nc_data <- nc_data_create("binary", n_studies, n_outcomes, n_datapoints, n_trt, NULL, data, ref_trt)

  return(list(nc_data = nc_data, mu = mu_tmp, delta = delta, d = d, x = x,
    Sigma_M = Sigma_M, Gamma = Gamma))
}
