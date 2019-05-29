#' Data simulation for the 'netcopula' package.
#'
#' Function for simulating data to use with the 'netcopula' package. The function
#' simulates data for the case of two arms studies, with the first half comparing
#' treatment 1 (reference) and 2 and the second half treatment 2 and 3.
#'
#' @param d Numeric matrix with dimension (#treatement - 1, #outcomes) containing
#' the true treatment effects used to simulate the data
#' @param sigma2_vec Numeric vector containing the true common between-study
#' variances.
#' @param rho_vec Numeric vector containing the true common between-study
#' correlations.
#' @param gamma_vec Numeric vector containing the true common Gaussian copula
#' correlation parameters.
#' @param n_studies A length-one numeric vector representing the number of
#' studies.
#' @param n_trt A length-one numeric vector representing the total number
#' of treatments involved in the analysis.
#' @param n_bin A length-one numeric vector representing the size of the
#' binomial distribution to sample the event numbers from.
#' @param n_outcomes A length-one numeric vector representing the number of
#' outcomes.
#' @param mu_sigma2 A length-one numeric vector representing the variance of
#' the mu parameters (study-specific baseline effects).
#' @param ref_trt A length-one integer vector providing the reference
#' treatment (which should not be confused with the baseline treatment for each
#' study).
#'
#' @return Returns a list with an instance of the \code{\link{nc_data-class}}
#' and the objects containing the parameters used to generate the data.
#'
#' @name nc_data_simulate
#' @rdname nc_data_simulate
#'
#' @export
#'
#' @examples
#' set.seed(101)
#' 
#' n_outcomes <- 3
#' n_trt <- 3
#' n_studies <- 30
#' n_bin <- rep(100, n_outcomes) # binomial distributions sizes
#' d_init <- array(data = c(0, -0.5, -1, 0, 0.8, 0.3), dim = c(n_trt - 1, n_outcomes))
#' sigma2_vec <- c(1, 1.6, 1.8) # standard deviations in Sigma_M
#' rho_vec <- rep(.5, n_outcomes)
#' gamma_vec <- c(.3, .1, .5) # components of copula association matrix
#' mu_sigma2_sim <- 1e-1
#' ref_trt <- 1
#' 
#' nc_data_sim <- nc_data_simulate(d_init, sigma2_vec, rho_vec, gamma_vec, n_studies,
#'   n_trt, n_bin, n_outcomes, mu_sigma2_sim, ref_trt)
#' summary(nc_data_sim$nc_data)
nc_data_simulate <- function(d, sigma2_vec, rho_vec, gamma_vec, n_studies, n_trt, n_bin, n_outcomes, mu_sigma2, ref_trt) {
  if (nrow(d) != (n_trt - 1)) {
    stop("the d matrix must have (n_trt - 1) rows")
  }
  if (ncol(d) != n_outcomes) {
    stop("the d matrix has must have n_outcomes columns")
  }
  if (length(sigma2_vec) != n_outcomes) {
    stop("incompatible sigma2_vec dimension")
  }
  if (length(rho_vec) > 1) {
    if (length(rho_vec) != (n_outcomes*(n_outcomes - 1)/2)) {
      stop("incompatible rho_vec dimension")
    }
    if (length(gamma_vec) != (n_outcomes*(n_outcomes - 1)/2)) {
      stop("incompatible gamma_vec dimension")
    }
  }

  n_datapoints <- 2*n_studies # i.e. 2 arms per study
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
  if(length(sigma2_vec) > 1) {
    Sigma_M <- diag(sigma2_vec)
  } else {
    Sigma_M <- sigma2_vec
  }
  Sigma_M <- Sigma_M %*% Corr_M %*% Sigma_M

  Gamma <- diag(n_outcomes)
  Gamma[lower.tri(Gamma)] <- Gamma[upper.tri(Gamma)] <- gamma_vec

  d_trt <- matrix(NA, nrow = n_trt, ncol = n_outcomes)
  trt_count <- 1
  for (q in seq_along(1:n_trt)[-ref_trt]) {
    d_trt[q, ] <- d[trt_count, ]
    trt_count <- trt_count + 1
  }
  d_trt[ref_trt, ] <- 0

  mu_mean <- rep(0, n_outcomes)
  mu_Sigma <- mu_sigma2*diag(n_outcomes)
  mu_tmp <- rmvn_arma(n_studies, mu_mean, mu_Sigma)
  mu <- matrix(0, nrow = n_datapoints, ncol = n_outcomes)
  mu[seq(1, n_datapoints, by = 2), ] <- mu_tmp     # this is not general
  mu[seq(1, n_datapoints, by = 2) + 1, ] <- mu_tmp # this is not general

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
  nc_data <- new("nc_data",
    type = "binary",
    n_study = n_studies,
    n_outcomes = n_outcomes,
    n_datapoints = n_datapoints,
    n_treatments = n_trt,
    study_id =  NULL,
    study_data =  data,
    ref_trt = ref_trt)

  return(list(nc_data = nc_data, mu = mu_tmp, delta = delta, d = d, x = x,
    Sigma_M = Sigma_M, Gamma = Gamma))
}

#' Missing data simulation for the 'netcopula' package.
#'
#' Function for adding missing values to a simulated data set. The function
#' currently generate missing values only according to the MCAR mechanism.
#'
#' @param nc_data_list List as generated by the \code{\link{nc_data_simulate}}
#' function.
#' @param n_miss A length-one numeric vector representing the number of missing
#' values.
#'
#' @return Returns a list with an instance of the \code{\link{nc_data-class}}
#' and the objects containing the parameters used to generate the data.
#'
#' @name nc_data_missing
#' @rdname nc_data_missing
#'
#' @export
#'
#' @examples
#' set.seed(101)
#' 
#' n_outcomes <- 3
#' n_trt <- 3
#' n_studies <- 30
#' n_bin <- rep(100, n_outcomes) # binomial distributions sizes
#' d_init <- array(data = c(0, -0.5, -1, 0, 0.8, 0.3), dim = c(n_trt - 1, n_outcomes))
#' sigma2_vec <- c(1, 1.6, 1.8) # standard deviations in Sigma_M
#' rho_vec <- rep(.5, n_outcomes)
#' gamma_vec <- c(.3, .1, .5) # components of copula association matrix
#' mu_sigma2_sim <- 1e-1
#' ref_trt <- 1
#' 
#' nc_data_sim <- nc_data_simulate(d_init, sigma2_vec, rho_vec, gamma_vec, n_studies,
#'   n_trt, n_bin, n_outcomes, mu_sigma2_sim, ref_trt)
#' nc_data_sim <- nc_data_missing(nc_data_sim, n_miss = 10)
#' summary(nc_data_sim$nc_data)
nc_data_missing <- function(nc_data_list, n_miss = 5) {
  if (!is.list(nc_data_list)) {
    stop("the provided data must be a list returned by the 'nc_data_simulate' function")
  }
  if (!all(names(nc_data_list) == c("nc_data", "mu", "delta", "d", "x", "Sigma_M", "Gamma"))) {
    stop("the provided data must be a list returned by the 'nc_data_simulate' function")
  }

  nc_data <- nc_data_list$nc_data
  n_study <- nc_data@n_study
  n_outcomes <- nc_data@n_outcomes
  n_treatments <- nc_data@n_treatments
  n_datapoints <- nc_data@n_datapoints
  study_data <- nc_data@study_data

  if (n_miss > n_datapoints/4) {
    stop("too many missing values selected; try reducing 'n_miss'")
  }

  # only MCAR mechanism is implemented
  # n_miss studies are randomly selected
  idx_study <- sample(x = 1:n_study, size = n_miss, replace = FALSE)
  # for each of these studies we randomly select an outcome
  idx_outcome <- sample(x = 1:n_outcomes, size = n_miss, replace = TRUE)
  # we then set to NA the corresponding outcome in the selected studies
  for (i in 1:n_miss) {
    study_data[study_data$study == idx_study[i], paste0("n", idx_outcome[i])] <- NA
    study_data[study_data$study == idx_study[i], paste0("y", idx_outcome[i])] <- NA
  }
  # we repeat the operation but only for the studies not yet selected
  studies_resid <- setdiff(1:n_study, idx_study)
  idx_study <- sample(x = studies_resid, size = n_miss, replace = FALSE)
  idx_outcome <- sample(x = 1:n_outcomes, size = n_miss, replace = TRUE)
  for (i in 1:n_miss) {
    study_data[study_data$study == idx_study[i], paste0("n", idx_outcome[i])] <- NA
    study_data[study_data$study == idx_study[i], paste0("y", idx_outcome[i])] <- NA
  }

  nc_data_list$nc_data@study_data <- study_data
  y <- study_data[, paste0("y", 1:n_outcomes)]
  nc_data_list$x[is.na(y)] <- NA

  return(nc_data_list)
}
