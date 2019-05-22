#' \code{netcopula}: a package for estimating copula based models for
#' multivariate network meta-analysis (NMA).
#'
#' The \code{netcopula} package provides methods for fitting copula based
#' models for network meta-analysis within a fully Bayesian framework.
#' Currently, it includes methods for binary data only.
#'
#' @docType package
#'
#' @name netcopula
#'
#' @useDynLib netcopula
#' @import methods
#' @import utils
#' @importFrom coda as.mcmc
#' @importFrom coda mcmc
#' @importFrom Rcpp evalCpp
NULL

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
nc_init <- function(data, prm_init, random_start = TRUE, init = NULL, prm_prior = NULL, prm_prop = NULL, tuning = NULL, adaptation = NULL) {
  if (!random_start & is.null(init)) {
    stop("either choose 'random_start = TRUE' or provide some starting values.")
  }

  n <- data@n_study
  M <- data@n_outcomes
  n_datapoints <- data@n_datapoints
  n_trt <- data@n_treatments
  n_y <- data@study_data[, (7):(7 + M - 1)]
  y <- data@study_data[, (7 + M):(6 + 2*M)]
  trt <- data@study_data[, "trt"]
  # logpost <- -Inf

  # while (is.na(logpost) | is.infinite(logpost)) {
    if (random_start & is.null(init)) {
      mu_mean <- rep(0, M)
      mu_Sigma <- (prm_init$mu_sigma2)*diag(M)
      d_mean <- rep(0, M)
      d_Sigma <- (prm_init$d_sigma2)*diag(M)

      # initialize mu (study-specific baseline effects)
      mu <- rmvn_arma(n, mu_mean, mu_Sigma)

      # initialize d (pooled treatment effects across trials)
      d <- rmvn_arma(n_trt, d_mean, d_Sigma)
      d[data@ref_trt, ] <- 0 # reference treatment

      # initialize Sigma_M (common between-study covariance structure)
      beta_sigma <- prm_init$beta_sigma
      Sigma_M <- rlogchol_arma(M, beta_sigma)

      # initialize delta (study-specific [random] treatment effects)
      delta <- matrix(NA, nrow = n_datapoints, ncol = M)
      trt_cum <- 1
      for (i in 1:n) {
        trt_study <- data@study_data[data@study_data[, "study"] == i, "trt"]
        base_study <- data@study_id[i, "baseline"]
        trt_no_base_study <- setdiff(trt_study, base_study)
        d_study <- d[trt_no_base_study, ] - d[base_study, ]
        d_study <- as.numeric(t(d_study))
        Sigma <- Sigma_block(Sigma_M/100, data@study_id[i, "narms"] - 1)
        delta_tmp <- rmvn_arma(1, d_study, Sigma)
        delta_tmp <- matrix(delta_tmp, nrow = length(trt_no_base_study), ncol = M, byrow = TRUE)
        delta[trt_cum:(trt_cum + length(trt_study) - 1), ] <- rbind(rep(0, M), delta_tmp)
        trt_cum <- trt_cum + length(trt_study)
      }

      # initialize Gamma_q's (Gaussian copula correlation structures)
      Gamma <- list()
      for (q in 1:prm_init$nGamma) {
        # Gamma[[q]] <- diag(M)
        if (M > 1) {
          Gamma[[q]] <- rlkj_arma(M, 1)
        } else {
          Gamma[[q]] <- as.matrix(1)
        }
      }

      # initialize x (latent variables)
      x <- matrix(NA, nrow = n_datapoints, ncol = M)
      for (i in 1:n_datapoints) {
        if (prm_init$nGamma == n_trt) {
          Gamma_study <- Gamma[[data@study_data[i, "trt"]]]
        } else {
          Gamma_study <- Gamma[[1]]
        }
        x[i, ] <- rmvn_arma(1, rep(0, M), as.matrix(Gamma_study))
      }

      # initialize D's (latent variables scale matrices)
      D <- list()
      for (q in 1:prm_init$nGamma) {
        if (prm_init$nGamma == n_trt) {
          if (M > 1) {
            D[[q]] <- as.matrix(diag(sqrt(1/colSums(x[data@study_data[, "trt"] == q, ]^2))))
          } else {
            D[[q]] <- as.matrix(sqrt(1/colSums(x[data@study_data[, "trt"] == q, ]^2)))
          }
        } else if (prm_init$nGamma == 1) {
          if (M > 1) {
            D[[q]] <- as.matrix(diag(sqrt(1/colSums(x^2))))
          } else {
            D[[q]] <- as.matrix(sqrt(1/colSums(x^2)))
          }
        }
      }

      x[is.na(y)] <- NA # sets NA in the x matrix
    } else if (!is.list(init)) {
      if (init == "univariate") {
        adaptation_uni <- adaptation
        adaptation_uni$maxiter <- 160
        start_vals <- init_univariate(data = data, burnin = 8000, nsim = 2000, nthin = 10, prm_prior = prm_prior, prm_prop = prm_prop, prm_init = prm_init, tuning = tuning, adaptation = adaptation_uni, verbose = FALSE)

        mu <- start_vals$mu
        delta <- start_vals$delta
        d <- start_vals$d
        Sigma_M <- start_vals$Sigma_M
        Gamma <- start_vals$Gamma
        D <- start_vals$D
        x <- start_vals$x
      }
    } else {
      mu <- init$mu
      delta <- init$delta
      d <- init$d
      Sigma_M <- init$Sigma_M
      Gamma <- init$Gamma
      D <- init$D
      x <- init$x
    }

  #   loglik <- nc_loglik(as.matrix(y), as.matrix(n_y), x, trt, mu, delta, Gamma)
  #   logprior <- nc_logprior(mu, sqrt(prm_init$mu_sigma2), d, sqrt(prm_init$d_sigma2), Sigma_M, prm_init$beta_sigma, ref_trt)
  #   logpost <- loglik + logprior
  # }

  return(list(mu = mu, delta = delta, d = d, Sigma_M = Sigma_M, Gamma = Gamma, D = D, x = x))
}

#' Estimation of copula based models for multivariate network meta-analysis.
#'
#' Estimation of copula based models for multivariate network meta-analysis.
#'
#' @param data Object of class \code{nc_data} containing the data to analyze.
#' @param burnin Length one integer vector providing the number of MCMC burnin
#' iterations.
#' @param nsim Length one integer vector providing the number of MCMC
#' iterations.
#' @param nthin Length one integer vector providing the thinning interval.
#' @param prm_prior Named list of parameters for the proposal distributions of
#' the \eqn{Z} latent positions and \eqn{\alpha} parameters; the elements'
#' names must be \code{"z"} and \code{"alpha"} and both need to be a single
#' number representing the standard deviations of the normal proposal
#' distributions for each \eqn{z_{ih}^{g}}{z_ih^g} and \eqn{\alpha_g}
#' respectively.
#' @param prm_prop Named list of hyperparameters for the prior distributions
#' of the \eqn{\sigma^2} and \eqn{\lambda} parameters; the elements' names
#' must be \code{"sigma2"} and \code{"lambda"} with the first one a two-sized
#' vector containing the \eqn{\sigma^2} prior hyperparameters and the second
#' element a vector providing the \eqn{\lambda} hyperparameters.
#' @param random_start Length-one logical vector; if \code{TRUE} (default), it
#' allows the algorithm to start from a random starting point.
#' @param init Named list providing user-defined starting values. If provided,
#' it overrides the value given for \code{random_start}.
#' @param tuning Named list of tuning parameters (see details).
#' @param adaptation Named list of adaptation parameters (see details).
#' @param verbose Length-one logical vector; if \code{TRUE}, messages
#' informing about the progress of the computations will be printed.
#'
#' @return \code{netcopula_old} returns an object of class \code{\link{nc_mcmc}}.
#' @export
#'
#' @examples
#' \dontrun{
#' # netcopula_old(...)
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
netcopula_old <- function(data, burnin = 10000, nsim = 5000, nthin = 1, prm_prior, prm_prop = NULL, prm_init, random_start = TRUE, init = NULL, tuning, adaptation, verbose = FALSE, mu_delta_mh = TRUE) {
  if (class(data) != "nc_data") {
    stop("the data argument must be a 'nc_data' object; see the help of the 'nc_data_create()' function.")
  }
  if (adaptation$maxiter > (burnin + nsim)/adaptation$every) {
    stop("the 'maxiter' adaptation parameter is too large.")
  }
  if (adaptation$miniter > adaptation$maxiter) {
    stop("the 'miniter' adaptation parameter must be smaller than or equal to 'maxiter'.")
  }
  if (adaptation$alpha < 0 | adaptation$alpha > 1 |
    adaptation$beta < 0 | adaptation$beta > 1 |
    adaptation$gamma < 0 | adaptation$gamma > 1 |
    adaptation$tar < 0 | adaptation$tar > 1) {
    stop("the 'alpha', 'beta', 'gamma' and 'tar' adaptation parameters must be in [0, 1].")
  }
  if (adaptation$tol <= 0) {
    stop("the 'tol' adaptation parameter must be positive.")
  }

  n <- slot(data, "n_study")
  M <- slot(data, "n_outcomes")
  n_datapoints <- slot(data, "n_datapoints")
  tot_iter <- burnin + nsim

  if (verbose) cat("Initialization of the algorithm...")

  if (is.null(prm_init)) {
    prm_init <- prm_prior
  }
  nc_start <- nc_init(data = data, prm_init = prm_init, random_start = random_start, init = init)

  if (verbose) cat("done!\n")

  # start iteration
  if (verbose) cat("Running the MCMC simulation...\n")

  if (mu_delta_mh) {
    res <- nc_mcmc_mh(
      data = data,
      init = nc_start,
      totiter = tot_iter,
      prior = prm_prior,
      prop = prm_prop,
      tuning = tuning,
      adapt = adaptation,
      verbose = verbose
    )
  } else {
    # TO DO: alla fine eliminare questa possibilita'
    res <- nc_mcmc_opt(
      data = data,
      init = nc_start,
      totiter = tot_iter,
      prior = prm_prior,
      prop = prm_prop,
      tuning = tuning,
      adapt = adaptation,
      verbose = verbose
    )
  }

  if (verbose) cat("done!\n")

  out <- new("nc_mcmc",
    mu = as.array(res[["mu"]]),
    delta = as.array(res[["delta"]]),
    d = as.array(res[["d"]]),
    Sigma = as.mcmc(res[["Sigma"]]),
    Gamma = as.array(res[["Gamma"]]),
    x = as.array(res[["x"]]),
    x_unadj = as.array(res[["x_unadj"]]),
    # a = as.array(res[["a"]]),
    # b = as.array(res[["b"]]),
    # D = as.array(res[["D"]]),
    # S_q = as.array(res[["S_q"]]),
    # Sigma_q_prop = as.array(res[["Sigma_q_prop"]]),
    # rho_mu = as.array(res[["rho_mu"]]),
    # cov_mu = as.array(res[["cov_mu"]]),
    # ar_mu_vec = as.array(res[["ar_mu_vec"]]),
    # mu_rate = as.array(res[["mu_rate"]]),
    # mu_prop = as.array(res[["mu_prop"]]),
    # tp_mu = as.array(res[["tp_mu"]]),
    # delta_prop = as.array(res[["delta_prop"]]),
    accept = as.list(res[["accept"]]),
    dens = list(loglik = as.numeric(res[["loglik"]]), logprior = as.numeric(res[["logprior"]]), logpost = as.numeric(res[["logpost"]])),
    control = list(burnin = burnin, nsim = nsim, nthin = nthin, prm_prop = prm_prop, prm_prior = prm_prior),
    dim = list(n_study = n, n_outcomes = M, n_datapoints = n_datapoints, n_treatments = slot(data, "n_treatments"), ref_trt = data@ref_trt),
    data = data,
    call = match.call()
  )

  return(out)
}

#' Estimation of copula based models for multivariate network meta-analysis.
#'
#' Estimation of copula based models for multivariate network meta-analysis.
#'
#' @param data Object of class \code{nc_data} containing the data to analyze.
#' @param burnin Length one integer vector providing the number of MCMC burnin
#' iterations.
#' @param nsim Length one integer vector providing the number of MCMC
#' iterations.
#' @param nthin Length one integer vector providing the thinning interval.
#' @param prm_prior Named list of parameters for the proposal distributions of
#' the \eqn{Z} latent positions and \eqn{\alpha} parameters; the elements'
#' names must be \code{"z"} and \code{"alpha"} and both need to be a single
#' number representing the standard deviations of the normal proposal
#' distributions for each \eqn{z_{ih}^{g}}{z_ih^g} and \eqn{\alpha_g}
#' respectively.
#' @param prm_prop Named list of hyperparameters for the prior distributions
#' of the \eqn{\sigma^2} and \eqn{\lambda} parameters; the elements' names
#' must be \code{"sigma2"} and \code{"lambda"} with the first one a two-sized
#' vector containing the \eqn{\sigma^2} prior hyperparameters and the second
#' element a vector providing the \eqn{\lambda} hyperparameters.
#' @param random_start Length-one logical vector; if \code{TRUE} (default), it
#' allows the algorithm to start from a random starting point.
#' @param init Named list providing user-defined starting values. If provided,
#' it overrides the value given for \code{random_start}.
#' @param tuning Named list of tuning parameters (see details).
#' @param adaptation Named list of adaptation parameters (see details).
#' @param verbose Length-one logical vector; if \code{TRUE}, messages
#' informing about the progress of the computations will be printed.
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
netcopula <- function(data, burnin = 10000, nsim = 5000, nthin = 1, prm_prior, prm_prop = NULL, prm_init, random_start = TRUE, init = NULL, tuning, adaptation, verbose = FALSE, mu_delta_mh = TRUE) {
  if (class(data) != "nc_data") {
    stop("the data argument must be a 'nc_data' object; see the help of the 'nc_data_create()' function.")
  }
  if (adaptation$maxiter > (burnin + nsim)/adaptation$every) {
    stop("the 'maxiter' adaptation parameter is too large.")
  }
  if (adaptation$miniter > adaptation$maxiter) {
    stop("the 'miniter' adaptation parameter must be smaller than or equal to 'maxiter'.")
  }
  if (adaptation$alpha < 0 | adaptation$alpha > 1 |
    adaptation$beta < 0 | adaptation$beta > 1 |
    adaptation$gamma < 0 | adaptation$gamma > 1 |
    adaptation$tar < 0 | adaptation$tar > 1) {
    stop("the 'alpha', 'beta', 'gamma' and 'tar' adaptation parameters must be in [0, 1].")
  }
  if (adaptation$tol <= 0) {
    stop("the 'tol' adaptation parameter must be positive.")
  }

  n <- slot(data, "n_study")
  M <- slot(data, "n_outcomes")
  n_datapoints <- slot(data, "n_datapoints")
  tot_iter <- burnin + nsim

  if (verbose) cat("Initialization of the algorithm...")

  if (is.null(prm_init)) {
    prm_init <- prm_prior
  }
  nc_start <- nc_init(data = data, prm_init = prm_init, random_start = random_start, init = init)

  if (verbose) cat("done!\n")

  # start iteration
  if (verbose) cat("Running the MCMC simulation...\n")

  if (mu_delta_mh) {
    res <- nc_mcmc_mh_new(
      data = data,
      init = nc_start,
      totiter = tot_iter,
      prior = prm_prior,
      prop = prm_prop,
      tuning = tuning,
      adapt = adaptation,
      verbose = verbose
    )
  } else {
    # TO DO: alla fine eliminare questa possibilita'
    res <- nc_mcmc_opt(
      data = data,
      init = nc_start,
      totiter = tot_iter,
      prior = prm_prior,
      prop = prm_prop,
      tuning = tuning,
      adapt = adaptation,
      verbose = verbose
    )
  }

  if (verbose) cat("done!\n")

  out <- new("nc_mcmc",
    mu = as.array(res[["mu"]]),
    delta = as.array(res[["delta"]]),
    d = as.array(res[["d"]]),
    Sigma = as.mcmc(res[["Sigma"]]),
    Gamma = as.array(res[["Gamma"]]),
    x = as.array(res[["x"]]),
    x_unadj = as.array(res[["x_unadj"]]),
    # a = as.array(res[["a"]]),
    # b = as.array(res[["b"]]),
    # D = as.array(res[["D"]]),
    # S_q = as.array(res[["S_q"]]),
    # Sigma_q_prop = as.array(res[["Sigma_q_prop"]]),
    # rho_mu = as.array(res[["rho_mu"]]),
    # cov_mu = as.array(res[["cov_mu"]]),
    # ar_mu_vec = as.array(res[["ar_mu_vec"]]),
    # mu_rate = as.array(res[["mu_rate"]]),
    # mu_prop = as.array(res[["mu_prop"]]),
    # tp_mu = as.array(res[["tp_mu"]]),
    # delta_prop = as.array(res[["delta_prop"]]),
    accept = as.list(res[["accept"]]),
    dens = list(loglik = as.numeric(res[["loglik"]]), logprior = as.numeric(res[["logprior"]]), logpost = as.numeric(res[["logpost"]])),
    control = list(burnin = burnin, nsim = nsim, nthin = nthin, prm_prop = prm_prop, prm_prior = prm_prior),
    dim = list(n_study = n, n_outcomes = M, n_datapoints = n_datapoints, n_treatments = slot(data, "n_treatments"), ref_trt = data@ref_trt),
    data = data,
    call = match.call()
  )

  return(out)
}

#' Estimation of copula based models for multivariate network meta-analysis.
#'
#' Estimation of copula based models for multivariate network meta-analysis.
#'
#' @param data Object of class \code{nc_data} containing the data to analyze.
#' @param burnin Length one integer vector providing the number of MCMC burnin
#' iterations.
#' @param nsim Length one integer vector providing the number of MCMC
#' iterations.
#' @param nthin Length one integer vector providing the thinning interval.
#' @param prm_prior Named list of parameters for the proposal distributions of
#' the \eqn{Z} latent positions and \eqn{\alpha} parameters; the elements'
#' names must be \code{"z"} and \code{"alpha"} and both need to be a single
#' number representing the standard deviations of the normal proposal
#' distributions for each \eqn{z_{ih}^{g}}{z_ih^g} and \eqn{\alpha_g}
#' respectively.
#' @param prm_prop Named list of hyperparameters for the prior distributions
#' of the \eqn{\sigma^2} and \eqn{\lambda} parameters; the elements' names
#' must be \code{"sigma2"} and \code{"lambda"} with the first one a two-sized
#' vector containing the \eqn{\sigma^2} prior hyperparameters and the second
#' element a vector providing the \eqn{\lambda} hyperparameters.
#' @param random_start Length-one logical vector; if \code{TRUE} (default), it
#' allows the algorithm to start from a random starting point.
#' @param init Named list providing user-defined starting values. If provided,
#' it overrides the value given for \code{random_start}.
#' @param tuning Named list of tuning parameters (see details).
#' @param adaptation Named list of adaptation parameters (see details).
#' @param verbose Length-one logical vector; if \code{TRUE}, messages
#' informing about the progress of the computations will be printed.
#'
#' @return \code{netcopula_new} returns an object of class \code{\link{nc_mcmc}}.
#' @export
#'
#' @examples
#' \dontrun{
#' # netcopula_new(...)
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
netcopula_new <- function(data, burnin = 10000, nsim = 5000, nthin = 1, prm_prior, prm_prop = NULL, prm_init, random_start = TRUE, init = NULL, tuning, adaptation, verbose = FALSE, mu_delta_mh = TRUE) {
  if (class(data) != "nc_data") {
    stop("the data argument must be a 'nc_data' object; see the help of the 'nc_data_create()' function.")
  }
  if (adaptation$maxiter > (burnin + nsim)/adaptation$every) {
    stop("the 'maxiter' adaptation parameter is too large.")
  }
  if (adaptation$miniter > adaptation$maxiter) {
    stop("the 'miniter' adaptation parameter must be smaller than or equal to 'maxiter'.")
  }
  if (adaptation$alpha < 0 | adaptation$alpha > 1 |
    adaptation$beta < 0 | adaptation$beta > 1 |
    adaptation$gamma < 0 | adaptation$gamma > 1 |
    adaptation$tar < 0 | adaptation$tar > 1) {
    stop("the 'alpha', 'beta', 'gamma' and 'tar' adaptation parameters must be in [0, 1].")
  }
  if (adaptation$tol <= 0) {
    stop("the 'tol' adaptation parameter must be positive.")
  }

  n <- slot(data, "n_study")
  M <- slot(data, "n_outcomes")
  n_datapoints <- slot(data, "n_datapoints")
  tot_iter <- burnin + nsim

  if (verbose) cat("Initialization of the algorithm...")

  if (is.null(prm_init)) {
    prm_init <- prm_prior
  }
  if (!is.list(init) & !is.null(init)) {
    if (init == "univariate") {
      nc_start <- nc_init(data = data, prm_init = prm_init, random_start = random_start, init = init, prm_prior = prm_prior, prm_prop = prm_prop, tuning = tuning, adaptation = adaptation)
    }
  } else {
    nc_start <- nc_init(data = data, prm_init = prm_init, random_start = random_start, init = init)
  }

  if (verbose) cat("done!\n")

  # start iteration
  if (verbose) cat("Running the MCMC simulation...\n")

  if (mu_delta_mh) {
    res <- nc_mcmc_mh_new2(
      data = data,
      init = nc_start,
      totiter = tot_iter,
      prior = prm_prior,
      prop = prm_prop,
      tuning = tuning,
      adapt = adaptation,
      verbose = verbose
    )
  } else {
    # TO DO: alla fine eliminare questa possibilita'
    res <- nc_mcmc_opt(
      data = data,
      init = nc_start,
      totiter = tot_iter,
      prior = prm_prior,
      prop = prm_prop,
      tuning = tuning,
      adapt = adaptation,
      verbose = verbose
    )
  }

  if (verbose) cat("done!\n")

  out <- new("nc_mcmc",
    mu = as.array(res[["mu"]]),
    delta = as.array(res[["delta"]]),
    d = as.array(res[["d"]]),
    Sigma = as.mcmc(res[["Sigma"]]),
    Gamma = as.array(res[["Gamma"]]),
    r = as.array(res[["r"]]),
    x = as.array(res[["x"]]),
    x_unadj = as.array(res[["x_unadj"]]),
    # a = as.array(res[["a"]]),
    # b = as.array(res[["b"]]),
    # D = as.array(res[["D"]]),
    # S_q = as.array(res[["S_q"]]),
    # Sigma_q_prop = as.array(res[["Sigma_q_prop"]]),
    # rho_mu = as.array(res[["rho_mu"]]),
    # cov_mu = as.array(res[["cov_mu"]]),
    # ar_mu_vec = as.array(res[["ar_mu_vec"]]),
    # mu_rate = as.array(res[["mu_rate"]]),
    # mu_prop = as.array(res[["mu_prop"]]),
    # tp_mu = as.array(res[["tp_mu"]]),
    # delta_prop = as.array(res[["delta_prop"]]),
    accept = as.list(res[["accept"]]),
    dens = list(loglik = as.numeric(res[["loglik"]]), logprior = as.numeric(res[["logprior"]]), logpost = as.numeric(res[["logpost"]])),
    control = list(burnin = burnin, nsim = nsim, nthin = nthin, prm_prop = prm_prop, prm_prior = prm_prior),
    dim = list(n_study = n, n_outcomes = M, n_datapoints = n_datapoints, n_treatments = slot(data, "n_treatments"), ref_trt = data@ref_trt),
    data = data,
    call = match.call()
  )

  return(out)
}

#' Estimation of copula based models for multivariate network meta-analysis.
#'
#' Estimation of copula based models for multivariate network meta-analysis.
#'
#' @param data Object of class \code{nc_data} containing the data to analyze.
#' @param burnin Length one integer vector providing the number of MCMC burnin
#' iterations.
#' @param nsim Length one integer vector providing the number of MCMC
#' iterations.
#' @param nthin Length one integer vector providing the thinning interval.
#' @param prm_prior Named list of parameters for the proposal distributions of
#' the \eqn{Z} latent positions and \eqn{\alpha} parameters; the elements'
#' names must be \code{"z"} and \code{"alpha"} and both need to be a single
#' number representing the standard deviations of the normal proposal
#' distributions for each \eqn{z_{ih}^{g}}{z_ih^g} and \eqn{\alpha_g}
#' respectively.
#' @param prm_prop Named list of hyperparameters for the prior distributions
#' of the \eqn{\sigma^2} and \eqn{\lambda} parameters; the elements' names
#' must be \code{"sigma2"} and \code{"lambda"} with the first one a two-sized
#' vector containing the \eqn{\sigma^2} prior hyperparameters and the second
#' element a vector providing the \eqn{\lambda} hyperparameters.
#' @param random_start Length-one logical vector; if \code{TRUE} (default), it
#' allows the algorithm to start from a random starting point.
#' @param init Named list providing user-defined starting values. If provided,
#' it overrides the value given for \code{random_start}.
#' @param tuning Named list of tuning parameters (see details).
#' @param adaptation Named list of adaptation parameters (see details).
#' @param verbose Length-one logical vector; if \code{TRUE}, messages
#' informing about the progress of the computations will be printed.
#'
#' @return \code{netcopula_uni} returns an object of class \code{\link{nc_mcmc}}.
#' @export
#'
#' @examples
#' \dontrun{
#' # netcopula_uni(...)
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
netcopula_uni <- function(data, burnin = 10000, nsim = 5000, nthin = 1, prm_prior, prm_prop = NULL, prm_init, random_start = TRUE, init = NULL, tuning, adaptation, verbose = FALSE, mu_delta_mh = TRUE) {
  if (class(data) != "nc_data") {
    stop("the data argument must be a 'nc_data' object; see the help of the 'nc_data_create()' function.")
  }
  if (adaptation$maxiter > (burnin + nsim)/adaptation$every) {
    stop("the 'maxiter' adaptation parameter is too large.")
  }
  if (adaptation$miniter > adaptation$maxiter) {
    stop("the 'miniter' adaptation parameter must be smaller than or equal to 'maxiter'.")
  }
  if (adaptation$alpha < 0 | adaptation$alpha > 1 |
    adaptation$beta < 0 | adaptation$beta > 1 |
    adaptation$gamma < 0 | adaptation$gamma > 1 |
    adaptation$tar < 0 | adaptation$tar > 1) {
    stop("the 'alpha', 'beta', 'gamma' and 'tar' adaptation parameters must be in [0, 1].")
  }
  if (adaptation$tol <= 0) {
    stop("the 'tol' adaptation parameter must be positive.")
  }

  n <- slot(data, "n_study")
  M <- slot(data, "n_outcomes")
  n_datapoints <- slot(data, "n_datapoints")
  tot_iter <- burnin + nsim

  if (verbose) cat("Initialization of the algorithm...")

  if (is.null(prm_init)) {
    prm_init <- prm_prior
  }
  nc_start <- nc_init(data = data, prm_init = prm_init, random_start = random_start, init = init)

  if (verbose) cat("done!\n")

  # start iteration
  if (verbose) cat("Running the MCMC simulation...\n")

  if (mu_delta_mh) {
    res <- nc_mcmc_mh_new2(
      data = data,
      init = nc_start,
      totiter = tot_iter,
      prior = prm_prior,
      prop = prm_prop,
      tuning = tuning,
      adapt = adaptation,
      verbose = verbose
    )
  } else {
    # TO DO: alla fine eliminare questa possibilita'
    res <- nc_mcmc_opt(
      data = data,
      init = nc_start,
      totiter = tot_iter,
      prior = prm_prior,
      prop = prm_prop,
      tuning = tuning,
      adapt = adaptation,
      verbose = verbose
    )
  }

  if (verbose) cat("done!\n")

  out <- new("nc_mcmc",
    mu = as.array(res[["mu"]]),
    delta = as.array(res[["delta"]]),
    d = as.array(res[["d"]]),
    Sigma = as.mcmc(res[["Sigma"]]),
    Gamma = as.array(res[["Gamma"]]),
    r = as.array(res[["r"]]),
    x = as.array(res[["x"]]),
    x_unadj = as.array(res[["x_unadj"]]),
    # a = as.array(res[["a"]]),
    # b = as.array(res[["b"]]),
    # D = as.array(res[["D"]]),
    # S_q = as.array(res[["S_q"]]),
    # Sigma_q_prop = as.array(res[["Sigma_q_prop"]]),
    # rho_mu = as.array(res[["rho_mu"]]),
    # cov_mu = as.array(res[["cov_mu"]]),
    # ar_mu_vec = as.array(res[["ar_mu_vec"]]),
    # mu_rate = as.array(res[["mu_rate"]]),
    # mu_prop = as.array(res[["mu_prop"]]),
    # tp_mu = as.array(res[["tp_mu"]]),
    # delta_prop = as.array(res[["delta_prop"]]),
    accept = as.list(res[["accept"]]),
    dens = list(loglik = as.numeric(res[["loglik"]]), logprior = as.numeric(res[["logprior"]]), logpost = as.numeric(res[["logpost"]])),
    control = list(burnin = burnin, nsim = nsim, nthin = nthin, prm_prop = prm_prop, prm_prior = prm_prior),
    dim = list(n_study = n, n_outcomes = M, n_datapoints = n_datapoints, n_treatments = slot(data, "n_treatments"), ref_trt = data@ref_trt),
    data = data,
    call = match.call()
  )

  return(out)
}

#' Estimation of copula based models for multivariate network meta-analysis.
#'
#' Estimation of copula based models for multivariate network meta-analysis.
#'
#' @param data Object of class \code{nc_data} containing the data to analyze.
#' @param burnin Number of MCMC burnin iterations.
#' @param nsim Number of MCMC iterations.
#' @param prm_prior Named list of parameters for the proposal distributions of
#' the \eqn{Z} latent positions and \eqn{\alpha} parameters; the elements'
#' names must be \code{"z"} and \code{"alpha"} and both need to be a single
#' number representing the standard deviations of the normal proposal
#' distributions for each \eqn{z_{ih}^{g}}{z_ih^g} and \eqn{\alpha_g}
#' respectively.
#' @param prm_prop Named list of hyperparameters for the prior distributions
#' of the \eqn{\sigma^2} and \eqn{\lambda} parameters; the elements' names
#' must be \code{"sigma2"} and \code{"lambda"} with the first one a two-sized
#' vector containing the \eqn{\sigma^2} prior hyperparameters and the second
#' element a vector providing the \eqn{\lambda} hyperparameters.
#' @param tuning Named list of tuning parameters.
#' @param random_start Length-one logical vector; if \code{TRUE} (default), it
#' allows the algorithm to start from a random starting point.
#' @param init Named list providing user-defined starting values. If provided,
#' it overrides the value given for \code{random_start}.
#' @param verbose Length-one logical vector; if \code{TRUE}, messages
#' informing about the progress of the computations will be printed.
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
as.mcmc.nc_mcmc <- function(res, latent_var = FALSE, ...) {
  if (!inherits(res, "nc_mcmc")) {
    stop("method as.mcmc.nc_mcmc() is only intended for nc_mcmc objects.")
  }

  moreargs <- list(...)
  start <- res@control$burnin + 1
  end <- res@control$burnin + res@control$nsim
  if (!is.null(res@control$nthin)) {
    thin <- res@control$nthin
  } else {
    if (length(moreargs) & !is.na(match("nthin", names(moreargs)))) {
      thin <- moreargs[[match("nthin", names(moreargs))]]
    } else {
      stop("to use as.mcmc.nc_mcmc() you must provide the nthin argument.")
    }
  }
  ns <- res@dim$n_study
  M <- res@dim$n_outcomes
  nt <- res@dim$n_treatments
  n <- res@dim$n_datapoints
  nGamma <- dim(res@Gamma)[[1]]
  nchains <- 1
  data_tmp <- try(is.null(res@data), silent = TRUE)
  if (class(data_tmp) != "try-error") {
    ref_trt <- res@dim$ref_trt
    data <- res@data
  } else {
    if (length(moreargs) & !is.na(match("data", names(moreargs)))) {
      data <- moreargs[[match("data", names(moreargs))]]
      ref_trt <- data@ref_trt
    } else {
      stop("to use as.mcmc.nc_mcmc() you must provide the data argument.")
    }
  }

  mu <- res@mu
  nm <- character(ns*M)
  for (j in 1:M) {
    for (i in 1:ns) {
     nm[i + ns*(j - 1)] <- paste0("mu[", i, ",", j, "]")
    }
  }
  attr(mu, "dim") <- c(dim(mu)[1]*dim(mu)[2], dim(mu)[3])
  mu <- t(mu)
  colnames(mu) <- nm

  delta <- res@delta
  nm <- character(n*M)
  for (j in 1:M) {
    for (i in 1:n) {
      nm[i + n*(j - 1)] <- paste0("delta[", i, ",", j, "]")
    }
  }
  attr(delta, "dim") <- c(dim(delta)[1]*dim(delta)[2], dim(delta)[3])
  delta <- t(delta)
  colnames(delta) <- nm
  delta <- delta[, rep(data@study_data$trt != data@study_data$baseline, M)]

  d <- res@d
  nm <- character(nt*M)
  for (j in 1:M) {
    for (i in 1:nt) {
      nm[i + nt*(j - 1)] <- paste0("d[", i, ",", j, "]")
    }
  }
  attr(d, "dim") <- c(dim(d)[1]*dim(d)[2], dim(d)[3])
  d <- t(d)
  colnames(d) <- nm
  trt_tmp <- setdiff(1:nt, ref_trt)
  d_tmp <- trt_tmp
  if (M > 1) {
    for (i in 1:(M - 1)) {
      d_tmp <- append(d_tmp, trt_tmp + nt*i)
    }
  }
  d <- d[, d_tmp]

  Sigma <- res@Sigma
  nm <- character(M*(M + 1)/2)
  count <- 1
  for (j in 1:M) {
    for (i in j:M) {
      nm[count] <- paste0("Sigma_M[", i, ",", j, "]")
      count <- count + 1
    }
  }
  colnames(Sigma) <- nm

  Gamma <- res@Gamma
  nelem <- ifelse(M > 1, M*(M - 1)/2, 1)
  nm <- character(nGamma*nelem)
  count <- 1
  for (h in 1:nGamma) {
    if (M > 1) {
      for (j in 1:(M - 1)) {
        for (i in (j + 1):M) {
          nm[count + nelem*(h - 1)] <- paste0("Gamma[", h, ",", i, ",", j, "]")
          count <- count + 1
        }
      }
    } else {
      nm[count + h - 1] <- paste0("Gamma[", h, ",", i, ",", j, "]")
    }
    count <- 1
  }
  attr(Gamma, "dim") <- c(dim(Gamma)[1]*dim(Gamma)[2], dim(Gamma)[3])
  Gamma <- t(Gamma)
  colnames(Gamma) <- nm

  if (.hasSlot(res, "S_q")) {
    if(!all(is.na(res@S_q))) {
      S_q <- res@S_q
      nelem <- M*(M + 1)/2
      nm <- character(nGamma*nelem)
      count <- 1
      for (h in 1:nGamma) {
        if (M > 1) {
          for (j in 1:M) {
            for (i in j:M) {
              nm[count + nelem*(h - 1)] <- paste0("S_q[", h, ",", i, ",", j, "]")
              count <- count + 1
            }
          }
        } else {
          nm[count + h - 1] <- paste0("S_q[", h, ",", i, ",", j, "]")
        }
        count <- 1
      }
      attr(S_q, "dim") <- c(dim(S_q)[1]*dim(S_q)[2], dim(S_q)[3])
      S_q <- t(S_q)
      colnames(S_q) <- nm

      Sigma_q_prop <- res@Sigma_q_prop
      nelem <- M*(M + 1)/2
      nm <- character(nGamma*nelem)
      count <- 1
      for (h in 1:nGamma) {
        if (M > 1) {
          for (j in 1:M) {
            for (i in j:M) {
              nm[count + nelem*(h - 1)] <- paste0("Sigma_q_prop[", h, ",", i, ",", j, "]")
              count <- count + 1
            }
          }
        } else {
          nm[count + h - 1] <- paste0("Sigma_q_prop[", h, ",", i, ",", j, "]")
        }
        count <- 1
      }
      attr(Sigma_q_prop, "dim") <- c(dim(Sigma_q_prop)[1]*dim(Sigma_q_prop)[2], dim(Sigma_q_prop)[3])
      Sigma_q_prop <- t(Sigma_q_prop)
      colnames(Sigma_q_prop) <- nm
    }
  }

  if (latent_var) {
    x <- res@x
    x_na <- as.logical(is.na(x[, , 1]))
    nm <- character(n*M)
    for (j in 1:M) {
      for (i in 1:n) {
        nm[i + n*(j - 1)] <- paste0("x[", i, ",", j, "]")
      }
    }
    attr(x, "dim") <- c(dim(x)[1]*dim(x)[2], dim(x)[3])
    x <- t(x)
    colnames(x) <- nm
    x <- x[, !x_na]

    x_unadj <- res@x_unadj
    x_na <- as.logical(is.na(x_unadj[, , 1]))
    nm <- character(n*M)
    for (j in 1:M) {
      for (i in 1:n) {
        nm[i + n*(j - 1)] <- paste0("x_unadj[", i, ",", j, "]")
      }
    }
    attr(x_unadj, "dim") <- c(dim(x_unadj)[1]*dim(x_unadj)[2], dim(x_unadj)[3])
    x_unadj <- t(x_unadj)
    colnames(x_unadj) <- nm
    x_unadj <- x_unadj[, !x_na]

    if (.hasSlot(res, "S_q")) {
      if (!all(is.na(res@S_q))) {
        tmp <- cbind(mu, delta, d, Sigma, Gamma, S_q, Sigma_q_prop, x, x_unadj)
      } else {
        tmp <- cbind(mu, delta, d, Sigma, Gamma, x, x_unadj)
      }
    } else {
      tmp <- cbind(mu, delta, d, Sigma, Gamma, x, x_unadj)
    }
  } else {
    if (.hasSlot(res, "S_q")) {
      if (!all(is.na(res@S_q))) {
        tmp <- cbind(mu, delta, d, Sigma, Gamma, S_q, Sigma_q_prop)
      } else {
        tmp <- cbind(mu, delta, d, Sigma, Gamma)
      }
    } else {
      tmp <- cbind(mu, delta, d, Sigma, Gamma)
    }
  }
  to_keep <- seq(start, end, by = thin)
  ord <- order(dimnames(tmp)[[2]])
  tmp <- tmp[to_keep, ]
  # tmp <- tmp[to_keep, ord]
  out <- mcmc(tmp, start = start, end = end, thin = thin)

  return(out)
}

#' Estimation of copula based models for multivariate network meta-analysis.
#'
#' Estimation of copula based models for multivariate network meta-analysis.
#'
#' @param data Object of class \code{nc_data} containing the data to analyze.
#' @param burnin Number of MCMC burnin iterations.
#' @param nsim Number of MCMC iterations.
#' @param prm_prior Named list of parameters for the proposal distributions of
#' the \eqn{Z} latent positions and \eqn{\alpha} parameters; the elements'
#' names must be \code{"z"} and \code{"alpha"} and both need to be a single
#' number representing the standard deviations of the normal proposal
#' distributions for each \eqn{z_{ih}^{g}}{z_ih^g} and \eqn{\alpha_g}
#' respectively.
#' @param prm_prop Named list of hyperparameters for the prior distributions
#' of the \eqn{\sigma^2} and \eqn{\lambda} parameters; the elements' names
#' must be \code{"sigma2"} and \code{"lambda"} with the first one a two-sized
#' vector containing the \eqn{\sigma^2} prior hyperparameters and the second
#' element a vector providing the \eqn{\lambda} hyperparameters.
#' @param tuning Named list of tuning parameters.
#' @param random_start Length-one logical vector; if \code{TRUE} (default), it
#' allows the algorithm to start from a random starting point.
#' @param init Named list providing user-defined starting values. If provided,
#' it overrides the value given for \code{random_start}.
#' @param verbose Length-one logical vector; if \code{TRUE}, messages
#' informing about the progress of the computations will be printed.
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
mu_logpost_func <- function(mu, varargs) {
  delta <- varargs$delta
  y <- varargs$y
  n <- varargs$n
  w <- varargs$w
  gamma <- varargs$gamma
  mu_sigma <- varargs$mu_sigma
  eps <- varargs$eps
  eps_ab <- varargs$eps_ab

  return(mu_logpost(mu, delta, y, n, w, gamma, mu_sigma, eps, eps_ab))
}

#' Estimation of copula based models for multivariate network meta-analysis.
#'
#' Estimation of copula based models for multivariate network meta-analysis.
#'
#' @param data Object of class \code{nc_data} containing the data to analyze.
#' @param burnin Number of MCMC burnin iterations.
#' @param nsim Number of MCMC iterations.
#' @param prm_prior Named list of parameters for the proposal distributions of
#' the \eqn{Z} latent positions and \eqn{\alpha} parameters; the elements'
#' names must be \code{"z"} and \code{"alpha"} and both need to be a single
#' number representing the standard deviations of the normal proposal
#' distributions for each \eqn{z_{ih}^{g}}{z_ih^g} and \eqn{\alpha_g}
#' respectively.
#' @param prm_prop Named list of hyperparameters for the prior distributions
#' of the \eqn{\sigma^2} and \eqn{\lambda} parameters; the elements' names
#' must be \code{"sigma2"} and \code{"lambda"} with the first one a two-sized
#' vector containing the \eqn{\sigma^2} prior hyperparameters and the second
#' element a vector providing the \eqn{\lambda} hyperparameters.
#' @param tuning Named list of tuning parameters.
#' @param random_start Length-one logical vector; if \code{TRUE} (default), it
#' allows the algorithm to start from a random starting point.
#' @param init Named list providing user-defined starting values. If provided,
#' it overrides the value given for \code{random_start}.
#' @param verbose Length-one logical vector; if \code{TRUE}, messages
#' informing about the progress of the computations will be printed.
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
delta_logpost_func <- function(delta, varargs) {
  mu <- varargs$mu
  tau <- varargs$tau
  eta <- varargs$eta
  y <- varargs$y
  n <- varargs$n
  w <- varargs$w
  gamma <- varargs$gamma
  eps <- varargs$eps
  eps_ab <- varargs$eps_ab

  return(delta_logpost(delta, mu, tau, eta, y, n, w, gamma, eps, eps_ab))
}

#' @export
mu_logpost_quad <- function(mu, varargs) {
  coef <- varargs$coef
  x <- c(1, mu, mu^2)

  return(sum(coef*x))
}

#' @export
delta_logpost_quad <- function(delta, varargs) {
  coef <- varargs$coef
  x <- c(1, delta, delta^2)

  return(sum(coef*x))
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
#' @return \code{init_univariate} returns an object of class \code{\link{nc_mcmc}}.
#' @export
#'
#' @examples
#' \dontrun{
#' # init_univariate(...)
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
init_univariate <- function(data, burnin, nsim, nthin, prm_prior, prm_prop, prm_init, tuning, adaptation, verbose) {
  M <- data@n_outcomes
  n_study <- data@n_study
  n_datapoints <- data@n_datapoints
  n_treatments <- data@n_treatments
  ref_trt <- data@ref_trt
  eps <- .001

  res_uni <- list()
  mu <- matrix(NA, nrow = n_study, ncol = M)
  delta <- x <- matrix(NA, nrow = n_datapoints, ncol = M)
  d <- matrix(NA, nrow = n_treatments, ncol = M)
  Sigma_M <- matrix(0, nrow = M, ncol = M)
  x_all <- array(NA, dim = c(n_datapoints, M, burnin + nsim))
  y <- matrix(NA, nrow = n_datapoints, ncol = M)

  data_tmp <- data@study_data
  for (m in 1:M) {
    data_m <- data_tmp[, c(2:6, 6 + m, 9 + m)]
    colnames(data_m)[c(6, 7)] <- c("n1", "y1")
    nc_data <- nc_data_create("binary", n_study, 1, n_datapoints, n_treatments, NULL, data_m, ref_trt)
    
    # MCMC simulation
    res_uni[[m]] <- netcopula_uni(nc_data, burnin = burnin, nsim = nsim, nthin = nthin, prm_prior, prm_prop, prm_init, tuning = tuning, adaptation = adaptation, verbose = verbose, mu_delta_mh = TRUE, random_start = TRUE, init = NULL)

    mu[, m] <- apply(drop(res_uni[[m]]@mu), 1, mean)
    delta[, m] <- apply(drop(res_uni[[m]]@delta), 1, mean)
    d[, m] <- apply(drop(res_uni[[m]]@d), 1, mean)
    Sigma_M[m, m] <- mean(drop(res_uni[[m]]@Sigma))

    mu_ikm <- res_uni[[m]]@mu
    delta_ikm <- drop(res_uni[[m]]@delta)
    mu_ikm <- apply(mu_ikm, 3, param_long, nc_data@study_id$narms, FALSE)
    theta_ikm <- mu_ikm + delta_ikm
    p_ikm <- apply(theta_ikm, 2, expit_rcpp)
    n_ik <- nc_data@study_data$n1
    y_ik <- nc_data@study_data$y1
    y[, m] <- y_ik
    phi_ikm <- x_ikm <- matrix(NA, nrow = n_datapoints, ncol = burnin + nsim)
    for (ik in 1:n_datapoints) {
      if (!is.na(y_ik[ik])) {
        phi_ikm[ik, ] <- pbinom(y_ik[ik], size = n_ik[ik], prob = p_ikm[ik, ])
        phi_ikm[ik, ][phi_ikm[ik, ] == 1] <- 1 - eps
        phi_ikm[ik, ][phi_ikm[ik, ] == 0] <- eps
        x_ikm[ik, ] <- qnorm(phi_ikm[ik, ])
      } else {
        x_ikm[ik, ] <- NA
      }
    }
    x_all[, m, ] <- x_ikm
  }
  Corr_M <- diag(n_outcomes)
  Corr_M[lower.tri(Corr_M)] <- Corr_M[upper.tri(Corr_M)] <- rep(0.5, M)
  Sigma_M <- Sigma_M^0.5 %*% Corr_M %*% Sigma_M^0.5

  x <- apply(x_all, c(1, 2), mean, na.rm = TRUE)
  x[is.na(y)] <- NA # sets NA in the x matrix

  Gamma <- D <- list()
  for (q in 1:prm_init$nGamma) {
    if (prm_init$nGamma == n_trt) {
      if (M > 1) {
        Gamma[[q]] <- cor(x[data@study_data[, "trt"] == q, ])
        D[[q]] <- as.matrix(diag(sqrt(1/colSums(x[data@study_data[, "trt"] == q, ]^2, na.rm = TRUE))))
      } else {
        Gamma[[q]] <- as.matrix(1)
        D[[q]] <- as.matrix(sqrt(1/colSums(x[data@study_data[, "trt"] == q, ]^2, na.rm = TRUE)))
      }
    } else if (prm_init$nGamma == 1) {
      if (M > 1) {
        Gamma[[q]] <- cor(x)
        D[[q]] <- as.matrix(diag(sqrt(1/colSums(x^2, na.rm = TRUE))))
      } else {
        Gamma[[q]] <- as.matrix(1)
        D[[q]] <- as.matrix(sqrt(1/colSums(x^2, na.rm = TRUE)))
      }
    }
  }

  return(list(mu = mu, delta = delta, d = d, Sigma_M = Sigma_M, Gamma = Gamma, D = D, x = x))
}
