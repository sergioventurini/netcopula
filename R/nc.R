#' Estimation of a Copula-Based Model for Multivariate Network Meta-Analysis.
#'
#' Estimation of a Copula-Based Model for Multivariate Network Meta-Analysis.
#'
#' @param data Object of class \code{nc_data} containing the data to analyze.
#' @param control A list of control parameters that affect the sampling
#'   but do not affect the posterior distribution. See
#'   \code{\link{nc_control}()} for more details.
#' @param prior A list containing the prior hyperparameters. See
#'   \code{\link{nc_prior}()} for more details.
#' @param prm_init An object of class \code{list} containing parameters used
#'   for initializing the simulation.
#' @param init Named \code{list} providing user-defined starting values.
#'   If provided, it overrides the value given for \code{random_start}.
#' @param tuning Named \code{list} of tuning parameters (see details).
#' @param adaptation Named \code{list} of adaptation parameters (see details).
#'
#' @return \code{netcopula} returns an object of class \code{\link{nc_mcmc}}.
#' @export
#'
#' @examples
#' \dontrun{
#' demo("example_homesafety", package = "netcopula")
#' }
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @seealso
#' \code{\link{nc_data-class}},
#' \code{\link{nc_mcmc}}
netcopula <- function(data, control = nc_control(), prior, prm_init, init = NULL, tuning, adaptation) {
  if (class(data) != "nc_data") {
    stop("the data argument must be a 'nc_data' object; type ?nc_data-class for more help.")
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

  control <- check_list_na(control, nc_control())
  if (!check_control(control))
    stop("the control list is not correct; see the documentation for more details.")

  if (is.null(prior)) {
    prior <- nc_prior()
  } else {
    prior <- check_list_na(prior, nc_prior())
  }
  if (!check_prior(prior)) {
    stop("the prior hyperparameter list is not correct; see the documentation for more details.")
  }

  n <- slot(data, "n_study")
  M <- slot(data, "n_outcomes")
  n_datapoints <- slot(data, "n_datapoints")
  nsim <- control[["nsim"]]
  burnin <- control[["burnin"]]
  thin <- control[["thin"]]
  tot_iter <- burnin + nsim
  random_start <- control[["random.start"]]
  store_burnin <- control[["store.burnin"]]
  verbose <- control[["verbose"]]
  prm_prop <- list(Gamma_update = control[["Gamma.update"]], eta_prop = control[["eta.prop"]],
    sigma_r_prop = control[["sigma.r.prop"]])
  prm_prior <- prior

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

  if (verbose) cat("done!\n")

  out <- new("nc_mcmc",
    mu = as.array(res[["mu"]]),
    delta = as.array(res[["delta"]]),
    d = as.array(res[["d"]]),
    Sigma = as.array(res[["Sigma"]]),
    Gamma = as.array(res[["Gamma"]]),
    r = as.array(res[["r"]]),
    x = as.array(res[["x"]]),
    x_unadj = as.array(res[["x_unadj"]]),
    accept = as.list(res[["accept"]]),
    dens = list(loglik = as.numeric(res[["loglik"]]), logprior = as.numeric(res[["logprior"]]), logpost = as.numeric(res[["logpost"]])),
    control = control,
    prior = prm_prior,
    dim = list(n_study = n, n_outcomes = M, n_datapoints = n_datapoints, n_treatments = slot(data, "n_treatments"), ref_trt = data@ref_trt),
    data = data,
    call = match.call()
  )

  return(out)
}

#' Estimation of a Copula-Based Model for Multivariate Network Meta-Analysis.
#'
#' Estimation of a Copula-Based Model for Multivariate Network Meta-Analysis.
#'
#' @param data Object of class \code{nc_data} containing the data to analyze.
#' @param control A list of control parameters that affect the sampling
#'   but do not affect the posterior distribution. See
#'   \code{\link{nc_control}()} for more details.
#' @param prior A list containing the prior hyperparameters. See
#'   \code{\link{nc_prior}()} for more details.
#' @param prm_init An object of class \code{list} containing parameters used
#'   for initializing the simulation.
#' @param init Named \code{list} providing user-defined starting values.
#'   If provided, it overrides the value given for \code{random_start}.
#' @param tuning Named \code{list} of tuning parameters (see details).
#' @param adaptation Named \code{list} of adaptation parameters (see details).
#'
#' @return \code{netcopula} returns an object of class \code{\link{nc_mcmc}}.
#' @export
#'
#' @examples
#' \dontrun{
#' demo("example_homesafety", package = "netcopula")
#' }
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @seealso
#' \code{\link{nc_data-class}},
#' \code{\link{nc_mcmc}}
netcopula_uni <- function(data, control = nc_control(), prior, prm_init, init = NULL, tuning, adaptation) {
  if (class(data) != "nc_data") {
    stop("the data argument must be a 'nc_data' object; type ?nc_data-class for more help.")
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

  control <- check_list_na(control, nc_control())
  if (!check_control(control))
    stop("the control list is not correct; see the documentation for more details.")

  if (is.null(prior)) {
    prior <- nc_prior()
  } else {
    prior <- check_list_na(prior, nc_prior())
  }
  if (!check_prior(prior)) {
    stop("the prior hyperparameter list is not correct; see the documentation for more details.")
  }

  n <- slot(data, "n_study")
  M <- slot(data, "n_outcomes")
  n_datapoints <- slot(data, "n_datapoints")
  nsim <- control[["nsim"]]
  burnin <- control[["burnin"]]
  thin <- control[["thin"]]
  tot_iter <- burnin + nsim
  random_start <- control[["random.start"]]
  store_burnin <- control[["store.burnin"]]
  verbose <- control[["verbose"]]
  prm_prop <- list(Gamma_update = control[["Gamma.update"]], eta_prop = control[["eta.prop"]],
    sigma_r_prop = control[["sigma.r.prop"]])
  prm_prior <- prior

  if (verbose) cat("Initialization of the algorithm...")

  if (is.null(prm_init)) {
    prm_init <- prm_prior
  }
  nc_start <- nc_init(data = data, prm_init = prm_init, random_start = random_start, init = init)

  if (verbose) cat("done!\n")

  # start iteration
  if (verbose) cat("Running the MCMC simulation...\n")

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

  if (verbose) cat("done!\n")

  out <- new("nc_mcmc",
    mu = as.array(res[["mu"]]),
    delta = as.array(res[["delta"]]),
    d = as.array(res[["d"]]),
    Sigma = as.array(res[["Sigma"]]),
    Gamma = as.array(res[["Gamma"]]),
    r = as.array(res[["r"]]),
    x = as.array(res[["x"]]),
    x_unadj = as.array(res[["x_unadj"]]),
    accept = as.list(res[["accept"]]),
    dens = list(loglik = as.numeric(res[["loglik"]]), logprior = as.numeric(res[["logprior"]]), logpost = as.numeric(res[["logpost"]])),
    control = control,
    prior = prm_prior,
    dim = list(n_study = n, n_outcomes = M, n_datapoints = n_datapoints, n_treatments = slot(data, "n_treatments"), ref_trt = data@ref_trt),
    data = data,
    call = match.call()
  )

  return(out)
}
