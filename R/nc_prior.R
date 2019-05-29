#' Auxiliary Function for Setting Model Priors.
#' 
#' @description{
#' \code{nc_prior()} is an auxiliary function as user interface for
#'   \code{netcopula()} fitting. Typically only used when calling the
#'   \code{netcopula()} function. It is used to set prior hyperparameters.
#' 
#' \code{prior_nc()} is an alias for \code{nc_prior()}.
#' 
#' \code{check_prior()} is an auxiliary function that verifies the
#'   correctness of the prior hyperparameters provided before a model is
#'   fitted with \code{\link{netcopula}()}.
#' }
#' 
#' @param mu_sigma2 A length-one numeric vector representing the study-specific
#'   baseline effects prior variance.
#'   contain two numeric vectors, namely \code{a} and \code{b}.
#' @param d_sigma2 A length-one numeric vector representing the pooled (across
#'    trials) of treatment effects prior variance.
#' @param beta_sigma A length-one numeric vector representing the prior variance
#'   of the Cholesky factorization factors for the between-study variance.
#' @param nGamma A length-one numeric vector representing the number of Gaussian
#'   copula correlation matrices.
#' @param eta_prior A length-one numeric vector representing the Gaussian
#'   copula correlation matrices prior parameter.
#' @param prior A named list of prior hyperparameters.
#' @return A list with the prior hyperparameters as components.
#' @author Sergio Venturini \email{sergio.venturini@@unibocconi.it}
#' @seealso \code{\link{netcopula}()}
#' @keywords NMA
#' @export
nc_prior <- function(mu_sigma2 = 10^3, d_sigma2 = 10^3, beta_sigma = 10^(-1/8), nGamma = 1,
  eta_prior = 1) {
  prior <- list()
  for (arg in names(formals(sys.function())))
    prior[[arg]] <- get(arg)
  prior
}

#' @rdname nc_prior
#' @export
prior_nc <- nc_prior


#' @rdname nc_prior
#' @export
check_prior <- function(prior) {
  prior_ok <- TRUE

  if (!is.list(prior)) {
    prior_ok <- FALSE
    return(prior_ok)
  }

  if (any(prior[["mu_sigma2"]] <= 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["d_sigma2"]] <= 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["beta_sigma"]] <= 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["nGamma"]] != 1)) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["eta_prior"]] <= 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }

  return(prior_ok)
}
