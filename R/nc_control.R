#' Auxiliary Function for Controlling Model Fitting in the \code{netcopula } Package
#' 
#' @description{
#' \code{nc_control()} is an auxiliary function as user interface for
#'   \code{netcopula()} fitting. Typically only used when calling the \code{netcopula()}
#'   function. It is used to set parameters that affect the sampling but do
#'   not affect the posterior distribution.
#' 
#' \code{control_nc()} is an alias for \code{nc_control()}.
#' 
#' \code{check_control()} is an auxiliary function that verifies the
#'   correctness of the controls provided before a netcopula is fitted with
#'   \code{\link{netcopula}()}.
#' }
#'
#' @param nsim A length-one numeric vector for the number of draws to be taken
#'   from the posterior distribution.
#' @param burnin A length-one numeric vector for the number of initial MCMC
#'   iterations (usually to be discarded).
#' @param thin A length-one numeric vector for the number of iterations between
#'   consecutive draws.
#' @param Gamma.update A length-one character vector providing the algorithm to
#'   use for updating the Gaussian copula correlation matrices.
#' @param eta.prop A length-one numeric vector providing the parameter of the
#'   proposal distribution for the Gaussian copula correlation matrices.
#' @param sigma.r.prop A length-one numeric vector providing the parameter of the
#'   proposal distribution for the Gaussian copula correlation matrices (used only
#'   when Gamma.update is set to 'DanaherSmith').
#' @param random.start A length-one logical vector. If \code{TRUE} the starting
#'   values are drawn randomly, otherwise.
#' @param store.burnin A logical scalar. If \code{TRUE}, the samples from the
#'   burnin are also stored and returned.
#' @param verbose A logical scalar. If \code{TRUE}, causes information to be
#'   printed out about the progress of the fitting.
#' @param control A list of control options.
#'
#' @return A named list with the control options as components.
#' @export
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @seealso \code{\link{netcopula}()}
#' @keywords NMA
nc_control <- function(nsim = 50000,
                       burnin = 100000,
                       thin = 1,
                       Gamma.update = "PX-RPMH",
                       eta.prop = 0.5,
                       sigma.r.prop = 0.01,
                       random.start = TRUE,
                       store.burnin = TRUE,
                       verbose = FALSE){
  control <- list()
  for (arg in names(formals(sys.function())))
    control[[arg]] <- get(arg)
  control
}

#' @rdname nc_control
#' @export
control_nc <- nc_control

#' @rdname nc_control
#' @export
check_control <- function(control) {
  control_ok <- TRUE

  if (!is.list(control)) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["nsim"]] < 1) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["burnin"]] < 1) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["thin"]] < 1) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["Gamma.update"]] != "PX-RPMH" && control[["Gamma.update"]] != "IMH" && 
    control[["Gamma.update"]] != "DanaherSmith" && control[["Gamma.update"]] != "JointR") {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["eta.prop"]] < 0) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["sigma.r.prop"]] < 0) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!is.logical(control[["random.start"]])) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!is.logical(control[["store.burnin"]])) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!is.logical(control[["verbose"]])) {
    control_ok <- FALSE
    return(control_ok)
  }

  return(control_ok)
}
