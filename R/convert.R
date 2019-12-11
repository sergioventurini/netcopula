#' Conversion of an \code{nc_mcmc} object to an object of class \code{mcmc}.
#' 
#' \code{nc_mcmc_to_mcmc} converts an object of class \code{nc_mcmc}
#'   to an object with class \code{\link{mcmc}} from the \code{coda}
#'   package.
#' 
#' @param res An object of type \code{nc_mcmc}.
#' @param include.burnin A logical scalar. If \code{TRUE} the burnin
#'   iterations (if available) are not removed.
#' @param verbose A logical scalar. If \code{TRUE} prints additional
#'   warnings during the conversion.
#' @return An object of type \code{mcmc}.
#' @seealso
#'   \code{\link{netcopula}()} for for fitting a copula-based multivariate
#'     network meta-analysis model;
#'   \code{\link{nc_mcmc-class}};
#'   \code{\link[coda]{mcmc}}.
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#' @export
nc_mcmc_to_mcmc <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  n_study <- res@dim[["n_study"]]
  n_outcomes <- res@dim[["n_outcomes"]]
  n_datapoints <- res@dim[["n_datapoints"]]
  n_treatments <- res@dim[["n_treatments"]]
  ref_trt <- res@dim[["ref_trt"]]
  nGamma <- res@prior[["nGamma"]]

  if (store.burnin) {
    if (include.burnin) {
      tokeep <- seq(1, totiter, by = thin)
      tokeep <- 1:length(tokeep)
    } else {
      todrop <- seq(1, burnin, by = thin)
      tokeep <- seq(1, totiter, by = thin)
      tokeep <- (length(todrop) + 1):length(tokeep)
    }
  } else {
    if (verbose && include.burnin)
      cat("warning: burnin iterations not shown because the 'store.burnin' option was set to FALSE.")
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

  theta <- matrix(NA, nrow = length(tokeep), ncol = ((n_study + n_datapoints + n_treatments + (n_outcomes + 1)/2 +
      2*nGamma*(n_outcomes - 1)/2 + 2*n_datapoints)*n_outcomes + 3))
  theta_nm <- character(ncol(theta))
  count_z <- 1
  for (m in 1:n_outcomes) {
    for (i in 1:n_study) {
      theta[, count_z] <- res@mu[i, m, tokeep]
      theta_nm[count_z] <- paste0("mu[", i, ", ", m, "]")
      count_z <- count_z + 1
    }
  }
  for (m in 1:n_outcomes) {
    for (i in 1:n_datapoints) {
      theta[, count_z] <- res@delta[i, m, tokeep]
      theta_nm[count_z] <- paste0("delta[", i, ", ", m, "]")
      count_z <- count_z + 1
    }
  }
  for (m in 1:n_outcomes) {
    for (i in 1:n_treatments) {
      theta[, count_z] <- res@d[i, m, tokeep]
      theta_nm[count_z] <- paste0("d[", i, ", ", m, "]")
      count_z <- count_z + 1
    }
  }
  for (j in 1:(n_outcomes*(n_outcomes + 1)/2)) {
    theta[, count_z] <- res@Sigma[tokeep, j]
    theta_nm[count_z] <- paste0("Sigma[", j, "]")
    count_z <- count_z + 1
  }
  for (g in 1:nGamma) {
    for (j in 1:(n_outcomes*(n_outcomes - 1)/2)) {
      theta[, count_z] <- res@Gamma[g, j, tokeep]
      theta_nm[count_z] <- paste0("Gamma[", g, ", ", j, "]")
      count_z <- count_z + 1
    }
  }
  for (g in 1:nGamma) {
    for (j in 1:(n_outcomes*(n_outcomes - 1)/2)) {
      theta[, count_z] <- res@r[g, j, tokeep]
      theta_nm[count_z] <- paste0("r[", g, ", ", j, "]")
      count_z <- count_z + 1
    }
  }
  for (m in 1:n_outcomes) {
    for (i in 1:n_datapoints) {
      theta[, count_z] <- res@x[i, m, tokeep]
      theta_nm[count_z] <- paste0("x[", i, ", ", m, "]")
      count_z <- count_z + 1
    }
  }
  for (m in 1:n_outcomes) {
    for (i in 1:n_datapoints) {
      theta[, count_z] <- res@x_unadj[i, m, tokeep]
      theta_nm[count_z] <- paste0("x_unadj[", i, ", ", m, "]")
      count_z <- count_z + 1
    }
  }

  theta[, (n_study + n_datapoints + n_treatments + (n_outcomes + 1)/2 +
      2*nGamma*(n_outcomes - 1)/2 + 2*n_datapoints)*n_outcomes + 1] <- res@dens$loglik[tokeep]
  theta_nm[(n_study + n_datapoints + n_treatments + (n_outcomes + 1)/2 +
      2*nGamma*(n_outcomes - 1)/2 + 2*n_datapoints)*n_outcomes + 1] <- paste0("loglik")
  theta[, (n_study + n_datapoints + n_treatments + (n_outcomes + 1)/2 +
      2*nGamma*(n_outcomes - 1)/2 + 2*n_datapoints)*n_outcomes + 2] <- res@dens$logprior[tokeep]
  theta_nm[(n_study + n_datapoints + n_treatments + (n_outcomes + 1)/2 +
      2*nGamma*(n_outcomes - 1)/2 + 2*n_datapoints)*n_outcomes + 2] <- paste0("logprior")
  theta[, (n_study + n_datapoints + n_treatments + (n_outcomes + 1)/2 +
      2*nGamma*(n_outcomes - 1)/2 + 2*n_datapoints)*n_outcomes + 3] <- res@dens$logpost[tokeep]
  theta_nm[(n_study + n_datapoints + n_treatments + (n_outcomes + 1)/2 +
      2*nGamma*(n_outcomes - 1)/2 + 2*n_datapoints)*n_outcomes + 3] <- paste0("logpost")
  
  colnames(theta) <- theta_nm

  for (j in ncol(theta):1) {
    if (all(is.na(theta[, j]))) {
      theta <- theta[, -j]
    }
  }

  if (store.burnin) {
    if (include.burnin) {
      out <- coda::mcmc(theta, start = 1, end = totiter, thin = thin)
    } else {
      out <- coda::mcmc(theta, start = (burnin + 1), end = totiter, thin = thin)
    }
  } else {
    out <- coda::mcmc(theta, start = (burnin + 1), end = totiter, thin = thin)
  }

  return(out)
}

#' Conversion of an \code{nc_mcmc} object to an object of class \code{list}.
#' 
#' \code{nc_mcmc_to_mcmc} converts an object of class \code{nc_mcmc}
#'   to an object with class \code{list}.
#' 
#' @param res An object of type \code{nc_mcmc}.
#' @param include.burnin A logical scalar. If \code{TRUE} the burnin
#'   iterations (if available) are not removed.
#' @param verbose A logical scalar. If \code{TRUE} prints additional
#'   warnings during the conversion.
#' @return An object of type \code{mcmc}.
#' @seealso
#'   \code{\link{netcopula}()} for for fitting a copula-based multivariate
#'     network meta-analysis model;
#'   \code{\link{nc_mcmc-class}}.
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#' @export
nc_mcmc_to_list <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  n_study <- res@dim[["n_study"]]
  n_outcomes <- res@dim[["n_outcomes"]]
  n_datapoints <- res@dim[["n_datapoints"]]
  n_treatments <- res@dim[["n_treatments"]]
  ref_trt <- res@dim[["ref_trt"]]
  nGamma <- res@prior[["nGamma"]]

  if (store.burnin) {
    if (include.burnin) {
      tokeep <- seq(1, totiter, by = thin)
      tokeep <- 1:length(tokeep)
    } else {
      todrop <- seq(1, burnin, by = thin)
      tokeep <- seq(1, totiter, by = thin)
      tokeep <- (length(todrop) + 1):length(tokeep)
    }
  } else {
    if (verbose && include.burnin)
      cat("warning: burnin iterations not shown because the 'store.burnin' option was set to FALSE.")
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

  theta <- matrix(NA, nrow = length(tokeep), ncol = ((n_study + n_datapoints + n_treatments + (n_outcomes + 1)/2 +
      2*nGamma*(n_outcomes - 1)/2 + 2*n_datapoints)*n_outcomes + 3))
  theta_nm <- character(ncol(theta))
  count_z <- 1
  for (m in 1:n_outcomes) {
    for (i in 1:n_study) {
      theta[, count_z] <- res@mu[i, m, tokeep]
      theta_nm[count_z] <- paste0("mu[", i, ", ", m, "]")
      count_z <- count_z + 1
    }
  }
  for (m in 1:n_outcomes) {
    for (i in 1:n_datapoints) {
      theta[, count_z] <- res@delta[i, m, tokeep]
      theta_nm[count_z] <- paste0("delta[", i, ", ", m, "]")
      count_z <- count_z + 1
    }
  }
  for (m in 1:n_outcomes) {
    for (i in 1:n_treatments) {
      theta[, count_z] <- res@d[i, m, tokeep]
      theta_nm[count_z] <- paste0("d[", i, ", ", m, "]")
      count_z <- count_z + 1
    }
  }
  for (j in 1:(n_outcomes*(n_outcomes + 1)/2)) {
    theta[, count_z] <- res@Sigma[tokeep, j]
    theta_nm[count_z] <- paste0("Sigma[", j, "]")
    count_z <- count_z + 1
  }
  for (g in 1:nGamma) {
    for (j in 1:(n_outcomes*(n_outcomes - 1)/2)) {
      theta[, count_z] <- res@Gamma[g, j, tokeep]
      theta_nm[count_z] <- paste0("Gamma[", g, ", ", j, "]")
      count_z <- count_z + 1
    }
  }
  for (g in 1:nGamma) {
    for (j in 1:(n_outcomes*(n_outcomes - 1)/2)) {
      theta[, count_z] <- res@r[g, j, tokeep]
      theta_nm[count_z] <- paste0("r[", g, ", ", j, "]")
      count_z <- count_z + 1
    }
  }
  for (m in 1:n_outcomes) {
    for (i in 1:n_datapoints) {
      theta[, count_z] <- res@x[i, m, tokeep]
      theta_nm[count_z] <- paste0("x[", i, ", ", m, "]")
      count_z <- count_z + 1
    }
  }
  for (m in 1:n_outcomes) {
    for (i in 1:n_datapoints) {
      theta[, count_z] <- res@x_unadj[i, m, tokeep]
      theta_nm[count_z] <- paste0("x_unadj[", i, ", ", m, "]")
      count_z <- count_z + 1
    }
  }

  theta[, (n_study + n_datapoints + n_treatments + (n_outcomes + 1)/2 +
      2*nGamma*(n_outcomes - 1)/2 + 2*n_datapoints)*n_outcomes + 1] <- res@dens$loglik[tokeep]
  theta_nm[(n_study + n_datapoints + n_treatments + (n_outcomes + 1)/2 +
      2*nGamma*(n_outcomes - 1)/2 + 2*n_datapoints)*n_outcomes + 1] <- paste0("loglik")
  theta[, (n_study + n_datapoints + n_treatments + (n_outcomes + 1)/2 +
      2*nGamma*(n_outcomes - 1)/2 + 2*n_datapoints)*n_outcomes + 2] <- res@dens$logprior[tokeep]
  theta_nm[(n_study + n_datapoints + n_treatments + (n_outcomes + 1)/2 +
      2*nGamma*(n_outcomes - 1)/2 + 2*n_datapoints)*n_outcomes + 2] <- paste0("logprior")
  theta[, (n_study + n_datapoints + n_treatments + (n_outcomes + 1)/2 +
      2*nGamma*(n_outcomes - 1)/2 + 2*n_datapoints)*n_outcomes + 3] <- res@dens$logpost[tokeep]
  theta_nm[(n_study + n_datapoints + n_treatments + (n_outcomes + 1)/2 +
      2*nGamma*(n_outcomes - 1)/2 + 2*n_datapoints)*n_outcomes + 3] <- paste0("logpost")
  
  colnames(theta) <- theta_nm

  for (j in ncol(theta):1) {
    if (all(is.na(theta[, j]))) {
      theta <- theta[, -j]
    }
  }

  return(list(theta))
}
