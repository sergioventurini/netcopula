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
#' betabinexch <- function(theta, varargs) {
#'   eta <- exp(theta[1])/(1 + exp(theta[1]))
#'   K <- exp(theta[2])
#'   data <- varargs[[1]]
#'   y <- data[, 1]
#'   n <- data[, 2]
#'   N <- length(y)
#'   logf <- function(y, n, K, eta) lbeta(K * eta + y, K * (1 - eta) + n - y) -
#'   lbeta(K * eta, K * (1 - eta))
#'   val <- sum(logf(y, n, K, eta))
#'   val <- val + theta[2] - 2 * log(1 + exp(theta[2]))
#'   return(val)
#' }
#' data(cancermortality, package = "LearnBayes")
#' laplace(betabinexch, c(-7, 6), list(cancermortality))
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
laplace <- function(logpost, mode, ...) {
    options(warn = -1)
    
    fit <- optim(mode, logpost, gr = NULL, ..., hessian = TRUE, control = list(fnscale = -1))
    
    options(warn = 0)

    mode <- fit$par
    h <- -solve(fit$hessian)
    p <- length(mode)
    int <- p/2 * log(2 * pi) + 0.5 * log(det(h)) + logpost(mode, ...)
    
    out <- list(mode = mode, var = h, int = int, converge = fit$convergence ==  0)
    
    return(out)
}

#' @export
betabinexch <- function(theta, varargs) {
  eta <- exp(theta[1])/(1 + exp(theta[1]))
  K <- exp(theta[2])
  data <- varargs[[1]]
  y <- data[, 1]
  n <- data[, 2]
  N <- length(y)
  logf <- function(y, n, K, eta) lbeta(K * eta + y, K * (1 - eta) + n - y) -
  lbeta(K * eta, K * (1 - eta))
  val <- sum(logf(y, n, K, eta))
  val <- val + theta[2] - 2 * log(1 + exp(theta[2]))
  return(val)
}

#' @export
betabinT <- function(theta, datapar) {
  tpar <- datapar$par
  d <- betabinexch(theta, datapar) - dmvt_arma(t(as.matrix(theta)), tpar$m, tpar$var, tpar$df, TRUE)
  return(d)
}

#' @export
groupeddatapost <- function(theta, varargs) {
  dj <- function(f, int.lo, int.hi, mu, sigma) {
    f * log(pnorm(int.hi, mu, sigma) - pnorm(int.lo, mu, sigma))
  }
  data <- varargs[[1]]
  mu <- theta[1]
  sigma <- exp(theta[2])
  sum(dj(data$f, data$int.lo, data$int.hi, mu, sigma))
}

#' @export
transplantpost <- function(theta, varargs) {
  data <- varargs[[1]]
  x <- data[, 1]
  y <- data[, 3]
  t <- data[, 2]
  d <- data[, 4]
  tau <- exp(theta[1])
  lambda <- exp(theta[2])
  p <- exp(theta[3])
  xnt <- x[t == 0]
  dnt <- d[t == 0]
  z <- x[t == 1]
  y <- y[t == 1]
  dt <- d[t == 1]
  logf <- function(xnt, dnt, lambda, p) {
    (dnt == 0) * (p * log(lambda) + log(p) - (p + 1) * log(lambda + xnt)) + (dnt == 1) * p * log(lambda/(lambda + xnt))
  }
  logg <- function(z, y, tau, lambda, p) {
    (dt == 0) * (p * log(lambda) + log(p * tau) - (p + 1) * log(lambda + y + tau * z)) + (dt == 1) * p * log(lambda/(lambda + y + tau * z))
  }
  val <- sum(logf(xnt, dnt, lambda, p)) + sum(logg(z, y, tau, lambda, p))
  val <- val + theta[1] + theta[2] + theta[3]
  return(val)
}

#' @export
weibullregpost <- function(theta, varargs) {
  logf <- function(t, c, x, sigma, mu, beta) {
    z = (log(t) - mu - x %*% beta)/sigma
    f = 1/sigma * exp(z - exp(z))
    S = exp(-exp(z))
    c * log(f) + (1 - c) * log(S)
  }
  data <- varargs[[1]]
  k <- dim(data)[2]
  p <- k - 2
  t <- data[, 1]
  c <- data[, 2]
  X <- data[, 3:k]
  sigma <- exp(theta[1])
  mu <- theta[2]
  beta <- array(theta[3:k], c(p, 1))
  return(sum(logf(t, c, X, sigma, mu, beta)))
}

#' @export
my_plot <- function(logf, limits, data, npoints = 50, type = "contour", ...) {
    if (type == "contour" | type == "persp") {
      logf_tmp <- function(theta, data) {
          if (is.matrix(theta) == TRUE) {
              val <- matrix(0, c(dim(theta)[1], 1))
              for (j in 1:dim(theta)[1]) {
                val[j] <- logf(theta[j, ], data)
              }
          }
          else val = logf(theta, data)
          return(val)
      }
      x0 <- seq(limits[1], limits[2], len = npoints)
      y0 <- seq(limits[3], limits[4], len = npoints)
      X <- outer(x0, rep(1, npoints))
      Y <- outer(rep(1, npoints), y0)
      n2 <- npoints^2
      Z <- logf_tmp(cbind(X[1:n2], Y[1:n2]), data)
      Z <- matrix(Z, c(npoints, npoints))
      if (type == "contour") {
        contour(x0, y0, Z, lwd = 2, ...)
      } else if (type == "persp") {
        persp(x0, y0, Z, ...)
      }
    } else if (type == "scatter") {
      logf_tmp <- function(theta, data) {
          val <- numeric(length(theta))
          for (i in 1:length(theta)) {
            val[i] <- logf(theta[i], data)
          }
          return(val)
      }
      x0 <- seq(limits[1], limits[2], len = npoints)
      Y <- logf_tmp(x0, data)
      plot(x0, Y, type = "l", ...)
    }
}

#' @export
build_output <- function(a) {
  nparams <- if (length(dim(a)) < 3) 1 else dim(a)[length(dim(a))]
  output <- matrix(NA, ncol = 7, nrow = nparams)
  if (length(dim(a)) == 2) {
    a <- array(a, c(dim(a), 1))
  }
  for (i in 1:nparams) {
    ai <- a[, , i, drop = FALSE]
    quantiles <- quantile(as.vector(ai), probs = c(.025, .25, .5, .75, .975))
    output[i, ] <- c(mean(ai), sd(as.vector(ai)), quantiles)
  }
  dimnames(output) <- list(dimnames(a)[[3]], c("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%"))
  return(output)
}

#' @export
sigma2corr <- function(sigma) {
  nn <- length(sigma)
  M <- (-1 + sqrt(1 + 8*nn))/2
  sigma_mat <- matrix(NA, nrow = M, ncol = M)
  sigma_mat[lower.tri(sigma_mat, diag = TRUE)] <- sigma
  sigma_mat[upper.tri(sigma_mat)] <- sigma_mat[lower.tri(sigma_mat)]
  corr <- cov2cor(sigma_mat)

  return(corr[lower.tri(corr)])
}

#' @export
add_corr <- function(res_mcmc) {
  Sigma_M <- res_mcmc[, grep("Sigma_M", dimnames(res_mcmc)[[2]])]
  Corr_M <- t(apply(Sigma_M, 1, sigma2corr))
  if (nrow(Corr_M) == 1) {
    Corr_M <- t(Corr_M)
  }
  nn <- dim(Sigma_M)[2]
  M <- (-1 + sqrt(1 + 8*nn))/2
  nm <- gsub("Sigma", "Corr", dimnames(Sigma_M)[[2]][-cumsum(1:M)])
  dimnames(Corr_M) <- list(NULL, nm)
  new_res_mcmc <- as.mcmc(cbind(res_mcmc, Corr_M))

  return(new_res_mcmc)
}

#' @export
col.alpha <- function(acol, alpha = 0.2) {
    acol <- col2rgb(acol)
    acol <- rgb(acol[1]/255, acol[2]/255, acol[3]/255, alpha)
    acol
}

#' @export
plot_est <- function(res_sum, true_values, param) {
  idx <- grep(paste0(param, "["), dimnames(res_sum$quantiles)[[1]], fixed = TRUE)
  plot(res_sum$statistics[idx, 1], true_values, ylim = c(min(as.numeric(res_sum$quantiles[idx, c(1, 5)]), true_values), max(as.numeric(res_sum$quantiles[idx, c(1, 5)]), true_values)), xlab = "posterior mean estimate", ylab = "true values", main = param, type = "n")
  for (i in idx) {
    segments(res_sum$statistics[i, 1], res_sum$quantiles[i, 1], res_sum$statistics[i, 1], res_sum$quantiles[i, 5], col = "gray", lwd = 1)
  }
  abline(h = 0, lty = 2, lwd = .2, col = 1)
  abline(0, 1, lty = 2, lwd = .2, col = 1)
  points(res_sum$statistics[idx, 1], res_sum$statistics[idx, 1], pch = 20, col = "lightblue")
  points(res_sum$statistics[idx, 1], true_values, pch = 20, col = "orange")
}
