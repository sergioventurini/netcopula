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

  control <- nc_control()
  control[["burnin"]] <- burnin
  control[["nsim"]] <- nsim
  control[["thin"]] <- nthin
  control[["verbose"]] <- verbose
  control[["Gamma.update"]] <- prm_prop$Gamma.update
  control[["eta.prop"]] <- prm_prop$eta.prop
  control[["sigma.r.prop"]] <- prm_prop$sigma.r.prop

  data_tmp <- data@study_data
  for (m in 1:M) {
    data_m <- data_tmp[, c(2:6, 6 + m, 9 + m)]
    colnames(data_m)[c(6, 7)] <- c("n1", "y1")
    nc_data <- new("nc_data", "binary", n_study, 1, n_datapoints, n_treatments, NULL, data_m, ref_trt)
    
    # MCMC simulation
    res_uni[[m]] <- netcopula_uni(nc_data, control, prm_prior, prm_init, tuning = tuning, adaptation = adaptation,
      init = NULL)

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

#' Function to compute the starting values before fitting a models with the
#'   \code{netcopula} package.
#'
#' \code{nc_init()} is the function that computes the initialization values for
#'   fitting a model with the \code{netcopula} package.
#'
#' @param data An object of class \code{link{nc_data-class}} containing the data
#'   to analyze.
#' @param prm_init An object of class \code{list} containing the parameters to
#'   use for initialization.
#' @param random_start A length-one logical vector. If \code{TRUE} the starting
#'   values are drawn randomly, otherwise.
#' @param init Either a length-one character vector or a named \code{list}
#'   \code{list} providing user-defined starting values.
#' @param prm_prior An object of class \code{list} containing the prior
#'   hyperparameters.
#' @param prm_prop An object of class \code{list} containing the proposal
#'   parameters.
#' @param tuning An object of class \code{list} containing tuning parameters
#'   used in the MCMC simulation.
#' @param adaptation An object of class \code{list} containing the parameters
#'   used for adaptation of the MCMC simulation.
#' @return A named \code{list} with the following items:
#'   \describe{
#'     \item{\code{mu}: }{array of latent coordinates starting values}
#'     \item{\code{delta}: }{numeric vector of initial cluster memberships}
#'     \item{\code{d}: }{numeric vector of initial cluster sizes}
#'     \item{\code{Sigma_M}: }{numeric vector of alpha starting values}
#'     \item{\code{Gamma}: }{numeric vector of eta starting values}
#'     \item{\code{D}: }{numeric vector of sigma2 starting values}
#'     \item{\code{x}: }{numeric vector of lambda starting values}
#'   }
#' @author Sergio Venturini \email{sergio.venturini@@unibocconi.it}
#' @seealso \code{\link{netcopula}()} for fitting a model with the
#'   \code{netcopula} package.
#' @export
nc_init <- function(data, prm_init, random_start, init = NULL, prm_prior = NULL, prm_prop = NULL,
  tuning = NULL, adaptation = NULL) {
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

  return(list(mu = mu, delta = delta, d = d, Sigma_M = Sigma_M, Gamma = Gamma, D = D, x = x))
}
