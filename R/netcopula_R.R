#' @export
netcopula_R <- function(data, burnin = 10000, nsim = 5000, nthin = 1, prm_prior, prm_prop = NULL, prm_init, random_start = TRUE, init = NULL, tuning, adaptation, verbose = FALSE, mu_delta_mh = TRUE) {
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
  mu <- nc_start$mu
  delta <- nc_start$delta
  d <- nc_start$d
  Sigma <- nc_start$Sigma
  Gamma <- nc_start$Gamma

  if (verbose) cat("done!\n")

  # start iteration
  if (verbose) cat("Running the MCMC simulation...\n")

  if (mu_delta_mh) {
    res <- nc_mcmc_mh_R(
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

#' @export
nc_mcmc_mh_R <- function(data, init, totiter, prior, prop, tuning, adapt, verbose) {
  niter <- 1
  adapt_k <- 0
  n_datapoints <- data@n_datapoints
  M <- data@n_outcomes
  n_trt <- data@n_treatments
  n_study <- data@n_study
  nn <- M*(M + 1)/2

  eps <- tuning$eps
  eps_ab <- tuning$eps_ab

  print_every <- 500
  every <- adapt$every
  maxiter <- adapt$maxiter

  alpha <- adapt$alpha
  beta <- adapt$beta
  gamma <- adapt$gamma
  tar <- adapt$tar
  tol <- adapt$tol

  study_data <- data@study_data
  study_id <- data@study_id
  cols_y <- cols_n <- numeric(M)
  for (i in 1:M) {
    cols_y[i] <- 6 + M + (i - 1)
    cols_n[i] <- 6 + (i - 1)
  }
  Gamma <- init$Gamma
  trt <- study_data$trt
  ref_trt <- data@ref_trt
  narms <- study_id$narms
  narms_study <- study_data$narms
  study_index <- study_data$study

  mu <- init$mu
  delta <- init$delta
  x <- init$x
  Sigma_M <- init$Sigma_M
  d <- init$d
  delta_arma <- delta
  delta_arma_prop <- delta
  y <- df_nm(study_data, cols_y)
  n_data <- df_nm(study_data, cols_n)
  n_imp <- n_imputed(n_data)
  x_imp <- y_imp <- x_imp_arma <- x_star_arma <- matrix(0, n_datapoints, M)

  D_list <- init$D
  D <- D_inv <- S_q <- D_prop <- matrix(0, M, M)
  D_tmp <- numeric(M)
  x2t_split <- vector(mode = "list", length = n_trt)

  nGamma <- prior$nGamma
  Gamma_update <- prop$Gamma_update
  eta_prior <- prior$eta_prior
  eta_prop <- prop$eta_prop
  Gamma_q_prop <- D_prop_inv <- Gamma_q_curr <- Sigma_q_prop <- matrix(0, M, M)
  accept_Gamma_q <- numeric(nGamma)
  D_chain <- array(0, dim = c(nGamma, M, totiter))
  S_q_chain <- Sigma_q_prop_chain <- array(0, dim = c(nGamma, nn, totiter))
  Gamma_chain <- array(0, dim = c(nGamma, M*(M - 1)/2, totiter))
  ran_unif <- Gamma_rate <- mu_rate <- delta_rate <- d_rate <- Sigma_M_rate <- Gamma_A <- Gamma_B <- target_Gamma_prop <- target_Gamma_curr <- 0

  Gamma_k <- Gamma_k_m <- matrix(0, M, M)
  Gamma_k_m_m <- matrix(0, 1, M)
  a_ikm <- b_ikm <- w_ikm <- gamma_ikm <- theta_ikm <- p_ikm <- 0
  a_chain <- b_chain <- x_chain <- array(0, dim = c(n_datapoints, M, totiter))

  mu_sigma <- sqrt(prior$mu_sigma2)

  mode_mu <- mode_delta <- numeric(1)
  var_mu <- var_delta <- matrix(0, 1, 1)
  Sigma_M_m_m <- delta_ik_m <- d_1 <- d_k <- d_1k <- d_1k_m <- matrix(0, 1, M)
  Sigma_M_m <- matrix(0, M, M)

  baseline <- study_data$baseline
  baseline_id <- study_id$baseline
  mu_long <- param_long(mu, narms, FALSE)

  accept_mu <- 0
  rho_mu <- matrix(2.38/sqrt(M*n_study), n_study, M)
  cov_mu <- matrix(1, n_study, M)
  cov_mu_mat <- matrix(0, 1, 1)
  mu_mu_vec <- numeric(1)
  mu_mu <- matrix(0, n_study, M)
  theta_mu <- numeric(every)
  ar_mu_vec <- array(0, dim = c(n_study, M, totiter))
  ar_mu <- array(0, dim = c(2, n_study, M))
  ar_mu_sub <- numeric(2)
  prop_mu <- vector(mode = "list", length = 4)
  cov_mu_chain <- rho_mu_chain <- mu_prop_chain <- mu_rate_chain <- tp_mu_chain <- mu_chain <- array(0, dim = c(n_study, M, totiter))
  accept_delta <- 0
  rho_delta <- matrix(2.38/sqrt(M*n_datapoints), n_datapoints, M)
  cov_delta <- matrix(1, n_datapoints, M)
  cov_delta_mat <- matrix(0, 1, 1)
  mu_delta_vec <- numeric(1)
  mu_delta <- matrix(0, n_datapoints, M)
  theta_delta <- numeric(every)
  ar_delta_vec <- array(0, dim = c(n_datapoints, M, totiter))
  ar_delta <- array(0, dim = c(2, n_datapoints, M))
  ar_delta_sub <- numeric(2)
  prop_delta <- vector(mode = "list", length = 4)
  delta_chain <- delta_prop_chain <- cov_delta_chain <- rho_delta_chain <- array(0, dim = c(n_datapoints, M, totiter))

  d_chain <- array(0, dim = c(n_trt, M, totiter))

  Sigma_M_chain <- beta_chain <- matrix(0, totiter, nn)

  accept_d <- 0
  rho_d <- 2.38/sqrt(M*(n_trt - 1))
  sd_multplier <- 1
  cov_d <- diag(M*(n_trt - 1))
  d_curr <- numeric(M*(n_trt - 1))
  d_prop <- matrix(0, 1, M*(n_trt - 1))
  d_curr_ref <- d_prop_ref <- matrix(0, n_trt, M)
  d_sigma <- sqrt(prior$d_sigma2)
  mu_d <- numeric(M*(n_trt - 1))
  theta_d <- matrix(0, n_trt, M, every)
  theta_d_reshaped <- matrix(0, every, M*(n_trt - 1))
  ar_d_vec <- numeric(totiter)
  ar_d <- numeric(2)
  prop_d <- vector(mode = "list", length = 4)

  accept_Sigma_M <- 0
  rho_beta <- 2.38/sqrt(nn)
  cov_beta <- diag(nn)
  Sigma_M_curr <- Sigma_M_prop <- numeric(M)
  beta_curr <- beta_prop <- numeric(nn)
  sigma_r <- prior$beta_sigma
  mu_beta <- numeric(nn)
  theta_beta <- matrix(0, every, nn)
  ar_Sigma_M_vec <- numeric(totiter)
  ar_Sigma_M <- numeric(2)
  prop_beta <- vector(mode = "list", length = 4)

  loglik <- logprior <- logpost <- numeric(totiter)

  tmp <- numeric(every)

  while ((niter - 1) < totiter) {
    # if (niter > 1) {
    #   x <- x_unadj_chain[, , niter - 1]
    # }

    # imputing x and y variables
    # x_imp <- x_imputed(x, Gamma, trt)
    x_imp <- x # CHECK THAT THIS PART IS OK WHEN IMPUTING!!!
    x_imp_arma <- x_imp
    y_imp <- y_imputed(y, x_imp, narms, mu, delta, n_imp)

    # updating mu (study-specific baseline effects) and delta (study-specific
    # [random] treatment effects)
    for (ik in 1:n_datapoints) {
      if (nGamma > 1) {
        Gamma_k <- Gamma[[trt[ik]]]
      } else {
        Gamma_k <- Gamma[[1]]
      }
      d_1 <- d[baseline[ik], ]
      d_k <- d[trt[ik], ]
      d_1k <- d_k - d_1
      if (trt[ik] == baseline[ik]) {
        # update mu (study-specific baseline effect)
        for (m in 1:M) {
          Gamma_k_m_m <- Gamma_k[m, -m]
          Gamma_k_m <- Gamma_k[-m, -m]
          x_ik_m <- x_imp_arma[ik, -m]
          w_ikm <- as.numeric(t(Gamma_k_m_m) %*% solve(Gamma_k_m) %*% x_ik_m)
          gamma_ikm <- sqrt(as.numeric(Gamma_k[m, m] - t(Gamma_k_m_m) %*% solve(Gamma_k_m) %*% Gamma_k_m_m))

          mu_curr <- mu_long[ik, m]
          mode_mu <- mu_curr
          var_mu[1, 1] <- rho_mu[study_index[ik], m]^2*cov_mu[study_index[ik], m]
          mu_prop <- as.numeric(rmvt_arma(1, mode_mu, var_mu, 7))
          mu_prop_chain[study_index[ik], m, niter] <- mu_prop
          rho_mu_chain[study_index[ik], m, niter] <- rho_mu[study_index[ik], m]
          cov_mu_chain[study_index[ik], m, niter] <- cov_mu[study_index[ik], m]

          target_mu_curr <- mu_logpost(mu_curr, delta[ik, m], as.double(y_imp[ik, m]), as.double(n_imp[ik, m]), w_ikm, gamma_ikm, mu_sigma, eps, eps_ab)
          target_mu_prop <- mu_logpost(mu_prop, delta[ik, m], as.double(y_imp[ik, m]), as.double(n_imp[ik, m]), w_ikm, gamma_ikm, mu_sigma, eps, eps_ab)
          tp_mu_chain[study_index[ik], m, niter] <- target_mu_prop

          mu_A <- target_mu_prop - target_mu_curr
          mu_B <- 0
          mu_rate <- exp(mu_A + mu_B)
          mu_rate_chain[study_index[ik], m, niter] <- mu_rate
          ran_unif <- runif(1)
          if (!is.na(target_mu_prop) & !is.na(target_mu_curr)) {
            if (ran_unif < mu_rate) {
              mu_long[ik, m] <- mu_prop
              mu <- param_wide(mu_long, narms, trt, baseline)
              mu_long <- param_long(mu, narms, FALSE) # is this necessary?
              accept_mu <- accept_mu + 1
            }
          } else {
            mu_rate <- 0
          }
          ar_mu_vec[study_index[ik], m, niter] <- min(mu_rate, 1, na.rm = TRUE)
        }
      } else {
        # update delta (study-specific [random] treatment effects)
        for (m in 1:M) {
          Gamma_k_m_m <- Gamma_k[m, -m]
          Gamma_k_m <- Gamma_k[-m, -m]
          x_ik_m <- x_imp_arma[ik, -m]
          w_ikm <- as.numeric(t(Gamma_k_m_m) %*% solve(Gamma_k_m) %*% x_ik_m)
          gamma_ikm <- sqrt(as.numeric(Gamma_k[m, m] - t(Gamma_k_m_m) %*% solve(Gamma_k_m) %*% Gamma_k_m_m))

          Sigma_M_m_m <- Sigma_M[m, -m]
          Sigma_M_m <- Sigma_M[-m, -m]
          delta_ik_m <- delta_arma[ik, -m]
          d_1k_m <- d_1k[-m]
          tau_ikm <- d_1k[m] + as.numeric(t(Sigma_M_m_m) %*% solve(Sigma_M_m) %*% (delta_ik_m - d_1k_m))
          eta_ikm <- sqrt(as.numeric(Sigma_M[m, m] - t(Sigma_M_m_m) %*% solve(Sigma_M_m) %*% Sigma_M_m_m))

          delta_curr <- delta[ik, m]
          mode_delta <- delta_curr
          var_delta[1, 1] <- rho_delta[ik, m]^2*cov_delta[ik, m]
          delta_prop <- as.numeric(rmvt_arma(1, mode_delta, var_delta, 7))
          delta_prop_chain[ik, m, niter] <- delta_prop
          rho_delta_chain[ik, m, niter] <- rho_delta[ik, m]
          cov_delta_chain[ik, m, niter] <- cov_delta[ik, m]

          target_delta_curr <- delta_logpost(delta_curr, mu_long[ik, m], tau_ikm, eta_ikm, as.double(y_imp[ik, m]), as.double(n_imp[ik, m]), w_ikm, gamma_ikm, eps, eps_ab)
          target_delta_prop <- delta_logpost(delta_prop, mu_long[ik, m], tau_ikm, eta_ikm, as.double(y_imp[ik, m]), as.double(n_imp[ik, m]), w_ikm, gamma_ikm, eps, eps_ab)

          delta_A <- target_delta_prop - target_delta_curr
          delta_B <- 0
          delta_rate <- exp(delta_A + delta_B)
          ran_unif <- runif(1)
          if (!is.na(target_delta_prop) & !is.na(target_delta_curr)) {
            if (ran_unif < delta_rate) {
              delta[ik, m] <- delta_prop
              accept_delta <- accept_delta + 1
            }
          } else {
            delta_rate <- 0
          }
          ar_delta_vec[ik, m, niter] <- min(delta_rate, 1, na.rm = TRUE)
        }
      }
    }
    mu_chain[, , niter] <- mu
    delta_chain[, , niter] <- delta
    delta_arma <- delta

    # updating x (latent variables)
    for (ik in 1:n_datapoints) {
      if (nGamma > 1) {
        Gamma_k <- Gamma(trt[ik])
      } else {
        Gamma_k <- Gamma[[1]]
      }
      for (m in 1:M) {
        theta_ikm <- mu_long[ik, m] + delta[ik, m]
        p_ikm <- expit_rcpp(theta_ikm)
        a_ikm <- pbinom(y_imp[ik, m] - 1, n_imp[ik, m], p_ikm, 1, 1) # log
        b_ikm <- pbinom(y_imp[ik, m], n_imp[ik, m], p_ikm, 1, 1) # log
        a_chain[ik, m, niter] <- a_ikm
        b_chain[ik, m, niter] <- b_ikm
        Gamma_k_m_m <- Gamma_k[m, -m]
        Gamma_k_m <- Gamma_k[-m, -m]
        x_ik_m <- x_imp_arma[ik, -m]
        w_ikm <- as.numeric(t(Gamma_k_m_m) %*% solve(Gamma_k_m) %*% x_ik_m)
        gamma_ikm <- sqrt(as.numeric(Gamma_k[m, m] - t(Gamma_k_m_m) %*% solve(Gamma_k_m) %*% Gamma_k_m_m))
        x[ik, m] <- as.numeric(rtruncnorm_rcpp(1, qnorm(a_ikm, 0, 1, 1, 1), qnorm(b_ikm, 0, 1, 1, 1), w_ikm, gamma_ikm))
        x_imp[ik, m] <- x[ik, m]
        x_imp_arma[ik, m] <- x[ik, m]
      }
    }
    x_chain[, , niter] <- x
    # x_unadj_chain[, , niter] <- x

    # updating Gamma (correlation of latent variables)
    if (Gamma_update == "PX-RPMH") {
      # if (niter > 1) {
      #   x <- x_adj_chain[, , niter - 1]
      #   x_imp <- x # CHECK THAT THIS PART IS OK WHEN IMPUTING!!!
      #   x_imp_arma <- x_imp
      # }

      if (nGamma > 1) {
        x2t_split <- split_nm(x_imp, trt)
        for (q in 1:n_trt) {
          x_q <- x2t_split[[q]]
          n_q <- nrow(x_q)
          D <- D_list[[q]]
          # D_inv <- solve(D)
          for (k in 1:n_q) {
            S_q <- S_q + D %*% t(x_q[k, , drop = FALSE]) %*% x_q[k, , drop = FALSE] %*% D
          }
          if (!is_positive_definite(S_q, 1423)) {
            Rprintf("S_q - nGamma > 1\n")
            S_q <- make_positive_definite(S_q)
          }
          S_q_chain[q, , niter] <- diag_tri(S_q)
          # in the following, degrees of freedom are set to (n_q > M) ? n_q : M)
          # otherwise, when n_q < M the IW routine returns an error
          # (indeterminate system)
          Sigma_q_prop <- rinvwish_arma(ifelse(n_q > M, n_q, M), S_q)
          Sigma_q_prop_chain[q, , niter] <- diag_tri(Sigma_q_prop)
          D_prop <- diag(sqrt(diag(Sigma_q_prop)))
          D_prop_inv <- solve(D_prop)
          Gamma_q_prop <- cov2cor(Sigma_q_prop)

          Gamma_q_curr <- Gamma[[q]]
          Gamma_rate <- 0.5*(M + 1)*(log(det(Gamma_q_prop)) - log(det(Gamma_q_curr)))
          ran_unif <- runif(1)
          if (ran_unif < exp(Gamma_rate)) {
            Gamma_q_curr <- Gamma_q_prop
            D <- D_prop
            accept_Gamma_q[q] <- accept_Gamma_q[q] + 1
          }
          Gamma[[q]] <- Gamma_q_curr
          D_list[[q]] <- D

          S_q <- matrix(0, M, M)
          D_chain[q, , niter] <- arma::diagvec(D)
        }
      } else {
        D <- diag(sqrt(1/colSums(x_imp_arma^2)))
        x_star_arma <- x_imp_arma %*% D
        S_q <- crossprod(x_star_arma)
        S_q_chain[1, , niter] <- diag_tri(S_q)
        Sigma_q_prop <- rinvwish_arma(ifelse(n_datapoints > M, n_datapoints, M), S_q)
        Sigma_q_prop_chain[1, , niter] <- diag_tri(Sigma_q_prop)
        D_prop <- diag(sqrt(diag(Sigma_q_prop)))
        Gamma_q_prop <- cov2cor(Sigma_q_prop)

        Gamma_q_curr <- Gamma[[1]]
        Gamma_rate <- 0.5*(M + 1)*(log(det(Gamma_q_prop)) - log(det(Gamma_q_curr)))
        ran_unif <- runif(1)
        if (ran_unif < exp(Gamma_rate)) {
          # adjust latent variables x according to the new scales in 'D_prop'
          D_prop_inv <- solve(D_prop)
          x_imp_arma <- x_imp <- x <- x_star_arma %*% D_prop_inv

          D <- D_prop
          Gamma_q_curr <- Gamma_q_prop
          accept_Gamma_q[1] <- accept_Gamma_q[1] + 1
        }
        Gamma[[1]] <- Gamma_q_curr

        S_q <- matrix(0, M, M)
      }
    } else if (Gamma_update == "IMH") {
      if (nGamma > 1) {
        x2t_split <- split_nm(x_imp, trt)
        for (q in 1:n_trt) {
          x_q <- x2t_split[[q]]
          Gamma_q_curr <- Gamma[[q]]
          Gamma_q_prop <- rlkj_arma(M, eta_prop)

          target_Gamma_curr <- Gamma_logpost(Gamma_q_curr, x_q, eta_prior)
          target_Gamma_prop <- Gamma_logpost(Gamma_q_prop, x_q, eta_prior)
          Gamma_A <- target_Gamma_prop - target_Gamma_curr
          Gamma_B <- dlkj_arma(Gamma_q_curr, eta_prop, TRUE) - dlkj_arma(Gamma_q_prop, eta_prop, TRUE)
          Gamma_rate <- exp(Gamma_A + Gamma_B)
          ran_unif <- runif(1)
          if (ran_unif < exp(Gamma_rate)) {
            Gamma_q_curr <- Gamma_q_prop
            accept_Gamma_q[q] <- accept_Gamma_q[q] + 1
          }
          Gamma[[q]] <- Gamma_q_curr
        }
      } else {
        Gamma_q_curr <- Gamma[[1]]
        Gamma_q_prop <- rlkj_arma(M, eta_prop)

        target_Gamma_curr <- Gamma_logpost(Gamma_q_curr, x_imp_arma, eta_prior)
        target_Gamma_prop <- Gamma_logpost(Gamma_q_prop, x_imp_arma, eta_prior)
        Gamma_A <- target_Gamma_prop - target_Gamma_curr
        Gamma_B <- dlkj_arma(Gamma_q_curr, eta_prop, TRUE) - dlkj_arma(Gamma_q_prop, eta_prop, TRUE)
        Gamma_rate <- exp(Gamma_A + Gamma_B)
        ran_unif <- runif(1)
        if (ran_unif < exp(Gamma_rate)) {
          Gamma_q_curr <- Gamma_q_prop
          accept_Gamma_q[1] <- accept_Gamma_q[1] + 1
        }
        Gamma[[1]] <- Gamma_q_curr
      }
    } else if (Gamma_update == "Talhouketal") {
      if (nGamma > 1) {
        print("STILL TO DO!")
      } else {
        Gamma_tmp <- Gamma[[1]]
        Gamma_inv <- solve(Gamma_tmp)
        for (m in 1:M) {
          D_tmp[m] <- rinvgamma_rcpp(1, (M + 1)/2.0, Gamma_inv[m, m]/2)
        }
        D <- diag(sqrt(D_tmp))
        W <- x_imp_arma %*% D
        S <- t(W) %*% W + diag(M)
        Gamma_new <- rinvwish_arma(2 + n_datapoints, S)
        D <- diag(sqrt(diag(Gamma_new)))
        Gamma[[1]] <- cov2cor(Gamma_new)
      }
    } else {
      stop("the specified update method for the Gamma parameter is not available. use either 'IMH' or 'PX-RPMH'.")
    }
    D_chain[1, , niter] <- diag(D)
    Gamma_chain[, , niter] <- list_mat(Gamma)
    x_chain[, , niter] <- x
    # x_adj_chain[, , niter] <- x

    # # updating d (pooled treatment effects across trials)
    d_curr <- mat_to_vec(d, TRUE, ref_trt)
    d_curr_ref <- d
    d_prop <- rmvt_arma(1, d_curr, (sd_multplier*rho_d)^2*cov_d, 7)
    d_prop_ref <- vec_to_mat(as.numeric(d_prop), M, TRUE, ref_trt)
    target_d_curr <- d_logpost(d_curr_ref, delta_arma, Sigma_M, trt, baseline, narms_study, d_sigma, ref_trt)
    target_d_prop <- d_logpost(d_prop_ref, delta_arma, Sigma_M, trt, baseline, narms_study, d_sigma, ref_trt)
    d_A <- target_d_prop - target_d_curr
    d_B <- 0
    d_rate <- exp(d_A + d_B)
    ran_unif <- runif(1)
    if (ran_unif < d_rate) {
      d <- d_prop_ref
      accept_d <- accept_d + 1
    }
    ar_d_vec[niter] <- min(d_rate, 1, na.rm = TRUE)
    d_chain[, , niter] <- d

    # # updating Sigma_M (common between-study covariance structure)
    Sigma_M_curr <- Sigma_M
    beta_curr <- Sigma_M_to_beta(Sigma_M_curr)
    beta_prop <- as.numeric(rmvt_arma(1, beta_curr, (sd_multplier*rho_beta)^2*cov_beta, 7))
    Sigma_M_prop <- beta_to_Sigma_M(beta_prop, M)
    target_Sigma_M_prop <- Sigma_M_logpost(d, delta_arma, Sigma_M_prop, trt, baseline, narms_study, sigma_r)
    target_Sigma_M_curr <- Sigma_M_logpost(d, delta_arma, Sigma_M_curr, trt, baseline, narms_study, sigma_r)
    Sigma_M_A <- target_Sigma_M_prop - target_Sigma_M_curr
    Sigma_M_B <- 0
    Sigma_M_rate <- exp(Sigma_M_A + Sigma_M_B)
    ran_unif <- runif(1)
    if (ran_unif < Sigma_M_rate) {
      Sigma_M <- Sigma_M_prop
      accept_Sigma_M <- accept_Sigma_M + 1
    }
    ar_Sigma_M_vec[niter] <- min(Sigma_M_rate, 1, na.rm = TRUE)
    Sigma_M_chain[niter, ] <- diag_tri(Sigma_M)
    beta_chain[niter, ] <- Sigma_M_to_beta(Sigma_M)

    # adaptation step
    if (((adapt_k + 1) <= maxiter) & ((niter %% every) == 0)) {
      if (verbose) {
        print(paste0("      --> performing adaptation [step ", adapt_k + 1, "/", maxiter, "]"))
      }

      # mu proposal parameters
      for (s in 1:n_study) {
        for (m in 1:M) {
          if (abs(ar_mu[1, s, m] - tar) > tol) {
            theta_mu <- as.matrix(mu_chain[s, m, (niter - every + 1):niter])
            if (adapt_k == 0) {
              mu_mu[s, m] <- mean(theta_mu)
            }
            mu_mu_vec[1] <- mu_mu[s, m]
            tmp <- ar_mu_vec[s, m, (niter - every + 1):niter]
            ar_mu[2, s, m] <- as.numeric(mean(tmp))
            ar_mu_sub <- ar_mu[1:2, s, m]
            cov_mu_mat[1, 1] <- cov_mu[s, m]
            prop_mu <- rwmh_adapt_R(theta_mu, mu_mu_vec, rho_mu[s, m], cov_mu_mat, ar_mu_sub, alpha, beta, gamma, tar, adapt_k, FALSE, 1)
            rho_mu[s, m] <- prop_mu$rho
            cov_mu_mat <- prop_mu$covariance
            if (!is_positive_definite(cov_mu_mat, 1791)) {
              print("cov_mu_mat")
              cov_mu_mat <- make_positive_definite(cov_mu_mat)
            }
            cov_mu[s, m] <- cov_mu_mat[1, 1]
            mu_mu_vec <- prop_mu$mu
            mu_mu[s, m] <- mu_mu_vec[1]
            ar_mu[1, s, m] <- prop_mu$ar
          }
        }
      }

      # delta proposal parameters
      for (ik in 1:n_datapoints) {
        if (trt[ik] != baseline[ik]) {
          for (m in 1:M) {
            if (abs(ar_delta[1, ik, m] - tar) > tol) {
              theta_delta <- as.matrix(delta_chain[ik, m, (niter - every + 1):niter])
              if (adapt_k == 0) {
                mu_delta[ik, m] <- mean(theta_delta)
              }
              mu_delta_vec[1] <- mu_delta[ik, m]
              tmp <- ar_delta_vec[ik, m, (niter - every + 1):niter]
              ar_delta[2, ik, m] <- as.numeric(mean(tmp))
              ar_delta_sub[1] <- ar_delta[1, ik, m]
              ar_delta_sub[2] <- ar_delta[2, ik, m]
              cov_delta_mat[1, 1] <- cov_delta[ik, m]
              prop_delta <- rwmh_adapt_R(theta_delta, mu_delta_vec, rho_delta[ik, m], cov_delta_mat, ar_delta_sub, alpha, beta, gamma, tar, adapt_k, FALSE, 2)
              rho_delta[ik, m] <- prop_delta$rho
              cov_delta_mat <- prop_delta$covariance
              if (!is_positive_definite(cov_delta_mat, 1820)) {
                print("cov_delta_mat")
                cov_delta_mat <- make_positive_definite(cov_delta_mat)
              }
              cov_delta[ik, m] <- cov_delta_mat[1, 1]
              mu_delta_vec <- prop_delta$mu
              mu_delta[ik, m] <- mu_delta_vec[1]
              ar_delta[1, ik, m] <- prop_delta$ar
            }
          }
        }
      }

      # d proposal parameters
      if (abs(ar_d[1] - tar) > tol) {
        theta_d <- d_chain[, , (niter - every + 1):niter, drop = FALSE]
        theta_d_reshaped <- cube_to_mat(theta_d, TRUE, ref_trt)
        if (adapt_k == 0) {
          mu_d <- as.numeric(colMeans(theta_d_reshaped))
        }
        ar_d[2] <- as.numeric(mean(ar_d_vec[(niter - every + 1):niter]))
        prop_d <- rwmh_adapt_R(theta_d_reshaped, mu_d, rho_d, cov_d, ar_d, alpha, beta, gamma, tar, adapt_k, FALSE, 3)
        rho_d <- prop_d$rho
        cov_d <- prop_d$covariance
        if (!is_positive_definite(cov_d, 1843)) {
          print("cov_d")
          cov_d <- make_positive_definite(cov_d)
        }
        mu_d <- prop_d$mu
        ar_d[1] <- prop_d$ar
      }

      # log Cholesky betas (Sigma_M) proposal parameters
      if (abs(ar_Sigma_M[1] - tar) > tol) {
        theta_beta <- as.matrix(beta_chain[(niter - every + 1):niter, ])
        if (adapt_k == 0) {
          mu_beta <- as.numeric(colMeans(theta_beta))
        }
        ar_Sigma_M[2] <- as.numeric(mean(ar_Sigma_M_vec[(niter - every + 1):niter]))
        prop_beta <- rwmh_adapt_R(theta_beta, mu_beta, rho_beta, cov_beta, ar_Sigma_M, alpha, beta, gamma, tar, adapt_k, FALSE, 4)
        rho_beta <- prop_beta$rho
        cov_beta <- prop_beta$covariance
        if (!is_positive_definite(cov_beta, 1861)) {
          print("cov_beta")
          cov_beta <- make_positive_definite(cov_beta)
        }
        mu_beta <- prop_beta$mu
        ar_Sigma_M[1] <- prop_beta$ar
      }

      adapt_k <- adapt_k + 1
    }

    # calculate the loglikelihood, logprior and logposterior
    # loglik[niter] <- nc_loglik(y_imp, n_imp, x_imp, trt, mu, delta, Gamma)
    # logprior[niter] <- nc_logprior(mu, mu_sigma, d, d_sigma, Sigma_M, sigma_r, ref_trt)

    # print the information
    if (((niter %% print_every) == 0) & verbose) {
      print(paste0("   iter. ", niter, "/", totiter, " ==> d: ", format(accept_d/(as.double(totiter)), digits = 3), " - Gamma: ", format(mean(accept_Gamma_q)/(as.double(totiter)), digits = 3), " - mu: ", format(accept_mu/(as.double(totiter*n_study*M)), digits = 3), " - delta: ", format(accept_delta/(as.double(totiter*(n_datapoints - n_study)*M)), digits = 3), " - Sigma_M: ", format(accept_Sigma_M/(as.double(totiter)), digits = 3)))
    }

    niter <- niter + 1
  }

  accept <- list("d" = accept_d/totiter,
                 "mu" = accept_mu/(totiter*n_study*M),
                 "delta" = accept_delta/(totiter*(n_datapoints - n_study)*M),
                 "Gamma" = accept_Gamma_q/totiter,
                 "Sigma" = accept_Sigma_M/totiter)

  return(list("mu" = mu_chain,
              "delta" = delta_chain,
              "d" = d_chain,
              "Sigma" = Sigma_M_chain,
              "Gamma" = Gamma_chain,
              "x" = x_chain,
              # "x" = x_adj_chain,
              "x_unadj" = x_chain,
              # "x_unadj" = x_unadj_chain,
              # "a" = a_chain,
              # "b" = b_chain,
              # "D" = D_chain,
              # "S_q" = S_q_chain,
              # "Sigma_q_prop" = Sigma_q_prop_chain))
              "loglik" = loglik,
              "logprior" = logprior,
              "logpost" = (loglik + logprior),
              "accept" = accept))
              # "rho_mu" = rho_mu_chain,
              # "cov_mu" = cov_mu_chain,
              # "ar_mu_vec" = ar_mu_vec,
              # "mu_rate" = mu_rate_chain,
              # "mu_prop" = mu_prop_chain,
              # "tp_mu" = tp_mu_chain,
              # "delta_prop" = delta_prop_chain))
}

#' @export
rwmh_adapt_R <- function(theta, mu, rho, cov, ar, alpha, beta, gamma, tar, k, iter_cols, what) {
  ar_new <- (1 - alpha)*ar[1] + alpha*ar[2]
  beta_k <- beta/((k + 1)^gamma)
  rho_new <- rho*exp(beta_k*(qnorm(ar_new/2) - qnorm(tar/2)))
  # if (what == 2 & rho_new > 2) browser()
  if (iter_cols) {
    theta_k <- theta
  } else {
    theta_k <- t(theta)
  }
  niter <- ncol(theta_k)
  theta_hat <- as.numeric(rowMeans(theta_k))
  if (k > 0) {
    mu_new <- mu + beta_k*(theta_hat - mu)
  } else {
    mu_new <- mu
  }
  theta_k_c <- theta_k - theta_hat
  cov_hat <- (theta_k_c %*% t(theta_k_c))/niter
  cov_new <- (1 - beta_k)*cov + beta_k*cov_hat
  if ((what == 1 | what == 2)  & (cov_new[1, 1] <= 0)) {
    cov_new <- 1e-1
  }

  return(list("rho" = rho_new,
              "covariance" = cov_new,
              "mu" = mu_new,
              "ar" = ar_new))
}
