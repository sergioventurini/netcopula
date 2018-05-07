## ----simulate0, echo=FALSE, results='hide', fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
library(netcopula)

path_to_save <- "~/dev/netcopula/demo/tmp"
rng <- 1809
set.seed(rng)

# data generation parameters
N <- 1 # number of simulated data sets
n_outcomes <- 3
n_trt <- 3
n_studies <- 30
n_bin <- rep(100, n_outcomes) # binomial distributions sizes
d_init <- array(data = c(0, -2, -1, 0, 3, 2), dim = c(n_trt - 1, n_outcomes))
sigma_vec <- c(1, 1.6, 1.8) # variances in Sigma_M
rho_vec <- rep(.5, n_outcomes)
# gamma_vec <- c(.3, .7, 0) # components of copula association matrix
# gamma_vec <- c(.3, .1, .2) # components of copula association matrix
gamma_vec <- c(.3, .1, .5) # components of copula association matrix
# gamma_vec <- c(-.9, -.9, .99) # components of copula association matrix
# gamma_vec <- c(0, 0, 0) # components of copula association matrix
# gamma_vec <- c(.5, .5, .5) # components of copula association matrix
ref_trt <- 1

# inizialization parameters
nGamma <- 1 # nc_data@n_treatments
prm_init <- list(mu_sigma2 = 1e-1, d_sigma2 = 1e-1, beta_sigma = 10^(-1/8), nGamma = nGamma, df = 1e2, eta_prior = 1)

# prior parameters
prm_prior <- list(mu_sigma2 = 1e-3, d_sigma2 = 1e3, beta_sigma = 10^(-1/8), nGamma = nGamma, df = 1e2, eta_prior = 1)

# MCMC settings
burnin <- 10000
nsim <- 10000
nthin <- 10 #100
eps <- eps_ab <- 1e-12
xmax <- 1000
xmin <- -xmax
tuning <- list(eps = eps, eps_ab = eps_ab, xmin = xmin, xmax = xmax)
every <- 50
maxiter <- min(floor((burnin + nsim)/every), floor(burnin/every))
adaptation <- list(every = every, maxiter = maxiter, miniter = 0, alpha = .75, beta = .8, gamma = 0, tar = 0.234, tol = .01)
# prm_prop <- list(Gamma_update = "Talhouketal", eta_prop = .5)
prm_prop <- list(Gamma_update = "PX-RPMH", eta_prop = .5)
# prm_prop <- list(Gamma_update = "IMH", eta_prop = .5)
mu_sigma2_sim <- 1e-9
n_datapoints <- 2*n_studies

res_Gamma <- matrix(NA, nrow = N, ncol = n_outcomes)
res_x <- array(NA, dim = c(N, n_datapoints, n_outcomes))
for (i in 1:N) {
  # data creation
  nc_data_tmp <- nc_data_simulate(d_init, sigma_vec, rho_vec, gamma_vec, n_studies, n_trt, n_bin, n_outcomes, mu_sigma2_sim, ref_trt)
  # debugonce(nc_data_missing)
  # nc_data_tmp <- nc_data_missing(nc_data_tmp, NULL, "MCAR", 5)
  
  nc_data <- nc_data_tmp$nc_data

  # true values
  Gamma_tmp <- nc_data_tmp$Gamma
  mu <- nc_data_tmp$mu
  delta <- nc_data_tmp$delta
  d <- rbind(0, nc_data_tmp$d)
  Sigma_M <- nc_data_tmp$Sigma_M
  x <- nc_data_tmp$x
  D <- diag(1/colSums(x^2))
  true_vals <- list(mu = mu, delta = delta, d = d, Sigma_M = Sigma_M, Gamma = list(Gamma_tmp), D = list(D), x = x)

  # starting values
  # start_vals <- true_vals
  start_vals <- nc_init(nc_data, prm_init)
  # start_vals$x <- x
  # start_vals$Gamma <- list(Gamma_tmp)
  # start_vals$Gamma <- list(diag(n_outcomes))
  # start_vals$mu <- mu
  # start_vals$delta <- delta
  # start_vals$d <- d
  # start_vals$Sigma_M <- Sigma_M

  # MCMC simulation
  # set.seed(rng)
  system.time(res <- netcopula(nc_data, burnin = burnin, nsim = nsim, nthin = nthin, prm_prior, prm_prop, prm_init, tuning = tuning, adaptation = adaptation, verbose = TRUE, mu_delta_mh = TRUE, random_start = FALSE, init = start_vals))

  res_Gamma[i, ] <- apply(res@Gamma[1, , ], 1, mean)
  res_x[i, , ] <- apply(res@x, c(1, 2), mean)

  # save("nc_data_tmp", "res", file = file.path(path_to_save, paste0("examples3_", i, ".RData")))
}

post_mean_Gamma <- apply(res@Gamma[1, , (burnin + 1):(nsim + burnin)], 1, mean)
# post_mean_Gamma <- mean(res@Gamma[1, , (burnin + 1):(nsim + burnin)])
post_mean_Gamma
Gamma_est <- diag(n_outcomes)
Gamma_est[lower.tri(Gamma_est)] <- Gamma_est[upper.tri(Gamma_est)] <- post_mean_Gamma
eigen(Gamma_tmp)$values
eigen(Gamma_est)$values
# library(ppcor)
# pcor(x)$estimate
# x_est <- apply(res@x_unadj[, , (burnin + 1):(nsim + burnin)], c(1, 2), mean)
# pcor(x_est)$estimate

plot(x = res, to_plot = "Gamma", type = "trace")
plot(x = res, to_plot = "Sigma_M", type = "trace")
plot(x = res, to_plot = "Corr_M", type = "trace")
plot(x = res, to_plot = "mu", type = "trace")
plot(x = res, to_plot = "delta", type = "trace")
plot(x = res, to_plot = "d", type = "trace")
plot(x = res, to_plot = "x_unadj", latent_var = TRUE, type = "trace")

res_mcmc <- as.mcmc(res, latent_var = TRUE)
out <- summary(res_mcmc)

corr_x <- apply(res@x_unadj[ , , (burnin + 1):(nsim + burnin)], 3, cor)
matrix(rowMeans(corr_x), nrow = n_outcomes, ncol = n_outcomes)

plot_est(out, Gamma_tmp[lower.tri(Gamma_tmp)], "Gamma")
plot_est(out, mu, "mu")
plot_est(out, delta[seq(2, 2*n_studies, by = 2), ], "delta")
plot_est(out, d[-1, ], "d")
plot_est(out, Sigma_M[lower.tri(Sigma_M, diag = TRUE)], "Sigma_M")
plot_est(out, x, "x_unadj")

## ----simulate5, echo=FALSE, results='hide', fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
library(netcopula)

path_to_save <- "~/dev/netcopula/demo/tmp"
rng <- 1809
set.seed(rng)

# data generation parameters
N <- 1 # number of simulated data sets
n_outcomes <- 3
n_trt <- 3
n_studies <- 30
n_bin <- rep(100, n_outcomes) # binomial distributions sizes
d_init <- array(data = c(0, -2, -1, 0, 3, 2), dim = c(n_trt - 1, n_outcomes))
sigma_vec <- c(1, 1.6, 1.8) # variances in Sigma_M
rho_vec <- rep(.5, n_outcomes)
# gamma_vec <- c(.3, .7, 0) # components of copula association matrix
# gamma_vec <- c(.3, .1, .2) # components of copula association matrix
gamma_vec <- c(.3, .1, .5) # components of copula association matrix
# gamma_vec <- c(-.9, -.9, .99) # components of copula association matrix
# gamma_vec <- c(0, 0, 0) # components of copula association matrix
# gamma_vec <- c(.5, .5, .5) # components of copula association matrix
ref_trt <- 1

# inizialization parameters
nGamma <- 1 # nc_data@n_treatments
prm_init <- list(mu_sigma2 = 1e-1, d_sigma2 = 1e-1, beta_sigma = 10^(-1/8), nGamma = nGamma, df = 1e2, eta_prior = 1)

# prior parameters
prm_prior <- list(mu_sigma2 = 1e-3, d_sigma2 = 1e3, beta_sigma = 10^(-1/8), nGamma = nGamma, df = 1e2, eta_prior = 1)

# MCMC settings
burnin <- 10000
nsim <- 10000
nthin <- 10 #100
eps <- eps_ab <- 1e-12
xmax <- 1000
xmin <- -xmax
tuning <- list(eps = eps, eps_ab = eps_ab, xmin = xmin, xmax = xmax)
every <- 50
maxiter <- min(floor((burnin + nsim)/every), floor(burnin/every))
adaptation <- list(every = every, maxiter = maxiter, miniter = 0, alpha = .75, beta = .8, gamma = 0, tar = 0.234, tol = .01)
# prm_prop <- list(Gamma_update = "Talhouketal", eta_prop = .5)
prm_prop <- list(Gamma_update = "PX-RPMH", eta_prop = .5)
# prm_prop <- list(Gamma_update = "IMH", eta_prop = .5)
mu_sigma2_sim <- 1e-9
n_datapoints <- 2*n_studies

res_Gamma <- matrix(NA, nrow = N, ncol = n_outcomes)
res_x <- array(NA, dim = c(N, n_datapoints, n_outcomes))
for (i in 1:N) {
  # data creation
  nc_data_tmp <- nc_data_simulate(d_init, sigma_vec, rho_vec, gamma_vec, n_studies, n_trt, n_bin, n_outcomes, mu_sigma2_sim, ref_trt)
  # debugonce(nc_data_missing)
  nc_data_tmp <- nc_data_missing(nc_data_tmp, NULL, "MCAR", 5)
  
  nc_data <- nc_data_tmp$nc_data

  # true values
  Gamma_tmp <- nc_data_tmp$Gamma
  mu <- nc_data_tmp$mu
  delta <- nc_data_tmp$delta
  d <- rbind(0, nc_data_tmp$d)
  Sigma_M <- nc_data_tmp$Sigma_M
  x <- nc_data_tmp$x
  D <- diag(1/colSums(x^2))
  true_vals <- list(mu = mu, delta = delta, d = d, Sigma_M = Sigma_M, Gamma = list(Gamma_tmp), D = list(D), x = x)

  # starting values
  # start_vals <- true_vals
  start_vals <- nc_init(nc_data, prm_init)
  # start_vals$x <- x
  # start_vals$Gamma <- list(Gamma_tmp)
  # start_vals$Gamma <- list(diag(n_outcomes))
  # start_vals$mu <- mu
  # start_vals$delta <- delta
  # start_vals$d <- d
  # start_vals$Sigma_M <- Sigma_M

  # MCMC simulation
  # set.seed(rng)
  system.time(res <- netcopula(nc_data, burnin = burnin, nsim = nsim, nthin = nthin, prm_prior, prm_prop, prm_init, tuning = tuning, adaptation = adaptation, verbose = TRUE, mu_delta_mh = TRUE, random_start = FALSE, init = start_vals))

  res_Gamma[i, ] <- apply(res@Gamma[1, , ], 1, mean)
  res_x[i, , ] <- apply(res@x, c(1, 2), mean)

  # save("nc_data_tmp", "res", file = file.path(path_to_save, paste0("examples3_", i, ".RData")))
}

post_mean_Gamma <- apply(res@Gamma[1, , (burnin + 1):(nsim + burnin)], 1, mean)
# post_mean_Gamma <- mean(res@Gamma[1, , (burnin + 1):(nsim + burnin)])
post_mean_Gamma
Gamma_est <- diag(n_outcomes)
Gamma_est[lower.tri(Gamma_est)] <- Gamma_est[upper.tri(Gamma_est)] <- post_mean_Gamma
eigen(Gamma_tmp)$values
eigen(Gamma_est)$values
# library(ppcor)
# pcor(x)$estimate
# x_est <- apply(res@x_unadj[, , (burnin + 1):(nsim + burnin)], c(1, 2), mean)
# pcor(x_est)$estimate

plot(x = res, to_plot = "Gamma", type = "trace")
plot(x = res, to_plot = "Sigma_M", type = "trace")
plot(x = res, to_plot = "Corr_M", type = "trace")
plot(x = res, to_plot = "mu", type = "trace")
plot(x = res, to_plot = "delta", type = "trace")
plot(x = res, to_plot = "d", type = "trace")
plot(x = res, to_plot = "x_unadj", latent_var = TRUE, type = "trace")

res_mcmc <- as.mcmc(res, latent_var = TRUE)
out <- summary(res_mcmc)

corr_x <- apply(res@x_unadj[ , , (burnin + 1):(nsim + burnin)], 3, cor)
matrix(rowMeans(corr_x), nrow = n_outcomes, ncol = n_outcomes)

plot_est(out, Gamma_tmp[lower.tri(Gamma_tmp)], "Gamma")
plot_est(out, mu, "mu")
plot_est(out, delta[seq(2, 2*n_studies, by = 2), ], "delta")
plot_est(out, d[-1, ], "d")
plot_est(out, Sigma_M[lower.tri(Sigma_M, diag = TRUE)], "Sigma_M")
plot_est(out, x, "x_unadj")

## ----simulate10, echo=FALSE, results='hide', fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
library(netcopula)

path_to_save <- "~/dev/netcopula/demo/tmp"
rng <- 1809
set.seed(rng)

# data generation parameters
N <- 1 # number of simulated data sets
n_outcomes <- 3
n_trt <- 3
n_studies <- 30
n_bin <- rep(100, n_outcomes) # binomial distributions sizes
d_init <- array(data = c(0, -2, -1, 0, 3, 2), dim = c(n_trt - 1, n_outcomes))
sigma_vec <- c(1, 1.6, 1.8) # variances in Sigma_M
rho_vec <- rep(.5, n_outcomes)
# gamma_vec <- c(.3, .7, 0) # components of copula association matrix
# gamma_vec <- c(.3, .1, .2) # components of copula association matrix
gamma_vec <- c(.3, .1, .5) # components of copula association matrix
# gamma_vec <- c(-.9, -.9, .99) # components of copula association matrix
# gamma_vec <- c(0, 0, 0) # components of copula association matrix
# gamma_vec <- c(.5, .5, .5) # components of copula association matrix
ref_trt <- 1

# inizialization parameters
nGamma <- 1 # nc_data@n_treatments
prm_init <- list(mu_sigma2 = 1e-1, d_sigma2 = 1e-1, beta_sigma = 10^(-1/8), nGamma = nGamma, df = 1e2, eta_prior = 1)

# prior parameters
prm_prior <- list(mu_sigma2 = 1e-3, d_sigma2 = 1e3, beta_sigma = 10^(-1/8), nGamma = nGamma, df = 1e2, eta_prior = 1)

# MCMC settings
burnin <- 10000
nsim <- 10000
nthin <- 10 #100
eps <- eps_ab <- 1e-12
xmax <- 1000
xmin <- -xmax
tuning <- list(eps = eps, eps_ab = eps_ab, xmin = xmin, xmax = xmax)
every <- 50
maxiter <- min(floor((burnin + nsim)/every), floor(burnin/every))
adaptation <- list(every = every, maxiter = maxiter, miniter = 0, alpha = .75, beta = .8, gamma = 0, tar = 0.234, tol = .01)
# prm_prop <- list(Gamma_update = "Talhouketal", eta_prop = .5)
prm_prop <- list(Gamma_update = "PX-RPMH", eta_prop = .5)
# prm_prop <- list(Gamma_update = "IMH", eta_prop = .5)
mu_sigma2_sim <- 1e-9
n_datapoints <- 2*n_studies

res_Gamma <- matrix(NA, nrow = N, ncol = n_outcomes)
res_x <- array(NA, dim = c(N, n_datapoints, n_outcomes))
for (i in 1:N) {
  # data creation
  nc_data_tmp <- nc_data_simulate(d_init, sigma_vec, rho_vec, gamma_vec, n_studies, n_trt, n_bin, n_outcomes, mu_sigma2_sim, ref_trt)
  # debugonce(nc_data_missing)
  nc_data_tmp <- nc_data_missing(nc_data_tmp, NULL, "MCAR", 10)
  
  nc_data <- nc_data_tmp$nc_data

  # true values
  Gamma_tmp <- nc_data_tmp$Gamma
  mu <- nc_data_tmp$mu
  delta <- nc_data_tmp$delta
  d <- rbind(0, nc_data_tmp$d)
  Sigma_M <- nc_data_tmp$Sigma_M
  x <- nc_data_tmp$x
  D <- diag(1/colSums(x^2))
  true_vals <- list(mu = mu, delta = delta, d = d, Sigma_M = Sigma_M, Gamma = list(Gamma_tmp), D = list(D), x = x)

  # starting values
  # start_vals <- true_vals
  start_vals <- nc_init(nc_data, prm_init)
  # start_vals$x <- x
  # start_vals$Gamma <- list(Gamma_tmp)
  # start_vals$Gamma <- list(diag(n_outcomes))
  # start_vals$mu <- mu
  # start_vals$delta <- delta
  # start_vals$d <- d
  # start_vals$Sigma_M <- Sigma_M

  # MCMC simulation
  # set.seed(rng)
  system.time(res <- netcopula(nc_data, burnin = burnin, nsim = nsim, nthin = nthin, prm_prior, prm_prop, prm_init, tuning = tuning, adaptation = adaptation, verbose = TRUE, mu_delta_mh = TRUE, random_start = FALSE, init = start_vals))

  res_Gamma[i, ] <- apply(res@Gamma[1, , ], 1, mean)
  res_x[i, , ] <- apply(res@x, c(1, 2), mean)

  # save("nc_data_tmp", "res", file = file.path(path_to_save, paste0("examples3_", i, ".RData")))
}

post_mean_Gamma <- apply(res@Gamma[1, , (burnin + 1):(nsim + burnin)], 1, mean)
# post_mean_Gamma <- mean(res@Gamma[1, , (burnin + 1):(nsim + burnin)])
post_mean_Gamma
Gamma_est <- diag(n_outcomes)
Gamma_est[lower.tri(Gamma_est)] <- Gamma_est[upper.tri(Gamma_est)] <- post_mean_Gamma
eigen(Gamma_tmp)$values
eigen(Gamma_est)$values
# library(ppcor)
# pcor(x)$estimate
# x_est <- apply(res@x_unadj[, , (burnin + 1):(nsim + burnin)], c(1, 2), mean)
# pcor(x_est)$estimate

plot(x = res, to_plot = "Gamma", type = "trace")
plot(x = res, to_plot = "Sigma_M", type = "trace")
plot(x = res, to_plot = "Corr_M", type = "trace")
plot(x = res, to_plot = "mu", type = "trace")
plot(x = res, to_plot = "delta", type = "trace")
plot(x = res, to_plot = "d", type = "trace")
plot(x = res, to_plot = "x_unadj", latent_var = TRUE, type = "trace")

res_mcmc <- as.mcmc(res, latent_var = TRUE)
out <- summary(res_mcmc)

corr_x <- apply(res@x_unadj[ , , (burnin + 1):(nsim + burnin)], 3, cor)
matrix(rowMeans(corr_x), nrow = n_outcomes, ncol = n_outcomes)

plot_est(out, Gamma_tmp[lower.tri(Gamma_tmp)], "Gamma")
plot_est(out, mu, "mu")
plot_est(out, delta[seq(2, 2*n_studies, by = 2), ], "delta")
plot_est(out, d[-1, ], "d")
plot_est(out, Sigma_M[lower.tri(Sigma_M, diag = TRUE)], "Sigma_M")
plot_est(out, x, "x_unadj")

## ----simulate15, echo=FALSE, results='hide', fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
library(netcopula)

path_to_save <- "~/dev/netcopula/demo/tmp"
rng <- 1809
set.seed(rng)

# data generation parameters
N <- 1 # number of simulated data sets
n_outcomes <- 3
n_trt <- 3
n_studies <- 30
n_bin <- rep(100, n_outcomes) # binomial distributions sizes
d_init <- array(data = c(0, -2, -1, 0, 3, 2), dim = c(n_trt - 1, n_outcomes))
sigma_vec <- c(1, 1.6, 1.8) # variances in Sigma_M
rho_vec <- rep(.5, n_outcomes)
# gamma_vec <- c(.3, .7, 0) # components of copula association matrix
# gamma_vec <- c(.3, .1, .2) # components of copula association matrix
gamma_vec <- c(.3, .1, .5) # components of copula association matrix
# gamma_vec <- c(-.9, -.9, .99) # components of copula association matrix
# gamma_vec <- c(0, 0, 0) # components of copula association matrix
# gamma_vec <- c(.5, .5, .5) # components of copula association matrix
ref_trt <- 1

# inizialization parameters
nGamma <- 1 # nc_data@n_treatments
prm_init <- list(mu_sigma2 = 1e-1, d_sigma2 = 1e-1, beta_sigma = 10^(-1/8), nGamma = nGamma, df = 1e2, eta_prior = 1)

# prior parameters
prm_prior <- list(mu_sigma2 = 1e-3, d_sigma2 = 1e3, beta_sigma = 10^(-1/8), nGamma = nGamma, df = 1e2, eta_prior = 1)

# MCMC settings
burnin <- 10000
nsim <- 10000
nthin <- 10 #100
eps <- eps_ab <- 1e-12
xmax <- 1000
xmin <- -xmax
tuning <- list(eps = eps, eps_ab = eps_ab, xmin = xmin, xmax = xmax)
every <- 50
maxiter <- min(floor((burnin + nsim)/every), floor(burnin/every))
adaptation <- list(every = every, maxiter = maxiter, miniter = 0, alpha = .75, beta = .8, gamma = 0, tar = 0.234, tol = .01)
# prm_prop <- list(Gamma_update = "Talhouketal", eta_prop = .5)
prm_prop <- list(Gamma_update = "PX-RPMH", eta_prop = .5)
# prm_prop <- list(Gamma_update = "IMH", eta_prop = .5)
mu_sigma2_sim <- 1e-9
n_datapoints <- 2*n_studies

res_Gamma <- matrix(NA, nrow = N, ncol = n_outcomes)
res_x <- array(NA, dim = c(N, n_datapoints, n_outcomes))
for (i in 1:N) {
  # data creation
  nc_data_tmp <- nc_data_simulate(d_init, sigma_vec, rho_vec, gamma_vec, n_studies, n_trt, n_bin, n_outcomes, mu_sigma2_sim, ref_trt)
  # debugonce(nc_data_missing)
  nc_data_tmp <- nc_data_missing(nc_data_tmp, NULL, "MCAR", 15)
  
  nc_data <- nc_data_tmp$nc_data

  # true values
  Gamma_tmp <- nc_data_tmp$Gamma
  mu <- nc_data_tmp$mu
  delta <- nc_data_tmp$delta
  d <- rbind(0, nc_data_tmp$d)
  Sigma_M <- nc_data_tmp$Sigma_M
  x <- nc_data_tmp$x
  D <- diag(1/colSums(x^2))
  true_vals <- list(mu = mu, delta = delta, d = d, Sigma_M = Sigma_M, Gamma = list(Gamma_tmp), D = list(D), x = x)

  # starting values
  # start_vals <- true_vals
  start_vals <- nc_init(nc_data, prm_init)
  # start_vals$x <- x
  # start_vals$Gamma <- list(Gamma_tmp)
  # start_vals$Gamma <- list(diag(n_outcomes))
  # start_vals$mu <- mu
  # start_vals$delta <- delta
  # start_vals$d <- d
  # start_vals$Sigma_M <- Sigma_M

  # MCMC simulation
  # set.seed(rng)
  system.time(res <- netcopula(nc_data, burnin = burnin, nsim = nsim, nthin = nthin, prm_prior, prm_prop, prm_init, tuning = tuning, adaptation = adaptation, verbose = TRUE, mu_delta_mh = TRUE, random_start = FALSE, init = start_vals))

  res_Gamma[i, ] <- apply(res@Gamma[1, , ], 1, mean)
  res_x[i, , ] <- apply(res@x, c(1, 2), mean)

  # save("nc_data_tmp", "res", file = file.path(path_to_save, paste0("examples3_", i, ".RData")))
}

post_mean_Gamma <- apply(res@Gamma[1, , (burnin + 1):(nsim + burnin)], 1, mean)
# post_mean_Gamma <- mean(res@Gamma[1, , (burnin + 1):(nsim + burnin)])
post_mean_Gamma
Gamma_est <- diag(n_outcomes)
Gamma_est[lower.tri(Gamma_est)] <- Gamma_est[upper.tri(Gamma_est)] <- post_mean_Gamma
eigen(Gamma_tmp)$values
eigen(Gamma_est)$values
# library(ppcor)
# pcor(x)$estimate
# x_est <- apply(res@x_unadj[, , (burnin + 1):(nsim + burnin)], c(1, 2), mean)
# pcor(x_est)$estimate

plot(x = res, to_plot = "Gamma", type = "trace")
plot(x = res, to_plot = "Sigma_M", type = "trace")
plot(x = res, to_plot = "Corr_M", type = "trace")
plot(x = res, to_plot = "mu", type = "trace")
plot(x = res, to_plot = "delta", type = "trace")
plot(x = res, to_plot = "d", type = "trace")
plot(x = res, to_plot = "x_unadj", latent_var = TRUE, type = "trace")

res_mcmc <- as.mcmc(res, latent_var = TRUE)
out <- summary(res_mcmc)

corr_x <- apply(res@x_unadj[ , , (burnin + 1):(nsim + burnin)], 3, cor)
matrix(rowMeans(corr_x), nrow = n_outcomes, ncol = n_outcomes)

plot_est(out, Gamma_tmp[lower.tri(Gamma_tmp)], "Gamma")
plot_est(out, mu, "mu")
plot_est(out, delta[seq(2, 2*n_studies, by = 2), ], "delta")
plot_est(out, d[-1, ], "d")
plot_est(out, Sigma_M[lower.tri(Sigma_M, diag = TRUE)], "Sigma_M")
plot_est(out, x, "x_unadj")

