## WARNING: the following simulation study will take a long time! ##

library(netcopula)

# replace the following path with your own to save the simulation results
path_to_save <- "~/dev/netcopula/demo"
set.seed(1406)

# data generation parameters
N <- 200 # number of simulated data sets
n_outcomes <- 3
n_trt <- 3
n_studies <- 30
n_bin <- rep(100, n_outcomes) # binomial distributions sizes
d_init <- array(data = c(0, -0.5, -1, 0, 0.8, 0.3), dim = c(n_trt - 1, n_outcomes))
sigma_vec <- c(1, 1.6, 1.8) # standard deviations in Sigma_M
rho_vec <- rep(.5, n_outcomes)
gamma_vec <- c(.3, .1, .5) # components of copula association matrix
ref_trt <- 1

# inizialization parameters
nGamma <- 1 # nc_data@n_treatments
prm_init <- list(mu_sigma2 = 1e-1, d_sigma2 = 1e-1, beta_sigma = 10^(-1/8),
  nGamma = nGamma, eta_prior = 1)

# prior parameters
prm_prior <- list(mu_sigma2 = 10^3, d_sigma2 = 10^3,
  beta_sigma = 10^(-1/8), nGamma = nGamma, eta_prior = 1)

# MCMC settings
burnin <- 50000
nsim <- 10000
nthin <- 1
prm_prior <- list(Gamma.update = "PX-RPMH", eta = .5, sigma.r = .01)
control <- list(burnin = burnin, nsim = nsim,
  Gamma.update = prm_prior[["Gamma.update"]], eta.prop = prm_prior[["eta"]],
  sigma.r.prop = prm_prior[["sigma.r"]], random.start = TRUE, verbose = FALSE,
  thin = nthin, store.burnin = TRUE)
eps <- eps_ab <- 1e-12
xmax <- 1000
xmin <- -xmax
tuning <- list(eps = eps, eps_ab = eps_ab, xmin = xmin, xmax = xmax)
every <- 50
maxiter <- max(25, floor(burnin/every))
adaptation <- list(every = every, maxiter = maxiter, miniter = 5,
  alpha = .75, beta = .8, gamma = 0, tar = 0.234, tol = .01)

mu_sigma2_sim <- 1e-1
n_datapoints <- 2*n_studies

res_summary <- array(NA,
  dim = c(N, n_datapoints*n_outcomes + n_outcomes*(n_outcomes + 1)/2 + n_outcomes + n_outcomes*2, 8))
doagain <- TRUE
for (i in 1:N) {
  print(paste0("--- SIMULATION ", i, "/", N, " ---"))
  
  # data creation
  nc_data_tmp <- nc_data_simulate(d_init, sigma_vec, rho_vec, gamma_vec,
    n_studies, n_trt, n_bin, n_outcomes, mu_sigma2_sim, ref_trt)
  # nc_data_tmp <- nc_data_missing(nc_data_tmp, n_miss = 10)
  
  nc_data <- nc_data_tmp$nc_data
  
  # true values
  Gamma_tmp <- nc_data_tmp$Gamma
  mu <- nc_data_tmp$mu
  delta <- nc_data_tmp$delta
  d <- rbind(0, nc_data_tmp$d)
  Sigma_M <- nc_data_tmp$Sigma_M
  x <- nc_data_tmp$x
  D <- diag(1/colSums(x^2))
  true_vals <- list(mu = mu, delta = delta, d = d, Sigma_M = Sigma_M,
    Gamma = list(Gamma_tmp), D = list(D), x = x)
  
  # starting values
  start_vals <- nc_init(nc_data, prm_init)
  
  # MCMC simulation
  while (doagain) {
    res <- try(netcopula(nc_data, control, prm_prior, prm_init, init = start_vals,
      tuning = tuning, adaptation = adaptation), silent = TRUE)
    if (class(res) == "try-error") {
      doagain <- TRUE
      print(" --> Do it again! :(")
    } else {
      doagain <- FALSE
    }
  }

  # res_summary[i, , ] <- summary(res, latent = FALSE, print = FALSE)
  # dimnames(res_summary)[2:3] <- dimnames(summary(res, latent = FALSE, print = FALSE))
  
  doagain <- TRUE

  # uncomment the following lines to save the simulation results
  # save("nc_data_tmp", "res", file = file.path(path_to_save, paste0("sim_study_", i, ".RData")))
  # save.image(file = file.path(path_to_save, "sim_study.RData"))
}

# plotting results
summ <- apply(res_summary, c(2, 3), mean)

param <- "d"
idx <- grep(paste0(param, "["), dimnames(summ)[[1]], fixed = TRUE)
d_true <- as.numeric(d_init)
bias_d <- summ[idx, 1]
bias_d <- bias_d - d_true
cov_d_tmp <- res_summary[, , c(4, 8)]
cov_d <- numeric(length(d_true))
for (i in 1:length(d_true)) {
  cov_d[i] <- sum(cov_d_tmp[, i , 1] <= d_true[i] & cov_d_tmp[, i , 2] >= d_true[i])/N
}
plot(bias_d, xlab = "Parameter", ylab = "Bias", xaxt = "n", cex.axis = .7,
  ylim = c(min(bias_d, 0.1), max(bias_d, 0)))
axis(side = 1, at = 1:length(d_true), labels = dimnames(summ)[[1]][idx])
abline(h = 0, lty = 2, col = gray(.9))
plot(cov_d, xlab = "Parameter", ylab = "Coverage", xaxt = "n", cex.axis = .7,
  ylim = c(0, 1))
axis(side = 1, at = 1:length(d_true), labels = dimnames(summ)[[1]][idx])
abline(h = .95, lty = 2, col = gray(.9))

param <- "Sigma_M"
idx <- grep(paste0(param, "["), dimnames(summ)[[1]], fixed = TRUE)
Sigma_M_true <- as.numeric(Sigma_M[lower.tri(Sigma_M, diag = TRUE)])
bias_Sigma_M <- summ[grep(paste0(param, "["), dimnames(summ)[[1]], fixed = TRUE), 1]
bias_Sigma_M <- bias_Sigma_M - Sigma_M_true
cov_Sigma_M_tmp <- res_summary[, grep(paste0(param, "["), dimnames(summ)[[1]], fixed = TRUE), c(4, 8)]
cov_Sigma_M <- numeric(length(Sigma_M_true))
for (i in 1:length(Sigma_M_true)) {
  cov_Sigma_M[i] <- sum(cov_Sigma_M_tmp[, i , 1] <= Sigma_M_true[i] & cov_Sigma_M_tmp[, i , 2] >= Sigma_M_true[i])/N
}
plot(bias_Sigma_M, xlab = "Parameter", ylab = "Bias", xaxt = "n", cex.axis = .7,
  ylim = c(min(bias_Sigma_M, 0.1), max(bias_Sigma_M, 0)))
axis(side = 1, at = 1:length(d_true), labels = dimnames(summ)[[1]][idx],
  cex.axis = .7)
abline(h = 0, lty = 2, col = gray(.9))
plot(cov_Sigma_M, xlab = "Parameter", ylab = "Bias", xaxt = "n", cex.axis = .7,
  ylim = c(0, 1))
axis(side = 1, at = 1:length(d_true), labels = dimnames(summ)[[1]][idx],
  cex.axis = .7)
abline(h = .95, lty = 2, col = gray(.9))
