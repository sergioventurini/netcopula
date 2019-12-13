## WARNING: the following simulation study will take a long time! ##

library(netcopula)

## replace the following path with your own to save the simulation results
path_to_save <- "~/dev/netcopula/demo"
set.seed(1406)

## data generation parameters
N <- 500 # number of simulated data sets
n_outcomes <- 3
n_trt <- 3
n_studies <- 30
n_bin <- rep(100, n_outcomes) # binomial distributions sizes
d_init <- array(data = c(0, -0.5, -1, 0, 0.8, 0.3), dim = c(n_trt - 1, n_outcomes))
sigma_vec <- c(1, 1.6, 1.8) # standard deviations in Sigma_M
rho_vec <- rep(.5, n_outcomes)
gamma_vec <- c(.3, .1, .5) # components of copula association matrix
ref_trt <- 1

## inizialization parameters
nGamma <- 1 # nc_data@n_treatments
prm_init <- list(mu_sigma2 = 1e-1, d_sigma2 = 1e-1, beta_sigma = 10^(-1/8),
  nGamma = nGamma, eta_prior = 1)

## prior parameters
prm_prior <- list(mu_sigma2 = 10^3, d_sigma2 = 10^3,
  beta_sigma = 10^(-1/8), nGamma = nGamma, eta_prior = 1)

## MCMC settings
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

## array to store the result summaries
res_summary <- array(NA,
  dim = c(N, n_outcomes*(n_outcomes + 1)/2 + n_trt*n_outcomes + n_outcomes, 2))
doagain <- TRUE

for (i in 1:N) {
  print(paste0("--- SIMULATION ", i, "/", N, " ---"))
  
  ## data creation
  nc_data_tmp <- nc_data_simulate(d_init, sigma_vec, rho_vec, gamma_vec,
    n_studies, n_trt, n_bin, n_outcomes, mu_sigma2_sim, ref_trt)
  # nc_data_tmp <- nc_data_missing(nc_data_tmp, n_miss = 10)
  
  nc_data <- nc_data_tmp$nc_data
  
  ## true values
  Gamma_tmp <- nc_data_tmp$Gamma
  mu <- nc_data_tmp$mu
  delta <- nc_data_tmp$delta
  d <- rbind(0, nc_data_tmp$d)
  Sigma_M <- nc_data_tmp$Sigma_M
  x <- nc_data_tmp$x
  D <- diag(1/colSums(x^2))
  true_vals <- list(mu = mu, delta = delta, d = d, Sigma_M = Sigma_M,
    Gamma = list(Gamma_tmp), D = list(D), x = x)
  
  ## starting values
  start_vals <- nc_init(nc_data, prm_init)
  
  ## MCMC simulation
  while (doagain) {
    res <- try(netcopula(nc_data, control, prm_prior, prm_init,
      init = start_vals, tuning = tuning, adaptation = adaptation),
        silent = TRUE)
    if (class(res) == "try-error") {
      doagain <- TRUE
      print(" --> Do it again! :(")
    } else {
      doagain <- FALSE
    }
  }

  res_summary_tmp <- summary(res, regex_pars = c("d", "Sigma", "Gamma"))
  res_summary[i, , 1] <- res_summary_tmp[[1]][, 1] # posterior means
  res_summary[i, , 2] <- res_summary_tmp[[2]][, 3] # posterior medians
  
  doagain <- TRUE

  ## uncomment the following lines to save the simulation results
  # save("nc_data_tmp", "res", file = file.path(path_to_save,
  #   paste0("sim_study_", i, ".RData")))
  # save.image(file = file.path(path_to_save, "sim_study.RData"))
}

dimnames(res_summary) <- list(NULL, rownames(res_summary_tmp[[1]]), NULL)
res_summary <- res_summary[, -c(1, 4, 7), ]

summ_pmean <- colMeans(res_summary[, , 1], na.rm = TRUE)
summ_pmed <- colMeans(res_summary[, , 2], na.rm = TRUE)

## plotting results
toplot <- t(res_summary[, , 2]) # plot estimated posterior medians
true_val <- c(as.numeric(d[-1, ]), Sigma_M[lower.tri(Sigma_M, diag = TRUE)],
  Gamma_tmp[lower.tri(Gamma_tmp)])
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
  "#D55E00", "#CC79A7")
toinst <- require(scales)
if (!toinst) install.packages("scales")
par(mar = c(5, 6.5, 1.5, 1.5) + 0.1)
plot(toplot[1, ], rep(1, N), xlim = c(min(toplot), max(toplot)),
  ylim = c(0, nrow(toplot) + 1), type = "n",
  xlab = "Posterior median", ylab = "", yaxt = "n")
for (j in 1:6) {
  segments(-10, j, 10, j, lty = 2, col = gray(.9))
  points(toplot[j, ], rep(j, N), pch = 20, col = alpha(cbPalette[2], alpha = .25))
  points(summ_pmed[j], j, pch = 18, col = cbPalette[7], cex = 1.75)
  points(true_val[j], j, pch = 17, col = cbPalette[4], cex = 1.25)
}
for (j in 7:12) {
  segments(-10, j, 10, j, lty = 2, col = gray(.9))
  points(toplot[j, ], rep(j, N), pch = 20, col = alpha(cbPalette[3], alpha = .25))
  points(summ_pmed[j], j, pch = 18, col = cbPalette[7], cex = 1.75)
  points(true_val[j], j, pch = 17, col = cbPalette[4], cex = 1.25)
}
for (j in 13:15) {
  segments(-10, j, 10, j, lty = 2, col = gray(.9))
  points(toplot[j, ], rep(j, N), pch = 20, col = alpha(cbPalette[1], alpha = .25))
  points(summ_pmed[j], j, pch = 18, col = cbPalette[7], cex = 1.75)
  points(true_val[j], j, pch = 17, col = cbPalette[4], cex = 1.25)
}
axis(side = 2, at = 1:nrow(toplot), labels = rownames(toplot), las = 2, cex.axis = .7)
title(ylab = "Parameter", line = 4.5)
legend(x = "topright", legend = c("Average", "True"),
  col = cbPalette[c(7, 4)], pch = c(18, 17), pt.cex = c(1.75, 1.25), horiz = TRUE)
