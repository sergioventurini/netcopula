## WARNING: the following analysis may take a long time to run! ##

library(netcopula)

# data loading and preparation
data(homesafety)
id_arms <- match(unique(homesafety[, "study"]), homesafety[, "study"])
id <- data.frame(studyid = unique(homesafety$studyid), study = 1:22, narms = homesafety[id_arms, "narms"], baseline = homesafety[id_arms, "baseline"])
M <- 3
ref_trt <- 1
nc_data <- new("nc_data",
  type = "binary",
  n_study = 22,
  n_outcomes = 3,
  n_datapoints = 45,
  n_treatments = 9,
  study_id = id,
  study_data = homesafety[, 2:12],
  ref_trt = ref_trt)
nGamma <- 1 # nc_data@n_treatments
summary(nc_data)

# inizialization parameters
rng <- 101
set.seed(rng)
prm.init <- list(mu_sigma2 = 10^-1, d_sigma2 = 10^-1, beta_sigma = 10^(-1/8),
  nGamma = nGamma, eta_prior = 1)

# prior parameters
prm.prior <- nc_prior(mu_sigma2 = 10^3, d_sigma2 = 10^3,
  beta_sigma = 10^(-1/8), nGamma = nGamma, eta_prior = 1)

# MCMC settings
burnin <- 250000
nsim <- 50000
nthin <- 100
prm.prop <- list(Gamma.update = "PX-RPMH", eta = .5, sigma.r = .01)
control <- list(burnin = burnin, nsim = nsim,
  Gamma.update = prm.prop[["Gamma.update"]], eta.prop = prm.prop[["eta"]],
  sigma.r.prop = prm.prop[["sigma.r"]], random.start = TRUE, verbose = TRUE,
  thin = nthin, store.burnin = TRUE)
eps <- eps_ab <- 1e-12
xmax <- 1000
xmin <- -xmax
tuning <- list(eps = eps, eps_ab = eps_ab, xmin = xmin, xmax = xmax)
every <- 50
maxiter <- max(25, floor(burnin/every))
adaptation <- list(every = every, maxiter = maxiter, miniter = 5,
  alpha = .75, beta = .8, gamma = 0, tar = 0.234, tol = .01)

# MCMC simulation
seed <- floor(runif(1, 1, 1000))
set.seed(seed)
res <- netcopula(nc_data, control, prm.prior, prm.init, tuning = tuning,
  adaptation = adaptation)

summary(res, include.burnin = TRUE, regex_pars = "d")
summary(res, include.burnin = TRUE, regex_pars = "Sigma")
summary(res, include.burnin = TRUE, regex_pars = "Gamma")

plot(x = res, regex_pars = "Gamma", what = "trace")
plot(x = res, regex_pars = "Sigma", what = "trace")
plot(x = res, regex_pars = "Sigma", what = "hist")
plot(x = res, regex_pars = "Sigma", what = "pairs")
res_mcmc <- nc_mcmc_to_mcmc(res)
pars <- sort(colnames(res_mcmc)[startsWith(colnames(res_mcmc), "d[")])[-(1:3)]
plot(x = res, pars = pars, what = "intervals")
plot(x = res, pars = pars[1:6], what = "trace")
