library(netcopula)

data(homesafety)

id_arms <- match(unique(homesafety[, "study"]), homesafety[, "study"])
id <- data.frame(studyid = unique(homesafety$studyid), study = 1:22, narms = homesafety[id_arms, "narms"], baseline = homesafety[id_arms, "baseline"])
ref_trt <- 1

out <- nc_data_create("binary", 22, 3, 45, 9, NULL, homesafety[, 2:12], NULL)
out
summary(out)

out <- nc_data_create("binary", 22, 3, 45, 9, id, homesafety[, 2:12], ref_trt)
out
summary(out)

###

M <- 3
u <- rep(0.1, M)
Gamma <- c(.3, .5, .2)
Gamma_matrix <- matrix(NA, nrow = M, ncol = M)
Gamma_matrix[upper.tri(Gamma_matrix)] <-
  Gamma_matrix[lower.tri(Gamma_matrix)] <- Gamma
diag(Gamma_matrix) <- 1
all(eigen(Gamma_matrix)$values > 0) # positive definiteness check
rcpp_T <- gausscopdens(u, Gamma_matrix, TRUE, FALSE)
rcpp_F <- gausscopdens(qnorm(u), Gamma_matrix, FALSE, FALSE)

gausscopdens_R <- function(u, Gamma) {
  x <- qnorm(u)
  out <- det(Gamma)^(-.5)*exp((t(x) %*% (diag(M) - solve(Gamma)) %*% x)/2)
  return(as.numeric(out))
}
r <- gausscopdens_R(u, Gamma_matrix)
rcpp_T - r
rcpp_F - r

M <- 2
Gamma <- .3
Gamma_matrix <- matrix(NA, nrow = M, ncol = M)
Gamma_matrix[upper.tri(Gamma_matrix)] <-
  Gamma_matrix[lower.tri(Gamma_matrix)] <- Gamma
diag(Gamma_matrix) <- 1
len <- 100
u1 <- u2 <- seq(.01, .99, length.out = len)
dens <- matrix(NA, nrow = len, ncol = len)
dens_R <- matrix(NA, nrow = len, ncol = len)
for (i in 1:len) {
  for (j in 1:len) {
    dens[i, j] <- gausscopdens(c(u1[i], u2[j]), Gamma_matrix, TRUE, FALSE)
    dens_R[i, j] <- gausscopdens_R(c(u1[i], u2[j]), Gamma_matrix)
  }
}
persp(u1, u2, dens, theta = 120, phi = 25)
persp(u1, u2, dens_R, theta = 120, phi = 25)
plot(as.numeric(dens_R - dens))
abline(h = 0, lty = 2, col = "gray")

###

prova <- new("nc_mcmc")
prova

###

set.seed(123)

# Covariance matrix and mean vector
sigma <- matrix(c(1, 0.9, -0.3, 0.9, 1, -0.4, -0.3, -0.4, 1), ncol = 3)
mu <- c(10, 5, -3)

# Comparison
set.seed(123)
x_arma <- rmvn_arma(1, mu, sigma)
set.seed(123)
x_mvtnorm <- mvtnorm::rmvnorm(1, mu, sigma)
set.seed(123)
x_MASS <- MASS::mvrnorm(1, mu, sigma)
x_arma
x_mvtnorm
x_MASS

# Benchmarking
n <- 1e+4
require(rbenchmark)
benchmark(mvtnorm::rmvnorm(n, mu, sigma),
  MASS::mvrnorm(n, mu, sigma),
  rmvn_arma(n, mu, sigma),
  columns = c("test", "replications", "relative", "elapsed"),
  order = "relative")

###

set.seed(123)
k <- 8
S <- bayesm::rwishart(10, diag(k))$IW
nu <- round(runif(1, min = (k + 1), max = 30))
n <- 1e+2
X <- array(NA, dim = c(n, k, k))
for (i in 1:n) X[i, , ] <- MCMCpack::riwish(nu, S)

MCMCpack::diwish(X[1, , ], nu , S)
dinvwish_arma(X[1, , ], nu, S, TRUE)
ldiw(X[1, , ], nu, S)

arma <- apply(X, 1, dinvwish_arma, nu, S, FALSE)
mcmcpack <- apply(X, 1, MCMCpack::diwish, nu, S)
all.equal(arma, mcmcpack)
summary((arma - mcmcpack)/mcmcpack)

res <- replicate(n, rinvwish_arma(nu, S), simplify = "array")

###

data(homesafety)
id_arms <- match(unique(homesafety[, "study"]), homesafety[, "study"])
id <- data.frame(studyid = unique(homesafety$studyid), study = 1:22, narms = homesafety[id_arms, "narms"], baseline = homesafety[id_arms, "baseline"])
ref_trt <- 1
nc_data <- nc_data_create("binary", 22, 3, 45, 9, id, homesafety[, 2:12], ref_trt)
n_datapoints <- nc_data@n_datapoints
n_study <- nc_data@n_study
M <- nc_data@n_outcomes
n_trt <- nc_data@n_treatments
y <- as.matrix(nc_data@study_data[, c("y1", "y2", "y3")])
n <- as.matrix(nc_data@study_data[, c("n1", "n2", "n3")])
baseline <- nc_data@study_data[, "baseline"]
narms <- nc_data@study_data[, c("narms")]
narms_s <- nc_data@study_id[, 3]
trt <- nc_data@study_data[, "trt"]

rng <- 1406
set.seed(rng)
prm_init <- list(mu_sigma2 = 10^-1, d_sigma2 = 10^-1, beta_sigma = 10^(-1/8), nGamma = nc_data@n_treatments, df = 10^2)
param <- nc_init(nc_data, prm_init, random_start = TRUE, init = NULL)

x <- param$x
Gamma_q <- param$Gamma
mu <- param$mu
mu <- param_long(mu, narms_s, FALSE)
delta <- param$delta
x_imp <- x_imputed(x, Gamma_q, trt)
n_imp <- n_imputed(as.matrix(nc_data@study_data[, paste0("n", 1:M)]))
y_imp <- y_imputed(y, x_imp, narms_s, mu, delta, n_imp)
ncloglik_rcpp <- nc_loglik(y_imp, n_imp, x_imp, trt, mu, delta, Gamma_q)

gausscopdens_R2 <- function(u, Gamma_q, is_u = FALSE, logd = TRUE) {
  if (is_u) {
    x <- qnorm(u)
  } else {
    x <- u
  }
  out <- -.5*(determinant(Gamma_q)$modulus - (t(x) %*% (diag(M) - solve(Gamma_q)) %*% x))
  return(as.numeric(ifelse(logd, out, exp(out))))
}
gcd_rcpp <- gcd_r <- numeric(nc_data@n_datapoints)
ind_a_b <- matrix(NA, nc_data@n_datapoints, M)
for (i in 1:nc_data@n_datapoints) {
  gcd_rcpp[i] <- gausscopdens(x_imp[i, ], Gamma_q[[trt[i]]], FALSE, TRUE)
  gcd_r[i] <- gausscopdens_R2(x_imp[i, ], Gamma_q[[trt[i]]], FALSE, TRUE)
  for (m in 1:M) {
    ind_a_b[i, m] <- indic_a_b(y_imp[i, m], n_imp[i, m], x_imp[i, m], mu[i, m], delta[i, m])
  }
}
ncloglik_r <- sum(gcd_r, na.rm = TRUE)
all.equal(gcd_rcpp, gcd_r)
all.equal(ncloglik_rcpp, ncloglik_r)

###

i <- 1
m <- 1
y_ikm <- y[i, m]
n_ik <- n_ik[i, m]
x_ikm <- x[i, m]
mu_im <- mu[i, m]
delta_ikm <- delta[i, m]
is_baseline <- TRUE
# indic(y_ikm, n_ik, x_ikm, mu_im, delta_ikm, is_baseline)

expit <- function(x) 1/(1 + exp(-x))
theta <- as.numeric(mu_im + delta_ikm*(1 - is_baseline))
pi <- expit(theta)
a_ikm <- as.numeric(pbinom(y_ikm - 1, n_ik, pi, 1, 0))
b_ikm <- as.numeric(pbinom(y_ikm, n_ik, pi, 1, 0))
q_ikm <- as.numeric(pnorm(x_ikm, 0, 1, 1, 0))
c(theta = theta, pi = pi, a_ikm = a_ikm, b_ikm = b_ikm, q_ikm = q_ikm)

###

rlkj_R <- function(K, eta) {
  CPCs <- numeric(K*(K - 1)/2)
  alpha <- eta + 0.5 * (K - 1)
  count <- 1
  for (i in 1:(K - 1)) {
    alpha <- alpha - 0.5
    for (j in (i + 1):K) {
      CPCs[count] <- 2*rbeta(1, alpha, alpha) - 1
      count <- count + 1
    }
  }

  temp <- numeric(K - 1)
  acc <- rep(1, K - 1)
  # Cholesky factor of correlation matrix
  L <- matrix(0, nrow = K, ncol = K)

  position <- 1
  pull <- (K - 1)

  L[1, 1] <- 1.0
  L[tail(1:K, pull), 1] <- temp <- CPCs[head(1:(K*(K - 1)/2), pull)]
  acc[tail(1:(K - 1), pull)] <- 1 - temp^2
  for (i in 2:(K - 1)) {
    position <- position + pull
    pull <- pull - 1
    temp <- CPCs[position:(position + pull - 1)]
    L[i, i] <- sqrt(acc[i - 1])
    L[tail(1:K, pull), i] <- temp * sqrt(acc[tail(1:(K - 1), pull)])
    acc[tail(1:(K - 1), pull)] <- acc[tail(1:(K - 1), pull)] * (1 - temp^2)
  }
  L[K, K] <- sqrt(acc[K - 1])

  return (L %*% t(L))
}

K <- 5
eta <- 1
rng <- round(runif(1, 1, 10000))
set.seed(rng)
rccp <- rlkj_arma(K, eta)
set.seed(rng)
r <- rlkj_R(K, eta)
all.equal(rccp, r)

###

library(truncnorm)

n <- 20
a <- 0
b <- Inf
mean <- 10
sd <- 5
rng <- round(runif(1, 1, 10000))
set.seed(rng)
rcpp <- rtruncnorm_rcpp(n, a, b, mean, sd)
set.seed(rng)
r <- rtruncnorm(n, a, b, mean, sd)
identical(rcpp, r)

a <- 0
b <- Inf
mean <- 10
sd <- 1
rng <- round(runif(1, 1, 10000))
set.seed(rng)
n <- 20
x <- rtruncnorm_rcpp(n, a, b, mean, sd)
rcpp <- dtruncnorm_rcpp(x, a, b, mean, sd)
r <- dtruncnorm(x, a, b, mean, sd)
identical(rcpp, r)

###

M <- 5
sigma_r <- 10^-(1/3)
rng <- round(runif(1, 1, 10000))
set.seed(rng)
Sigma <- rlogchol_arma(M, sigma_r)
all(eigen(Sigma)$values > 0)
Sigma_sqrt <- solve(diag(diag(Sigma)))^0.5
Sigma_sqrt %*% Sigma %*% Sigma_sqrt

###

dlogchol <- function(Sigma, sigma_r, logd = FALSE) {
  Sigma_inv <- solve(Sigma)
  R <- chol(Sigma_inv)
  diag(R) <- log(diag(R))
  beta <- R[upper.tri(R, diag = TRUE)]
  out <- sum(dnorm(beta, 0, sigma_r, log = TRUE))
  if (!logd) out <- exp(out)
  return(out)
}

M <- 5
sigma_r <- 10^-(1/3)
rng <- round(runif(1, 1, 10000))
set.seed(rng)
Sigma <- rlogchol_arma(M, sigma_r)
logd <- TRUE
rcpp <- dlogchol_arma(Sigma, sigma_r, logd)
r <- dlogchol(Sigma, sigma_r, logd)
all.equal(rcpp, r)

###

data(homesafety)
id_arms <- match(unique(homesafety[, "study"]), homesafety[, "study"])
id <- data.frame(studyid = unique(homesafety$studyid), study = 1:22, narms = homesafety[id_arms, "narms"], baseline = homesafety[id_arms, "baseline"])
M <- 3
ref_trt <- 1
nc_data <- nc_data_create("binary", 22, 3, 45, 9, id, homesafety[, 2:12], ref_trt)
y <- nc_data@study_data[, 10:12]
trt <- nc_data@study_data[, 4]
n_trt <- nc_data@n_treatments
n_out <- nc_data@n_outcomes
n_data <- nc_data@n_datapoints
eta <- 1
rng <- 1406
set.seed(rng)
Gamma <- lapply(rep(eta, n_trt), function(prm) rlkj_arma(M, prm))

rng <- 123
set.seed(rng)
prm_init <- list(mu_sigma2 = 10^-1, d_sigma2 = 10^-1, beta_sigma = 10^(-1/8), nGamma = nc_data@n_treatments, df = 10^2)
param <- nc_init(nc_data, prm_init, random_start = TRUE, init = NULL)
x <- param$x

# R implementation
mat_block_diag_R <- function(x, n) {
  size <- nrow(x)
  y <- matrix(0, nrow = size*n, ncol = size*n)
  for (q in 1:n) {
    y[(size*(q - 1) + 1):(size*q), (size*(q - 1) + 1):(size*q)] <- x
  }
  return(y)
}

x_imputed_R <- function(x, Gamma, trt) {
  n_data <- nrow(x)
  n_out <- ncol(x)
  s2t_split <- split_iv(1:n_data, trt)
  n_split <- sapply(s2t_split, length)
  x2t_split <- split_nm(x, trt)
  n_q <- lapply(x2t_split, nrow)
  x_q_vec <- lapply(x2t_split, function(x) as.numeric(t(x)))
  x_q_vec_obs <- lapply(x_q_vec, function(x) x[!is.na(x)])
  n_q_vec <- lapply(x_q_vec, length)
  n_q_obs <- lapply(x_q_vec_obs, length)
  n_q_miss <- mapply(function(x, y) x - y, n_q_vec, n_q_obs, SIMPLIFY = FALSE)
  x_q_vec_miss_i <- lapply(x_q_vec, function(x) which(is.na(x)))
  x_q_vec_obs_i <- lapply(x_q_vec, function(x) which(!is.na(x)))
  Gamma_q_tilde <- mapply(mat_block_diag_R, Gamma, n_q)
  Gamma_q_miss <- mapply(function(x, miss_i) x[miss_i, miss_i, drop = FALSE], Gamma_q_tilde, x_q_vec_miss_i)
  Gamma_q_obs <- mapply(function(x, obs_i) x[obs_i, obs_i, drop = FALSE], Gamma_q_tilde, x_q_vec_obs_i)
  Gamma_q_21 <- mapply(function(x, obs_i, miss_i) x[obs_i, miss_i, drop = FALSE], Gamma_q_tilde, x_q_vec_obs_i, x_q_vec_miss_i)
  Gamma_q_12 <- lapply(Gamma_q_21, t)
  mu <- mapply(function(x, y, z) x%*%solve(y)%*%z, Gamma_q_12, Gamma_q_obs, x_q_vec_obs)
  Omega <- mapply(function(x, y, z, w) x - y%*%solve(z)%*%w, Gamma_q_miss, Gamma_q_12, Gamma_q_obs, Gamma_q_21)
  x_q_vec_miss <- mapply(function(mu, sigma) as.numeric(rmvn_arma(1, mu, sigma)), mu, Omega)
  x_q_vec_imp <- x_q_vec
  for (q in 1:n_trt) {
    x_q_vec_imp[[q]][is.na(x_q_vec_imp[[q]])] <- x_q_vec_miss[[q]]
  }
  x_q_imp <- lapply(x_q_vec_imp, function(x) matrix(x, ncol = n_out, byrow = TRUE))
  x_imp <- x
  for (q in 1:n_trt) {
    x_imp[s2t_split[[q]], ] <- x_q_imp[[q]]
  }
  return(x_imp)
}

rng <- 1406
set.seed(rng)
x_imp_R <- x_imputed_R(x, Gamma, trt)

# Rcpp implementation
set.seed(rng)
x_imp <- x_imputed(x, Gamma, trt)

all.equal(x_imp_R, x_imp)
print(cbind(x, x_imp_R, x_imp), digits = 5)

# Simulated data
set.seed(rng)
sigma <- bayesm::rwishart(10, diag(n_out))$IW
means <- rnorm(n_out)
n <- 1e+3
X <- mvtnorm::rmvnorm(n_data, means, sigma)
p <- .2
for (i in 1:n_data) {
  for (j in 1:n_out) {
    if (rbinom(n = 1, size = 1, prob = p)) {
      X[i, j] <- NA
    }
  }
}

# R implementation
set.seed(rng)
x_imp_R <- x_imputed_R(X, Gamma, trt)

# Rcpp implementation
set.seed(rng)
x_imp <- x_imputed(X, Gamma, trt)

x_imp_R
x_imp
X
all.equal(x_imp_R, x_imp)

# Microbenchmarking
set.seed(rng)
require(microbenchmark)
microbenchmark(x_imputed_R(X, Gamma, trt), x_imputed(X, Gamma, trt), times = 1e+3)

# Benchmarking
require(rbenchmark)
benchmark(x_imputed_R(X, Gamma, trt), x_imputed(X, Gamma, trt),
  order = "relative", replications = 1e+3)[, 1:4]

###

data(homesafety)
id <- data.frame(studyid = unique(homesafety$studyid), study = 1:22)
nc_data <- nc_data_create("binary", 22, 3, 45, 9, id, homesafety[, 2:12])
y <- nc_data@study_data[, 10:12]
x <- as.matrix(y)
trt <- nc_data@study_data[, 4]
n_trt <- nc_data@n_treatments
n_out <- nc_data@n_outcomes
n_data <- nc_data@n_datapoints

x_split_1 <- lapply(split(as.data.frame(x), trt), as.matrix)
x_split_2 <- split_nm(x, trt)

# Benchmarking
require(rbenchmark)
benchmark(lapply(split(as.data.frame(x), trt), as.matrix),
  split_nm(x, trt),
  order = "relative", replications = 1e+3)[, 1:4]

###

library(netcopula)

par <- c(-7, 6)
data(cancermortality, package = "LearnBayes")
options <- list(trace = 0, fnscale = -1, parscale = c(1, 1), ndeps = c(0.001, 0.001), maxit = 500, abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, gamma = 2)
out_rcpp <- optim_rcpp(par, betabinexch, list(cancermortality), options, hessian = TRUE)
out_r <- optim(par, betabinexch, gr = NULL, list(cancermortality), control = options, hessian = TRUE)
identical(out_rcpp$par, out_r$par)
identical(out_rcpp$hessian, out_r$hessian)

# x2 <- function(x, varargs) {
#   a <- varargs[[1]]
#   b <- varargs[[2]]
#   a[1]*x[1]^2 + a[2]*x[2]^2 + b
# }
# a <- c(-1, -1)
# b <- 5
# optim_rcpp(par, x2, list(a, b), options)
# optim(par, x2, gr = NULL, list(a, b), control = options, hessian = TRUE)

out_rcpp <- laplace_rcpp(betabinexch, par, list(cancermortality), options)
out_r <- laplace(betabinexch, par, list(cancermortality))
identical(out_rcpp$mode, out_r$mode)
all.equal(out_rcpp$var, out_r$var)
identical(out_rcpp$int, out_r$int)

fit <- laplace_rcpp(betabinexch, par, list(cancermortality), options)
tpar <- list(m = fit$mode, var = 2*fit$var, df = 4)
datapar <- list(data = cancermortality, par = tpar)
start <- c(-6.9, 12.4)
laplace_rcpp(betabinT, start, datapar, options)
laplace(betabinT, start, datapar)

d <- list(int.lo = c(-Inf, seq(66, 74, by = 2)), int.hi = c(seq(66, 74, by= 2 ), Inf), f = c(14, 30, 49, 70, 33, 15))
start <- c(70, 1)
laplace_rcpp(groupeddatapost, start, list(d), options)
laplace(groupeddatapost, start, list(d))

data(stanfordheart, package = "LearnBayes")
start <- c(0, 3, -1)
npar <- length(start)
options <- list(trace = 0, fnscale = -1, parscale = rep(1, npar), ndeps = rep(0.001, npar), maxit = 500, abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, gamma = 2)
laplace_rcpp(transplantpost, start, list(stanfordheart), options)
laplace(transplantpost, start, list(stanfordheart))

data(chemotherapy, package = "LearnBayes")
start <- c(-.5, 9, .5, -.05)
attach(chemotherapy)
d <- cbind(time, status, treat - 1, age)
npar <- length(start)
options <- list(trace = 0, fnscale = -1, parscale = rep(1, npar), ndeps = rep(0.001, npar), maxit = 500, abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, gamma = 2)
laplace_rcpp(weibullregpost, start, list(d), options)
laplace(weibullregpost, start, list(d))
detach(chemotherapy)

###

expit <- function(x) 1/(1 + exp(-x))
mu <- 0
delta <- 0
tau <- 1
eta <- 2
y <- 12
n <- 20
w <- 1
gamma <- 1
mu_sigma <- 10
eps <- 0
eps_ab <- 1e-3
mu_delta_logpost(mu, delta, tau, eta, y, n, w, gamma, mu_sigma, eps, eps_ab)
theta <- mu + delta
p <- expit(theta)
a <- pbinom(y - 1, n, p, 1, 1)
b <- pbinom(y, n, p, 1, 1)
log(pnorm((qnorm(b, log.p = TRUE) - w)/gamma, 0.0, 1.0, 1, 0) - pnorm((qnorm(a, log.p = TRUE) - w)/gamma, 0.0, 1.0, 1, 0)) + dnorm(delta, tau, eta, 1) + dnorm(mu, 0, mu_sigma, 1)
prm <- list(tau = tau, eta = eta, y = y, n = n, w = w, gamma = gamma, mu_sigma = mu_sigma, eps = eps, eps_ab = eps_ab)
mu_delta_logpost_func(c(mu, delta), prm)
limits <- c(-2, 3, -2, 3)
my_plot(mu_delta_logpost_func, limits, prm, type = "contour", npoints = 100, nlevels = 30, xlab = "mu", ylab = "delta")
abline(v = optim_rcpp_example()$par[1], lty = 2)
abline(h = optim_rcpp_example()$par[2], lty = 2)
my_plot(mu_delta_logpost_func, limits, prm, type = "persp", theta = 20, phi = 20, r = 50, ticktype = "detailed", xlab = "mu", ylab = "delta", zlab = "log posterior")

par <- c(-7, 6)
args <- list(tau = tau, eta = eta, y = y, n = n, w = w, gamma = gamma, mu_sigma = mu_sigma, eps = eps, eps_ab = eps_ab)
options <- list(trace = 0, fnscale = -1, parscale = c(1, 1), ndeps = c(0.001, 0.001), maxit = 500, abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, gamma = 2)
optim_rcpp(par, mu_delta_logpost_func, args, options, hessian = TRUE)
optim(par, mu_delta_logpost_func, gr = NULL, args, control = options, hessian = TRUE)

laplace_rcpp(mu_delta_logpost_func, par, args, options)
laplace(mu_delta_logpost_func, par, args)

optim_rcpp_example()

###

tau <- -0.78061
eta <- 0.37616
y_tmp <- 57
n_tmp <- 73
w <- 14.38137
gamma <- 0.51407
mu_sigma <- sqrt(10^3)
eps <- 0 #.Machine$double.eps
eps_ab <- 1e-3

mu <- 0
delta <- 0
mu_delta_logpost(mu, delta, tau, eta, y_tmp, n_tmp, w, gamma, mu_sigma, eps, eps_ab)

args <- list(tau = tau, eta = eta, y = y_tmp, n = n_tmp, w = w, gamma = gamma, mu_sigma = mu_sigma, eps = eps, eps_ab = eps_ab)
par <- c(mu, delta)
options <- list(trace = 0, fnscale = -1, parscale = c(1, 1), ndeps = c(0.001, 0.001), maxit = 500, abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, gamma = 2)
# optim_rcpp(par, mu_delta_logpost_func, args, options, hessian = TRUE)
res <- laplace_rcpp(mu_delta_logpost_func, par, args, options)
limits <- c(-10, 18, -10, 18)
my_plot(mu_delta_logpost_func, limits, args, type = "contour", npoints = 100, nlevels = 30, xlab = "mu", ylab = "delta")
abline(v = res$mode[1], lty = 2)
abline(h = res$mode[2], lty = 2)

###

# expit <- function(x) 1/(1 + exp(-x))
# mu <- -9.8515598252e+00
# delta <- 1.1904749185e+00
# p <- expit(mu + delta)
# tau <- 0.0000000000e+00
# eta <- 9.9308943505e-01
# y <- 0
# n <- 57
# w <- -4.1077310348e+01
# gamma <- 5.3654103603e-01
# mu_sigma <- sqrt(prm_prior$mu_sigma2)
# eps <- tuning$eps
# eps_ab <- tuning$eps_ab
# # mu_delta_logpost(mu, delta, tau, eta, y, n, w, gamma, mu_sigma, eps, eps_ab)
#
# a <- pbinom(y - 1, n, p, 1, 0)
# b <- pbinom(y, n, p, 1, 0)
#
# if (a >= (1 - eps_ab)) {
#   a <- 1 - eps_ab
# } else if (a <= eps_ab) {
#   a <- eps_ab
# }
# if (b >= (1 - eps_ab)) {
#   b <- 1 - eps_ab
# } else if (b <= eps_ab) {
#   b <- eps_ab
# }
# if (a == b) {
#   a <- a - eps
# }
# phi_inv_b <- qnorm(b, 0.0, 1.0, 1, 0)
# phi_inv_a <- qnorm(a, 0.0, 1.0, 1, 0)
# tmp <- log(pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 1, 0) - pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 1, 0))
# if (!is.finite(tmp)) {
#   tmp <- log(pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 0, 0) - pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 0, 0))
# }
# if (!is.finite(tmp)) {
#   tmp_a <- pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 1, 1)
#   tmp_b <- pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 1, 1)
#   tmp <- log(exp(tmp_b - tmp_a) - 1) + tmp_a
# }
# if (!is.finite(tmp)) {
#   tmp_div <- 700
#   tmp_a <- pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 0, 1)
#   tmp_b <- pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 0, 1)
#   tmp_diff <- tmp_a - tmp_b
#   k <- ifelse(tmp_diff > tmp_div, tmp_diff - tmp_div, 0)
#   tmp <- k + log(exp(tmp_div) - exp(-k)) + tmp_b
# }
# tmp + dnorm(delta, tau, eta, 1) + dnorm(mu, 0, mu_sigma, 1)

###

set.seed(123)
n <- 10
p <- 3
x <- rnorm(n*p)
X <- matrix(x, nrow = n, ncol = p)
M <- matrix(0, nrow = n, ncol = p)
U <- diag(n)
V <- diag(p)
dmatvn_arma(X, M, U, V)

set.seed(123)
n <- 10
p <- 3
M <- matrix(0, nrow = n, ncol = p)
U <- diag(n)
V <- diag(p)
rmatvn_arma(M, U, V)
set.seed(123)
rmvn_arma(n, rep(0, p), diag(p))

###

set.seed(123)
n <- 10
p <- 3
x <- rnorm(n*p)
X <- matrix(x, nrow = n, ncol = p)
M <- matrix(0, nrow = n, ncol = p)
Sigma <- diag(n)
Omega <- diag(p)
df <- 5
dmatvt_arma(X, M, Sigma, Omega, df)

x <- t(as.matrix(as.numeric(X)))
mean  <- as.numeric(M)
sigma <- kronecker(Omega, Sigma)
dmvt_arma(x, mean, sigma, df)

set.seed(123)
n <- 10
p <- 3
x <- rnorm(n*p)
M <- matrix(0, nrow = n, ncol = p)
Sigma <- diag(n)
Omega <- diag(p)
df <- 5
X <- rmatvt_arma(M, Sigma, Omega, df)
var(X)

###

set.seed(123)
n <- 4
p <- 3
nu <- round(runif(1, min = (n + 1), max = 30))
S <- bayesm::rwishart(10, diag(p))$IW
Omega <- solve(rinvwish_arma(nu, S))
Sigma <- diag(n)
kronecker(Sigma, Omega)
df <- 7
m <- matrix(0, n, p)
d <- rmatvt_arma(m, Sigma, Omega, df)
var(d)
var(t(d))

###

set.seed(123)

totiter <- 100
d_sigma2 <- 10^-1
ref_trt <- 1
theta_d <- array(NA, dim = c(n_trt, M, totiter))
theta_d_reshaped <- matrix(NA, nrow = totiter, ncol = M*(n_trt - 1))
for (i in 1:totiter) {
  d <- rmvn_arma(n_trt, rep(0, M), d_sigma2*diag(M))
  d[ref_trt, ] <- 0
  theta_d[, , i] <- d
  theta_d_reshaped[i, ] <- as.numeric(d[-ref_trt, ])
}

tmp <- cube_to_mat(theta_d, TRUE, ref_trt)

mu <- colMeans(theta_d_reshaped)
rho <- 1
cov <- bayesm::rwishart(30, diag(M*(n_trt - 1)))$W
ar <- c(0.3, .8)
alpha <- .75
beta <- .8
gamma <- 0
tar <- 0.234
k <- 3
adpt <- rwmh_adapt(theta_d_reshaped, mu, rho, cov, ar, alpha, beta, gamma, tar, k, FALSE)

###

d_logpost_R <- function(d, delta, Sigma_M, trt, baseline, narms_study, d_sigma, ref_trt) {
  n_datapoints <- nc_data@n_datapoints
  tmp <- d_prior <- 0
  count_trt <- 1
  for (ik in 1:n_datapoints) {
    if (trt[ik] != baseline[ik]) {
      if (narms_study[ik] == 2) {
        d_base <- d[baseline[ik], ]
        d_trt <- d[trt[ik], ]
        d_k <- d_trt - d_base
        delta_k <- delta[ik, ]
        tmp <- tmp + as.numeric(dmvn_arma(t(as.matrix(delta_k)), d_k, Sigma_M, TRUE))
      } else {
        if (count_trt == 1) {
          d_base <- d[baseline[ik], ]
          d_trt <- d[trt[ik], ]
          d_k <- d_trt - d_base
          delta_k <- delta[ik, ]
          tmp <- tmp + as.numeric(dmvn_arma(t(as.matrix(delta_k)), d_k, Sigma_M, TRUE))
        } else {
          d_base <- d[baseline[ik], ]
          d_trt <- d[trt[ik], ]
          d_k <- d_trt - d_base
          for (k in 1:(count_trt - 1)) {
            d_base <- d[baseline[ik], ]
            d_trt <- d[trt[ik - k], ]
            d_k_tmp <- d_trt - d_base
            delta_k <- delta[ik - k, ]
            d_k <- d_k + (delta_k - d_k_tmp)/count_trt
          }
          Sigma_M_k <- Sigma_M*(count_trt + 1)/(2*count_trt)
          delta_k <- delta[ik, ]
          tmp <- tmp + as.numeric(dmvn_arma(t(as.matrix(delta_k)), d_k, Sigma_M_k, TRUE))
        }
        count_trt <- count_trt + 1
        if (count_trt == narms_study[ik]) {
          count_trt <- 1
        }
      }
      # print(tmp, digits = 8)
    }
  }

  for (k in 1:n_trt) {
    if (k != ref_trt) {
      for (m in 1:M) {
        d_prior <- d_prior + dnorm(d[k, m], 0, d_sigma, 1)
      }
    }
  }

  tmp + d_prior
}

library(netcopula)

data(homesafety)
id_arms <- match(unique(homesafety[, "study"]), homesafety[, "study"])
id <- data.frame(studyid = unique(homesafety$studyid), study = 1:22, narms = homesafety[id_arms, "narms"], baseline = homesafety[id_arms, "baseline"])
M <- 3
ref_trt <- 1
nc_data <- nc_data_create("binary", 22, 3, 45, 9, id, homesafety[, 2:12], ref_trt)
n_trt <- nc_data@n_treatments

rng <- 123
set.seed(rng)
prm_init <- list(mu_sigma2 = 10^-1, d_sigma2 = 10^-1, beta_sigma = 10^(-1/8), nGamma = nc_data@n_treatments, df = 10^2)
param <- nc_init(nc_data, prm_init, random_start = TRUE, init = NULL)
d <- param$d
delta <- param$delta
Sigma_M <- param$Sigma_M
trt <- nc_data@study_data[, "trt"]
baseline <- nc_data@study_data[, "baseline"]
narms_study <- nc_data@study_data[, "narms"]
d_sigma <- sqrt(prm_init$d_sigma2)
ref_trt <- nc_data@ref_trt
d_logpost(d, delta, Sigma_M, trt, baseline, narms_study, d_sigma, ref_trt)
d_logpost_R(d, delta, Sigma_M, trt, baseline, narms_study, d_sigma, ref_trt)

delta_tmp <- delta
delta_tmp[22, ] <- as.numeric(rmvn_arma(1, rep(0, M), Sigma_M))
trt_tmp <- trt
trt_tmp[22:23] <- c(7, 9)
baseline_tmp <- baseline
narms_study_tmp <- narms_study
narms_study_tmp[19:23] <- 5
d_logpost(d, delta_tmp, Sigma_M, trt_tmp, baseline_tmp, narms_study_tmp, d_sigma, ref_trt)
d_logpost_R(d, delta_tmp, Sigma_M, trt_tmp, baseline_tmp, narms_study_tmp, d_sigma, ref_trt)

###

library(corpcor)

# A <- matrix(c(0.643421821, -0.643422058, -0.643422058, 0.643422058), nrow = 2)
# A <- matrix(c(0.487058814, -0.487059067, -0.487059067, 0.487059067), nrow = 2)
# A <- matrix(c(29.359141667, -0.456688339, -0.456688339, 0.470286998), nrow = 2)
# A <- matrix(c(0.000000000, 0.000000000, 0.605415384, 0.000000000), nrow = 2)
# A <- matrix(c(2.6878931861, 3.0070342581, -1.0701050045, 3.0070341246, 3.3640693773, -1.1997664551, -1.0701049569, -1.1997664551, 4.4999602767), nrow = 3)
A <- matrix(c(2.8346e+03, -2.8139e+01, -1.4683e+03, -2.8139e+01, 2.9739e+03, 1.0831e+03, -1.4683e+03, 1.0831e+03, 2.4957e+03), nrow = 3)

is.positive.definite(A)
dA <- make.positive.definite(A)
is.positive.definite(dA)

is_positive_definite(A)
dA_arma <- make_positive_definite(A)
is_positive_definite(dA_arma)

all.equal(dA, dA_arma)

eval <- eigen(A, only.values = TRUE, symmetric = TRUE)$values
tol <- max(dim(A))*max(abs(eval))*.Machine$double.eps
sum(eval > tol)

###

num <- 1e+4
all_equal <- logical(num)
for (i in 1:num) {
  prm_init <- list(mu_sigma2 = 10^-1, d_sigma2 = 10^-1, beta_sigma = 10^(-1/8), nGamma = nGamma, df = 10^2)
  param <- nc_init(nc_data, prm_init, random_start = TRUE, init = NULL)

  delta <- param$delta
  d <- param$d
  Sigma_M <- param$Sigma_M
  trt_arms <- nc_data@study_data$trt
  baseline <- nc_data@study_id$baseline
  baseline_arms <- nc_data@study_data[, "baseline"]
  narms <- nc_data@study_id$narms
  narms_study <- nc_data@study_data[, "narms"]
  d_sigma <- sqrt(prm_init$d_sigma2)
  ref_trt <- nc_data@ref_trt

  # delta_logprior(delta, d, Sigma_M, trt_arms, baseline, narms)
  old <- d_logpost(d, delta, Sigma_M, trt_arms, baseline_arms, narms_study, d_sigma, ref_trt)
  new <- d_logpost_multi(d, delta, Sigma_M, trt_arms, baseline, narms, d_sigma, ref_trt)
  all_equal[i] <- all.equal(old, new)
}
all(all_equal)

###

library(netcopula)

# data loading and preparation
data(homesafety)
id_arms <- match(unique(homesafety[, "study"]), homesafety[, "study"])
id <- data.frame(studyid = unique(homesafety$studyid), study = 1:22, narms = homesafety[id_arms, "narms"], baseline = homesafety[id_arms, "baseline"])
M <- 3
ref_trt <- 1
nc_data <- nc_data_create("binary", 22, 3, 45, 9, id, homesafety[, 2:12], ref_trt)
nGamma <- 1 # nc_data@n_treatments
summary(nc_data)

# inizialization parameters
rng <- 1406
set.seed(rng)
prm_init <- list(mu_sigma2 = 10^-1, d_sigma2 = 10^-1, beta_sigma = 10^(-1/8), nGamma = nGamma, df = 10^2)
param <- nc_init(nc_data, prm_init, random_start = TRUE, init = NULL)

# prior parameters
prm_prior <- list(mu_sigma2 = 10^3, d_sigma2 = 10^3, beta_sigma = 10^(-1/8), nGamma = nGamma, df = 10^2)

# MCMC settings
burnin <- 10000
nsim <- 20000
eps <- eps_ab <- 1e-12
tuning <- list(eps = eps, eps_ab = eps_ab)
every <- 10
maxiter <- max(25, floor(burnin/every))
adaptation <- list(every = every, maxiter = maxiter, miniter = 5, alpha = .75, beta = .8, gamma = 0, tar = 0.234, tol = .01)
prm_prop <- list()

# MCMC simulation
seed <- floor(runif(1, 1, 1000))
seed <- 93
set.seed(seed)
res <- netcopula(nc_data, burnin = burnin, nsim = nsim, prm_prior, prm_prop, prm_init, tuning = tuning, adaptation = adaptation, verbose = TRUE)

## MCMC diagnostics
par(mfrow = c(4, 2), mar = c(2.5, 3, 2, 2))
plot(res@Gamma[1, 1, ], type = "l", main = "Gamma")
plot(as.numeric(res@Sigma[, 6]), type = "l", main = "Sigma_M")
plot(res@mu[2, 1, ], type = "l", main = "mu")
plot(res@delta[2, 1, ], type = "l", main = "delta")
plot(res@d[2, 1, ], type = "l", main = "d")
plot(res@x[10, 2, ], type = "l", main = "x")
plot(res@dens$loglik, type = "l", main = "loglik")

###

mu_delta_logpost_r <- function(mu, delta, tau, eta, y, n, w, gamma, mu_sigma, eps, eps_ab) {
  theta <- mu + delta
  p <- expit_rcpp(theta)

  a <- pbinom(y - 1, n, p, 1, 0)
  b <- pbinom(y, n, p, 1, 0)
  if (a > (1 - eps_ab)) {
    a <- 1 - eps_ab
  } else if (a < eps_ab) {
    a <- eps_ab
  }
  if (b > (1 - eps_ab)) {
    b <- 1 - eps_ab
  } else if (b < eps_ab) {
    b <- eps_ab
  }
  eps_a <- a/2.0
  eps_b <- (1.0 - b)/2.0
  if (a == b) {
    if (a == 0) {
      b <- b + eps
    } else if (b == 1) {
      a <- a - eps
    } else {
      a <- a - min(eps/2.0, eps_a)
      b <- b + min(eps/2.0, eps_b)
    }
  }
  phi_inv_a <- qnorm(a, 0.0, 1.0, 1, 0)
  phi_inv_b <- qnorm(b, 0.0, 1.0, 1, 0)
  tmp <- log(pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 1, 0) - pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 1, 0))
  if (!is.finite(tmp)) {
    tmp <- log(pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 0, 0) - pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 0, 0))
  }
  if (!is.finite(tmp)) {
    tmp_a <- pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 1, 1)
    tmp_b <- pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 1, 1)
    tmp <- log(exp(tmp_b - tmp_a) - 1) + tmp_a
  }
  if (!is.finite(tmp)) {
    tmp_a <- pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 0, 1)
    tmp_b <- pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 0, 1)
    tmp <- log(exp(tmp_a - tmp_b) - 1) + tmp_b
  }
  # if (!is.finite(tmp)) {
  #   tmp_div <- 700
  #   tmp_a <- pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 0, 1)
  #   tmp_b <- pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 0, 1)
  #   tmp_diff <- tmp_a - tmp_b
  #   k <- ifelse((tmp_diff > tmp_div), (tmp_diff - tmp_div), 0)
  #   tmp <- k + log(exp(tmp_div) - exp(-k)) + tmp_b
  # }

  if (!is.finite(tmp)) {
    return(NA)
  }

  lpost <- tmp + dnorm(delta, tau, eta, 1) + dnorm(mu, 0, mu_sigma, 1)

  return(lpost)
}

par <- c(0, 0)
# par <- c(7.2500213125, 0.00000000)
tau <- 0
eta <- 1.1546390022
y <- 200
n <- 469
w <- 0
gamma <- 1
# y <- 469
# n <- 469
# w <- -17.2256852690
# gamma <- 0.0007461466
mu_sigma <- sqrt(prm_prior$mu_sigma2)
eps <- eps_ab <- 1e-12

args <- list(tau = tau, eta = eta, y = y, n = n, w = w, gamma = gamma, mu_sigma = mu_sigma, eps = eps, eps_ab = eps_ab)
options <- list(trace = 0, fnscale = -1, parscale = c(1, 1), ndeps = c(0.001, 0.001), maxit = 500, abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, gamma = 2)

laplace_rcpp(mu_delta_logpost_func, par, args, options)

res <- optim_rcpp(par, mu_delta_logpost_func, args, options, hessian = TRUE)
optim(par, mu_delta_logpost_func, gr = NULL, args, control = options, hessian = TRUE)

limits <- c(-1, 1, -1, 1)
my_plot(mu_delta_logpost_func, limits, args, type = "persp", npoints = 50, xlab = "mu", ylab = "delta", theta = 50, phi = 10)
my_plot(mu_delta_logpost_func, limits, args, type = "contour", npoints = 100, xlab = "mu", ylab = "delta", nlevels = 10)
abline(v = res$par[1], lty = 2)
abline(h = res$par[2], lty = 2)

mu_delta_logpost(res$par[1], res$par[2], tau, eta, y, n, w, gamma, mu_sigma, eps, eps_ab)
mu_delta_logpost_r(res$par[1], res$par[2], tau, eta, y, n, w, gamma, mu_sigma, eps, eps_ab)

mu <- 0.0119746762
delta <- -2.5879706424
tau <- -2.5877457911
eta <- 0.6787593925
y <- 0
n <- 73
w <- 10.1843580495
gamma <- 0.2463625031
mu_sigma <- sqrt(prm_prior$mu_sigma2)
eps <- eps_ab <- 0
mu_delta_logpost(mu, delta, tau, eta, y, n, w, gamma, mu_sigma, eps, eps_ab)

args <- list(tau = tau, eta = eta, y = y, n = n, w = w, gamma = gamma, mu_sigma = mu_sigma, eps = eps, eps_ab = eps_ab)

res <- optim_rcpp(par, mu_delta_logpost_func, args, options, hessian = TRUE)
optim(par, mu_delta_logpost_func, gr = NULL, args, control = options, hessian = TRUE)

limits <- c(-100, 10, -20, 20)
my_plot(mu_delta_logpost_func, limits, args, type = "persp", npoints = 50, xlab = "mu", ylab = "delta", theta = 30, phi = 10)
my_plot(mu_delta_logpost_func, limits, args, type = "contour", npoints = 100, xlab = "mu", ylab = "delta", nlevels = 10)
abline(v = res$par[1], lty = 2)
abline(h = res$par[2], lty = 2)

mu_delta_r <- function(mu, delta, tau, eta, y, n, w, gamma, mu_sigma, eps, eps_ab) {
  theta <- mu + delta
  p <- expit_rcpp(theta)

  a <- pbinom(y - 1, n, p, 1, 0)
  b <- pbinom(y, n, p, 1, 0)
  if (a > (1 - eps_ab)) {
    a <- 1 - eps_ab
  } else if (a < eps_ab) {
    a <- eps_ab
  }
  if (b > (1 - eps_ab)) {
    b <- 1 - eps_ab
  } else if (b < eps_ab) {
    b <- eps_ab
  }
  eps_a <- a/2.0
  eps_b <- (1.0 - b)/2.0
  if (a == b) {
    if (a == 0) {
      b <- b + eps
    } else if (b == 1) {
      a <- a - eps
    } else {
      a <- a - min(eps/2.0, eps_a)
      b <- b + min(eps/2.0, eps_b)
    }
  }
  phi_inv_a <- qnorm(a, 0.0, 1.0, 1, 0)
  phi_inv_b <- qnorm(b, 0.0, 1.0, 1, 0)
  tmp <- pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 1, 0) - pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 1, 0)

  if (!is.finite(tmp)) {
    return(NA)
  }

  return(tmp)
}

mu_delta_func <- function(theta, varargs) {
  mu <- theta[1]
  delta <- theta[2]
  tau <- varargs$tau
  eta <- varargs$eta
  y <- varargs$y
  n <- varargs$n
  w <- varargs$w
  gamma <- varargs$gamma
  mu_sigma <- varargs$mu_sigma
  eps <- varargs$eps
  eps_ab <- varargs$eps_ab

  return(mu_delta_r(mu, delta, tau, eta, y, n, w, gamma, mu_sigma, eps, eps_ab))
}

par <- c(0, 0)
tau <- 0
eta <- 1
y <- 0
n <- 73
w <- 0
gamma <- 1
mu_sigma <- sqrt(prm_prior$mu_sigma2)
eps <- eps_ab <- 0

args <- list(tau = tau, eta = eta, y = y, n = n, w = w, gamma = gamma, mu_sigma = mu_sigma, eps = eps, eps_ab = eps_ab)
options <- list(trace = 0, fnscale = -1, parscale = c(1, 1), ndeps = c(0.001, 0.001), maxit = 500, abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, gamma = 2)

res <- optim_rcpp(par, mu_delta_func, args, options, hessian = TRUE)
optim(par, mu_delta_func, gr = NULL, args, control = options, hessian = TRUE)

limits <- c(-15, 1, -5, 10)
my_plot(mu_delta_func, limits, args, type = "persp", npoints = 50, xlab = "mu", ylab = "delta", theta = 30, phi = 10)
my_plot(mu_delta_func, limits, args, type = "contour", npoints = 100, xlab = "mu", ylab = "delta", nlevels = 10)

mu_delta_func(par, args)

delta <- 0
mu <- seq(-15, 10, by = .1)
mu_val <- numeric(length(mu))
for (i in 1:length(mu)) {
  par <- c(mu[i], delta)
  mu_val[i] <- mu_delta_func(par, args)
}
plot(mu, mu_val, type = "l")

mu <- -2
delta <- seq(-15, 10, by = .1)
delta_val <- numeric(length(delta))
for (i in 1:length(delta)) {
  par <- c(mu, delta[i])
  delta_val[i] <- mu_delta_func(par, args)
}
plot(delta, delta_val, type = "l")
which.max(delta_val)

###

start <- res@control$burnin + 1
end <- res@control$burnin + res@control$nsim
ns <- res@dim$n_study
M <- res@dim$n_outcomes
nt <- res@dim$n_treatments
n <- res@dim$n_datapoints
nGamma <- dim(res@Gamma)[[1]]

mu <- res@mu
nm <- character(ns*M)
for (j in 1:M) {
  for (i in 1:ns) {
   nm[i + ns*(j - 1)] <- paste0("mu[", i, ",", j, "]")
  }
}
attr(mu, "dim") <- c(dim(mu)[1]*dim(mu)[2], dim(mu)[3])
mu <- t(mu)
dimnames(mu)[[2]] <- nm

delta <- res@delta
nm <- character(n*M)
for (j in 1:M) {
  for (i in 1:n) {
    nm[i + n*(j - 1)] <- paste0("delta[", i, ",", j, "]")
  }
}
attr(delta, "dim") <- c(dim(delta)[1]*dim(delta)[2], dim(delta)[3])
delta <- t(delta)
dimnames(delta)[[2]] <- nm

d <- res@d
nm <- character(nt*M)
for (j in 1:M) {
  for (i in 1:nt) {
    nm[i + nt*(j - 1)] <- paste0("d[", i, ",", j, "]")
  }
}
attr(d, "dim") <- c(dim(d)[1]*dim(d)[2], dim(d)[3])
d <- t(d)
dimnames(d)[[2]] <- nm

Sigma <- res@Sigma
nm <- character(M*(M + 1)/2)
count <- 1
for (j in 1:M) {
  for (i in j:M) {
    nm[count] <- paste0("Sigma_M[", i, ",", j, "]")
    count <- count + 1
  }
}
dimnames(Sigma)[[2]] <- nm

Gamma <- res@Gamma
nelem <- M*(M - 1)/2
nm <- character(nGamma*nelem)
count <- 1
for (h in 1:nGamma) {
  for (j in 1:(M - 1)) {
    for (i in (j + 1):M) {
      nm[count + nelem*(h - 1)] <- paste0("Gamma[", h, "][", i, ",", j, "]")
      count <- count + 1
    }
  }
  count <- 1
}
attr(Gamma, "dim") <- c(dim(Gamma)[1]*dim(Gamma)[2], dim(Gamma)[3])
Gamma <- t(Gamma)
dimnames(Gamma)[[2]] <- nm

x <- res@x
nm <- character(n*M)
for (j in 1:M) {
  for (i in 1:n) {
    nm[i + n*(j - 1)] <- paste0("x[", i, ",", j, "]")
  }
}
attr(x, "dim") <- c(dim(x)[1]*dim(x)[2], dim(x)[3])
x <- t(x)
dimnames(x)[[2]] <- nm

tmp <- cbind(mu, delta, d, Sigma, Gamma)#, x)
# ord <- order(dimnames(tmp)[[2]])
# tmp <- tmp[, ord]
out <- mcmc(tmp, start = start, end = end, thin = thin)
out <- as.mcmc(out)
summary(out)

###

summary(res, digits_summary = 5, latent_var = FALSE, data = nc_data, nthin = 50)

library(mcmcplots)
library(coda)

res_mcmc <- as.mcmc(res, data = nc_data, nthin = 50, latent_var = TRUE)
out <- summary(res, digits_summary = 5, latent_var = FALSE, print = FALSE, data = nc_data, nthin = 50)
to_plot <- "d"
to_keep <- substr(rownames(out), 1, nchar(to_plot) + 1) == paste0(to_plot, "[")
tmp <- out[to_keep, ]
caterplot(res_mcmc, to_plot, reorder = TRUE, quantiles = list(outer = c(0.025, 0.975), inner = c(0.1, 0.9)))
caterpoints(tmp[sort.int(tmp[, 5], decreasing = TRUE, index.return = TRUE)$ix, 1], pch = "x", col = "darkorange")
title(paste0("parameters ", to_plot, " (data from Achana et al., 2014)"), cex.main = 1.25)

traplot(res_mcmc, parms = to_plot, mar = c(1.0, 1.5, .75, .15) + .1, col = NULL, lty = 1, plot.title = NULL, main = NULL, greek = FALSE, style = "plain", cex.main = .85, cex.axis = .7, tcl = -.3, xaxt = "n")

###

library(coda)

res_mcmc <- as.mcmc(res, data = nc_data, nthin = 50, latent_var = TRUE)
str(res_mcmc)
tmp <- summary(res_mcmc)
plot(res_mcmc, smooth = FALSE, ask = TRUE)
autocorr.plot(res_mcmc, ask = TRUE)
gelman.diag(res_mcmc, autoburnin = FALSE)
gelman.plot(res_mcmc, autoburnin = FALSE, ask = TRUE)
raftery.diag(res_mcmc)

###

res_mcmc <- add_corr(res_mcmc)
str(res_mcmc)
res_mcmc[3, ]

###

f <- function(x, a) (x - a)^2
optimize(f, c(0, 1), tol = 0.0001, a = 1/3)$minimum
f_2 <- function(x, varargs) {
  a <- varargs[[1]]
  -(x - a)^2
}
optimize_rcpp(f_2, 0, 1, tol = 0.0001, list(a = 1/3))

## See where the function is evaluated:
res_r <- optimize(function(x) x^2*(print(x) - 1), lower = 0, upper = 10)$minimum
res_rcpp <- optimize_rcpp(function(x, varargs = list()) -x^2*(print(x) - 1), 0, 10, .Machine$double.eps^0.25, list())
identical(res_r, res_rcpp)

## "wrong" solution with unlucky interval and piecewise constant f():
f  <- function(x) ifelse(x > -1, ifelse(x < 4, exp(-1/abs(x - 1)), 10), 10)
fp <- function(x) { print(x); f(x) }
fp_2 <- function(x, vargs = list()) { print(x); -f(x) }

plot(f, -2,5, ylim = 0:1, col = 2)
optimize(fp, lower = -4, upper = 20)$minimum   # doesn't see the minimum
optimize_rcpp(fp_2, -4, 20, .Machine$double.eps^0.25, list())   # doesn't see the minimum
optimize(fp, lower = -7, upper = 20)$minimum   # ok
optimize_rcpp(fp_2, -7, 20, .Machine$double.eps^0.25, list())   # ok

library(numDeriv)
f <- function(x, a) (x - a)^2
x <- optimize(f, c(0, 1), tol = 0.0001, a = 1/3)$minimum
hessian(f, x, a = 1/3)

f_2 <- function(x, varargs) {
  a <- varargs[[1]]
  -(x - a)^2
}
args <- list(a = 1/3)
options <- list(trace = 0, fnscale = -1, parscale = 1, ndeps = 0.001, maxit = 500, abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, gamma = 2)
xmax <- optimize_rcpp(f_2, 0, 1, tol = 0.0001, args)
optimhess_rcpp(xmax, f_2, args, options)[1, 1]

f  <- function(x) ifelse(x > -1, ifelse(x < 4, exp(-1/abs(x - 1)), 10), 10)
xmin <- optimize(f, lower = -7, upper = 20)$minimum
hessian(f, xmin)
fp_2 <- function(x, vargs = list()) -f(x)
xmax <- optimize_rcpp(fp_2, -7, 20, .Machine$double.eps^0.25, list())   # ok
optimhess_rcpp(xmax, fp_2, list(), options)[1, 1]
laplace_u_rcpp(fp_2, -7, 20, tol = 0.0001, list(), options)

###

mu_logpost_r <- function(mu, delta, y, n, w, gamma, mu_sigma, eps, eps_ab) {
  theta <- mu + delta
  p <- expit_rcpp(theta)

  a <- pbinom(y - 1, n, p, 1, 1)
  b <- pbinom(y, n, p, 1, 1)
  # if (a > (1 - eps_ab)) {
  #   a <- 1 - eps_ab
  # } else if (a < eps_ab) {
  #   a <- eps_ab
  # }
  # if (b > (1 - eps_ab)) {
  #   b <- 1 - eps_ab
  # } else if (b < eps_ab) {
  #   b <- eps_ab
  # }
  # eps_a <- a/2.0
  # eps_b <- (1.0 - b)/2.0
  # if (a == b) {
  #   if (a == 0) {
  #     b <- b + eps
  #   } else if (b == 1) {
  #     a <- a - eps
  #   } else {
  #     a <- a - min(eps/2.0, eps_a)
  #     b <- b + min(eps/2.0, eps_b)
  #   }
  # }
  phi_inv_a <- qnorm(a, 0.0, 1.0, 1, 1)
  phi_inv_b <- qnorm(b, 0.0, 1.0, 1, 1)
  tmp <- log(pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 1, 0) - pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 1, 0))
  if (!is.finite(tmp)) {
    tmp <- log(pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 0, 0) - pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 0, 0))
  }
  if (!is.finite(tmp)) {
    tmp_a = pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 1, 1)
    tmp_b = pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 1, 1)
    tmp = log(exp(tmp_b - tmp_a) - 1) + tmp_a
  }
  if (!is.finite(tmp)) {
    tmp_a = pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 0, 1)
    tmp_b = pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 0, 1)
    tmp = log(exp(tmp_a - tmp_b) - 1) + tmp_b
  }
  if (!is.finite(tmp)) {
    return(NA)
  }

  lpost <- tmp + dnorm(mu, 0, mu_sigma, 1)

  return(lpost)
}

delta_logpost_r <- function(delta, mu, tau, eta, y, n, w, gamma, eps, eps_ab) {
  theta <- mu + delta
  p <- expit_rcpp(theta)

  a <- pbinom(y - 1, n, p, 1, 1)
  b <- pbinom(y, n, p, 1, 1)
  # if (a > (1 - eps_ab)) {
  #   a <- 1 - eps_ab
  # } else if (a < eps_ab) {
  #   a <- eps_ab
  # }
  # if (b > (1 - eps_ab)) {
  #   b <- 1 - eps_ab
  # } else if (b < eps_ab) {
  #   b <- eps_ab
  # }
  # eps_a <- a/2.0
  # eps_b <- (1.0 - b)/2.0
  # if (a == b) {
  #   if (a == 0) {
  #     b <- b + eps
  #   } else if (b == 1) {
  #     a <- a - eps
  #   } else {
  #     a <- a - min(eps/2.0, eps_a)
  #     b <- b + min(eps/2.0, eps_b)
  #   }
  # }
  phi_inv_a <- qnorm(a, 0.0, 1.0, 1, 1)
  phi_inv_b <- qnorm(b, 0.0, 1.0, 1, 1)
  tmp <- log(pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 1, 0) - pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 1, 0))
  if (!is.finite(tmp)) {
    tmp <- log(pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 0, 0) - pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 0, 0))
  }
  if (!is.finite(tmp)) {
    tmp_a = pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 1, 1)
    tmp_b = pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 1, 1)
    tmp = log(exp(tmp_b - tmp_a) - 1) + tmp_a
  }
  if (!is.finite(tmp)) {
    tmp_a = pnorm((phi_inv_a - w)/gamma, 0.0, 1.0, 0, 1)
    tmp_b = pnorm((phi_inv_b - w)/gamma, 0.0, 1.0, 0, 1)
    tmp = log(exp(tmp_a - tmp_b) - 1) + tmp_b
  }
  if (!is.finite(tmp)) {
    return(NA)
  }

  lpost <- tmp + dnorm(delta, tau, eta, 1)

  return(lpost)
}

#---
mu <- 4.2322400926
delta <- -3.2951930489
tau <- -3.4524806642
eta <- 0.5908990878
y <- 160
n <- 163
w <- 6.5050731659
gamma <- 0.3718791088
eps <- eps_ab <- 1e-12
delta_logpost(delta, mu, tau, eta, y, n, w, gamma, eps, eps_ab)
delta_logpost_r(delta, mu, tau, eta, y, n, w, gamma, eps, eps_ab)

prm <- list(mu = mu, tau = tau, eta = eta, y = y, n = n, w = w, gamma = gamma, eps = eps, eps_ab = eps_ab)
xmin <- -1e+01
xmax <- 1e+01
delta_grid <- seq(xmin, xmax, by = .01)
delta_lp <- numeric(length(delta_grid))
for (i in 1:length(delta_grid)) {
  delta_lp[i] <- delta_logpost_func(delta_grid[i], prm)
}
plot(delta_grid, delta_lp, type = "l")
options <- list(trace = 0, fnscale = -1, parscale = 1, ndeps = 0.001, maxit = 500, abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, gamma = 2)
tol <- .Machine$double.eps^0.25
x_sol <- optimize_rcpp(delta_logpost_func, xmin, xmax, tol, prm)
abline(v = x_sol, lty = 2)
laplace_u_rcpp(delta_logpost_func, xmin, xmax, tol, prm, options)

-1/optimhess_rcpp(x_sol, delta_logpost_func, prm, options)

#---
mu = 35.1854870702
delta = 0.3908172659
tau = -0.1690032097
eta = 0.5937496633
y =  139
n =  139
w = 25.7475589281
gamma = 0.4806157901
eps <- eps_ab <- 1e-3

delta_logpost_r(delta, mu, tau, eta, y, n, w, gamma, eps, eps_ab)

prm <- list(mu = mu, tau = tau, eta = eta, y = y, n = n, w = w, gamma = gamma, eps = eps, eps_ab = eps_ab)
xmin <- -1e+01
xmax <- 1e+01
delta_grid <- seq(xmin, xmax, by = .01)
delta_lp <- numeric(length(delta_grid))
for (i in 1:length(delta_grid)) {
  delta_lp[i] <- delta_logpost_func(delta_grid[i], prm)
}
plot(delta_grid, delta_lp, type = "l")
options <- list(trace = 0, fnscale = -1, parscale = 1, ndeps = 0.001, maxit = 500, abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, gamma = 2)
tol <- .Machine$double.eps^0.25
x_sol <- optimize_rcpp(delta_logpost_func, xmin, xmax, tol, prm)
abline(v = x_sol, lty = 2)
laplace_u_rcpp(delta_logpost_func, xmin, xmax, tol, prm, options)

-1/optimhess_rcpp(x_sol, delta_logpost_func, prm, options)

#---
mu <- -21.2769010855
delta <- 0.0000000000
tau <- 0.6098463488
eta <- 0.7939631189
y <- 0
n <- 57
w <- 18.2797027337
gamma <- 0.2791747819
eps <- eps_ab <- 1e-12

delta_logpost_r(delta, mu, tau, eta, y, n, w, gamma, eps, eps_ab)
delta_logpost(delta, mu, tau, eta, y, n, w, gamma, eps, eps_ab)

prm_delta <- list(mu = mu, tau = tau, eta = eta, y = y, n = n, w = w, gamma = gamma, eps = eps, eps_ab = eps_ab)
xmin <- -1e+02
xmax <- 1e+01
delta_grid <- seq(xmin, xmax, by = .2)
delta_lp <- numeric(length(delta_grid))
for (i in 1:length(delta_grid)) {
  delta_lp[i] <- delta_logpost(delta_grid[i], mu, tau, eta, y, n, w, gamma, eps, eps_ab)
}
plot(delta_grid, delta_lp, type = "l")
options <- list(trace = 0, fnscale = -1, parscale = 1, ndeps = 0.001, maxit = 500, abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, gamma = 2)
tol <- .Machine$double.eps^0.25
x_sol <- optimize_rcpp(delta_logpost_func, xmin, xmax, tol, prm_delta)
abline(v = x_sol, lty = 2)
# laplace_u_rcpp(delta_logpost_func, xmin, xmax, tol, prm_delta, options)

y_val <- delta_lp
y_na <- is.na(y_val)
y_val <- y_val[!y_na]
x <- delta_grid
x <- x[!y_na]
y_lm <- lm(y_val ~ x + I(x^2), na.action = "na.omit")
y_hat <- predict(y_lm, newdata = data.frame(x = delta_grid))
lines(delta_grid, y_hat, col = "green", lty = 2)

coef_r <- as.numeric(coef(y_lm))
coef_rcpp <- as.numeric(ols_coef(xmin, xmax, prm_delta, TRUE))
all.equal(coef_r, coef_rcpp)

opt <- laplace_u_rcpp(delta_logpost_quad, xmin, xmax, tol, list(coef = coef_rcpp), options)
abline(v = opt$mode, lty = 2, col = "green")

#---
mu <- -36.2203837777
delta <- 0.0000000000
tau <- 0.8769834137
eta <- 0.7939912184
y <- 0
n <- 669
w <- 6.2582675848
gamma <- 0.4354416419

# mu <- -14.8403996815
# delta <- 0.0000000000
# tau <- 1.0721924294
# eta <- 0.7942526792
# y <- 4
# n <- 57
# w <- 19.4328250016
# gamma <- 0.2924673485

# mu <- -21.2769010855
# delta <- 0.0000000000
# tau <- 0.6098463488
# eta <- 0.7939631189
# y <- 0
# n <- 57
# w <- 18.2797027337
# gamma <- 0.2791747819
mu_sigma <- sqrt(prm_prior$mu_sigma2)
eps <- eps_ab <- 1e-3

mu_logpost_r(mu, delta, y, n, w, gamma, mu_sigma, eps, eps_ab)
mu_logpost(mu, delta, y, n, w, gamma, mu_sigma, eps, eps_ab)

prm_mu <- list(delta = delta, y = y, n = n, w = w, gamma = gamma, mu_sigma = mu_sigma, eps = eps, eps_ab = eps_ab)
xmin <- -1e+02
xmax <- -1e+01
mu_grid <- seq(xmin, xmax, by = .2)
mu_lp <- numeric(length(mu_grid))
for (i in 1:length(mu_grid)) {
  mu_lp[i] <- mu_logpost(mu_grid[i], delta, y, n, w, gamma, mu_sigma, eps, eps_ab)
}
plot(mu_grid, mu_lp, type = "l")
options <- list(trace = 0, fnscale = -1, parscale = 1, ndeps = 0.001, maxit = 500, abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, gamma = 2)
tol <- .Machine$double.eps^0.25
x_sol <- optimize_rcpp(mu_logpost_func, xmin, xmax, tol, prm_mu)
abline(v = x_sol, lty = 2)
# laplace_u_rcpp(mu_logpost_func, xmin, xmax, tol, prm_mu, options)

y_val <- mu_lp
y_na <- is.na(y_val)
y_val <- y_val[!y_na]
x <- mu_grid
x <- x[!y_na]
y_lm <- lm(y_val ~ x + I(x^2), na.action = "na.omit")
y_hat <- predict(y_lm, newdata = data.frame(x = mu_grid))
lines(mu_grid, y_hat, col = "green", lty = 2)

coef_r <- as.numeric(coef(y_lm))
coef_rcpp <- as.numeric(ols_coef(xmin, xmax, prm_mu, FALSE))
all.equal(coef_r, coef_rcpp)

opt <- laplace_u_rcpp(mu_logpost_quad, xmin, xmax, tol, list(coef = coef_rcpp), options)
abline(v = opt$mode, lty = 2, col = "green")

tmp <- ols_pred(coef_rcpp, mu)
abline(h = tmp, lty = 2, col = "green")
abline(v = mu, lty = 2, col = "green")

###

library(ggmcmc)

res_mcmc <- as.mcmc(res)
res_ggmcmc <- ggs(res_mcmc, family = "mu")
# str(res_ggmcmc)
ggmcmc(res_ggmcmc)

###

set.seed(1406)

eta <- 4
R <- rlkj_arma(K = 2, eta = eta)
dlkj_arma(R, eta, FALSE)
rethinking::dlkjcorr(R, eta, TRUE)

# plot density of correlation
n <- 1e+4
R <- array(NA, dim = c(n, 2, 2))
for (i in 1:n) {
  R[i, , ] <- rlkj_arma(K = 2, eta = 4)
}
hist(R[, 1, 2])

# visualize 3x3 matrix
set.seed(123)
R <- array(NA, dim = c(n, 3, 3))
for (i in 1:n) {
  R[i, , ] <- rlkj_arma(K = 3, eta = 20)
}
plot(R[, 1, 2], R[, 1, 3], col = col.alpha("black", 0.2), pch = 16)

set.seed(123)
R_2 <- rethinking::rlkjcorr(n, K = 3, eta = 20)
plot(R_2[, 1, 2], R_2[, 1, 3], col = col.alpha("black", 0.2), pch = 16)

###

r12 <- seq(-1, 1, .001)
d_r12 <- d_r12_2 <- numeric(length(r12))
eta <- 1
for (i in 1:length(r12)) {
  R <- matrix(c(1, rep(r12[i], 2), 1), nrow = 2, ncol = 2)
  d_r12[i] <- dlkj_arma(R, eta, FALSE)
  d_r12_2[i] <- rethinking::dlkjcorr(R, eta, FALSE)
}
plot(r12, d_r12, type = "l")
plot(r12, d_r12_2, type = "l")

###

load("~/dev/netcopula/demo/tmp/examples_10.RData")
res3 <- res
load("~/dev/netcopula/demo/tmp/examples_11.RData")
res3i <- res
load("~/dev/netcopula/demo/tmp/examples2_10.RData")
res2 <- res
load("~/dev/netcopula/demo/tmp/examples2_11.RData")
res2i <- res
rm(res)

res2i_mcmc <- as.mcmc(res2i)
d2i <- summary(exp(res2i_mcmc[, 1:8]))
res3i_mcmc <- as.mcmc(res3i)
d3i <- summary(exp(res3i_mcmc[, 1:24]))
nm <- match(dimnames(d2i$statistics)[[1]], dimnames(d3i$statistics)[[1]])
d3i_sub <- summary(exp(res3i_mcmc[, nm]))
plot(d2i$statistics[, "Mean"], d3i_sub$statistics[, "Mean"], xlim = c(0, 5), ylim = c(0, 5))
abline(0, 1, lty = 2, col = "gray")

res2_mcmc <- as.mcmc(res2)
d2 <- summary(exp(res2_mcmc[, 1:8]))
res3_mcmc <- as.mcmc(res3)
d3 <- summary(exp(res3_mcmc[, 1:24]))
nm <- match(dimnames(d2$statistics)[[1]], dimnames(d3$statistics)[[1]])
d3_sub <- summary(exp(res3_mcmc[, nm]))
plot(d2$statistics[, "Mean"], d3_sub$statistics[, "Mean"], xlim = c(0, 70), ylim = c(0, 70))
abline(0, 1, lty = 2, col = "gray")

plot(d2i$statistics[, "Mean"], d2$statistics[, "Mean"], xlim = c(0, 70), ylim = c(0, 70))
abline(0, 1, lty = 2, col = "gray")

plot(d3i$statistics[, "Mean"], d3$statistics[, "Mean"], xlim = c(0, 16), ylim = c(0, 16))
abline(0, 1, lty = 2, col = "gray")

###

n_outcomes <- n_trt <- 3
n_studies <- 30
n_bin <- rep(100, n_outcomes) # binomial distributions sizes
d <- array(data = c(0, -2, -1, 0, 3, 2), dim = c(2, 3))
sigma_vec <- c(1, 1.6, 1.8) # variances in Sigma_M
rho_vec <- c(0.5, 0.5, 0.5)
gamma_vec <- rep(0.9, n_outcomes) # components of copula association matrix
ref_trt <- 1
mu_sigma2 <- 10^-1

seed <- 2301
set.seed(seed)
data_all <- nc_data_simulate(d, sigma_vec, rho_vec, gamma_vec, n_studies, n_trt, n_bin, n_outcomes, mu_sigma2, ref_trt)
nc_data <- data_all$nc_data

###

library(netcopula)

a <- 2.5
b <- 4
mu <- 2
sigma <- 1
n <- 1e+6

set.seed(123)
rn1 <- rtruncnorm_rcpp(n, a, b, mu, sigma)

set.seed(123)
U <- runif(n)
Phi_a <- pnorm((a - mu)/sigma)
Phi_b <- pnorm((b - mu)/sigma)
rn2 <- qnorm(Phi_a + U*(Phi_b - Phi_a))*sigma + mu

summary(rn1)
summary(rn2)

hist(rn1, freq = FALSE, breaks = 100, border = rgb(0, 0, 0, alpha = .5))
hist(rn2, freq = FALSE, breaks = 100, add = TRUE, border = rgb(1, 0, 0, alpha = .4))

plot(sort(rn1), sort(rn2), pch = 20, cex = .7, col = rgb(0, 0, 0, alpha = .02))

###

set.seed(1406)
n_datapoints <- 2*n_studies
M <- n_outcomes
narms <- nc_data@study_id$narms
y_imp <- nc_data@study_data[, paste0("y", 1:M)]
n_imp <- nc_data@study_data[, paste0("n", 1:M)]
x_new <- matrix(NA, n_datapoints, M)
iter <- 2
max_iter <- 150
plot_it <- TRUE

while (iter <= max_iter) {
  x_iter <- res@x[, , iter - 1]
  D_iter <- diag(res@D[, , iter - 1])
  Gamma_iter <- matrix(1, M, M)
  Gamma_iter[lower.tri(Gamma_iter)] <- Gamma_iter[upper.tri(Gamma_iter)] <- res@Gamma[, , iter - 1]
  S_q_iter <- S_q_curr <- Sigma_q_prop_iter <- matrix(0, M, M)
  S_q_iter[lower.tri(S_q_iter, diag = TRUE)] <- res@S_q[, , iter - 1]
  S_q_iter[upper.tri(S_q_iter)] <- S_q_iter[lower.tri(S_q_iter)]
  S_q_curr[lower.tri(S_q_curr, diag = TRUE)] <- res@S_q[, , iter]
  S_q_curr[upper.tri(S_q_curr)] <- S_q_curr[lower.tri(S_q_curr)]
  Sigma_q_prop_iter[lower.tri(Sigma_q_prop_iter, diag = TRUE)] <- res@Sigma_q_prop[, , iter]
  Sigma_q_prop_iter[upper.tri(Sigma_q_prop_iter)] <- Sigma_q_prop_iter[lower.tri(Sigma_q_prop_iter)]
  S_q <- matrix(0, n_outcomes, n_outcomes)
  for (i in 1:nrow(x_iter)) {
    S_q <- S_q + D_iter %*% tcrossprod(x_iter[i, ]) %*% D_iter
  }
  if (!is_positive_definite(S_q, 1)) {
    S_q <- make_positive_definite(S_q)
  }
  par(mfrow = c(1, 3))
  if (plot_it) plot(S_q, S_q_curr, main = paste0("S_q - ", all.equal(S_q, S_q_curr))); abline(0, 1, lty = 2)
  Sigma_q_prop <- rinvwish_arma(n_datapoints, S_q)
  # Sigma_q_prop <- solve(rWishart(1, n_datapoints, solve(S_q))[, , 1])
  # Sigma_q_prop <- bayesm::rwishart(n_datapoints, solve(S_q))$IW
  if (plot_it) plot(Sigma_q_prop, Sigma_q_prop_iter, main = "Sigma_q_prop"); abline(0, 1, lty = 2)
  # D_prop <- diag(sqrt(diag(Sigma_q_prop)))
  # Gamma_q_prop <- cov2cor(Sigma_q_prop)
  # Gamma_rate <- 0.5*(M + 1)*(log(det(Gamma_q_prop)) - log(det(Gamma_iter)))
  # ran_unif <- runif(1, 0.0, 1.0)
  # accept_Gamma_q <- 0
  # if (ran_unif < exp(Gamma_rate)) {
  #   Gamma_iter <- Gamma_q_prop
  #   D_iter <- D_prop
  #   accept_Gamma_q <- accept_Gamma_q + 1
  # }

  mu_long_iter <- param_long(res@mu[, , iter], narms, FALSE)
  delta_iter <- res@delta[, , iter]
  Gamma_k <- Gamma_iter
  for (ik in 1:n_datapoints) {
    for (m in 1:M) {
      theta_ikm <- mu_long_iter[ik, m] + delta_iter[ik, m]
      p_ikm <- expit_rcpp(theta_ikm)
      a_ikm <- pbinom(y_imp[ik, m] - 1, n_imp[ik, m], p_ikm, 1, 1)
      b_ikm <- pbinom(y_imp[ik, m], n_imp[ik, m], p_ikm, 1, 1)
      Gamma_k_m_m <- Gamma_k[m, ][-m]
      Gamma_k_m <- Gamma_k[-m, -m]
      x_ik_m <- x_iter[ik, -m]
      w_ikm <- as.numeric(Gamma_k_m_m %*% solve(Gamma_k_m) %*% x_ik_m)
      gamma_ikm <- as.numeric(sqrt(Gamma_k[m, m] - Gamma_k_m_m %*% solve(Gamma_k_m) %*% Gamma_k_m_m))
      x_new[ik, m] <- rtruncnorm_rcpp(1, qnorm(a_ikm, 0.0, 1.0, 1, 1), qnorm(b_ikm, 0.0, 1.0, 1, 1), w_ikm, gamma_ikm)
    }
  }
  if (plot_it) plot(x_new, res@x[, , iter], main = "x"); abline(0, 1, lty = 2)

  iter <- iter + 1
}

# cov_x <- apply(res@x, 3, cov)

# par(mar = c(4, 4, 1.5, .15) + .1, mfrow = c(2, 3))
# for (ri in 1:M) {
#   for (ci in ri:M) {
#     cp <- apply(res@x, c(1, 3), function(x) tcrossprod(x)[ri, ci])
#     cp <- colSums(cp)
#     plot(cp, type = "l", ylab = paste0("crossprod[", ri, ", ", ci, "]"), xlab = "iteration")
#   }
# }

###

max_iter <- 10000
Gamma_prop <- array(NA, dim = c(n_outcomes, n_outcomes, max_iter))
iscor <- logical(max_iter)
for (iter in 1:max_iter) {
  x <- rmvn_arma(2*n_studies, rep(0, n_outcomes), Gamma_tmp)

  D <- 1/colSums(x^2)
  # D*colSums(x^2)
  x_star <- x %*% sqrt(diag(D))
  # colSums(x_star^2)

  # t(x_star) %*% x_star
  # crossprod(x_star)
  S_q <- matrix(0, n_outcomes, n_outcomes)
  for (i in 1:nrow(x)) {
    S_q <- S_q + tcrossprod(x_star[i, ])
  }
  Sigma_prop <- rinvwish_arma(n_datapoints, S_q)
  D_prop <- diag(sqrt(diag(Sigma_prop)))
  # Gamma_prop[, , iter] <- solve(D_prop) %*% Sigma_prop %*% solve(D_prop)
  Gamma_prop[, , iter] <- cov2cor(Sigma_prop)
  iscor[iter] <- is_correlation(Gamma_prop[, , iter])

  x_new <- x_star %*% solve(D_prop)
  # plot(x, x_new)
}

all(iscor)
ri <- 2
cj <- 1
Gamma_i_j <- Gamma_prop[ri, cj, ]
plot(Gamma_i_j, type = "l")
c(mean(Gamma_i_j), quantile(Gamma_i_j, probs = c(.025, .975)))

###

iter <- 1
x_unadj <- res@x_unadj[, , iter]
x_adj <- res@x[, , iter]
D <- diag(1/sqrt(colSums(x_unadj^2)))
x_star <- x_unadj %*% D
colSums(x_star^2)

S_q <- matrix(NA, n_outcomes, n_outcomes)
S_q[lower.tri(S_q, diag = TRUE)] <- res@S_q[, , iter]
S_q[upper.tri(S_q)] <- S_q[lower.tri(S_q)]
all.equal(S_q, crossprod(x_star))

Gamma <- matrix(1, n_outcomes, n_outcomes)
Gamma[lower.tri(Gamma)] <- Gamma[upper.tri(Gamma)] <- res@Gamma[, , iter]
Sigma_q_prop <- matrix(NA, n_outcomes, n_outcomes)
Sigma_q_prop[lower.tri(Sigma_q_prop, diag = TRUE)] <- res@Sigma_q_prop[, , iter]
Sigma_q_prop[upper.tri(Sigma_q_prop)] <- Sigma_q_prop[lower.tri(Sigma_q_prop)]
D_prop <- diag(sqrt(diag(Sigma_q_prop)))
Gamma_prop <- cov2cor(Sigma_q_prop)
all.equal(Gamma_prop, Gamma)

plot(as.numeric(x_unadj - x_adj))
plot(x_unadj, x_adj)
all.equal(x_adj, x_star %*% solve(D_prop))
cor(x_unadj)
cor(x_adj)
cor(x_star)

# plot(res@D[1, 3, ], type = "l")
# apply(res@D[1, , ], 1, summary)
# plot(res@D[1, 2, ], res@D[1, 3, ], pch = ".")

###

iter <- 1
n_datapoints <- 2*n_studies
M <- n_outcomes
Gamma_curr <- matrix(1, M, M)
Gamma_curr[lower.tri(Gamma_curr)] <- Gamma_curr[upper.tri(Gamma_curr)] <- res@Gamma[1, , iter]
x_unadj <- x
mu_long <- param_long(res@mu[, , iter], nc_data@study_id$narms, FALSE)
delta <- res@delta[, , iter]
y_imp <- nc_data@study_data[, paste0("y", 1:M)]
n_imp <- nc_data@study_data[, paste0("n", 1:M)]
max_iter <- 1000
x_all <- array(NA, dim = c(n_datapoints, M, max_iter))
for (i_iter in 1:max_iter) {
  print(paste0("iter ", i_iter, "/", max_iter))
  for (ik in 1:n_datapoints) {
    Gamma_k <- Gamma_curr
    for (m in 1:M) {
      theta_ikm <- mu_long[ik, m] + delta[ik, m]
      p_ikm <- expit_rcpp(theta_ikm)
      a_ikm <- pbinom(y_imp[ik, m] - 1, n_imp[ik, m], p_ikm, TRUE, TRUE)
      b_ikm <- pbinom(y_imp[ik, m], n_imp[ik, m], p_ikm, TRUE, TRUE)
      Gamma_k_m_m <- Gamma_k[m, -m]
      Gamma_k_m <- Gamma_k[-m, -m]
      x_ik_m <- x_unadj[ik, -m]
      w_ikm <- as.numeric(t(Gamma_k_m_m) %*% solve(Gamma_k_m) %*% x_ik_m)
      gamma_ikm <- sqrt(as.numeric(Gamma_k[m, m] - t(Gamma_k_m_m) %*% solve(Gamma_k_m) %*% Gamma_k_m_m))
      x_unadj[ik, m] <- rtruncnorm_rcpp(1, qnorm(a_ikm, 0.0, 1.0, TRUE, TRUE), qnorm(b_ikm, 0.0, 1.0, TRUE, TRUE), w_ikm, gamma_ikm)
    }
  }
  x_all[, , i_iter] <- x_unadj
}

corr_x <- apply(x_all, 3, cor)
matrix(rowMeans(corr_x), nrow = 3, ncol = 3)

###

set.seed(101)

M <- 3
eta <- 1
C <- rlkj_arma(M, eta)

x1 <- rmvn_arma(n, rep(0, 3), C)
x2 <- rmvn_arma(n, rep(0, 3), C)
gausscopdens_logpost(x1, C) - gausscopdens_logpost(x2, C)
r_logpost(C, x1) - r_logpost(C, x2)
