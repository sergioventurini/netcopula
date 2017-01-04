## ----global_setting, echo=TRUE, results='hide', message=FALSE------------
# input parameters
r_sim <- function(n, M, sigma2_r) {
  n_beta <- M*(M + 1)/2
  n_r <- M*(M - 1)/2
  Rho_arr <- array(NA, dim = c(n, n_r))
  for (i in 1:n) {
    beta <- rnorm(n_beta, mean = 0, sd = sqrt(sigma2_r))
    R <- matrix(0, nrow = M, ncol = M)
    R_diag <- upper.tri(R)
    diag(R_diag) <- TRUE
    R[R_diag] <- beta
    diag(R) <- exp(diag(R))
    R_inv <- solve(R)
    Sigma_M <- R_inv%*%t(R_inv)
    sd_M <- solve(diag(sqrt(diag(Sigma_M))))
    Rho_M <- sd_M%*%Sigma_M%*%sd_M
    Rho_arr[i, ] <- Rho_M[upper.tri(Rho_M)]
	}
	return(as.numeric(Rho_arr))
}

set.seed(1406)
n <- 1000
M <- 5
n_r <- M*(M - 1)/2
sigma2_r_exp <- sort(c(-3:1, -1/2:9))
sigma2_r <- 10^sigma2_r_exp

## ----prior, echo=FALSE, results="asis"-----------------------------------
res <- array(NA, dim = c(length(sigma2_r), n*n_r))
for (k in seq_along(sigma2_r)) {
  res[k, ] <- r_sim(n = n, M = M, sigma2_r = sigma2_r[k])
  hist(res[k, ], breaks = 30, xlab = "Correlations induced by Sigma_M", main = paste0("sigma2_r = 10^", sigma2_r_exp[k]), xlim = c(-1, 1))
  rug(res[k, ])
  cat("\\clearpage \n\n")
}

