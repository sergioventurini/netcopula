## ----global_options, include=FALSE---------------------------------------
# knitr::opts_chunk$set(warning=FALSE, message=FALSE, comment='')
path <- "~/dev/netcopula/demo/tmp"
gen_res <- function(to_plot) {
  to_keep <- substr(rownames(out), 1, nchar(to_plot) + 1) == paste0(to_plot, "[")
  tmp <- out[to_keep, ]
  caterplot(res_mcmc, to_plot, reorder = TRUE, quantiles = list(outer = c(0.025, 0.975), inner = c(0.1, 0.9)))
  caterpoints(tmp[sort.int(tmp[, 5], decreasing = TRUE, index.return = TRUE)$ix, 1], pch = "x", col = "darkorange")
  title(paste0("parameters ", to_plot, " (data from Achana et al., 2014)"), cex.main = 1.25)
  
  traplot(res_mcmc, parms = to_plot, mar = c(1.0, 1.5, .75, .15) + .1, col = NULL, lty = 1, plot.title = NULL, main = NULL, greek = FALSE, style = "plain", cex.main = .85, cex.axis = .7, tcl = -.3, xaxt = "n")
}

## ----loading, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
library(netcopula)
library(mcmcplots)
library(coda)

load(file = file.path(path, "examples_5.RData"))
if (!.hasSlot(res, "S_q")) {
  res@S_q <- array()
  res@Sigma_q_prop <- array()
  res@rho_mu <- array()
  res@cov_mu <- array()
  res@ar_mu_vec <- array()
  res@mu_rate <- array()
  res@mu_prop <- array()
  res@tp_mu <- array()
  res@delta_prop <- array()
  res@data <- nc_data
}

## ----Gamma_fix, echo=FALSE-----------------------------------------------
cor(na.omit(nc_data@study_data[, colnames(nc_data@study_data)[regexpr("y", colnames(nc_data@study_data)) == 1]]))

## ----loading_1, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
res_mcmc <- as.mcmc(res, data = nc_data, nthin = 50, latent_var = TRUE)
res_mcmc <- add_corr(res_mcmc)
out <- summary(res, digits_summary = 5, latent_var = TRUE, print = FALSE, data = nc_data, nthin = 50)

## ----Gamma_1, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("Gamma")

## ----Sigma_M_1, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("Sigma_M")

## ----Corr_M_1, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("Corr_M")

## ----mu_1, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("mu")

## ----delta_1, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("delta")

## ----d_1, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("d")

## ----x_1, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("x")

## ----loading_2, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
load(file = file.path(path, "examples_6.RData"))
if (!.hasSlot(res, "S_q")) {
  res@S_q <- array()
  res@Sigma_q_prop <- array()
  res@rho_mu <- array()
  res@cov_mu <- array()
  res@ar_mu_vec <- array()
  res@mu_rate <- array()
  res@mu_prop <- array()
  res@tp_mu <- array()
  res@delta_prop <- array()
  res@data <- nc_data
}
res_mcmc <- as.mcmc(res, data = nc_data, nthin = 50, latent_var = TRUE)
res_mcmc <- add_corr(res_mcmc)
out <- summary(res, digits_summary = 5, latent_var = TRUE, print = FALSE, data = nc_data, nthin = 50)

## ----Gamma_2, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("Gamma")

## ----Sigma_M_2, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("Sigma_M")

## ----Corr_M_2, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("Corr_M")

## ----mu_2, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("mu")

## ----delta_2, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("delta")

## ----d_2, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("d")

## ----x_2, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("x")

## ----loading_3, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
load(file = file.path(path, "examples_7.RData"))
if (!.hasSlot(res, "S_q")) {
  res@S_q <- array()
  res@Sigma_q_prop <- array()
  res@rho_mu <- array()
  res@cov_mu <- array()
  res@ar_mu_vec <- array()
  res@mu_rate <- array()
  res@mu_prop <- array()
  res@tp_mu <- array()
  res@delta_prop <- array()
  res@data <- nc_data
}
res_mcmc <- as.mcmc(res, data = nc_data, nthin = 50, latent_var = TRUE)
res_mcmc <- add_corr(res_mcmc)
out <- summary(res, digits_summary = 5, latent_var = TRUE, print = FALSE, data = nc_data, nthin = 50)

## ----Gamma_3, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("Gamma")

## ----Sigma_M_3, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("Sigma_M")

## ----Corr_M_3, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("Corr_M")

## ----mu_3, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("mu")

## ----delta_3, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("delta")

## ----d_3, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("d")

## ----x_3, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
gen_res("x")

