## ----global_options, include=FALSE---------------------------------------
# knitr::opts_chunk$set(warning=FALSE, message=FALSE, comment='')
path <- "~/dev/netcopula/demo/tmp"

## ----loading, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
library(netcopula)

load(file = file.path(path, "examples_9.RData"))
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

## ----Sigma_M_1, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE, results="asis"----
plot(x = res, to_plot = "Sigma_M", type = "trace")
cat("\\clearpage")
plot(x = res, to_plot = "Sigma_M", type = "caterpillar")

## ----Corr_M_1, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE, results="asis"----
plot(x = res, to_plot = "Corr_M", type = "trace")
cat("\\clearpage")
plot(x = res, to_plot = "Corr_M", type = "caterpillar")

## ----mu_1, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE, results="asis"----
plot(x = res, to_plot = "mu", type = "trace")
cat("\\clearpage")
plot(x = res, to_plot = "mu", type = "caterpillar")

## ----delta_1, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE, results="asis"----
plot(x = res, to_plot = "delta", type = "trace")
cat("\\clearpage")
plot(x = res, to_plot = "delta", type = "caterpillar")

## ----d_1, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE, results="asis"----
plot(x = res, to_plot = "d", type = "trace")
cat("\\clearpage")
plot(x = res, to_plot = "d", type = "caterpillar")

## ----x_1, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE, results="asis"----
plot(x = res, to_plot = "x", latent_var = TRUE, type = "trace")
cat("\\clearpage")
plot(x = res, to_plot = "x", latent_var = TRUE, type = "caterpillar")

