## ----global_options, include=FALSE---------------------------------------
# knitr::opts_chunk$set(warning=FALSE, message=FALSE, comment='')
path <- "~/dev/netcopula/demo/tmp"
library(netcopula)

## ----loading_1, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
load(file = file.path(path, "examples_11.RData"))

## ----Gamma_1, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE, results="asis"----
plot(x = res, to_plot = "Gamma", type = "trace")
cat("\\clearpage")
plot(x = res, to_plot = "Gamma", type = "caterpillar")

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

## ----loading_2, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
load(file = file.path(path, "examples2_11.RData"))

## ----Gamma_2, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE, results="asis"----
plot(x = res, to_plot = "Gamma", type = "trace")
cat("\\clearpage")
plot(x = res, to_plot = "Gamma", type = "caterpillar")

## ----Sigma_M_2, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE, results="asis"----
plot(x = res, to_plot = "Sigma_M", type = "trace")
cat("\\clearpage")
plot(x = res, to_plot = "Sigma_M", type = "caterpillar")

## ----Corr_M_2, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE, results="asis"----
plot(x = res, to_plot = "Corr_M", type = "trace")
cat("\\clearpage")
plot(x = res, to_plot = "Corr_M", type = "caterpillar")

## ----mu_2, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE, results="asis"----
plot(x = res, to_plot = "mu", type = "trace")
cat("\\clearpage")
plot(x = res, to_plot = "mu", type = "caterpillar")

## ----delta_2, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE, results="asis"----
plot(x = res, to_plot = "delta", type = "trace")
cat("\\clearpage")
plot(x = res, to_plot = "delta", type = "caterpillar")

## ----d_2, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE, results="asis"----
plot(x = res, to_plot = "d", type = "trace")
cat("\\clearpage")
plot(x = res, to_plot = "d", type = "caterpillar")

## ----x_2, echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE, results="asis"----
plot(x = res, to_plot = "x", latent_var = TRUE, type = "trace")
cat("\\clearpage")
plot(x = res, to_plot = "x", latent_var = TRUE, type = "caterpillar")

