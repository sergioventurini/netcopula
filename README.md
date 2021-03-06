# netcopula

<!-- badges: start -->

[![](http://cranlogs.r-pkg.org/badges/grand-total/netcopula?color=blue)](https://cran.r-project.org/package=netcopula)
[![Travis build status](https://travis-ci.org/sergioventurini/netcopula.svg?branch=master)](https://travis-ci.org/sergioventurini/netcopula)
[![CRAN
status](https://www.r-pkg.org/badges/version/netcopula)](https://cran.r-project.org/package=netcopula)
<!-- badges: end -->

## Overview

###### Current release: 0.1.3
###### R version required: at least 3.6.0
`R` package implementing a Bayesian copula-based model for multivariate network meta-analysis (NMA).

## Installation

Since the package requires some code to be compiled, you need a working C++
compiler. To get it:

- On Windows, install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
- On Mac, install Xcode from the app store.
- On Linux, `sudo apt-get install r-base-dev` or similar.

Then, the easiest way to get the package is to install it from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("sergioventurini/netcopula")
```

See the help pages of the `netcopula()` function for some examples
or have a look at the demos in the package.

