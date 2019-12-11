#' Copula-Based Model for Multivariate Network Meta-Analysis.
#'
#' @docType package
#'
#' @name netcopula-package
#' @aliases netcopula-pkg
#' @aliases netcopula-p
#'
#' @description
#' The \pkg{netcopula} package implements a Bayesian algorithm for fitting
#' a copula-based model for multivariate meta-analysis.
#'
#' For efficiency reasons, the core computations in the package are implemented
#' using the \code{C} programming language and the \pkg{RcppArmadillo} package.
#'
#' Plotting functionalities are imported from the nice \pkg{bayesplot} package.
#' Currently, the package includes methods for binary data only. In future
#' releases routines will be added specifically for continuous (i.e. normal),
#' and count data.
#'
#' @section \pkg{netcopula} classes:
#' The \pkg{netcopula} package defines the following new classes:
#' \itemize{
#'   \item{\code{\link{nc_data}}: }{defines the data to use in a netcopula model.}
#'   \item{\code{\link{nc_mcmc}}: }{defines the results of a netcopula analysis
#'     for a single MCMC chain.}
#' }
#' The package includes \code{print}, \code{summary} and \code{plot} methods
#'   for each one of these classes.
#'
#' @section Resources:
#' \itemize{
#'  \item{\strong{Bug reports}:}{
#'  If you have noticed a bug that needs to be fixed, please let us know at the
#'  \pkg{netcopula} issue tracker on GitHub:
#'
#'  \url{https://github.com/sergioventurini/netcopula/issues/}.
#'  }
#'  \item{\strong{General questions and help}:}{
#'  To ask a question about \pkg{netcopula} send and email to:
#'
#'  \email{sergio.venturini@unito.it}.
#' }
#' }
#'
#' @seealso \code{\link[bayesplot]{theme_default}} for the default ggplot theme
#'  used by \pkg{bayesplot}.
#' @seealso \code{\link[bayesplot]{bayesplot-colors}} to set or view the color
#'  scheme used for plotting with \pkg{bayesplot}.
#' @seealso \code{\link[ggplot2]{ggsave}} in \pkg{ggplot2} for saving plots.
#'
#' @references
#'   Graziani, R., Venturini, S. (2019), "Network Meta-Analysis for Adverse Events:
#'   An Exact (Bayesian) Approach Using Gaussian Copulas", Technical report.
#'
#' @useDynLib netcopula
#' @import bayesplot
#' @import coda
#' @import ggplot2
#' @import graphics
#' @import methods
#' @import stats
#' @import utils
NULL

#' Home safety data
#'
#' Data from a Cochrane systematic review of home safety education and
#' provision of safety equipment for injury prevention in children.
#'
#' @format A data frame with 45 rows and 12 variables. The variables are as
#' follows:
#'
#' \itemize{
#'   \item studyid: study identification number
#'   \item study: study progressive number in the sample
#'   \item arm: single arm of a specific study
#'   \item trt: intervention in the specific study arm; 1 means the
#'   control/baseline treatment
#'   \item baseline: baseline treatment for a specific study
#'   \item narms: total number of arms for a specific study
#'   \item n1, n2, n3: sample sizes for each study arm
#'   \item y1, y2, y3: outcomes (number of events) for each study arm
#' }
#'   There are 9 interventions corresponding to:
#' \itemize{
#'   \item treatment 1: usual care
#'   \item treatment 2: education
#'   \item treatment 3: education + provision of free/low cost equipment
#'   \item treatment 4: education + provision of free/low cost equipment +
#'   home safety inspection
#'   \item treatment 5: education + provision of free/low cost equipment +
#'   fitting of equipment
#'   \item treatment 6: education + home safety inspection
#'   \item treatment 7: education + provision of free/low cost equipment +
#'   home safety inspection + fitting of equipment
#'   \item treatment 8: education + home visit
#'   \item treatment 9: provision of free/low cost equipment
#' }
#'
#' @source \url{http://www.biomedcentral.com/1471-2288/14/92}
#' @name homesafety
"homesafety"

#' Alcohol dependence data
#'
#' Data from a A systematic review of pharmacological treatments for alcohol
#' dependence, published in two consecutive Cochrane reviews.
#'
#' @format A data frame with 87 rows and 12 variables. The variables are as
#' follows:
#'
#' \itemize{
#'   \item studyid: study identification number
#'   \item study: study progressive number in the sample
#'   \item arm: single arm of a specific study
#'   \item trt: intervention in the specific study arm; 1 means the
#'   control/baseline treatment
#'   \item baseline: baseline treatment for a specific study
#'   \item narms: total number of arms for a specific study
#'   \item n1, n2, n3: sample sizes for each study arm
#'   \item y1, y2, y3: outcomes (number of events) for each study arm
#' }
#'   There are 4 interventions corresponding to:
#' \itemize{
#'   \item treatment 1: placebo
#'   \item treatment 2: naltrexone
#'   \item treatment 3: acamprosate
#'   \item treatment 4: naltrexone + acamprosate
#' }
#'
#' @source Liu, Y., DeSantis, S. M., and Chen, Y. (2018) Bayesian mixed
#' treatment comparisons meta-analysis for correlated outcomes subject to
#' reporting bias. Journal of the Royal Statistical Society C (67) 127–144
#' @name alcoholdependence
"alcoholdependence"
