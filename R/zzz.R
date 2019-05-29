# Package-wide global variables
.netcopulaEnv <- new.env()

#' @importFrom utils globalVariables
#' @importFrom tools file_path_as_absolute
.onLoad <- function(lib, pkg){
  # needed to avoid annoying notes in R CMD CHECK
  # (see https://github.com/tidyverse/magrittr/issues/29)
  if (getRversion() >= "2.15.1") {
    # utils::globalVariables(c(".", "S", "cl", "italic", "lbl", "p_i", "p_j"))
  }
  .netcopulaEnv$path.to.me <- tools::file_path_as_absolute(lib)
  .netcopulaEnv$nlog.double.eps <- -log(.Machine[["double.eps"]])
  .netcopulaEnv$allowedfamilies <- c("binomial", "multinomial", "gaussian", "poisson")
  .netcopulaEnv$current_family <- "binomial"
}

.onAttach <- function(lib, pkg) {
  packageStartupMessage(sprintf("Package %s (%s) loaded.\nTo cite, type citation(\"%s\")",
    pkg, utils::packageDescription(pkg)$Version, pkg))
}

.onUnload <- function(lib) {
  library.dynam.unload("netcopula", lib)
}
