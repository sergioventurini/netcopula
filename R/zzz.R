.onAttach <- function(lib, pkg) {
  packageStartupMessage(sprintf("Package %s (%s) loaded.\nTo cite, type citation(\"%s\")",
    pkg, utils::packageDescription(pkg)$Version, pkg))
}

.onUnload <- function(lib) {
  library.dynam.unload("netcopula", lib)
}
