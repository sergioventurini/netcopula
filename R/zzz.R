.onAttach <- function(lib, pkg) {
	packageStartupMessage(sprintf("Package %s (%s) loaded.
To cite, type citation(\"%s\")", pkg, packageDescription(pkg)$Version, pkg))
}

.onUnload <- function(libpath) {
	library.dynam.unload("netcopula", libpath)
}
