.onLoad <- function(libname, pkgname) {
  options(error = NULL)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to BASTION!")
}

.onUnload <- function(libpath) {
  message("Unloading BASTION.")
}
