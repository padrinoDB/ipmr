#' @importFrom utils packageDescription

.onAttach <- function(libname = find.package("ipmr"), pkgname = "ipmr") {

    packageStartupMessage(
      "Welcome to `ipmr`! `browseVignettes('ipmr')` to get started."
    )

}
