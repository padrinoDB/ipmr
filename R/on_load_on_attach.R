#' @importFrom utils packageDescription

.onAttach <- function(libname = find.package("ipmr"), pkgname = "ipmr") {

  # Warn development users that package will soon be migrating to padrinoDB

  desc <- utils::packageDescription("ipmr")

  if(!is.null(desc$RemoteType)) {

    if(isTRUE(desc$RemoteType == "github" && desc$GithubUsername == "levisc8"))

    packageStartupMessage(
      "It looks like you've installed the development vesrion of ipmr from the",
      " repository 'levisc8/ipmr'.\n",
      "ipmr will soon be moving to the padrinoDB Github organization.\n",
      "In the future, you will need to install the development",
      " version from there.\nremotes::install_github('padrinoDB/ipmr')"
    )

  }

}
