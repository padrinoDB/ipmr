#' @title Domains and State Variables
#' @rdname doms-n-states
#'
#' @description Functions to assist with constructing and manipulating
#' the domains of state variables in \code{ipmr}.
#'
#' @param proto_ipm The model to add the domain information to.
#' @param ... Numeric vectors of length 3 where the name corresponds to the
#' state variable, the first entry is the lower bound of the domain, the second
#' is the upper bound of the domain, and the third entry is the number of
#' meshpoints.
#'
#' @return The \code{proto_ipm} with the domain information added.
#'
#' @examples
#'
#'
#' domains <- define_domains(dbh = c(0, 10, 50), ht = c(0, 20, 80))
#'
#' my_ipm <- make_ipm(pre_defined_proto, domain_list = domains)
#'
#'@export

define_domains <- function(proto_ipm, ...) {

  doms <- rlang::enexprs(...) %>%
    lapply(eval)

  .check_domain_inputs(doms)
  proto <- .match_domains_to_kernels(proto_ipm, doms)
  return(proto)

}

# add more checks as things come up
#' @noRd
.check_domain_inputs <- function(doms) {

  lns <- vapply(doms, length, integer(1))
  if(any(lns != 3)) {
    ind <- which(lns != 3)

    msg <- paste('The following entry or entries are not the right length: ',
                 paste(names(doms)[ind], collapse = ', '),
                 '. Inputs must be numeric vectors of length 3.', sep = "")

    stop(msg)
  }

  cls <- vapply(doms, class, character(1))

  if(any(!cls %in% c('integer', 'numeric'))) {
    stop('All inputs must either be integers or real numbers.')
  }

  invisible(TRUE)
}

# Should work for all cases of domains
.match_domains_to_kernels <- function(proto_ipm, domain_list) {
  dom_names <- names(domain_list)

  res <- proto_ipm$domain
  for(i in seq_along(dom_names)) {
    nm <- dom_names[i]

    temp <- purrr::map(res,
                       function(.x) {
                         X <- .x
                         nms <- names(X)
                         ind <- grepl(nm, nms)
                         X[ind] <- domain_list[i]
                         return(X)
                       })

    res <- temp
  }


  proto_ipm$domain <- res
  return(proto_ipm)

}
