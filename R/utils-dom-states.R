#' @rdname define_star
#'
#'
#' @export

define_domains <- function(proto_ipm, ...) {

  doms <- rlang::enquos(...) %>%
    lapply(FUN = function(x) rlang::eval_tidy(x))

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
                         X      <- .x
                         nms    <- names(X)
                         ind    <- grepl(nm, nms)
                         X[ind] <- domain_list[i]
                         return(X)

                       })

    res <- temp
  }

  proto_ipm$domain <- res

  return(proto_ipm)

}

#' @importFrom utils globalVariables

utils::globalVariables(c('.', "e1", "e2"), add = FALSE)
