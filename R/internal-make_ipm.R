# make_ipm internal helpers

#' @noRd
.generate_kernel_env <- function(parameters, domain_env) {

  kernel_env <- rlang::child_env(.parent = domain_env)

  env_bind(kernel_env,
           !!! parameters)

  return(kernel_env)
}

#' @noRd
.parse_vr_formulae <- function(text, kern_env) {

  vr_forms <- rlang::parse_exprs(text)
  vr_quos <- rlang::enquos(vr_forms)

  out <- lapply(vr_quos,
                function(x, to_set) rlang::quo_set_env(x, to_set),
                to_set = kern_env)

  return(out)
}


#' @noRd
#' @importFrom purrr flatten
.generate_domain_env <- function(domain_list) {

  dom_env <- rlang::new_environment()

  domain_list <- purrr::flatten(domain_list)

  domain_list <- domain_list[!duplicated(names(domain_list))]

  domains <- purrr::map(domain_list, function(x) .make_domain_seqs(x))

  rlang::env_bind(dom_env,
                  !!! domains)

  base::attach(dom_env, name = 'domain_env', warn.conflicts = FALSE)

  invisible(dom_env)
}

.make_domain_seqs <- function(dom_vec) {
  if(all(!is.na(dom_vec))) {
    out <- seq(dom_vec[1], dom_vec[2], length.out = dom_vec[3])
    return(out)
  } else {
    return(NULL)
  }


}

#' @noRd
.check_pop_state <- function(proto_ipm) {

  #DEFINE ME
}

#' @noRd
.check_env_params <- function(proto_ipm) {

  #DEFINE ME
}

#' @noRd
