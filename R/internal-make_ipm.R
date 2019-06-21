# make_ipm internal helpers

#'@noRd
.make_k_simple <- function(k_row, proto_ipm, sub_kern_list, domain_env) {
  params <- k_row$params[[1]]
  formula <- params$formula

  # set up the environment and bind the subkernels to it
  k_env <- rlang::child_env(.parent = domain_env,
                            !!! sub_kern_list)

  k_form <- ipmr:::.parse_vr_formulae(formula,
                                      k_env)
  names(k_form) <- k_row$kernel_id

  rlang::env_bind_lazy(k_env,
                       !!! k_form,
                       .eval_env = k_env)

  out <- rlang::env_get_list(k_env, nms = k_row$kernel_id)

  return(out)

}

#' @noRd
.generate_kernel_env <- function(parameters, domain_env) {

  kernel_env <- rlang::child_env(.parent = domain_env)

  rlang::env_bind(kernel_env,
                  !!! parameters)

  return(kernel_env)
}

#' @noRd
.parse_vr_formulae <- function(text, kern_env) {

  # parse the text expressions and then convert to list of depth 1.
  # This is critical as otherwise, env_bind_lazy bind a list containing the
  # expression rather than the expression itself!

  vr_forms <- lapply(text, function(x) rlang::parse_exprs(x)) %>%
    purrr::flatten()

  # convert to quosures and set the environment for evaluation
  out <-  lapply(vr_forms, function(x, to_set){
    temp <- rlang::enquo(x)
    rlang::quo_set_env(temp, to_set)
  },
  to_set = kern_env)

  return(out)
}


#' @noRd
#' @importFrom purrr flatten
.generate_domain_env <- function(domain_list) {

  # Inherits from whatever is 2nd on search path. all loaded functions/packges
  # should still be findable, but objects in the global environment should not
  # be to prevent overscoping!

  dom_env <- new.env(parent = as.environment(search()[2]))

  domain_list <- purrr::flatten(domain_list)

  domain_list <- domain_list[!duplicated(names(domain_list))]

  bounds <- purrr::map(domain_list, function(x) ipmr:::.make_domain_seqs(x))

  mids <- purrr::map(bounds, .f = function(x) {
    l <- length(x) - 1
    out_domain <- 0.5 * (x[1:l] + x[2:(l + 1)])
    return(out_domain)

  })

  for(i in seq(1, length(mids), by = 2)){

    domain_grid <- expand.grid(mids[[i]],
                               mids[[i + 1]])

    names(domain_grid) <- c(names(mids)[i], names(mids)[i + 1])

    rlang::env_bind(dom_env,
                    !!! domain_grid)

    sv <- strsplit(names(domain_list)[i], '_[0-9]')[[1]][1]

    nm <- paste0('cell_size_', sv, sep = "")

    h <- domain_grid[2, 1] - domain_grid[1, 1]

    assign(nm, h, envir = dom_env)

  }


  invisible(dom_env)
}

.make_domain_seqs <- function(dom_vec) {
  if(all(!is.na(dom_vec))) {
    out <- seq(dom_vec[1], dom_vec[2], length.out = dom_vec[3] + 1)
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
