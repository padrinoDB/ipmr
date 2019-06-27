# make_ipm internal helpers

#' @noRd
.generate_kernel_env <- function(parameters,
                                 domain_env,
                                 param_tree) {

  kernel_env <- rlang::child_env(.parent = domain_env)
  rlang::env_bind(kernel_env,
                  !!! parameters)

  kern_quos <- .parse_vr_formulae(param_tree$vr_text,
                                  kernel_env)

  # Bind the vital rate expressions so the initial discretization can take
  # place
  rlang::env_bind_lazy(kernel_env,
                       !!! kern_quos,
                       .eval_env = kernel_env)

  invisible(kernel_env)
}

#' @noRd
.parse_vr_formulae <- function(text, kernel_env) {

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
  to_set = kernel_env)

  return(out)
}


#' @noRd
#' @importFrom purrr flatten map_dbl
.generate_domain_env <- function(domain_list, usr_funs) {

  # Inherits from whatever is 2nd on search path. all loaded functions/packges
  # should still be findable, but objects in the global environment should not
  # be to prevent overscoping!

  dom_env <- new.env(parent = as.environment(search()[2]))

  domain_list <- purrr::flatten(domain_list)

  domain_list <- domain_list[!duplicated(names(domain_list))]

  bounds <- purrr::map(domain_list, function(x) .make_domain_seqs(x))

  n_mesh_p <- purrr::map_dbl(domain_list, ~.x[3])

  names(n_mesh_p) <- paste('n_', names(domain_list), sep = "")

  rlang::env_bind(dom_env,
                  !!! n_mesh_p)

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

  rlang::env_bind(dom_env,
                  !!! usr_funs)


  invisible(dom_env)
}

.extract_kernel_from_eval_env <- function(kernel_env,
                                          kernel_id,
                                          sub_kernel_list,
                                          family,
                                          pos) {

  sub_kernel_list[[pos]] <- rlang::env_get(kernel_env, kernel_id)
  names(sub_kernel_list)[pos] <- kernel_id
  class(sub_kernel_list[[pos]]) <- family

  return(sub_kernel_list)
}

.make_domain_seqs <- function(dom_vec) {
  if(all(!is.na(dom_vec))) {
    out <- seq(dom_vec[1], dom_vec[2], length.out = dom_vec[3] + 1)
    return(out)
  } else {
    return(NULL)
  }


}

.make_kern_seq <- function(proto, kernels, iterations, env_seq) {

  if(is.null(env_seq)) {
    seq_type <- 'int_generated'
  } else {

    test <- is.matrix(env_seq)

    if(test) {
      seq_type <- 'mc_mat'

    } else {
      seq_type <- 'usr_specified'
    }
  }

  out <- switch(seq_type,
                'int_generated' = .make_internal_seq(kernels, iterations),
                'mc_mat' = .make_markov_seq(proto, kernels, env_seq, iterations),
                'usr_specified' = .make_usr_seq(kernels, env_seq, iterations))

  return(out)

}

.make_internal_seq <- function(kernels, iterations) {

  n_kerns <- length(kernels)

  out <- round(runif(iterations, min = 1, max = n_kerns))

  return(out)

}

.make_usr_seq <- function(kernels, env_seq, iterations) {

  int_test <- vapply(env_seq, function(x) is.integer(x), logical(1))

  if(!all(int_test)) {
    stop("All values in 'env_seq' must be integers.")
  }

  max_test <- max(env_seq)

  if(max_test > length(kernels)) {
    stop('Maximum value of env_seq cannot exceed the number of kernels.')
  }

  if(length(env_seq) > iterations) {
    warning("'length(env_seq)' is greater than requested 'iterations'.",
            " Simulation will only run for as many 'iterations'.")
  }

  return(env_seq)
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
