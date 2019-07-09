# make_ipm internal helpers

.make_sub_kernel <- function(proto, env_list, return_envs = FALSE) {
  out <- list()
  master_env <- env_list$master_env

  for(i in seq_len(dim(proto)[1])) {
    param_tree <- proto$params[[i]]

    kern_env <- .generate_kernel_env(param_tree$params,
                                     master_env,
                                     param_tree)

    kern_form <- .parse_vr_formulae(param_tree$formula,
                                    kern_env)

    names(kern_form) <- proto$kernel_id[i]

    if(proto$evict[i]) {
      kern_env <- .correct_eviction(proto$evict_fun[[i]],
                                    kern_env)
    }

    rlang::env_bind_lazy(kern_env,
                         !!! kern_form,
                         .eval_env = kern_env)

    out <- .extract_kernel_from_eval_env(kern_env,
                                         proto$kernel_id[i],
                                         out,
                                         proto$params[[i]]$family,
                                         pos = i)

    if(return_envs) {
      env_list <- purrr::splice(env_list, list(kern_env))
      names(env_list)[(i + 1)] <- proto$kernel_id[i]
    }

  } # end sub-kernel construction

  res <- list(sub_kernels = out, env_list = env_list)

  return(res)

}

#' @noRd
# makes sub-kernels, but ensures that stochastic parameters are sampled from
# their respective distributions one time for each iteration. One of many reasons
# kernel resampling is preferred when a viable alternative (at least within
# the context of ipmr).

.make_sub_kernel_lazy <- function(proto, master_env, return_envs = FALSE) {

  env_state_funs <- lapply(
    proto$env_state,
    function(x, master_env) {
      temp <- x$env_quos
      out <- lapply(temp,
                    function(x, master_env) {
                      rlang::quo_set_env(x,
                                         master_env)
                    },
                    master_env = master_env)
      return(out)
    },
    master_env = master_env)

  env_state_funs <- env_state_funs[!duplicated(names(env_state_funs))]

  master_env <- .eval_env_exprs(master_env, env_state_funs)

  env_list <- list(master_env = master_env)

  sys <- .make_sub_kernel(proto, env_list, return_envs = return_envs)

  out <- list(ipm_system = sys,
              master_env = master_env)
  return(out)
}

#' @noRd
# makes sure the expressions for each stochastic parameter are evaluated
# only one time per iteration of the whole model

.eval_env_exprs <- function(master_env, env_funs) {

  nms <- names(env_funs)

  for(i in unique(nms)) {

    assign(i, rlang::eval_tidy(env_funs[[i]]), envir = master_env)

  }

  return(master_env)

}

#'
.prep_param_resamp_output <- function(others, k_row, proto_ipm) {

  out <- list(iterators = list(),
              sub_kernels = list(),
              env_list = list(),
              env_seq = NA_character_, # placeholder
              pop_state = others$pop_state[[1]], # placeholder,
              proto_ipm = proto_ipm)

  env_vars <- lapply(proto_ipm$env_state, function(x) x$env_quos)

  env_var_nms <- lapply(env_vars, names) %>%
    unlist() %>%
    unique()

  n_env_vars <- length(env_var_nms)

  env_holder <- matrix(0,
                       nrow = 1,
                       ncol = n_env_vars,
                       dimnames = list(c(NA),
                                       c(env_var_nms)))

  out$env_seq <- env_holder

  return(out)
}

#' @noRd

.update_param_resamp_output <- function(sub_kernels,
                                        iterator,
                                        pop_vec,
                                        data_envs = NA_character_,
                                        master_env,
                                        output) {

  pop_state_nms <- names(output$pop_state)

  pop_state_temp <- rlang::env_get_list(master_env,
                                        pop_state_nms,
                                        default = NA_real_)

  env_vars <- dimnames(output$env_seq)[[2]]

  env_temp <- rlang::env_get_list(master_env,
                                  env_vars,
                                  default = NA_real_) %>%
    unlist()

  if(!is.na(data_envs)) {
    out$sub_kernel_envs <- purrr::splice(out$sub_kernel_envs, data_envs)
  }

  output$env_seq <- rbind(output$env_seq, env_temp)

  output$sub_kernels <- purrr::splice(output$sub_kernels, sub_kernels)
  output$iterators <- purrr::splice(output$iterators, iterator)

  # This MUST BE GENERALIZED for multiple state variables!!!!!!!!!!!!!!!!!
  output$pop_state <- cbind(output$pop_state, pop_state_temp)

  return(output)

}

#' @noRd
# Generates evaluation environment for a sub-kernel. Assumes that all parameters
# name/value pairs have been generated. Inherits from master_env so that it can
# access domains, stochastic parameters, and population states.

.generate_kernel_env <- function(parameters,
                                 master_env,
                                 param_tree) {

  kernel_env <- rlang::child_env(.parent = master_env)
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

parse_k_formulae <- function(text, kernel_env) {

  text <- lapply(text, function(x)
    gsub('n_', 'pop_state_', x)
  )

  text <- lapply(text, function(x) {
    temp <- strsplit(x, split = '[=]')

    if(length(temp[[1]]) == 1) { # defined by iteration matrix only
      out <- temp[[1]][1]
    } else {                     # defined by iteration matrix and pop_vector
      out <- temp[[1]][2]
    }
    return(out)
  })

  out <- ipmr:::.parse_vr_formulae(text, kernel_env)

  return(out)
}

#' @noRd
#' @importFrom purrr flatten map_dbl
#'
# Rename to master_env or something like that - this doesn't strictly hold
# domain information anymore

.generate_master_env <- function(domain_list, usr_funs) {

  # Inherits from whatever is 2nd on search path. all loaded functions/packges
  # should still be findable, but objects in the global environment should not
  # be to prevent overscoping!

  master_env <- new.env(parent = as.environment(search()[2]))

  domain_list <- purrr::flatten(domain_list)

  domain_list <- domain_list[!duplicated(names(domain_list))]

  bounds <- purrr::map(domain_list, function(x) .make_domain_seqs(x))

  n_mesh_p <- purrr::map_dbl(domain_list, ~.x[3])

  names(n_mesh_p) <- paste('n_', names(domain_list), sep = "")

  rlang::env_bind(master_env,
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

    rlang::env_bind(master_env,
                    !!! domain_grid)

    sv <- strsplit(names(domain_list)[i], '_[0-9]')[[1]][1]

    nm <- paste0('cell_size_', sv, sep = "")

    h <- domain_grid[2, 1] - domain_grid[1, 1]

    assign(nm, h, envir = master_env)

  }

  rlang::env_bind(master_env,
                  !!! usr_funs)


  invisible(master_env)
}

#' @noRd
# Generates sequences for the domain of midpoint rule integration (and midpoint
# rule only!!!!!!!!!). This must be generalized for handling trapezoid, g-l, etc.

.make_domain_seqs <- function(dom_vec) {
  if(all(!is.na(dom_vec))) {
    out <- seq(dom_vec[1], dom_vec[2], length.out = dom_vec[3] + 1)
    return(out)
  } else {
    return(NULL)
  }


}



#' @noRd
# Pulls out the evaluated sub-kernels from their evaluation environment
# and splices them into a list to hold them.

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

#' @noRd
# Generates a sequence of indices to sample kernels during stochastic
# kernel resampling procedures.

.make_kern_seq <- function(proto, kernels, iterations, kernel_seq) {


  if(is.null(kernel_seq)) {

    seq_type   <- 'internal_generated'

  } else {

    test       <- is.matrix(kernel_seq)

    if(test) {

      seq_type <- 'markov_chain_mat'

    } else {

      seq_type <- 'usr_specified'
    }
  }

  out <- switch(seq_type,
                'internal_generated' = .make_internal_seq(kernels,
                                                          iterations),
                'markov_chain_mat'   = .make_markov_seq(proto,
                                                        kernels,
                                                        kernel_seq,
                                                        iterations),
                'usr_specified'      = .make_usr_seq(kernels,
                                                     kernel_seq,
                                                     iterations))

  return(out)

}

#' @importFrom stats runif
.make_internal_seq <- function(kernels, iterations) {

  n_kerns <- length(kernels)

  out     <- sample.int(seq(1, n_kerns, by =  1),
                        size = iterations,
                        replace = TRUE)

  return(out)

}

.make_usr_seq <- function(kernels, kernel_seq, iterations) {

  int_test <- vapply(kernel_seq, function(x) is.integer(x), logical(1))

  if(!all(int_test)) {
    stop("All values in 'kernel_seq' must be integers.")
  }

  max_test <- max(kernel_seq)

  if(max_test > length(kernels)) {
    stop("Maximum value of 'kernel_seq' cannot exceed the number of kernels.")
  }

  if(length(kernel_seq) > iterations) {
    warning("'length(kernel_seq)' is greater than requested 'iterations'.",
            " Simulation will only run for as many 'iterations'.")
  }

  return(kernel_seq)
}

#' @noRd
.check_ipm_definition <- function(proto_ipm ,iterate) {

  .check_pop_state(proto_ipm)
  .check_env_state(proto_ipm)

  ipm_type <- class(proto_ipm)[1]

  if(grepl('_param|dd', ipm_type) & !iterate) {
    stop("Stochastic, parameter resampled and density dependent models must be\n",
         "iterated! Set 'iterate' to 'TRUE' and re-run.")
  }

  # probably want to add more here -------


  invisible(TRUE)

}

#' @noRd
.check_pop_state <- function(proto_ipm) {

  # ipm type is always first in class(proto)
  ipm_type <- class(proto_ipm)[1]

  pop_state <- unlist(proto_ipm$pop_state) %>%
    unique()

  state_vars <- unlist(proto_ipm$state_var) %>%
    unique()

  if(!all(names(pop_state) %in% names(state_vars))) {
    stop("Names of state variables do not match names of 'pop_state'!")
  }

  # density dependent IPMs must have pop_state defined!
  if(grepl('_dd_', ipm_type)) {
    if(all(is.na(pop_state))) {
      stop("Density dependent IPMs must have 'pop_state' defined!\n",
           "See '?define_pop_state()' for more details.")
    }
  }

  invisible(TRUE)
}

#' @noRd
.check_env_state <- function(proto_ipm) {

  # ipm type is always first in class(proto)
  ipm_type <- class(proto_ipm)[1]

  env_state <- unlist(proto_ipm$env_state, recursive = FALSE)

  if(grepl('_param', ipm_type)) {
    if(all(is.na(env_state))) {
      stop("Stochastic parameter-resampled IPMs must have 'env_state' defined!\n",
           "See '?define_env_state()' for more details.")
    }

  }

  invisible(TRUE)

}

.bind_all_exprs <- function(pop_state = NA_real_,
                            env_state = NA_real_,
                            env_to_bind) {

  if(!all(is.na(pop_state))) {
    if(rlang::is_list(pop_state)) {
      pop_state <- .flatten_to_depth(pop_state, 1)
    }
  }

  if(!all(is.na(env_state))) {
    if(rlang::is_list(env_state)) {
      env_state <- .flatten_to_depth(env_state, 1)
    }
  }

  temp <- purrr::splice(pop_state, env_state)

  to_bind <- .drop_duplicated_names_and_splice(temp)

  rlang::env_bind_lazy(env_to_bind,
                       !!! to_bind,
                       .eval_env = env_to_bind)


  return(env_to_bind)

}





