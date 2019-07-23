# make_ipm internal helpers

#' @noRd

.make_sub_kernel <- function(proto, env_list, return_envs = FALSE) {
  out <- list()
  master_env <- env_list$master_env

  for(i in seq_len(dim(proto)[1])) {
    param_tree <- proto$params[[i]]

    kern_env         <- .generate_kernel_env(param_tree$params,
                                             master_env,
                                             param_tree)

    kern_text        <- .append_dz_to_kern_form(param_tree$formula,
                                                proto,
                                                i)

    kern_form        <- .parse_vr_formulae(kern_text,
                                           kern_env)

    names(kern_form) <- proto$kernel_id[i]

    if(proto$evict[i]) {
      kern_env       <- .correct_eviction(proto$evict_fun[[i]],
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
      env_list                 <- purrr::splice(env_list, list(kern_env))
      names(env_list)[(i + 1)] <- proto$kernel_id[i]
    }

  } # end sub-kernel construction

  res <- list(sub_kernels = out, env_list = env_list)

  return(res)

}

#' @noRd

.append_dz_to_kern_form <- function(kern_text, proto, id) {

  sv <- names(proto$domain[[id]])
  sv <- gsub('_[0-9]', "", sv)

  if(all(!is.na(sv))){
    use_var <- sv[1]
  } else {
    ind <- which(!is.na(sv))
    use_var <- sv[ind]
  }

  out <- paste(kern_text, ' * cell_size_', use_var, sep = "")

  return(out)
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

      if(rlang::is_quosure(temp[[1]]) || rlang::is_quosures(temp[[1]])) {
        out <- lapply(temp,
                      function(x, master_env) {
                        rlang::quo_set_env(x,
                                           master_env)
                      },
                      master_env = master_env)
      } else {
        out <- NULL
      }

      return(out)
    },
    master_env = master_env)

  nms <- lapply(env_state_funs, names) %>% unlist()

  ind <- duplicated(nms)

  env_state_funs <- env_state_funs[!ind]

  master_env     <- .bind_env_exprs(master_env, env_state_funs)

  env_list       <- list(master_env = master_env)

  sys            <- .make_sub_kernel(proto,
                                     env_list,
                                     return_envs = return_envs)

  out            <- list(ipm_system = sys,
                         master_env = master_env)

  return(out)
}

#' @noRd
# makes sure the expressions for each stochastic parameter are evaluated
# only one time per iteration of the whole model

.bind_env_exprs <- function(master_env, env_funs) {

  nms <- lapply(env_funs, names) %>% unlist()

  for(i in seq_len(length(nms))) {

    ass_nm <- nms[i]

    assign(ass_nm, rlang::eval_tidy(env_funs[[i]][[1]]), envir = master_env)

  }

  return(master_env)

}

#' @noRd

.prep_di_output <- function(others, k_row, proto_ipm, iterations) {

  out <- list(iterators   = list(),
              sub_kernels = list(),
              env_list    = list(),
              env_seq     = NA_real_, # placeholder
              pop_state   = NA_real_,
              proto_ipm   = proto_ipm)


  out$pop_state <- .init_pop_state_list(others, iterations)

  return(out)
}

#' @noRd

.init_pop_state_list <- function(others,
                                 iterations) {

  pop_state <- others$pop_state[[1]]

  out <- list()
  # If pop_vec is specified, then initialize an array to hold the output
  if(!rlang::is_empty(pop_state)){
    if(rlang::is_list(pop_state)) {

      # multiple states

      for(i in seq_along(pop_state)) {

        dim_pop_out <- ifelse(is.matrix(pop_state[[i]]),
                              dim(pop_state[[i]]),
                              length(pop_state[[i]]))

        # Need to work out exactly how to know the indexing procedure here -
        # higher dimensional kernels will have time as the 3rd, 4th, or 5th
        # dimension (so trippy!) and normal bracket notation won't necessarily
        # work without some awful if{...}else{} sequence. Right now, this will
        # only work for distinct continuous state vars

        pop_out            <- array(NA_real_, dim = c(dim_pop_out, iterations + 1))



        pop_out[ , 1]      <- pop_state[[1]]

        out[[i]]           <- pop_out
      }

      names(out)           <- gsub("_t", "", names(pop_state))

    } else {

      # single continuous state
      dim_pop_out <- ifelse(is.matrix(pop_state),
                            dim(pop_state),
                            length(pop_state))

      # Need to work out exactly how to know the indexing procedure here -
      # higher dimensional kernels will have time as the 3rd, 4th, or 5th
      # dimension (so trippy!) and normal bracket notation won't necessarily
      # work without some awful if{...}else{} sequence. Right now, this will
      # only work for single continuous state vars

      pop_out            <- array(0, dim = c(dim_pop_out, iterations + 1))



      pop_out[ , 1]      <- eval(others$pop_state[[1]])

      out[[1]]           <- pop_out
      names(out)         <- gsub('_t', '', names(pop_state))
    }
  }

  return(out)
}

#' @noRd

.update_param_resamp_output <- function(sub_kernels,
                                        ipm_system,
                                        data_envs = NA_character_,
                                        master_env,
                                        output,
                                        tot_iterations,
                                        current_iteration) {

  # Determine if env_state is comprised of functions. If so, get whatever
  # they returned for that iteration. If not, grab the constants (I think this
  # is more useful for troubleshooting than anyone actually using it - if the
  # environment isn't varying, then they shouldn't be using this method anyway).
  if(!rlang::is_empty(names(output$proto_ipm$env_state[[1]]$env_quos))) {

    env_vars <- names(output$proto_ipm$env_state[[1]]$env_quos)

  } else {

    env_vars <- names(output$proto_ipm$env_state[[1]]$constants)
  }

  env_temp   <- rlang::env_get_list(master_env,
                                    env_vars,
                                    default = NA_real_,
                                    inherit = FALSE) %>%
    unlist()

  if(current_iteration == 1) {

    output$env_seq <- matrix(NA_real_,
                             nrow = tot_iterations,
                             ncol = length(env_temp),
                             byrow = TRUE,
                             dimnames = list(c(NULL),
                                             c(names(env_temp))))

  }

  output$env_seq[current_iteration, ] <-  env_temp

  # On to the rest of the output

  if(!all(is.na(data_envs))) {
    out$sub_kernel_envs <- purrr::splice(out$sub_kernel_envs, data_envs)
  }

  # Determine who's a kernel and who's a pop_vector!
  ipm_system <- .flatten_to_depth(ipm_system, 1)

  kern_ind   <- lapply(ipm_system, function(x) dim(x)[2] != 1) %>%
    unlist()

  iterator   <- ipm_system[kern_ind]


  # make names a bit prettier to help distinguish between iterations
  names(sub_kernels) <- paste(names(sub_kernels),
                               current_iteration,
                               sep = "_")

  names(iterator)    <- paste(names(iterator), current_iteration, sep = "_")

  output$sub_kernels <- purrr::splice(output$sub_kernels, sub_kernels)
  output$iterators   <- purrr::splice(output$iterators, iterator)

  ps_ind             <- lapply(ipm_system, function(x) dim(x)[2] == 1) %>%
    unlist()

  if(sum(ps_ind) > 0){
    output$pop_state <- .update_pop_state(output$pop_state,
                                          ipm_system[ps_ind],
                                          current_iteration)
  }


  return(output)

}

#' @noRd

.update_pop_state <- function(pop_history, pop_out_t_1, iteration) {

  pop_out_t_1 <- unlist(pop_out_t_1)

  if(is.list(pop_history)) {
    for(i in seq_along(pop_history)) {

      pop_history[[i]][ , (iteration + 1)] <- pop_out_t_1
    }
  } else if(is.matrix(pop_history)) {

    pop_history[ , (iteration + 1)]        <- pop_out_t_1
  }

  return(pop_history)

}

#' @noRd

.update_master_env <- function(pop_state, master_env, iteration) {

  for(i in seq_along(pop_state)) {
    nm <- paste0(names(pop_state), '_t', sep = "")
    nm <- gsub('n_', 'pop_state_', nm)

    assign(nm, pop_state[[i]][, (iteration + 1)], master_env)
  }

  return(master_env)
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

  kern_quos  <- .parse_vr_formulae(param_tree$vr_text,
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

  out      <-  lapply(vr_forms, function(x, to_set) {

    temp   <- rlang::enquo(x)

    rlang::quo_set_env(temp, to_set)
  },
  to_set = kernel_env)

  return(out)
}

.parse_k_formulae <- function(text, kernel_env) {

  text <- lapply(text, function(x) {
    gsub('n_', 'pop_state_', x)
    }
  )

  text   <- lapply(text, function(x) {

    temp <- strsplit(x, split = '[=]')

    if(length(temp[[1]]) == 1) { # defined by iteration matrix only

      out <- temp[[1]][1]

    } else {                     # defined by iteration matrix and pop_vector

      out <- temp[[1]][2]

    }

    return(out)
  })

  names(text) <- gsub('n_', 'pop_state_', names(text))
  out         <- .parse_vr_formulae(text, kernel_env)

  return(out)
}

#' @noRd
#' @importFrom purrr flatten map_dbl
#'
# Rename to master_env or something like that - this doesn't strictly hold
# domain information anymore

.make_master_env <- function(domain_list, usr_funs) {

  # Parent is whatever is 2nd on search path. all loaded functions/packges
  # should still be findable, but objects in the global environment should not
  # be to prevent overscoping!

  master_env      <- new.env(parent = as.environment(search()[2]))

  domain_list     <- purrr::flatten(domain_list)

  domain_list     <- domain_list[!duplicated(names(domain_list))]

  bounds          <- purrr::map(domain_list, function(x) .make_domain_seqs(x))

  n_mesh_p        <- purrr::map_dbl(domain_list, ~.x[3])

  names(n_mesh_p) <- paste('n_', names(domain_list), sep = "")

  rlang::env_bind(master_env,
                  !!! n_mesh_p)

  mids <- purrr::map(bounds, .f = function(x) {
    l <- length(x) - 1
    out_domain <- 0.5 * (x[1:l] + x[2:(l + 1)])
    return(out_domain)

  })

  for(i in seq(1, length(mids), by = 2)){

    domain_grid        <- expand.grid(mids[[i]],
                                      mids[[i + 1]])

    names(domain_grid) <- c(names(mids)[i], names(mids)[i + 1])

    rlang::env_bind(master_env,
                    !!! domain_grid)

    sv <- strsplit(names(domain_list)[i], '_[0-9]')[[1]][1]

    nm <- paste0('cell_size_', sv, sep = "")

    h  <- domain_grid[2, 1] - domain_grid[1, 1]

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

  sub_kernel_list[[pos]]        <- rlang::env_get(kernel_env, kernel_id)
  names(sub_kernel_list)[pos]   <- kernel_id
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
#' @noRd

.make_internal_seq <- function(kernels, iterations) {

  n_kerns <- length(kernels)

  out     <- sample.int(n = n_kerns,
                        size = iterations,
                        replace = TRUE)

  return(out)

}

#' @noRd

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
  ipm_type   <- class(proto_ipm)[1]

  pop_state  <- unlist(proto_ipm$pop_state) %>%
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
  ipm_type  <- class(proto_ipm)[1]

  env_state <- unlist(proto_ipm$env_state, recursive = FALSE)

  if(grepl('_param', ipm_type)) {
    if(all(is.na(env_state))) {
      stop("Stochastic parameter-resampled IPMs must have 'env_state' defined!\n",
           "See '?define_env_state()' for more details.")
    }

  }

  invisible(TRUE)

}

#' @noRd

.bind_all_constants <- function(pop_state = NA_real_,
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

  temp    <- purrr::splice(pop_state, env_state)

  to_bind <- .drop_duplicated_names_and_splice(temp)

  rlang::env_bind(env_to_bind,
                  !!! to_bind)


  return(env_to_bind)

}

#' @noRd

.iterate_kerns_simple <- function(iterators, iterations, kern_seq, pop_state) {

  if(rlang::is_quosure(pop_state)) {
    pop_state <- rlang::eval_tidy(pop_state)
  } else {
    pop_state <- pop_state[[1]][ , 1]
  }

  pop_holder  <- array(NA_real_, dim = c(length(pop_state), (iterations + 1)))

  pop_holder[ , 1] <- n_t <-  pop_state

  for(i in seq_len(iterations)) {

    .check_n_t(n_t)

    k_selector <- kern_seq[i]

    n_t_1      <- right_mult(iterators[[k_selector]], n_t)

    pop_holder[ , (i + 1)] <- n_t <- n_t_1

  }

  return(list(pop_state = pop_holder))
}

#' @noRd

.check_n_t <- function(n_t) {

  if(any(n_t < 0)) {
    stop('some elements of the population vector are less than 0!')
  }

  if(!any(is.finite(n_t))) {
    stop("some elements of the population vector have become infinite!")
  }

  invisible(TRUE)

}

#' @noRd
# Returns a list with entries others and k_row with hier_effs split out
# Checks ipm definition
.initialize_kernels <- function(proto_ipm, iterate) {

  # checks pop_state, env_state, domain definitions
  .check_ipm_definition(proto_ipm, iterate)


  # Split out K from others so it isn't evaluated until we're ready. If it
  # isn't there, then proceed as usual

  K_row    <- which(grepl("K|^n_.*?_t", proto_ipm$kernel_id))

  if(length(K_row) > 0) {

    k_row  <- proto_ipm[K_row, ]
    others <- proto_ipm[-c(K_row), ]

  } else {

    others <- proto_ipm
    k_row <- NA_character_

  }

  # If vital rates are fit with a hierarchical model of any kind,
  # then split those out into their respective years/plots/what-have-you

  if(any(others$has_hier_effs) | any(k_row$has_hier_effs)) {

    others <- .split_hier_effs(others)
    k_row  <- .split_hier_effs(k_row)

  }

  out <- list(others = others,
              k_row  = k_row)

  return(out)

}

