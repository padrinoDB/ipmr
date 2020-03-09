# make_ipm internal helpers

#' @noRd

.make_sub_kernel_general <- function(proto, env_list, return_envs = FALSE) {

  out <- list()
  master_env <- env_list$master_env

  for(i in seq_len(dim(proto)[1])) {

    if(proto$evict[i]) {
      proto[i, ] <- .correct_eviction(proto[i, ])
    }

    param_tree <- proto$params[[i]]

    kern_env         <- .generate_kernel_env(param_tree$params,
                                             master_env,
                                             param_tree)

    kern_text        <- param_tree$formula

    kern_form        <- .parse_vr_formulae(kern_text,
                                           kern_env)

    names(kern_form) <- proto$kernel_id[i]

    rlang::env_bind_lazy(kern_env,
                         !!! kern_form,
                         .eval_env = kern_env)

    temp     <- .extract_kernel_from_eval_env(kern_env,
                                              proto$kernel_id[i],
                                              out,
                                              proto$params[[i]]$family,
                                              pos = i)

    out[[i]] <- .fun_to_iteration_mat(temp[[i]],
                                      state_var_start = names(proto$domain[[i]])[1],
                                      state_var_end   = names(proto$domain[[i]])[2],
                                      master_env      = master_env)

    names(out)[i] <- proto$kernel_id[i]

    if(return_envs) {
      env_list                 <- purrr::splice(env_list, list(kern_env))
      names(env_list)[(i + 1)] <- proto$kernel_id[i]
    }

  } # end sub-kernel construction

  res <- list(sub_kernels = out, env_list = env_list)

  return(res)
}

.make_sub_kernel_general_lazy <- function(proto, master_env, return_envs = FALSE) {

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

  sys            <- .make_sub_kernel_general(proto,
                                             env_list,
                                             return_envs = return_envs)

  out            <- list(ipm_system = sys,
                         master_env = master_env)

  return(out)
}

#' @noRd

.make_sub_kernel_simple <- function(proto, env_list, return_envs = FALSE) {

  out <- list()
  master_env <- env_list$master_env

  for(i in seq_len(dim(proto)[1])) {

    if(proto$evict[i]) {
      proto[i, ] <- .correct_eviction(proto[i, ])
    }

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

    rlang::env_bind_lazy(kern_env,
                         !!! kern_form,
                         .eval_env = kern_env)

    temp     <- .extract_kernel_from_eval_env(kern_env,
                                              proto$kernel_id[i],
                                              out,
                                              proto$params[[i]]$family,
                                              pos = i)

    out[[i]] <- .fun_to_iteration_mat(temp[[i]],
                                      state_var_start = names(proto$domain[[i]])[1],
                                      state_var_end   = names(proto$domain[[i]])[2],
                                      master_env      = master_env)

    names(out)[i] <- proto$kernel_id[i]

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

  out <- paste(kern_text, ' * d_', unique(sv), sep = "")

  return(out)
}

#' @noRd
# makes sub-kernels, but ensures that stochastic parameters are sampled from
# their respective distributions one time for each iteration. One of many reasons
# kernel resampling is preferred when a viable alternative (at least within
# the context of ipmr).

.make_sub_kernel_simple_lazy <- function(proto, master_env, return_envs = FALSE,
                                         dd = 'n') {

  out <- switch(dd,
                'n' = .make_sub_kernel_simple_lazy_di(proto,
                                                      master_env,
                                                      return_envs),
                'y' = .make_sub_kernel_simple_lazy_dd(proto,
                                                      master_env,
                                                      return_envs))

  return(out)

}

.make_sub_kernel_simple_lazy_dd <- function(proto,
                                            master_env,
                                            return_envs) {


  out <- .make_sub_kernel_simple(proto,
                                 list(master_env = master_env),
                                 return_envs = return_envs)

  return(out)

}

.make_sub_kernel_simple_lazy_di <- function(proto, master_env, return_envs = FALSE) {

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

  sys            <- .make_sub_kernel_simple(proto,
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

  # Flatten and drop duplicates. Duplication is likely to occur in simple_*
  # ipms because every kernel will have the same population state associated
  # with it. General IPMs with hierarchical effects, on the other hand, will need
  # different population states for almost every kernel row.

  pop_state <- .flatten_to_depth(others$pop_state, 1L)

  pop_state <- pop_state[!duplicated(names(pop_state))]

  out <- list()
  # If pop_vec is specified, then initialize an array to hold the output
  if(!rlang::is_empty(pop_state)){
    if(rlang::is_list(pop_state)) {

      # multiple states

      for(i in seq_along(pop_state)) {

        if(rlang::is_quosure(pop_state[[i]])) {
          pop_state[[i]] <- rlang::eval_tidy(pop_state[[i]])
        }

        dim_pop_out <- ifelse(is.matrix(pop_state[[i]]),
                              dim(pop_state[[i]]),
                              length(pop_state[[i]]))

        # Need to work out exactly how to know the indexing procedure here -
        # higher dimensional kernels will have time as the 3rd, 4th, or 5th
        # dimension (so trippy!) and normal bracket notation won't necessarily
        # work without some awful if{...}else{} sequence. Right now, this will
        # only work for distinct continuous state vars

        pop_out            <- array(NA_real_, dim = c(dim_pop_out, iterations + 1))



        pop_out[ , 1]      <- pop_state[[i]]

        out[[i]]           <- pop_out
      }

      names(out)           <- gsub("_t", "", names(pop_state))

    }
  }

  return(out)
}

#' @noRd

.update_param_simple_output <- function(sub_kernels,
                                        ipm_system,
                                        data_envs = NA_character_,
                                        master_env,
                                        output,
                                        tot_iterations,
                                        current_iteration) {

  # Updates env_seq and data_environments part of output. env, perhaps confusingly,
  # refers to environment in both the programming and the biological sense

  output <- .update_env_output(output,
                               master_env        = master_env,
                               data_envs         = data_envs,
                               tot_iterations    = tot_iterations,
                               current_iteration = current_iteration)

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

.update_param_dd_output <- function(sub_kernels,
                                    ipm_system,
                                    data_envs = NA_character_,
                                    master_env,
                                    output,
                                    tot_iterations,
                                    current_iteration) {

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

.update_param_general_output <- function(sub_kernels,
                                         pop_state,
                                         data_envs,
                                         master_env,
                                         output,
                                         tot_iterations,
                                         current_iteration) {

  # Updates env_seq and data_environments part of output. env, perhaps confusingly,
  # refers to environment in both the programming and the biological sense

  output <- .update_env_output(output            = output,
                               master_env        = master_env,
                               data_envs         = data_envs,
                               tot_iterations    = tot_iterations,
                               current_iteration = current_iteration)

  # The rest differs from update_param_simple in that we don't really have
  # "iterators" - we just use the formulae in K to relate sub_kernels to pop_state
  # and solve the system numerically. Thus, iterator slot is NA_real_

  names(sub_kernels) <- paste(names(sub_kernels),
                              current_iteration,
                              sep = "_")

  output$sub_kernels <- purrr::splice(output$sub_kernels, sub_kernels)

  output$pop_state   <- pop_state

  return(output)


}

# Modifies output in place - updates environmental parameter sequence
# AND the data_envs list slot. "env" refers to both programming environments
# and to biological environments - might need to improve terminology for long
# term maintenance.

.update_env_output <- function(output,
                               master_env,
                               data_envs,
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

    names(data_envs)       <- paste(names(data_envs),
                                    current_iteration,
                                    sep = "_")

    output$sub_kernel_envs <- purrr::splice(output$sub_kernel_envs, data_envs)

  }

  return(output)
}

#' @noRd

.update_pop_state <- function(pop_history, pop_out_t_1, iteration) {



  if(is.list(pop_history)) {
    for(i in seq_along(pop_history)) {

      pop_history[[i]][ , (iteration + 1)] <- pop_out_t_1[[i]]

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
    nm <- gsub('^n_', 'pop_state_', nm)

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

  # This is very shaky - requires code to have whitespace, and will fail
  # when users don't put spaces between variables. I think that's bad practice
  # in general, but I also think the whole n_ -> pop_state naming business
  # needs a rethink. This is something that will happen once the rest of the
  # di_* methods are written, so that tests can be re-written in bulk

  text <- lapply(text, function(x) {
    gsub(' n_', ' pop_state_', x, perl = TRUE)
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

  names(text) <- gsub('^n_', 'pop_state_', names(text))
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

  # Create helper vars for user-facing formula writing
  Ls              <- purrr::map_dbl(domain_list, ~.x[1])
  Us              <- purrr::map_dbl(domain_list, ~.x[2])
  n_mesh_p        <- purrr::map_dbl(domain_list, ~.x[3])

  # Generate unique names

  names(Ls)       <- paste('L_', names(domain_list), sep = "")
  names(Us)       <- paste('U_', names(domain_list), sep = "")
  names(n_mesh_p) <- paste('n_', names(domain_list), sep = "")

  rlang::env_bind(master_env,
                  !!! Ls,
                  !!! Us,
                  !!! n_mesh_p)

  # Generate midpoints for integration mesh

  mids <- purrr::map(bounds, .f = function(x) {

    l <- length(x) - 1
    out_domain <- 0.5 * (x[1:l] + x[2:(l + 1)])
    return(out_domain)

  })

  # For general IPMs, we also need indices to extract the correct vectors
  # from the evaluated kernels. For CC and DD, these are just empty vectors because
  # we want the complete result. However, DC and CD will still generate vectors the
  # same length as CC, even though we only want the first row for CD and first column
  # for DC. CD is easy, it's just a sequence 1:n_mesh_p. DC is (0:n_mesh_p * n_mesh_p) + 1
  # to get the first value in each row (e.g. the first column of the iteration matrix).
  # Append names of the state variable to each so that we can access them later.
  # We do not want to do this for continuous states that don't exist - e.g. n_mesh_p = NA
  # This throws an error, so reduce the list to !is.na() entries.

  dc_cd_nms     <- names(domain_list)[!grepl('_not_applicable', names(domain_list))]
  n_mesh_p_cont <- n_mesh_p[!is.na(n_mesh_p)]

  for(i in seq_along(dc_cd_nms)) {

    cd_nm <- paste('cd_ind_', dc_cd_nms[i], sep = "")
    dc_nm <- paste('dc_ind_', dc_cd_nms[i], sep = "")

    to_bind <- rlang::list2(!! cd_nm := seq(1,
                                            n_mesh_p_cont[[i]],
                                            by = 1),
                            !! dc_nm := (seq(0,
                                             n_mesh_p_cont[[i]] - 1,
                                             by = 1) * n_mesh_p_cont[[i]]) + 1)

    rlang::env_bind(master_env,
                    !!! to_bind)

  }

  # Loop over the different domains. Each continuous variable will have two,
  # hence the "by = 2". This only applies for midpoint rule IPMs, others will need
  # different weights, etc.

  cont_svs <- strsplit(names(domain_list), '_[0-9]') %>% unlist()
  cont_svs <- cont_svs[!grepl('not_applicable', cont_svs)] %>%
    unique()

  for(i in cont_svs) {

    mid_ind <- grepl(i, names(mids))

    domain_grid        <- expand.grid(mids[mid_ind])

    names(domain_grid) <- c(names(mids)[mid_ind])

    rlang::env_bind(master_env,
                    !!! domain_grid)

    nm <- paste0('d_', i, sep = "")

    h  <- domain_grid[2, 1] - domain_grid[1, 1]

    assign(nm, h, envir = master_env)

  }

  # Add in user specified functions. These need to be in the master_env
  # so all kernels can access them during evaluation

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

    seq_type   <- 'NA'

  } else {

    test       <- is.matrix(kernel_seq)

    if(test) {

      seq_type <- 'markov_chain_mat'

    } else {

      seq_type <- 'usr_specified'
    }
  }

  out <- switch(seq_type,
                'NA'                 = NULL,
                'markov_chain_mat'   = .make_markov_seq(proto,
                                                        kernels,
                                                        kernel_seq,
                                                        iterations),
                'usr_specified'      = .make_usr_seq(kernels,
                                                     kernel_seq,
                                                     iterations))

  return(out)

}

.make_markov_seq <- function(proto,
                             kernels,
                             kernel_seq,
                             iterations) {

  stop('markov chain environmental sequences not yet supported',
       call. = FALSE)
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

  if(is.integer(kernel_seq)) {

    kernel_seq <- as.character(kernel_seq)

  }

  # Make sure everything in kernel_seq appears is actually an option

  nms_test <- logical(length(unique(kernel_seq)))

  for(i in seq_along(unique(kernel_seq))) {

    nms_test[i] <- any(grepl(kernel_seq[i], names(kernels)))

  }

  if(! all(nms_test)) {
    stop("Not all values of 'kern_seq' are present in kernel names. Please ",
         "check the model definition.")
  }

  if(length(kernel_seq) > iterations) {
    warning("'length(kernel_seq)' is greater than requested 'iterations'.",
            " Simulation will only run for as many 'iterations'.")
  }

  return(kernel_seq)
}

#' @noRd
.check_ipm_definition <- function(proto_ipm, iterate) {

  .check_pop_state(proto_ipm, iterate)
  .check_env_state(proto_ipm)

  ipm_type <- class(proto_ipm)[1]

  if(grepl('general|_param|dd', ipm_type) & !iterate) {
    stop("Stochastic, parameter resampled and density dependent models must be\n",
         "iterated! Set 'iterate' to 'TRUE' and re-run.")
  }


  invisible(TRUE)

}

#' @noRd

.check_pop_state <- function(proto_ipm, iterate) {

  # ipm type is always first in class(proto)
  ipm_type   <- class(proto_ipm)[1]

  pop_state  <- unlist(proto_ipm$pop_state) %>%
    unique()

  if(any(is.na(pop_state)) && iterate) {
    stop("'iterate = TRUE' but 'pop_state' is not defined!",
         call. = FALSE)
  }

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

.bind_all_constants <- function(env_state = NA_real_,
                                env_to_bind) {


  if(!all(is.na(env_state))) {
    if(rlang::is_list(env_state)) {
      env_state <- .flatten_to_depth(env_state, 1)

      to_bind <- .drop_duplicated_names_and_splice(env_state)

      rlang::env_bind(env_to_bind,
                      !!! to_bind)

    }
  }




  return(env_to_bind)

}

#' @noRd

.iterate_kerns_simple <- function(iterators,
                                  sub_kern_list,
                                  iterations,
                                  current_iteration,
                                  kern_seq,
                                  pop_state,
                                  master_env,
                                  proto_ipm,
                                  k_row) {

  if(rlang::is_quosure(pop_state)) {
    pop_state <- rlang::eval_tidy(pop_state)
  }


  for(i in seq_len(iterations)) {

    # dd_* uses iterate_kerns with "iterations = 1" and a loop
    # in the outer layer, whereas di_* uses this loop for iterating
    # through time. Thus, we need a switch to ensure that the index
    # for inserting the pop_state at t+1 is correctly drawn.

    iteration_ind <- ifelse(is.na(current_iteration),
                            i,
                            current_iteration)

    # Select kernels from kern_seq, or just use all of them if it's
    # a deterministic simulation. Similarly, we need to subset the k_rows
    # object so that .eval_general_det doesn't loop over ones which include
    # expressions for other levels of 'kern_seq'

    if(!is.null(kern_seq)) {

      if(is.character(kern_seq)) {

        # Almost identical to integer kern_seq method, but need to return
        # a character value from vapply for exact matching. Earlier version
        # used fuzzy matching which I think is too risky.

        kern_ind <- vapply(names(sub_kern_list),
                           function(x) strsplit(x, '_')[[1]][2],
                           character(1L))
        kern_ind <- which(kern_ind == kern_seq[iteration_ind])
        k_row_ind <- which(grepl(kern_seq[iteration_ind], k_row$kernel_id))


        use_kerns <- sub_kern_list[kern_ind]
        use_k     <- k_row[k_row_ind, ]

      } else {

        kern_ind <- vapply(names(sub_kern_list),
                           function(x) as.integer(strsplit(x, '_')[[1]][2]),
                           integer(1L)) %>%
          which(. == kern_seq[iteration_ind])

        k_row_ind <- which(grepl(kern_seq[iteration_ind], k_row$kernel_id))


        use_kerns <- sub_kern_list[kern_ind]

        use_k     <- k_row[k_row_ind, ]
      }

    } else {

      use_kerns <- sub_kern_list
      use_k     <- k_row

    }

    # Need to return an iterated pop_state object

    pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                      proto_ipm     = proto_ipm,
                                      sub_kern_list = use_kerns,
                                      pop_state     = pop_state,
                                      master_env    = master_env)

    # make names match pop_state and then reorder the list for easy
    # insertion
    names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

    pop_list_t_1        <- pop_list_t_1[names(pop_state)]

    pop_state           <- purrr::map2(.x = pop_state,
                                       .y = pop_list_t_1,
                                       .f = function(.x, .y, iteration) {

                                         .x[ , (iteration + 1)] <- .y

                                         return(.x)
                                       },
                                       iteration = iteration_ind
    )

    # Now, update the names so that we can bind the new n_*_t to
    # master_env so that the next iteration is evaluated correctly

    names(pop_list_t_1) <- names(pop_list_t_1) %>%
      paste(., '_t', sep = "")


    rlang::env_bind(.env = master_env,
                    !!! pop_list_t_1)

  }

  return(pop_state)
}

#' @noRd
#' @importFrom purrr map2
.iterate_kerns_general <- function(k_row,
                                   proto_ipm,
                                   sub_kern_list,
                                   iterations,
                                   current_iteration,
                                   kern_seq,
                                   pop_state,
                                   master_env) {

  # If it's an expression, make sure it's evaluated.
  # It can only be an expression or a numeric vector/matrix, otherwise
  # we get an error earlier on, so there shouldn't be an "else". This should
  # almost never be triggered though, as .init_pop_state_list *should* be
  # generating matrices/vectors before hand.

  if(rlang::is_quosures(pop_state) || rlang::is_quosure(pop_state)) {
    pop_state <- lapply(pop_state, rlang::eval_tidy)
  }

  for(i in seq_len(iterations)) {

    # stoch_param uses iterate_kerns with "iterations = 1" and a loop
    # in the outer layer, whereas stoch_kern uses this loop for iterating
    # through time. Thus, we need a switch to ensure that the index
    # for inserting the pop_state at t+1 is correctly drawn.

    iteration_ind <- ifelse(is.na(current_iteration),
                            i,
                            current_iteration)

    # Select kernels from kern_seq, or just use all of them if it's
    # a deterministic simulation. Similarly, we need to subset the k_rows
    # object so that .eval_general_det doesn't loop over ones which include
    # expressions for other levels of 'kern_seq'

    if(!is.null(kern_seq)) {

      if(is.character(kern_seq)) {

        use_kerns <- sub_kern_list[grepl(kern_seq[i], names(sub_kern_list))]

        use_k     <- k_row[grepl(kern_seq[i], k_row$kernel_id), ]

      } else {

        use_kerns <- sub_kern_list[kern_seq[i]]

        use_k     <- k_row[kern_seq[i] , ]
      }

    } else {

      use_kerns <- sub_kern_list
      use_k     <- k_row

    }

    pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                      proto_ipm     = proto_ipm,
                                      sub_kern_list = use_kerns,
                                      pop_state     = pop_state,
                                      master_env    = master_env)


    # make names match pop_state and then reorder the list for easy
    # insertion
    names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

    pop_list_t_1        <- pop_list_t_1[names(pop_state)]

    pop_state           <- purrr::map2(.x        = pop_state,
                                       .y        = pop_list_t_1,
                                       .f        = function(.x, .y, iteration) {

                                         .x[ , (iteration + 1)] <- .y

                                         return(.x)
                                       },
                                       iteration = iteration_ind
    )

    # Now, update the names so that we can bind the new n_*_t to
    # master_env so that the next iteration is evaluated correctly

    names(pop_list_t_1) <- names(pop_list_t_1) %>%
      paste(., '_t', sep = "")


    rlang::env_bind(.env = master_env,
                    !!! pop_list_t_1)

  }

  return(pop_state)

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


  # Split out K from Fothers so it isn't evaluated until we're ready. If it
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

#' @noRd
# Returns a list with entries others and k_row with hier_effs split out
# Checks ipm definition. slightly modified so that it can manipulate
# vr_text as well.

.initialize_kernels_dd <- function(proto_ipm, iterate) {

  # checks pop_state, env_state, domain definitions
  .check_ipm_definition(proto_ipm, iterate)

  # modifies vr_text so that it density dependent parts are also gsub'd

  proto_ipm <- .sub_dd_terms(proto_ipm)

  # Split out K from Fothers so it isn't evaluated until we're ready. If it
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

.sub_dd_terms <- function(proto) {

  possible_svs <- unlist(proto$state_var) %>%
    unique()

  for(i in seq_along(proto$params)) {

    kern_row <- proto$params[[i]]

    for(j in seq_along(possible_svs)) {

      test    <- paste('n_', possible_svs[j], '_t', sep = "")
      replace <- paste('pop_state_', possible_svs[j], '_t', sep = "")

      kern_row$vr_text <- lapply(kern_row$vr_text,
                                 function(vr_fun,
                                          test,
                                          to_sub) {

                                   gsub(test, to_sub, vr_fun)
                                 },
                                 test   = test,
                                 to_sub = replace)

      kern_row$formula <- lapply(kern_row$formula,
                                 function(vr_fun,
                                          test,
                                          to_sub) {
                                   gsub(test, to_sub, vr_fun)
                                 },
                                 test   = test,
                                 to_sub = replace)
    }



    proto$params[[i]] <- kern_row

  }

  return(proto)
}

.fun_to_iteration_mat <- function(fun,
                                  state_var_start,
                                  state_var_end,
                                  master_env) {

  # store class of function to append to result. It gets stripped out
  # by `[` I think.

  fun_cls <- class(fun)[1]

  test_vec <- c(state_var_start, state_var_end)

  ind_lgl  <- vapply(test_vec,
                     function(x) grepl('_not_applicable', x),
                     logical(1))

  # If there is no starting domain, its a discrete to continuous, and will
  # produce a column vector (or scalar for DD). Otherwise, get number of
  # meshpoints

  if(state_var_start == 'start_not_applicable' ||
     is.null(state_var_start) ||
     is.na(state_var_start)) {

    n_mesh_p_start <- 1

    ind <- paste('master_env$dc_ind_', state_var_end, sep = "")
    ind <- rlang::parse_expr(ind)

  } else {

    mesh_p_start_text <- paste('n_', state_var_start, sep = "")
    n_mesh_p_start    <- master_env[[mesh_p_start_text]]

  }

  # Same as above but for domain end and it makes a row vector (or scalar)

  if(state_var_end == 'end_not_applicable' ||
     is.na(state_var_end) ||
     is.null(state_var_end)) {

    n_mesh_p_end <- 1

    ind <- paste('master_env$cd_ind_', state_var_start, sep = "")
    ind <- rlang::parse_expr(ind)

  } else {

    mesh_p_end_text   <- paste('n_', state_var_end, sep = "")

    n_mesh_p_end      <- master_env[[mesh_p_end_text]]

  }

  # ind_lgl will be all false for CC and all true for DD, so if either is not true
  # then we need to use ind to subset the values contained in fun. Then, we can
  # actually make an iteration matrix

  if(!all(ind_lgl) && !all(!ind_lgl)) {
    fun <- fun[eval(ind)]
  }

  out <- matrix(fun, nrow = n_mesh_p_end, ncol = n_mesh_p_start, byrow = TRUE)

  class(out) <- c(fun_cls, class(out))

  return(out)

}

#' @noRd
# Helper to ensure returned objects have the right ipmr classes associated with
# them.

set_ipmr_classes <- function(to_set, cls = NULL) {

  # Sets ipmr_matrix as default. Used to make sure lightweight S3 system for plotting
  # works correctly!
  if(! is.null(cls)) {

    set_as <- paste('ipmr_', cls, sep = "")

  } else {

    set_as <- 'ipmr_matrix'

  }

  out <- lapply(
    to_set,
    function(x, set_as) {

      class(x) <- c(set_as, class(x))

      return(x)
    },
    set_as = set_as
  )

  return(out)

}


#' @noRd
# Function to prepare env_list, kern_seq, and pop_state. ifelse() returns
# values that are the same length as the input, so a logical(1L) input
# returns a single entry. This isn't ideal because pop_state, env_list,
# and kern_seq can all be lists of length > 1.

.prep_other_output <- function(env_list,
                               kern_seq,
                               pop_state,
                               return_all,
                               iterate) {

  out <- list()

  if(return_all) {
    out$env_ret <- env_list
  } else {
    out$env_ret <- NA_character_
  }

  if(iterate) {
    out$env_seq_ret <- kern_seq
    out$pop_ret     <- pop_state
  } else {
    out$env_seq_ret <- NA_integer_
    out$pop_ret     <- NA_real_
  }

  return(out)
}


#' @noRd
# Attaches usr_funs to the proto_ipm. These are duplicated for every kernel in
# case a user subsets out some kernels in subsequent re-implementations
# ( or something)

.append_usr_funs_to_proto <- function(proto, usr_funs) {

  # Generates a list of length dim(proto)[1] consisting of duplicates of
  # usr_funs

  proto$usr_funs <- I(
    lapply(
      seq_len(
        dim(proto)[1]
      ),
      function(x, y) y,
      y = usr_funs
    )
  )

  return(proto)

}

