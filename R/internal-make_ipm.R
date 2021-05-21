# make_ipm internal helpers

#' @noRd

.make_sub_kernel_general <- function(proto, env_list, return_envs = FALSE) {

  out <- list()
  main_env <- env_list$main_env

  for(i in seq_len(dim(proto)[1])) {

    if(proto$evict[i]) {
      proto[i, ] <- .correct_eviction(proto[i, ])
    }

    param_tree <- proto$params[[i]]

    kern_env         <- .make_kernel_env(param_tree$params,
                                         main_env,
                                         proto[i, ])

    kern_text        <- param_tree$formula

    kern_form        <- .parse_vr_formulae(kern_text,
                                           kern_env,
                                           proto[i, ],
                                           main_env)

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
                                      main_env        = main_env,
                                      kern_name       = proto$kernel_id[i])

    names(out)[i] <- proto$kernel_id[i]

    if(return_envs) {
      env_list                 <- c(env_list, list(kern_env))
      names(env_list)[(i + 1)] <- proto$kernel_id[i]
    }

  } # end sub-kernel construction

  res <- list(sub_kernels = out, env_list = env_list)

  return(res)
}

.make_sub_kernel_general_lazy <- function(proto, main_env, return_envs = FALSE) {

  env_state_funs <- lapply(
    proto$env_state,
    function(x, main_env) {

      temp <- x$env_quos

      if(rlang::is_quosure(temp[[1]]) || rlang::is_quosures(temp[[1]])) {
        out <- lapply(temp,
                      function(x, main_env) {
                        rlang::quo_set_env(x,
                                           main_env)
                      },
                      main_env = main_env)
      } else {
        out <- NULL
      }

      return(out)
    },
    main_env = main_env)

  nms <- lapply(env_state_funs, names) %>% unlist()

  ind <- duplicated(nms)

  env_state_funs <- env_state_funs[!ind]

  main_env     <- .bind_env_exprs(main_env, env_state_funs)

  env_list       <- list(main_env = main_env)

  sys            <- .make_sub_kernel_general(proto,
                                             env_list,
                                             return_envs = return_envs)

  out            <- list(ipm_system = sys,
                         main_env = main_env)

  return(out)
}

#' @noRd

.make_sub_kernel_simple <- function(proto, env_list, return_envs = FALSE) {

  out <- list()
  main_env <- env_list$main_env

  for(i in seq_len(dim(proto)[1])) {

    if(proto$evict[i]) {
      proto[i, ] <- .correct_eviction(proto[i, ])
    }


    param_tree <- proto$params[[i]]
    integrate  <- param_tree$integrate

    kern_env         <- .make_kernel_env(param_tree$params,
                                         main_env,
                                         proto[i, ])

    kern_text        <- .append_dz_to_kern_form(param_tree$formula,
                                                proto,
                                                i,
                                                integrate)

    kern_form        <- .parse_vr_formulae(kern_text,
                                           kern_env,
                                           proto[i, ],
                                           main_env)

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
                                      main_env      = main_env,
                                      kern_name       = proto$kernel_id[i])

    names(out)[i] <- proto$kernel_id[i]

    if(return_envs) {
      env_list                 <- c(env_list, list(kern_env))
      names(env_list)[(i + 1)] <- proto$kernel_id[i]
    }

  } # end sub-kernel construction

  res <- list(sub_kernels = out, env_list = env_list)

  return(res)

}

#' @noRd

.append_dz_to_kern_form <- function(kern_text, proto, id,
                                    integrate) {

  # If the user has specified not to integrate, then this step is not
  # necessary

  if(!integrate) return(kern_text)
  # If discrete_extrema is used, then the d_z has already been appended
  # somewhere on that kernel. Thus, just return the kernel text

  if(proto$evict[id]) {

    quo_l <- .flatten_to_depth(proto$evict_fun[[id]], 1L)
    if(rlang::call_name(quo_l[[1]]) == "discrete_extrema") {
      return(kern_text)
    }
  }

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

.make_sub_kernel_simple_lazy <- function(proto, main_env, return_envs = FALSE,
                                         dd = 'n') {

  out <- switch(dd,
                'n' = .make_sub_kernel_simple_lazy_di(proto,
                                                      main_env,
                                                      return_envs),
                'y' = .make_sub_kernel_simple_lazy_dd(proto,
                                                      main_env,
                                                      return_envs))

  return(out)

}

.make_sub_kernel_simple_lazy_dd <- function(proto,
                                            main_env,
                                            return_envs) {


  out <- .make_sub_kernel_simple(proto,
                                 list(main_env = main_env),
                                 return_envs = return_envs)

  return(out)

}

.make_sub_kernel_simple_lazy_di <- function(proto, main_env, return_envs = FALSE) {

  env_state_funs <- lapply(
    proto$env_state,
    function(x, main_env) {

      temp <- x$env_quos

      if(rlang::is_quosure(temp[[1]]) || rlang::is_quosures(temp[[1]])) {
        out <- lapply(temp,
                      function(x, main_env) {
                        rlang::quo_set_env(x,
                                           main_env)
                      },
                      main_env = main_env)
      } else {
        out <- NULL
      }

      return(out)
    },
    main_env = main_env)

  nms <- lapply(env_state_funs, names) %>% unlist()

  ind <- duplicated(nms)

  env_state_funs <- env_state_funs[!ind]

  main_env       <- .bind_env_exprs(main_env, env_state_funs)

  env_list       <- list(main_env = main_env)

  sys            <- .make_sub_kernel_simple(proto,
                                            env_list,
                                            return_envs = return_envs)

  out            <- list(ipm_system = sys,
                         main_env = main_env)

  return(out)
}

#' @noRd
# makes sure the expressions for each stochastic parameter are evaluated
# only one time per iteration of the whole model. Creates data in 2 formats:
#
# 1. Individual values of each parameter that are bound to a correspoding symbol
# in the main environment. This means that users can reference each variable
# by NAME without using left hand side of the env_state expression. For example
# in a vital rate expression, env_params$g_r_yr becomes g_r_yr, no env_params$.
#
# 2. a list named by the left hand side of the env_state expression that contains
# all of the values it creates, also named. This is so that .update_env_output
# can grab that list, unlist it, and stick it into a matrix. Matching all of those
# things by names provided by the user would probably get a bit more convoluted
# and be error prone.

.bind_env_exprs <- function(main_env, env_funs) {

  nms <- lapply(env_funs, names) %>% unlist()
  env_funs <- .flatten_to_depth(env_funs, 1L)

  for(i in seq_along(nms)) {

    # This does the binding so that values are accessible by the names
    # the user gives them.

    temp <- rlang::eval_tidy(env_funs[[i]])

    rlang::env_bind(main_env, !!! temp)

    # This creates a list containing the same values so that .update_env_output
    # can find them to store the env_seq data to return to the user.

    ass_nm <- nms[i]

    env_param_list <- rlang::list2(!! ass_nm := temp) %>%
      .flatten_to_depth(1L)

    # This handles binding values when expressions are supplied rather
    # than functions. This is always the case in PADRINO, and may be the case
    # for interactive use as well.

    if(all(names(env_funs)[i] %in% names(env_param_list))){

      rlang::env_bind(main_env, !!! env_param_list)

    } else {

      assign(ass_nm, env_param_list, envir = main_env)

    }
  }

  return(main_env)

}

#' @noRd

.prep_di_output <- function(others, k_row, proto_ipm, iterations, normal) {

  # all placeholders.
  out <- list(iterators   = list(),
              sub_kernels = list(),
              env_list    = list(),
              env_seq     = NA_real_,
              pop_state   = NA_real_,
              proto_ipm   = proto_ipm)

  out$pop_state <- .init_pop_state_list(others, iterations, normal)

  return(out)
}

#' @noRd

.init_pop_state_list <- function(others,
                                 iterations,
                                 normal) {

  # Flatten and drop duplicates. Duplication is likely to occur in simple_*
  # ipms because every kernel will have the same population state associated
  # with it. General IPMs with parameter sets, on the other hand, will need
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

        pop_out            <- array(NA_real_, dim = c(dim_pop_out,
                                                      iterations + 1))



        pop_out[ , 1]      <- pop_state[[i]]

        out[[i]]           <- pop_out
      }

      names(out)           <- gsub("_t", "", names(pop_state))

    }
  }

  # Rescale population vector to 1 if requested.

  if(normal) {

    out <- .norm_pop_size(out, time_step = 1L)

  }

  out$lambda <- matrix(NA_real_, nrow = 1, ncol = iterations + 1)

  return(out)
}

.pop_size <- function(pop, time_step) {

  # Remove lambda slot - this shouldn't be counted or scaled!

  pop <- pop[!grepl("lambda", names(pop))]

  temp <- vapply(pop,
                 function(x, time_step) sum(x[ , time_step]),
                 numeric(1L),
                 time_step = time_step)

  out <- Reduce('+', temp, init = 0)

  return(out)

}

#' @noRd

.norm_pop_size <- function(pop_list, time_step) {

  # compute total population size and re-scale population vectors
  # by that. Do not re-scale lambda. .pop_size knows to remove this
  # column, so no need to worry about that in step 1.

  tot_size <- .pop_size(pop_list, time_step)
  lams     <- pop_list[grepl("lambda", names(pop_list))]
  pop      <- pop_list[! grepl("lambda", names(pop_list))]

  out <- lapply(
    pop,
    function(x, pop_size, time_step) {

      x[ , time_step] <- x[ , time_step] / pop_size
      return(x)

    },
    pop_size = tot_size,
    time_step = time_step
  )

  # instert lambdas back into the pop_state object
  out$lambda <- lams
  return(out)

}

#' @noRd

.update_param_output <- function(sub_kernels,
                                        pop_state,
                                        data_envs = NA_character_,
                                        main_env,
                                        output,
                                        tot_iterations,
                                        current_iteration) {

  # Updates env_seq and data_environments part of output. env, perhaps confusingly,
  # refers to environment in both the programming and the biological sense

  output <- .update_env_output(output,
                               main_env          = main_env,
                               data_envs         = data_envs,
                               tot_iterations    = tot_iterations,
                               current_iteration = current_iteration)


  output$sub_kernels <- c(output$sub_kernels, sub_kernels)

  output$pop_state   <- pop_state


  return(output)

}



# Modifies output in place - updates environmental parameter sequence
# AND the data_envs list slot. "env" refers to both programming environments
# and to biological environments - might need to improve terminology for long
# term maintenance.

.update_env_output <- function(output,
                               main_env,
                               data_envs,
                               tot_iterations,
                               current_iteration) {

  # First, figure out if this is a stoch_param model. If so, then
  # determine if env_state is comprised of functions. If so, get whatever
  # they returned for that iteration. If not, grab the constants (I think this
  # is more useful for troubleshooting than anyone actually using it - if the
  # environment isn't varying, then they shouldn't be using this method anyway).

  if(any(grepl("stoch_param", class(output$proto_ipm)))){
    if(!rlang::is_empty(names(output$proto_ipm$env_state[[1]]$env_quos))) {

      env_vars <- names(output$proto_ipm$env_state[[1]]$env_quos)

    } else {

      env_vars <- names(output$proto_ipm$env_state[[1]]$constants)
    }


    env_temp   <- rlang::env_get_list(main_env,
                                      env_vars,
                                      default = NA_real_,
                                      inherit = FALSE) %>%
      unlist()

    if(current_iteration == 1) {

      env_var_nms <- names(env_temp) %>%
        strsplit('\\.') %>%
        vapply(function(x) x[length(x)], character(1L))


      output$env_seq <- matrix(NA_real_,
                               nrow = tot_iterations,
                               ncol = length(env_temp),
                               byrow = TRUE,
                               dimnames = list(c(NULL),
                                               c(env_var_nms)))

    }

    output$env_seq[current_iteration, ] <-  env_temp
  }

  # On to the rest of the output

  if(!all(is.na(data_envs))) {

    output$sub_kernel_envs <- c(output$sub_kernel_envs, data_envs)

  }

  return(output)
}


#' @noRd
# Generates evaluation environment for a sub-kernel. Assumes that all parameters
# name/value pairs have been generated. Inherits from main_env so that it can
# access domains, stochastic parameters, and population states.

.make_kernel_env <- function(parameters,
                             main_env,
                             proto) {

  param_tree <- proto$params

  kernel_env <- rlang::child_env(.parent = main_env)
  rlang::env_bind(kernel_env,
                  !!! parameters)

  kern_quos  <- .parse_vr_formulae(param_tree[[1]]$vr_text,
                                   kernel_env,
                                   proto,
                                   main_env)

  # Bind the vital rate expressions so the initial discretization can take
  # place
  rlang::env_bind_lazy(kernel_env,
                       !!! kern_quos,
                       .eval_env = kernel_env)

  invisible(kernel_env)
}

#' @noRd
.parse_vr_formulae <- function(text,
                               kernel_env,
                               proto,
                               main_env) {

  # Check for calls to sum() in text. These need to be modified to
  # divide by the number of meshpoints over the domain of the function IF there
  # is only one domain for the function (which is what the user will want, not
  # the sum of the entire result of expand.grid(z, z1)).

  test <- vapply(text, .check_sum_calls, logical(1L))

  if(any(test)) {

    text[test] <- lapply(text[test],
                   function(x, proto_ipm, main_env) {
                     .prep_sum_calls(x,
                                     proto_ipm = proto_ipm,
                                     main_env  = main_env)
                   },
                   proto_ipm = proto,
                   main_env  = main_env)

  }

  # parse the text expressions and then convert to list of depth 1.
  # This is critical as otherwise, env_bind_lazy bind a list containing the
  # expression rather than the expression itself!

  vr_forms <- lapply(text, function(x) rlang::parse_exprs(x)) %>%
    .flatten_to_depth(1L)

  # convert to quosures and set the environment for evaluation

  out      <-  lapply(vr_forms, function(x, to_set) {

    temp   <- rlang::enquo(x)

    rlang::quo_set_env(temp, to_set)
  },
  to_set = kernel_env)

  return(out)
}

.parse_k_formulae <- function(text, kernel_env, proto, main_env) {

  # Modify forms to make sense w/ n->pop_state renaming

  text <- lapply(text, function(x) {
    gsub(' n_', ' pop_state_', x, perl = TRUE)
  }
  )

  names(text) <- gsub('^n_', 'pop_state_', names(text))
  out         <- .parse_vr_formulae(text, kernel_env, proto, main_env)

  return(out)
}

#' @noRd
#' @importFrom purrr flatten map_dbl
#'
# Rename to main_env or something like that - this doesn't strictly hold
# domain information anymore

.make_main_env <- function(domain_list, usr_funs, age_size) {

  # Parent is whatever is 2nd on search path. all loaded functions/packges
  # should still be findable, but objects in the global environment should not
  # be to prevent overscoping!

  main_env    <- new.env(parent = as.environment(search()[2]))
  domain_list <- .flatten_to_depth(domain_list, 1L)
  domain_list <- domain_list[!duplicated(names(domain_list))]

  # Filter out some common troubles w/ age X size models
  if(age_size) {
    names(domain_list) <- gsub("_(age)|_(max)", "", names(domain_list))
    names(domain_list) <- gsub("_[0-9]", "", names(domain_list))
    domain_list        <- domain_list[!duplicated(names(domain_list))] %>%
      c(., .) %>%
      setNames(paste(names(.), c(1:2), sep = "_"))

  }



  rm_ind      <- vapply(domain_list,
                        function(x) any(is.na(x)),
                        logical(1L))
  domain_list <- domain_list[ ! rm_ind ]

  # This next bit guarantees that d_1 comes before d_2 every time.
  # testing partially par_setarchical models produced a bug where d_2
  # can occur before d_1 in the names of this list in one very specific, but
  # potentially not-uncommon case where the first non-NA domain name in domain_list
  # d_2

  nm_order        <- sort(names(domain_list))
  domain_list     <- domain_list[nm_order]
  domain_list     <- domain_list[!is.na(names(domain_list))]

  # Resume as before

  bounds          <- purrr::map(domain_list, function(x) .make_domain_seqs(x))

  # Create helper vars for user-facing formula writing
  Ls              <- purrr::map_dbl(domain_list, ~.x[1])
  Us              <- purrr::map_dbl(domain_list, ~.x[2])
  n_mesh_p        <- purrr::map_dbl(domain_list, ~.x[3])

  # Generate unique names

  names(Ls)       <- paste('L_', names(domain_list), sep = "")
  names(Us)       <- paste('U_', names(domain_list), sep = "")
  names(n_mesh_p) <- paste('n_', names(domain_list), sep = "")

  rlang::env_bind(main_env,
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

    rlang::env_bind(main_env,
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

    rlang::env_bind(main_env,
                    !!! domain_grid)

    nm <- paste0('d_', i, sep = "")

    h  <- domain_grid[2, 1] - domain_grid[1, 1]

    assign(nm, h, envir = main_env)

  }

  # Add in user specified functions. These need to be in the main_env
  # so all kernels can access them during evaluation

  rlang::env_bind(main_env,
                  !!! usr_funs)

  invisible(main_env)
}

#' @noRd
# Generates sequences for the domain of midpoint rule integration (and midpoint
# rule only!!!!!!!!!). This must be generalized for handling bin to bin, cdf, etc.

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

  out <- rlang::env_get(kernel_env, kernel_id)

  # Checks for negative/NA entries

  out <- .valid_it_mat(out, kernel_id)

  sub_kernel_list[[pos]]        <- out
  names(sub_kernel_list)[pos]   <- kernel_id
  class(sub_kernel_list[[pos]]) <- family



  return(sub_kernel_list)
}

#' @noRd
# Generates a sequence of indices to sample kernels during stochastic
# kernel resampling procedures.

.make_kern_seq <- function(proto, iterations, kernel_seq) {


  if(is.null(kernel_seq) || all(is.na(kernel_seq))) {

    seq_type   <- 'NA'

  } else {

    test       <- is.matrix(kernel_seq)

    if(test) {

      seq_type <- 'markov_chain_mat'

    } else if(all(kernel_seq == 'internal')) {

      seq_type <- "internal"

    } else {

      seq_type <- 'usr_specified'
    }
  }

  out <- switch(seq_type,
                'NA'                 = NULL,
                'markov_chain_mat'   = .make_markov_seq(proto,
                                                        kernel_seq,
                                                        iterations),
                'usr_specified'      = .make_usr_seq(proto,
                                                     kernel_seq,
                                                     iterations),
                'internal'           = .make_internal_seq(proto, iterations))

  return(out)

}

# Generates a random sequence from a uniform distribution. This only gets called
# if the user doesn't specify one on their own in _stoch_kern methods

.make_internal_seq <- function(proto, iterations) {

  par_sets <- proto$par_set_indices[proto$uses_par_sets]
  par_sets <- par_sets[!duplicated(par_sets)]

  opts <- .make_par_set_indices(par_sets)

  out  <- sample(opts, size = iterations, replace = TRUE)

  return(out)
}

.make_markov_seq <- function(proto,
                             kernel_seq,
                             iterations) {

  stop('markov chain environmental sequences not yet supported',
       call. = FALSE)
}

#' @noRd

.make_usr_seq <- function(proto, kernel_seq, iterations) {

  if(is.integer(kernel_seq)) {

    kernel_seq <- as.character(kernel_seq)

  }

  # Make sure everything in kernel_seq appears is actually an option

  nms_test <- logical(length(unique(kernel_seq)))

  pos_ids  <- proto[proto$uses_par_sets, ]

  for(i in seq_along(unique(kernel_seq))) {

    nms_test[i] <- any(grepl(kernel_seq[i], pos_ids))

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

  if(grepl('_param|dd', ipm_type) & !iterate) {
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

      # Heuristic check to see if we're binding scalars or if we need to
      # pass, say a data frame or something.

      len_test <- all(vapply(env_state, function(x) length(x) == 1, logical(1L)))

      if(all(len_test)) {

        to_bind <- .flatten_to_depth(env_state, 1) %>%
          .[!duplicated(names(.))]


      } else {
        to_bind <- env_state
      }


      rlang::env_bind(env_to_bind,
                      !!! to_bind)

    }
  }




  return(env_to_bind)

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

.initialize_kernels <- function(proto_ipm, iterate, ...) {

  UseMethod(".initialize_kernels")

}

#' @noRd
# Returns a list with entries others and k_row with par_sets split out
# Checks ipm definition

.initialize_kernels.default <- function(proto_ipm, iterate, iter_dir) {

  # checks pop_state, env_state, domain definitions
  .check_ipm_definition(proto_ipm, iterate)

  # Experimental function - automatically generates the k_row object.
  # is meant to take the place of define_k.

  k_row <- .init_iteration(proto_ipm, iterate, direction = iter_dir)

  # If vital rates are fit with a par_setarchical model of any kind,
  # then split those out into their respective years/plots/what-have-you

  if(any(proto_ipm$uses_par_sets)) {

    others <- .split_par_sets(proto_ipm)
    k_row  <- .split_par_sets(k_row)

  } else {

    others <- proto_ipm

  }

  for(i in seq_len(nrow(others))) {

    names(others$domain[[i]]) <- paste(names(others$domain[[i]]),
                                       1:2,
                                       sep = "_")

  }

  out <- list(others = others,
              k_row  = k_row)

  return(out)

}

#' @noRd

.initialize_kernels.age_x_size <- function(proto_ipm, iterate, iter_dir) {

  # checks pop_state, env_state, domain definitions
  .check_ipm_definition(proto_ipm, iterate)

  # Split out K from others so it isn't evaluated until we're ready. If it
  # isn't there, then proceed as usual

  k_row    <- .init_iteration(proto_ipm, iterate, direction = iter_dir)

  if(iterate){
    k_row  <- .make_age_k_row(k_row)

    if(any(k_row$uses_par_sets)) {

      k_row  <- .split_par_sets(k_row)

    }
  }

  others <- .split_sub_kern_ages(proto_ipm)

  if(any(others$uses_par_sets)) {

    others <- .split_par_sets(others)

  }

  # Finally, we want to remove the state_var_0 and state_var_age domain names.
  # These are going to mess things up later on, and aren't really needed (unless
  # someone really wants different domains for the same state across different
  # ages...)

  for(i in seq_len(nrow(others))) {

    names(others$domain[[i]]) <- gsub("_(age)", "", names(others$domain[[i]]))
    names(others$domain[[i]]) <- gsub("_(0)", "", names(others$domain[[i]]))
    names(others$domain[[i]]) <- paste(names(others$domain[[i]]),
                                       1:2,
                                       sep = "_")

  }

  out    <- list(others = others,
                 k_row  = k_row)

  return(out)

}

#' @noRd

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

#' @noRd

.fun_to_iteration_mat <- function(fun,
                                  state_var_start,
                                  state_var_end,
                                  main_env,
                                  kern_name) {

  # store class of function to append to result. It gets stripped out
  # by `[` I think.

  fun_cls <- class(fun)[1]

  start_cd <- substr(fun_cls, 1, 1)
  end_cd   <- substr(fun_cls, 2, 2)

  # Remove "max" generated by age_size models in state names.
  state_var_start <- gsub("(_max_)", "_", state_var_start)
  state_var_end   <- gsub("(_max_)", "_", state_var_end)

  # If there is no starting domain, its a discrete to continuous, and will
  # produce a column vector (or scalar for DD). Otherwise, get number of
  # meshpoints

  if(start_cd == "D") {

    n_mesh_p_start <- 1

    ind <- paste('main_env$dc_ind_', state_var_end, sep = "")
    ind <- rlang::parse_expr(ind)

  } else {

    mesh_p_start_text <- paste('n_', state_var_start, sep = "")
    n_mesh_p_start    <- main_env[[mesh_p_start_text]]

  }

  # Same as above but for domain end and it makes a row vector (or scalar)

  if(end_cd == "D") {

    n_mesh_p_end <- 1

    ind <- paste('main_env$cd_ind_', state_var_start, sep = "")
    ind <- rlang::parse_expr(ind)

  } else {

    mesh_p_end_text   <- paste('n_', state_var_end, sep = "")

    n_mesh_p_end      <- main_env[[mesh_p_end_text]]

  }

  # ind_lgl will be all false for CC and all true for DD, so if either is not true
  # then we need to use ind to subset the values contained in fun. Then, we can
  # actually make an iteration matrix

  if(fun_cls %in% c("DC", "CD")) {
    fun <- fun[eval(ind)]
  }

  out <- matrix(fun, nrow = n_mesh_p_end, ncol = n_mesh_p_start, byrow = TRUE)

  if(fun_cls %in% c("DC", "CD")) {

    out <- .check_vec(out, fun_cls, kern_name)

  }


  class(out) <- c(fun_cls, class(out))

  return(out)

}

#' @noRd

.valid_it_mat <- function(mat, kern_name) {

  # Finally, we need to check for floating point errors that generate
  # entries slightly less than 0, and correct those.

  if(any(mat < 0)) {

    min_0 <- min(mat[mat < 0])
    max_0 <- max(mat[mat < 0])

    if(isTRUE(all.equal(min_0, 0, tolerance = 1e-15)) &&
       isTRUE(all.equal(max_0, 0, tolerance = 1e-15))) {

      mat[mat < 0] <- 0

    } else {

      msg <- paste("Negative numbers greater than expected due",
                   " to  floating point error generated building: ",
                   kern_name,". Double check model parameteriztion.",
                   sep = "")

      stop(msg, call. = FALSE)

    }

  }


  if(any(is.na(mat))) {

    msg <- paste("NAs detected in kernel: ", kern_name,
                 ". Check kernel definition and model parameters.",
                 sep = "")

    stop(msg, call. = FALSE)

  }

  return(mat)

}

#' @noRd
# Utility that will check for the edge case when a DC or CD transition is *not*
# size dependent. When this happens, the kernel is a vector of constants, but
# NA's get inserted after the first entry by default by matrix(fun,...).

.check_vec <- function(it_mat, cls, kern_name) {

  if(any(is.na(it_mat))) {

    # check how many NAs we have. In the case specified above, the iteration
    # matrix should only have 1 non-NA value. In that case, it gets inserted
    # into all NA valued matrix elements.

    test_na <- sum(is.na(it_mat))

    if((test_na + 1) == max(dim(it_mat))) {

      it_mat[is.na(it_mat)] <- it_mat[!is.na(it_mat)]

    } else {

      msg <- paste("NAs detected in kernel: ", kern_name,
                   ". Check kernel definition and model parameters.",
                   sep = "")

      stop(msg, call. = FALSE)

    }

  }

  return(it_mat)

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
                               return_main_env,
                               return_all_envs,
                               iterate) {

  out <- list()

  if(return_all_envs) {
    out$env_ret <- env_list
  } else if(return_main_env) {

    out$env_ret <- list(main_env = env_list$main_env)
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

#' @noRd
# macro to insert updated pop_state names into vital rate expressions.
# differs from prep_dd_k_exprs in that formula is a scalar here,
# and probably a list there.

.prep_dd_vr_exprs <- function(proto) {

  possible_states <- unique(names(pop_state(proto))) %>%
    gsub(pattern = "pop_state_", replacement = "", x = .)

  state_nms <- paste("n_", possible_states, "_t", sep = "")
  to_sub    <- paste("as.vector(pop_state_", possible_states, "_t)", sep = "")

  if(any(proto$uses_par_sets)) {

    par_sets <- proto$par_set_indices %>%
      .flatten_to_depth(1L) %>%
      .[!duplicated(names(.))]

    levs <- .make_par_set_indices(par_sets)
    nms  <- names(par_sets)[!names(par_sets) %in% "to_drop"] %>%
      paste(collapse = "_")

    new_stn <- character()
    new_sub <- character()

    for(i in seq_along(levs)) {
      new_sub <- c(new_sub, gsub(nms, levs[i], to_sub))
      new_stn <- c(new_stn, gsub(nms, levs[i], state_nms))
    }

    to_sub    <- unique(new_sub)
    state_nms <- unique(new_stn)
  }

  out <- proto$params

  for(i in seq_along(state_nms)){

    out <- lapply(out,
                  function(x, target, to_sub) {

      x$formula <- gsub(target, to_sub, x$formula)
      x$vr_text <- lapply(x$vr_text, function(y, target, to_sub) {
        gsub(target, to_sub, y)
      },
      target = target,
      to_sub = to_sub)

      return(x)
    },
    target = state_nms[i],
    to_sub = to_sub[i])
  }

  proto$params <- out
  return(proto)

}

#' @noRd
# Basically the same as above, but preserves list names in ... used by
# define_k (i.e. the lapply(x$formula,...) and ignores x$vr_text)

.prep_dd_k_exprs <- function(proto) {

  possible_states <- unique(unlist(proto$state_var))

  state_nms <- paste("n_", possible_states, "_t", sep = "")
  to_sub    <- paste("as.vector(pop_state_", possible_states, "_t)", sep = "")

  for(i in seq_along(state_nms)){

    out <- lapply(proto$params, function(x, target, to_sub) {

      x$formula <- lapply(x$formula, function(y, target, to_sub){
        gsub(target, to_sub, y)
      },
      target = target,
      to_sub = to_sub)


      return(x)
    },
    target = state_nms[i],
    to_sub = to_sub[i])
  }

  proto$params <- out
  return(proto)

}


#' @noRd
# Evaluates the K exprs with their population vectors. Not prefaced with .make_*
# because for general IPMs, we don't really "make" a K, just iterate it with
# the sub kernels and use the population vectors. I think it may simplify the API
# eventually to revert this back to make_* naming though

.eval_general_det <- function(k_row,
                              proto_ipm,
                              sub_kern_list,
                              pop_state,
                              main_env) {

  pop_list <- list()

  n_ks <- dim(k_row)[1]

  for(i in seq_len(n_ks)) {

    id         <- k_row$kernel_id[i]

    param_tree <- k_row$params[[i]]

    # Part one of .generate_kernel env, but we don't want to bind kern_quos
    # because they do not exist yet!

    kern_env <- rlang::child_env(.parent = main_env,
                                 !!! param_tree$params)

    rlang::env_bind_lazy(kern_env,
                         !!! sub_kern_list,
                         .eval_env = kern_env)

    kern_form <- .parse_k_formulae(param_tree$formula,
                                   kern_env,
                                   k_row,
                                   main_env)

    pull_name <- ifelse(names(kern_form) == id, id, names(kern_form))

    rlang::env_bind_lazy(kern_env,
                         !!! kern_form,
                         .eval_env = kern_env)

    pop_list[[i]]        <- rlang::env_get_list(kern_env, nms = pull_name)

    names(pop_list[[i]]) <- pull_name

  }

  pop_list <- .flatten_to_depth(pop_list, 1L)

  return(pop_list)
}

#' @noRd

.add_pop_state_to_main_env <- function(pop_state, main_env) {


  if(!all(is.na(pop_state))) {
    if(rlang::is_list(pop_state)) {
      pop_state <- .flatten_to_depth(pop_state, 1)
    }

    # We don't want to add lambda to the main env. I don't think
    # it'd cause any problems, as I don't we don't have any other objects
    # referencing that value further downstream. This is a safety precaution.

    pop_state <- pop_state[names(pop_state != 'lambda')]

    # Turn pop_states for continuous vars into column vectors. discrete
    # vars get a 1x1 matrix

    for(i in seq_along(pop_state)) {

      # We have our matrix. Next, we need to create the n_*t helper variable.
      # This is always initialized as the first column of the size x time population
      # state matrix.

      nm <- paste0(names(pop_state)[i], '_t', sep = "")
      time_t <- pop_state[[i]][ , 1]
      assign(nm, time_t, envir = main_env)

      # Finally, we need to bind the complete pop_state_list with size x time
      # matrices so these can be updated at each iteration (if iteration is requested)
      # by .iterate_kerns_* functions. Not using env_bind because that will
      # sometimes return a list of zaps invisibly when used with assignment at
      # the next level up (I don't think it'd happen here since main_env is
      # explicitly returned, but just being safe!)

      assign(names(pop_state)[i], pop_state[[i]], envir = main_env)

    }
  }

  return(main_env)

}

