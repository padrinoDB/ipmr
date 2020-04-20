# .make_k_*() implementations and other helpers - probably want to split out
# some of these functions into separate files after all are implemented!


#'@noRd

.make_k_simple <- function(k_rows, proto_ipm, sub_kern_list, master_env) {

  params  <- k_rows$params[[1]]
  formula <- params$formula

  # set up the environment and bind the subkernels to it
  k_env <- rlang::child_env(.parent = master_env,
                            !!! sub_kern_list)

  k_form        <- .parse_vr_formulae(formula,
                                      k_env)

  rlang::env_bind_lazy(k_env,
                       !!! k_form,
                       .eval_env = k_env)

  out           <- rlang::env_get_list(k_env, nms = k_rows$kernel_id)

  return(out)

}

#' @noRd

.make_k_kern_samp <- function(k_rows, proto_ipm, sub_kernel_list, master_env) {

  # result storage
  k_list <- list()

  for(i in seq_len(dim(k_rows)[1])) {

    k_id       <- k_rows$kernel_id[i]

    to_bind    <- .get_sub_kernels_for_k(k_id, sub_kernel_list)

    param_tree <- k_rows$params[[i]]

    # Part one of .generate_kernel env, but we don't want to bind kern_quos
    # because they do not exist yet!

    kern_env <- rlang::child_env(.parent = master_env,
                                 !!! param_tree$parameters)

    rlang::env_bind_lazy(kern_env,
                         !!! to_bind,
                         .eval_env = kern_env)

    kern_form        <- .parse_k_formulae(param_tree$formula,
                                          kern_env)

    kern_form <- kern_form[k_id]

    rlang::env_bind_lazy(kern_env,
                         !!! kern_form,
                         .eval_env = kern_env)

    temp             <- rlang::env_get_list(kern_env, nms = k_id)

    k_list[[i]]      <- temp

    names(k_list)[i] <- k_id
  }

  k_list             <- .flatten_to_depth(k_list, 1)

  return(k_list)

}

#' @noRd

.make_k_param_samp <- function(k_rows,
                               sub_kernel_list,
                               master_env) {

  # result storage
  k_list <- list()

  for(i in seq_len(dim(k_rows)[1])) {

    k_id       <- k_rows$kernel_id[i]

    to_bind    <- .get_sub_kernels_for_k(k_id, sub_kernel_list)

    param_tree <- k_rows$params[[i]]

    # Part one of .generate_kernel env, but we don't want to bind kern_quos
    # because they do not exist yet!

    kern_env <- rlang::child_env(.parent = master_env,
                                 !!! param_tree$parameters)

    rlang::env_bind_lazy(kern_env,
                         !!! to_bind,
                         .eval_env = kern_env)

    kern_form <- .parse_k_formulae(param_tree$formula,
                                   kern_env)

    # Gets both the iteration kernel and the population state vector (if applicable)

    pull_name <- ifelse(names(kern_form) == k_id, k_id, names(kern_form))

    rlang::env_bind_lazy(kern_env,
                         !!! kern_form,
                         .eval_env = kern_env)

    temp               <- rlang::env_get_list(kern_env, nms = pull_name)

    k_list[[i]]        <- temp %>%
      .flatten_to_depth(1L)

    # set up dummy names, then use the square test to figure out who's
    # the iteration kernel and who's the pop vector. This only works
    # for simple IPMs (e.g. 1 state, iterator MUST be a square), but
    # then again, it only needs to work for them

    names(k_list[[i]])      <- rep('placeholder', length(k_list))

    k_ind <- vapply(temp, function(x) dim(x)[1] == dim(x)[2], logical(1L))

    names(k_list[[i]])[k_ind]  <- 'iterators'
    names(k_list[[i]])[!k_ind] <- pull_name[!k_ind]

  }

  k_list <- .flatten_to_depth(k_list, 1L)

  return(k_list)


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
                              master_env) {

  pop_list <- list()

  n_ks <- dim(k_row)[1]

  for(i in seq_len(n_ks)) {

    id         <- k_row$kernel_id[i]

    to_bind    <- .get_sub_kernels_for_k(id, sub_kern_list)

    param_tree <- k_row$params[[i]]

    # Part one of .generate_kernel env, but we don't want to bind kern_quos
    # because they do not exist yet!

    kern_env <- rlang::child_env(.parent = master_env,
                                 !!! param_tree$params)

    rlang::env_bind_lazy(kern_env,
                         !!! to_bind,
                         .eval_env = kern_env)

    kern_form <- .parse_k_formulae(param_tree$formula,
                                   kern_env)

    pull_name <- ifelse(names(kern_form) == id, id, names(kern_form))

    rlang::env_bind_lazy(kern_env,
                         !!! kern_form,
                         .eval_env = kern_env)

    pop_list[[i]]        <- rlang::env_get_list(kern_env, nms = pull_name)

    names(pop_list[[i]]) <- pull_name

  }

  pop_list <- .flatten_to_depth(pop_list, 1)

  return(pop_list)
}

#' @noRd

.add_pop_state_to_master_env <- function(pop_state, master_env) {


  if(!all(is.na(pop_state))) {
    if(rlang::is_list(pop_state)) {
      pop_state <- .flatten_to_depth(pop_state, 1)
    }

    # We don't want to add lambda to the master env. I don't think
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
      assign(nm, time_t, envir = master_env)

      # Finally, we need to bind the complete pop_state_list with size x time
      # matrices so these can be updated at each iteration (if iteration is requested)
      # by .iterate_kerns_* functions. Not using env_bind because that will
      # sometimes return a list of zaps invisibly when used with assignment at
      # the next level up (I don't think it'd happen here since master_env is
      # explicitly returned, but just being safe!)

      assign(names(pop_state)[i], pop_state[[i]], envir = master_env)

    }
  }

  return(master_env)

}

#' @noRd

.get_sub_kernels_for_k <- function(k_id, sub_kernel_list) {

  sk_names <- names(sub_kernel_list)

  # Presuming that kernels will not have multi-part names separated by underscores,
  # and that dropping the first part will get the suffix info for matching...

  suffix <- strsplit(k_id, '_')[[1]][-1] %>%
    paste(collapse = "_")

  sk_ind <- grepl(suffix, sk_names)

  out    <- sub_kernel_list[sk_ind]

  return(out)
}
