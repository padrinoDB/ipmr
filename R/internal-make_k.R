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
  names(k_form) <- k_rows$kernel_id

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

    names(kern_form) <- k_id


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

    k_list[[i]]        <- temp

    names(k_list[[i]]) <- pull_name

  }

  return(k_list)


}

.make_k_general_det <- function(k_row,
                                proto_ipm,
                                sub_kern_list,
                                master_env) {

  k_list <- list()


  n_ks <- sum(.flatten_to_depth(k_row$params, 1) %>%
                grepl("formula", names(.)))

  for(i in seq_len) {

    # FILL ME !
  }

}

#' @noRd

.add_pop_state_to_master_env <- function(pop_state, master_env) {

  if(!rlang::is_empty(pop_state)) {
    for(i in seq_along(pop_state)) {

      nm <- paste0(names(pop_state)[i], '_t', sep = "")
      assign(nm, pop_state[[i]][ , 1], envir = master_env)
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
