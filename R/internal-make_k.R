# .make_k_*() implementations


#'@noRd
.make_k_simple <- function(k_rows, proto_ipm, sub_kern_list, domain_env) {
  params <- k_rows$params[[1]]
  formula <- params$formula

  # set up the environment and bind the subkernels to it
  k_env <- rlang::child_env(.parent = domain_env,
                            !!! sub_kern_list)

  k_form <- .parse_vr_formulae(formula,
                               k_env)
  names(k_form) <- k_rows$kernel_id

  rlang::env_bind_lazy(k_env,
                       !!! k_form,
                       .eval_env = k_env)

  out <- rlang::env_get_list(k_env, nms = k_rows$kernel_id)

  return(out)

}

.make_k_kern_samp <- function(k_rows, proto_ipm, sub_kernel_list, domain_env) {

  # result storage
  k_list <- list()

  for(i in seq_len(dim(k_rows)[1])) {

    k_id <- k_rows$kernel_id[i]

    to_bind <- .get_sub_kernels_for_k(k_id, sub_kernel_list)

    param_tree <- k_rows$params[[i]]

    # Part one of .generate_kernel env, but we don't want to bind kern_quos
    # because they do not exist yet!

    kern_env <- rlang::child_env(.parent = domain_env,
                                 !!! param_tree$parameters)

    rlang::env_bind_lazy(kern_env,
                         !!! to_bind,
                         .eval_env = kern_env)

    kern_form <- .parse_vr_formulae(param_tree$formula,
                                    kern_env)

    names(kern_form) <- k_id

    # does anyone ever correct for eviction at the k stage? not even sure one could...
    if(k_rows$evict[i] & rlang::is_quosure(k_rows$evict_fun[[i]])) {

      kern_env <- .correct_eviction(k_rows$evict_fun[[i]][[1]],
                                    kern_env)
    }

    rlang::env_bind_lazy(kern_env,
                         !!! kern_form,
                         .eval_env = kern_env)

    temp <- rlang::env_get_list(kern_env, nms = k_rows$kernel_id[i])

    k_list[[i]] <- temp[[1]]
    names(k_list)[i] <- k_id

  }

  return(k_list)

}


.get_sub_kernels_for_k <- function(k_id, sub_kernel_list) {

  sk_names <- names(sub_kernel_list)

  # Presuming that kernels will not have multi-part names separated by underscores,
  # and that dropping the first part will get the suffix info for matching...

  suffix <- strsplit(k_id, '_')[[1]][-1] %>%
    paste(collapse = "_")

  sk_ind <- grepl(suffix, sk_names)

  out <- sub_kernel_list[sk_ind]

  return(out)
}
