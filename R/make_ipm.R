

make_ipm <- function(proto_ipm, ...) {
  UseMethod('make_ipm')
}


make_ipm.general_di_det <- function(proto_ipm) {

  # Split out K from others so it isn't evaluated until we're ready. If it
  # isn't there, then proceed as usual

  K_row <- which(grepl("K", proto_ipm$kernel_id))

  if(length(k_row) > 0) {
    k_row <- proto_ipm[K_row, ]
    sub_kernels <- proto_ipm[-c(K_row), ]
  } else {
    sub_kernels <- proto_ipm
  }

  # If vital rates are fit with a hierarchical model of any kind,
  # then split those out into their respective years/plots/what-have-you
  # BE SURE TO WRITE VIGNETTE ON THIS SYNTAX

  if(.has_hier_effs(others) | .has_hier_effs(k_row)) {
    others <- .split_hier_effs(others)
    k_row <- .split_hier_effs(k_row)
  }

  # Initialize the domain_environment so these values can all be found at
  # evaluation time

  domain_env <- .generate_domain_env(others$domain)

  # Loop over the kernels for evaluation
  sub_kern_list <- list()

  for(i in seq_len(dim(others)[1])) {

    param_tree <- others$params[[i]]

    # kern_env inherits from domain_env so that those variables are
    # findable at evaluation time

    kern_env <- .generate_kernel_env(param_tree$params, domain_env)

    kern_quos <- .parse_all_formulae(param_tree$formula, param_tree$vr_text)

    rlang::env_bind_lazy(kern_env,
                         !!! kern_quos,
                         .eval_env = kern_env)


    sub_kern_list[[i]] <- rlang::env_get(kern_env, "formula")
    names(sub_kern_list)[i] <- others$kernel_id[i]
  }



}
