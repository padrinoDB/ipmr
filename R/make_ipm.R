#' @rdname make_ipm
#' @title Methods to implement an IPM
#'
#' @description The \code{make_ipm.*} methods convert a \code{proto_ipm} into a
#' set of discretized kernels and population vectors. Methods have different
#' requirements, so carefully read the parameter documentation and the
#' \code{vignette("implementation", package = "ipmr")}.
#'
#' @param proto_ipm The proto_ipm object you wish to implement. This should be the
#' output of \code{add_kernel}, \code{add_K}, or the \code{define_*} functions.
#' @param ... Other arguments passed to methods
#' @param return_all A logical indicating whether to return the environments that
#' the kernel expressions are evaluated in. This is useful for developer debugging and not
#' much else.
#' @param domain_list An optional list of new domain information to implement
#' the IPM with.
#' @param usr_funs An optional list of user-specified functions that are passed
#' on to the evaluation environments. This can help make vital rate expressions
#' more concise and expressive. Names in this list should exactly match the names
#' of the function calls in the \code{...} or \code{formula}.
#' @param iterate A logical indicating whether or not iterate the model during or just
#' return the iteration kernels.
#' @param iterations If \code{iterate} is \code{TRUE}, then the number of iterations
#' to simulate.
#' @param k_seq The sequence of kernels to use during the iterations. This can either
#' be a vector of integers corresponding to kernel indices, a character vector corresponding
#' to kernel names, a Markov chain matrix with transition probabilities between
#' given states (NOT YET IMPLEMENTED), or empty. If empty, a random sequence will
#' be generated internally from a uniform distribution.
#'
#' @return The \code{make_ipm.*det} methods will always return a list of length 4
#' containing the following components:
#'
#' \itemize{
#'   \item{\strong{iterators}}{: iteration kernel(s) (if specified by \code{add_K}),
#'                             otherwise contains \code{NA}.}
#'   \item{\strong{sub_kernels}}{: the sub_kernels specified in \code{add_kernel}.}
#'   \item{\strong{pop_state}}{: population vectors stored as an instance of the
#'                              \code{pop_state} class.}
#'   \item{\strong{proto_ipm}}{: the \code{proto_ipm} object used to implement
#'                              the model.}
#' }
#'  The \code{make_ipm.*stoch} methods will always return a list of length 5
#' containing the following components:
#'
#' \itemize{
#'   \item{\strong{iterators}}{: iteration kernel(s) (if specified by \code{add_K}),
#'                             otherwise contains \code{NA}.}
#'   \item{\strong{sub_kernels}}{: the sub_kernels specified in \code{add_kernel}.}
#'   \item{\strong{pop_state}}{: population vectors stored as an instance of the
#'                              \code{pop_state} class.}
#'   \item{\strong{env_state}}{: a matrix with dimension \code{n_iterations} X 1 of
#'                              kernel indices indicating the order
#'                              in which kernels are to be/were resampled OR
#'                              a matrix with as many columns as stochastic parameters
#'                              \code{n_iterations} rows.}
#'   \item{\strong{proto_ipm}}{: the \code{proto_ipm} object used to implement
#'                              the model.}
#' }
#'
#'
#'
#'
#'
#' @author Sam Levin
#'
#' @export

make_ipm <- function(proto_ipm, return_all = FALSE, ...) {
  UseMethod('make_ipm')
}


#' @rdname make_ipm
#' @inheritParams make_ipm
#'
#' @importFrom methods hasArg
#' @export

make_ipm.simple_di_det <- function(proto_ipm,
                                   return_all = FALSE,
                                   domain_list = NULL,
                                   usr_funs = list(),
                                   ...) {

  # Split out K from others so it isn't evaluated until we're ready. If it
  # isn't there, then proceed as usual

  K_row <- which(grepl("K", proto_ipm$kernel_id))

  if(length(K_row) > 0) {
    k_row <- proto_ipm[K_row, ]
    others <- proto_ipm[-c(K_row), ]
  } else {
    others <- proto_ipm
  }

  # Initialize the domain_environment so these values can all be found at
  # evaluation time
  if(is.null(domain_list)){
    domain_env <- .generate_domain_env(others$domain, usr_funs)
  } else {
    domain_env <- .generate_domain_env(domain_list, usr_funs)
  }
  # Loop over the kernels for evaluation
  sub_kern_list <- list()
  env_list <- list(dom_env = domain_env)

  for(i in seq_len(dim(others)[1])) {

    param_tree <- others$params[[i]]

    # kern_env inherits from domain_env so that those variables are
    # findable at evaluation time

    kern_env <- .generate_kernel_env(param_tree$params,
                                     domain_env,
                                     param_tree)

    kern_form <- .parse_vr_formulae(param_tree$formula,
                                    kern_env)
    names(kern_form) <- others$kernel_id[i]


    if(others$evict[i] &
       rlang::is_quosure(others$evict_fun[[i]][[1]])) {

      # modifies the kernel
      kern_env <- .correct_eviction(others$evict_fun[[i]][[1]],
                                    kern_env)
    }

    rlang::env_bind_lazy(kern_env,
                         !!! kern_form,
                         .eval_env = kern_env)

    sub_kern_list <- .extract_kernel_from_eval_env(kern_env,
                                                   others$kernel_id[i],
                                                   sub_kern_list,
                                                   others$params[[i]]$family,
                                                   pos = i)

    if(return_all) {
      env_list <- purrr::splice(env_list, list(kern_env))
      names(env_list)[(i + 1)] <- others$kernel_id[i]
    }

  } # End sub-kernel construction

  if(length(K_row) > 0) {
    iterators <- .make_k_simple(k_row, proto_ipm, sub_kern_list, domain_env)
  } else {
    iterators <- NA_real_
  }

  out <- list(iterators = iterators,
              sub_kernels = sub_kern_list,
              data_envs = ifelse(return_all,
                                 env_list,
                                 NA),
              pop_state = ifelse(methods::hasArg(pop_state),
                                 pop_state,
                                 NA),
              proto_ipm = proto_ipm)

  return(out)

}

#' @rdname make_ipm
#'
#' @export
make_ipm.simple_di_stoch_kern <- function(proto_ipm,
                                          return_all = FALSE,
                                          domain_list = NULL,
                                          iterate = FALSE,
                                          iterations = 50,
                                          k_seq = NULL,
                                          usr_funs = list(),
                                          ...) {

  # Split out K from others so it isn't evaluated until we're ready. If it
  # isn't there, then proceed as usual

  K_row <- which(grepl("K", proto_ipm$kernel_id))

  if(length(K_row) > 0) {
    k_row  <- proto_ipm[K_row, ]
    others <- proto_ipm[-c(K_row), ]
  } else {
    others <- proto_ipm
  }

  # If vital rates are fit with a hierarchical model of any kind,
  # then split those out into their respective years/plots/what-have-you
  # BE SURE TO WRITE VIGNETTE ON THIS SYNTAX  ONCE IMPLEMENTED

  if(any(others$has_hier_effs) | any(k_row$has_hier_effs)) {
    others <- .split_hier_effs(others)
    k_row  <- .split_hier_effs(k_row)
  }

  # Initialize the domain_environment so these values can all be found at
  # evaluation time
  if(is.null(domain_list)){
    domain_env <- .generate_domain_env(others$domain, usr_funs)
  } else {
    domain_env <- .generate_domain_env(domain_list, usr_funs)
  }
  # Loop over the kernels for evaluation
  sub_kern_list <- list()
  env_list <- list(dom_env = domain_env)

  for(i in seq_len(dim(others)[1])) {

    param_tree <- others$params[[i]]

    # kern_env inherits from domain_env so that those variables are
    # findable at evaluation time

    kern_env <- .generate_kernel_env(param_tree$params,
                                     domain_env,
                                     param_tree)

    kern_form <- .parse_vr_formulae(param_tree$formula,
                                    kern_env)
    names(kern_form) <- others$kernel_id[i]

    if(others$evict[i] &
       rlang::is_quosure(others$evict_fun[[i]])) {

      # modifies the kernel
      kern_env <- .correct_eviction(others$evict_fun[[i]],
                                    kern_env)
    }

    rlang::env_bind_lazy(kern_env,
                         !!! kern_form,
                         .eval_env = kern_env)

    sub_kern_list <- .extract_kernel_from_eval_env(kern_env,
                                                   others$kernel_id[i],
                                                   sub_kern_list,
                                                   others$params[[i]]$family,
                                                   pos = i)
    if(return_all) {
      env_list[[i + 1]] <- kern_env
      names(env_list)[i + 1] <- others$kernel_id[i]
    }
  }

  iterators <- .make_k_kern_samp(k_row, proto_ipm, sub_kern_list, domain_env)

  kern_seq <- .make_kern_seq(proto_ipm, iterators, iterations, k_seq)

  if(iterate) {

    ## DEFINE ME!!
    out <- .iterate_kerns(iterators,
                          iterations,
                          kern_seq,
                          pop_state)

  } else {
    out <- list(iterators = iterators,
                sub_kernels = sub_kern_list,
                env_list = ifelse(return_all, env_list, NA),
                env_seq = kern_seq,
                pop_state = ifelse(methods::hasArg(pop_state),
                                   pop_state,
                                   NA),
                proto_ipm = proto_ipm)
  }
  return(out)
}

make_ipm.simple_di_stoch_param <- function(proto_ipm,
                                           return_all = FALSE,
                                           domain_list = NULL,
                                           iterate = FALSE,
                                           iterations = 50,
                                           usr_funs = list(),
                                           ...) {

  # Split out K from others so it isn't evaluated until we're ready. If it
  # isn't there, then proceed as usual

  K_row <- which(grepl("K", proto_ipm$kernel_id))

  if(length(K_row) > 0) {
    k_row  <- proto_ipm[K_row, ]
    others <- proto_ipm[-c(K_row), ]
  } else {
    others <- proto_ipm
  }

  # If vital rates are fit with a hierarchical model of any kind,
  # then split those out into their respective years/plots/what-have-you
  # BE SURE TO WRITE VIGNETTE ON THIS SYNTAX  ONCE IMPLEMENTED

  if(any(others$has_hier_effs) | any(k_row$has_hier_effs)) {
    others <- .split_hier_effs(others)
    k_row  <- .split_hier_effs(k_row)
  }

}

make_ipm.general_di_det <- function(proto_ipm, ...) {

  # DEFINE ME
}

make_ipm.general_di_stoch_kern <- function(proto_ipm, ...) {

  # DEFINE ME
}

make_ipm.general_di_stoch_param <- function(proto_ipm, ...) {

  # DEFINE ME
}

# Density dependent methods----------

make_ipm.simple_dd_det <- function(proto_ipm, ...) {

  # DEFINE ME
}

make_ipm.simple_dd_stoch_kern <- function(proto_ipm, ...) {

  # DEFINE ME
}

make_ipm.simple_dd_stoch_param <- function(proto_ipm, ...) {

  # DEFINE ME
}

make_ipm.general_dd_det <- function(proto_ipm, ...) {

  # DEFINE ME
}

make_ipm.general_dd_stoch_kern <- function(proto_ipm, ...) {

  # DEFINE ME
}

make_ipm.general_dd_stoch_param <- function(proto_ipm, ...) {

  # DEFINE ME
}
