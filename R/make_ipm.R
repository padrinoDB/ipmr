#' @title Methods to implement an IPM
#' @rdname make_ipm
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
#' the kernel expressions are evaluated in. This is useful for developer
#' debugging and not much else.
#' @param domain_list An optional list of new domain information to implement
#' the IPM with.
#' @param usr_funs An optional list of user-specified functions that are passed
#' on to the evaluation environments. This can help make vital rate expressions
#' more concise and expressive. Names in this list should exactly match the names
#' of the function calls in the \code{...} or \code{formula}.
#' @param iterate A logical indicating whether or not iterate the model during
#' or just return the iteration kernels. For density dependent (\code{dd}) and/or
#' stochastic parameter resampled models (\code{stoch_param}), this should always
#' be \code{TRUE}.
#' @param iterations If \code{iterate} is \code{TRUE}, then the number of iterations
#' to simulate.
#' @param kernel_seq The sequence of kernels to use during the iterations.
#' This can either be a vector of integers corresponding to kernel indices,
#' a character vector corresponding to kernel names, a Markov chain matrix with
#' transition probabilities between given states (NOT YET IMPLEMENTED), or empty.
#' If empty, a random sequence will be generated internally from a uniform
#' distribution.
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
                                   iterate = FALSE,
                                   ...) {
  # checks pop_state, env_state, domain_definitions
  .check_ipm_definition(proto_ipm, iterate)

  # Split out K from others so it isn't evaluated until we're ready. If it
  # isn't there, then proceed as usual

  K_row <- which(grepl("K", proto_ipm$kernel_id))

  if(length(K_row) > 0) {
    k_row <- proto_ipm[K_row, ]
    others <- proto_ipm[-c(K_row), ]
  } else {
    others <- proto_ipm
  }

  # Initialize the master_environment so these values can all be found at
  # evaluation time
  if(is.null(domain_list)){
    master_env <- .generate_master_env(others$domain, usr_funs)
  } else {
    master_env <- .generate_master_env(domain_list, usr_funs)
  }

  # construct the kernels from their function defintions
  env_list <- list(master_env = master_env)

  all_sub_kerns <- .make_sub_kernel(others,
                                    env_list,
                                    return_envs = return_all)

  sub_kern_list <- all_sub_kerns$sub_kernels


  if(length(K_row) > 0) {
    iterators <- .make_k_simple(k_row, proto_ipm, sub_kern_list, master_env)
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
                                          kernel_seq = NULL,
                                          usr_funs = list,
                                          ...) {


  # checks pop_state, env_state, domain_definitions
  .check_ipm_definition(proto_ipm, iterate)

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

  # Initialize the master_environment so these values can all be found at
  # evaluation time
  if(is.null(domain_list)){
    master_env <- .generate_master_env(others$domain, usr_funs)
  } else {
    master_env <- .generate_master_env(domain_list, usr_funs)
  }

  # construct the kernels from their function defintions
  env_list <- list(master_env = master_env)

  all_sub_kerns <- .make_sub_kernel(others,
                                    env_list,
                                    return_envs = return_all)

  sub_kern_list <- all_sub_kerns$sub_kernels

  # build up the iteration kernels from their sub-kernels
  iterators     <- .make_k_kern_samp(k_row,
                                     proto_ipm,
                                     sub_kern_list,
                                     master_env)

  kern_seq      <-  .make_kern_seq(proto_ipm,
                                   iterators,
                                   iterations,
                                   kernel_seq)

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

#' @inheritParams make_ipm
#' @rdname make_ipm
#'
#' @export
make_ipm.simple_di_stoch_param <- function(proto_ipm,
                                           return_all = FALSE,
                                           domain_list = NULL,
                                           iterate = TRUE,
                                           iterations = 50,
                                           usr_funs = list(),
                                           ...) {


  # checks pop_state, env_state, domain definitions
  .check_ipm_definition(proto_ipm, iterate)


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

  if(any(others$has_hier_effs) | any(k_row$has_hier_effs)) {
    others <- .split_hier_effs(others)
    k_row  <- .split_hier_effs(k_row)
  }

  # Initialize the master_environment so these values can all be found at
  # evaluation time
  if(is.null(domain_list)){
    master_env <- .generate_master_env(others$domain, usr_funs)
  } else {
    master_env <- .generate_master_env(domain_list, usr_funs)
  }

  # Bind env_exprs, constants, and pop_vectors to master_env so that
  # we can always find them and avoid that miserable repitition

  master_env <- .bind_all_exprs(pop_state = others$pop_state[[1]],
                                env_state = others$env_state[[1]],
                                env_to_bind = master_env)

  env_list <- list(master_env = master_env)

  out <- .prep_param_resamp_output(others, k_row, proto_ipm)

  for(i in seq_len(iterations)) {

    # Lazy variant makes sure that whatever functions that generate parameter
    # values stochastically are only evaluated 1 time per iteration! This is so
    # multiple parameters meant to come from a joint distribution really come from
    # the joint distribution!

    sys <- .make_sub_kernel_lazy(others,
                                 master_env,
                                 return_envs = return_all)

    sub_kernels <- sys$ipm_system$sub_kernels

    sys_i <- .make_k_param_samp(sub_kernels, k_row, master_env)


    out <- .update_param_resamp_output(sub_kernels,
                                       sys_i$iterator,
                                       sys_i$pop_vec,
                                       ifelse(return_all,
                                              sys$data_envs,
                                              NA_character_),
                                       master_env,
                                       out)

    master_env <- .update_master_env(sys$pop_vec, master_env)

  }


  return(out)


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
