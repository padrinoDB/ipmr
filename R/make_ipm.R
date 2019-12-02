#' @title Methods to implement an IPM
#' @rdname make_ipm
#'
#' @description The \code{make_ipm.*} methods convert a \code{proto_ipm} into a
#' set of discretized kernels and population vectors. Methods have different
#' requirements, so be sure to read the parameter documentation. \code{
#' vignettes('Introduction to ipmr', 'ipmr')} also contains helpful information.
#'
#' @param proto_ipm The proto_ipm object you wish to implement. This should be the
#' output of \code{define_kernel}, \code{define_k}, or the \code{define_*} functions.
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
#' @param iterate A logical indicating whether or not iterate the model before exiting
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
#' @return
#'  The \code{make_ipm.*} methods will always return a list of length 6
#' containing the following components:
#'
#' \itemize{
#'   \item{\strong{iterators}}{: iteration kernel(s) (if specified by \code{define_k}),
#'                             otherwise contains \code{NA}.}
#'   \item{\strong{sub_kernels}}{: the sub_kernels specified in \code{define_kernel}.}
#'   \item{\strong{env_list}}{: a list containing the evaluation environments of
#'                            kernel. This will always be empty unless \code{
#'                            return_all} is \code{TRUE}. Mostly here for developer
#'                            debugging.}
#'   \item{\strong{env_seq}}{: a matrix with dimension \code{n_iterations} X 1 of
#'                              kernel indices indicating the order
#'                              in which kernels are to be/were resampled OR
#'                              a matrix with as many columns as stochastic parameters
#'                              \code{n_iterations} rows.}
#'   \item{\strong{pop_state}}{: population vectors stored as an instance of the
#'                              stored as a list of arrays. The first dimension
#'                              of each array corresponds to the state variable,
#'                              and the second dimsension corresponds to time
#'                              steps.}
#'   \item{\strong{proto_ipm}}{: the \code{proto_ipm} object used to implement
#'                              the model.}
#' }
#'
#' In addition to the list class, each object will have the class from
#' \code{model_class} defined in \code{init_ipm} plus \code{'_ipm'}. This is to
#' facilitate \code{print}, \code{plot}, and \code{lambda} methods. For example,
#' a \code{'simple_di_stoch_kern'} model will have the class \code{'simple_di_stoch_kern_ipm'}
#' once it has been implemented using \code{make_ipm}.
#'
#'
#'
#' @author Sam Levin
#'
#' @export

make_ipm <- function(proto_ipm,
                     return_all = FALSE,
                     usr_funs = list(),
                     ...) {

  UseMethod('make_ipm')

}


#' @rdname make_ipm
#'
#' @importFrom methods hasArg
#' @export

make_ipm.simple_di_det <- function(proto_ipm,
                                   return_all = FALSE,
                                   usr_funs = list(),
                                   ...,
                                   domain_list = NULL,
                                   iterate = FALSE,
                                   iterations = 50
) {

  # checks pop_state, env_state, domain_definitions

  proto_list <- .initialize_kernels(proto_ipm, iterate)

  others <- proto_list$others
  k_row  <- proto_list$k_row   # k_row is either an NA or a proto_ipm object

  # Initialize the master_environment so these values can all be found at
  # evaluation time

  if(is.null(domain_list)){
    master_env <- .make_master_env(others$domain, usr_funs)
  } else {
    master_env <- .make_master_env(domain_list, usr_funs)
  }

  # construct the kernels from their function defintions
  env_list      <- list(master_env = master_env)

  all_sub_kerns <- .make_sub_kernel_simple(others,
                                           env_list,
                                           return_envs = return_all)

  sub_kern_list <- all_sub_kerns$sub_kernels

  if(!is.na(k_row)[1]) {

    iterators <- .make_k_simple(k_row, proto_ipm, sub_kern_list, master_env)

  } else {

    iterators <- NA_real_

  }

  if(iterate) {

    kern_seq  <- rep(1, iterations)

    pop_state <- .iterate_kerns_simple(iterators,
                                       iterations,
                                       kern_seq,
                                       pop_state)
  }

  if(return_all) {
    out_ret <- all_sub_kerns$env_list
  } else {
    out_ret <- NA_character_
  }

  if(iterate) {
    out_seq <- kern_seq
  } else {
    out_seq <- NA_integer_
  }

  if(!all(is.na(unlist(proto_ipm$pop_state)))) {
    out_pop <- pop_state
  } else {
    out_pop <- NA_real_
  }

  iterators     <- set_ipmr_classes(iterators)
  sub_kern_list <- set_ipmr_classes(sub_kern_list)

  out <- list(iterators   = iterators,
              sub_kernels = sub_kern_list,
              env_list    = out_ret,
              env_seq     = out_seq,
              pop_state   = out_pop,
              proto_ipm   = proto_ipm)

  class(out) <- c('simple_di_det_ipm', 'list')

  return(out)

}

#' @rdname make_ipm
#'
#' @export
make_ipm.simple_di_stoch_kern <- function(proto_ipm,
                                          return_all  = FALSE,
                                          usr_funs    = list(),
                                          ...,
                                          domain_list = NULL,
                                          iterate     = FALSE,
                                          iterations  = 50,
                                          kernel_seq  = NULL) {


  # checks pop_state, env_state, domain_definitions

  proto_list <- .initialize_kernels(proto_ipm, iterate)

  others <- proto_list$others
  k_row  <- proto_list$k_row

  # Initialize the master_environment so these values can all be found at
  # evaluation time

  if(is.null(domain_list)){
    master_env <- .make_master_env(others$domain, usr_funs)
  } else {
    master_env <- .make_master_env(domain_list, usr_funs)
  }

  # construct the kernels from their function defintions

  env_list      <- list(master_env = master_env)

  all_sub_kerns <- .make_sub_kernel_simple(others,
                                           env_list,
                                           return_envs = return_all)

  sub_kern_list <- all_sub_kerns$sub_kernels
  env_list      <- all_sub_kerns$env_list

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

    init_pop_state <- .init_pop_state_list(others, iterations)

    pop_state      <- .iterate_kerns_simple(iterators,
                                            iterations,
                                            kern_seq,
                                            init_pop_state)

  } else {

    pop_state <- NA_real_

  }

  if(return_all) {
    out_ret <- env_list
  } else {
    out_ret <- NA_character_
  }

  iterators     <- set_ipmr_classes(iterators)
  sub_kern_list <- set_ipmr_classes(sub_kern_list)

  out <- list(iterators   = iterators,
              sub_kernels = sub_kern_list,
              env_list    = out_ret,
              env_seq     = kern_seq,
              pop_state   = pop_state,
              proto_ipm   = proto_ipm)


  class(out) <- c('simple_di_stoch_kern_ipm', 'list')

  return(out)
}

#' @rdname make_ipm
#'
#' @export
make_ipm.simple_di_stoch_param <- function(proto_ipm,
                                           return_all  = FALSE,
                                           usr_funs    = list(),
                                           ...,
                                           domain_list = NULL,
                                           iterate     = TRUE,
                                           iterations  = 50) {

  proto_list <- .initialize_kernels(proto_ipm, iterate)

  others <- proto_list$others
  k_row  <- proto_list$k_row

  # Initialize the master_environment so these values can all be found at
  # evaluation time

  if(is.null(domain_list)) {
    master_env <- .make_master_env(others$domain, usr_funs)
  } else {
    master_env <- .make_master_env(domain_list, usr_funs)
  }

  # Bind env_exprs, constants, and pop_vectors to master_env so that
  # we can always find them and avoid that miserable repitition

  master_env <- .bind_all_constants(pop_state   = others$pop_state[[1]],
                                    env_state   = others$env_state[[1]]$constants,
                                    env_to_bind = master_env)

  out        <- .prep_di_output(others, k_row, proto_ipm, iterations)

  # initialize the pop_state vectors in master_env so they can be found
  # at evaluation time

  master_env <- .add_pop_state_to_master_env(out$pop_state,
                                             master_env)

  # list to hold the possibly returned evaluation environments

  env_list <- list(master_env = master_env)

  for(i in seq_len(iterations)) {

    # Lazy variant makes sure that whatever functions that generate parameter
    # values stochastically are only evaluated 1 time per iteration! This is so
    # multiple parameters meant to come from a joint distribution really come from
    # the joint distribution!

    sys         <- .make_sub_kernel_lazy(others,
                                         master_env,
                                         return_envs = return_all)

    sub_kernels <- sys$ipm_system$sub_kernels

    sys_i       <- .make_k_param_samp(k_row,
                                      sub_kernels,
                                      master_env)


    out         <- .update_param_resamp_output(sub_kernels,
                                               sys_i,
                                               ifelse(return_all,
                                                      sys$data_envs,
                                                      NA_character_),
                                               master_env,
                                               out,
                                               iterations,
                                               i)

    # turn current pop_state_t_1 into pop_state_t in master_env so next computation
    # can occur
    master_env <- .update_master_env(out$pop_state,
                                     master_env,
                                     i)

  }

  out$iterators   <- set_ipmr_classes(out$iterators)
  out$sub_kernels <- set_ipmr_classes(out$sub_kernels)

  if(return_all) out$data_envs <- purrr::splice(out$data_envs, env_list)

  class(out) <- c('simple_di_stoch_param_ipm', 'list')

  return(out)


}


#' @rdname make_ipm
#'
#' @export

make_ipm.general_di_det <- function(proto_ipm,
                                    return_all  = FALSE,
                                    usr_funs    = list(),
                                    ...,
                                    domain_list = NULL,
                                    iterate     = TRUE,
                                    iterations  = 50) {

  # initialize others + k_row

  proto_list <- .initialize_kernels(proto_ipm, iterate)

  others <- proto_list$others
  k_row  <- proto_list$k_row

  if(is.null(domain_list)) {
    master_env <- .make_master_env(others$domain, usr_funs)
  } else {
    master_env <- .make_master_env(domain_list, usr_funs)
  }

  # Bind env_exprs, constants, and pop_vectors to master_env so that
  # we can always find them and avoid that miserable repitition

  master_env <- .bind_all_constants(pop_state   = others$pop_state[[1]],
                                    env_state   = others$env_state[[1]],
                                    env_to_bind = master_env)


  temp       <- .prep_di_output(others, k_row, proto_ipm, iterations)

  # Thus far, I think general_* methods will have to use iteration for lambdas
  # as I'm not sure I want to work out the correct cbind(rbind(...)) rules for
  # creating a mega-matrix. So throw an error if pop_state isn't defined
  if(all(is.na(proto_ipm$pop_state[[1]]))) {

    stop("All general_* IPMs must have a 'pop_state' defined.",
         "See ?define_pop_state() for more details." )

  } else {

    master_env  <- .add_pop_state_to_master_env(temp$pop_state, master_env)

  }

  env_list      <- list(master_env = master_env)

  all_sub_kerns <- .make_sub_kernel_general(others,
                                            env_list,
                                            return_envs = return_all)

  sub_kern_list <- all_sub_kerns$sub_kernels
  env_list      <- all_sub_kerns$env_list

  # Things switch up here from the simple_* versions. Rather than construct a mega-K
  # through rbind + cbinding, we just iterate the population vector with the
  # sub kernels. This means I don't have to overload all of the arithmetic
  # operators and keeps things simpler internally. It also means there's no
  # need to implement sparse kernels/matrices for things like an age x size IPM.

  if(iterate) {

    kern_seq  <- rep(1, iterations)

    pop_state <- .iterate_kerns_general(k_row,
                                        proto_ipm,
                                        sub_kern_list,
                                        iterations,
                                        kern_seq,
                                        temp$pop_state,
                                        master_env)
  }


  sub_kern_list <- set_ipmr_classes(sub_kern_list)

  if(return_all) {
    env_ret <- env_list
  } else {
    env_ret <- NA_character_
  }

  if(iterate) {
    env_seq_ret <- kern_seq
    pop_ret     <- pop_state
  } else {
    env_seq_ret <- NA_integer_
    pop_ret     <- temp$pop_state
  }

  out <- list(
    iterators   = NA_character_,
    sub_kernels = sub_kern_list,
    env_list    = env_ret,
    env_seq     = env_seq_ret,
    pop_state   = pop_ret,
    proto_ipm   = proto_ipm
  )

  class(out) <- c('general_di_det_ipm', 'list')

  return(out)

}

make_ipm.general_di_stoch_kern <- function(proto_ipm,
                                           return_all = FALSE,
                                           usr_funs = list(),
                                           ...) {

  # DEFINE ME
}

make_ipm.general_di_stoch_param <- function(proto_ipm,
                                            return_all = FALSE,
                                            usr_funs = list(),
                                            ...) {

  # DEFINE ME
}

# Density dependent methods----------

make_ipm.simple_dd_det <- function(proto_ipm,
                                   return_all = FALSE,
                                   usr_funs = list(),
                                   ...) {

  # DEFINE ME
}

make_ipm.simple_dd_stoch_kern <- function(proto_ipm,
                                          return_all = FALSE,
                                          usr_funs = list(),
                                          ...) {

  # DEFINE ME
}

make_ipm.simple_dd_stoch_param <- function(proto_ipm,
                                           return_all = FALSE,
                                           usr_funs = list(),
                                           ...) {

  # DEFINE ME
}

make_ipm.general_dd_det <- function(proto_ipm,
                                    return_all = FALSE,
                                    usr_funs = list(),
                                    ...) {

  # DEFINE ME
}

make_ipm.general_dd_stoch_kern <- function(proto_ipm,
                                           return_all = FALSE,
                                           usr_funs = list(),
                                           ...) {

  # DEFINE ME
}

make_ipm.general_dd_stoch_param <- function(proto_ipm,
                                            return_all = FALSE,
                                            usr_funs = list(),
                                            ...) {

  # DEFINE ME
}
