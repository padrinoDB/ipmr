#' @title Methods to implement an IPM
#' @rdname make_ipm
#'
#' @description The \code{make_ipm.*} methods convert a \code{proto_ipm} into a
#' set of discretized kernels and population vectors. Methods have different
#' requirements, so be sure to read the parameter documentation. \code{
#' vignette('ipmr-introduction', 'ipmr')} also contains helpful information.
#'
#' @param proto_ipm The proto_ipm. This should be the
#' output of \code{define_kernel}, \code{define_k}, or the \code{define_*} functions.
#' @param ... Other arguments passed to methods
#' @param return_all A logical indicating whether to return the environments that
#' the kernel expressions are evaluated in. This is useful for developer
#' debugging and not much else.
#' @param domain_list An optional list of new domain information to implement
#' the IPM with.
#' @param usr_funs An optional list of user-specified functions that are passed
#' on to the  model building process. This can help make vital rate expressions
#' more concise and expressive. Names in this list should exactly match the names
#' of the function calls in the \code{...} or \code{formula}.
#' @param iterate A logical indicating whether or not iterate the model before exiting
#' or just return the iteration kernels. Only applies to density independent, deterministic
#' models.
#' @param iterations If \code{iterate} is \code{TRUE}, then the number of iterations
#' to simulate.
#' @param normalize_pop_size A logical indicating whether to re-scale the population
#' vector to sum before each iteration. Default is \code{TRUE} for \code{*_di_*}
#' methods and \code{FALSE} for \code{*_dd_*} methods.
#' @param kernel_seq For \code{*_stoch_kern} methods, the sequence of kernels
#' to use during the simulation process. It should have the same number of entries
#' as the number of \code{iterations}.
#' This can either be a vector of integers corresponding to kernel names (e.g.
#' kernels for different years - \code{2011:2018}),
#' a character vector corresponding to kernel names (e.g. kernels from different
#' sites - \code{'a', 'b', 'c'}), a Markov chain matrix with
#' transition probabilities between given states (NOT YET IMPLEMENTED), or empty.
#' If it is empty, \code{make_ipm} will try to generate a sequence internally using
#' a random selection of the \code{levels_hier_effs} defined in \code{define_kernel}.
#' @param report_progress A logical indicating whether or not to periodically
#' report progress for a stochastic simulation. Does not apply to deterministic
#' methods. Default is \code{FALSE}.
#'
#'
#' @return
#'  The \code{make_ipm.*} methods will always return a list of length 6
#' containing the following components:
#'
#' \itemize{
#'   \item{\strong{iterators}}{: iteration kernel(s) (if specified by \code{define_k})
#'                             for \code{simple_*}. otherwise contains \code{NA}.}
#'   \item{\strong{sub_kernels}}{: a list of sub_kernels specified in \code{define_kernel}.}
#'   \item{\strong{env_list}}{: a list containing the evaluation environments of
#'                            kernel. This will always be empty unless \code{
#'                            return_all} is \code{TRUE}. Mostly here for developer
#'                            debugging.}
#'   \item{\strong{env_seq}}{: a matrix with dimension \code{iterations} X 1 of
#'                              kernel indices indicating the order
#'                              in which kernels are to be/were resampled OR
#'                              a matrix with as many columns as stochastic parameters
#'                              and \code{n_iterations} rows.}
#'   \item{\strong{pop_state}}{: population vectors
#'                              stored as a list of arrays. The first dimension
#'                              of each array corresponds to the state variable distribution,
#'                              and the second dimsension corresponds to time
#'                              steps.}
#'   \item{\strong{proto_ipm}}{: the \code{proto_ipm} object used to implement
#'                              the model.}
#' }
#'
#' In addition to the list class, each object will have the class from
#' \code{model_class} defined in \code{init_ipm} plus \code{'_ipm'}. This is to
#' facilitate \code{print}, \code{plot}, and \code{lambda} methods. For example,
#' a \code{'simple_di_stoch_kern'} model will have the class
#' \code{'simple_di_stoch_kern_ipm'} once it has been implemented using
#' \code{make_ipm}.
#'
#' @details When \code{kernel_seq} is a character vector, names
#' are matched using \code{grepl()}.
#' When it is an integer vector, the vector
#' is first checked to make sure they are all present in kernel names. The model
#' procedure will stop if thay are not all present.
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
#' @export

make_ipm.simple_di_det <- function(proto_ipm,
                                   return_all = FALSE,
                                   usr_funs = list(),
                                   ...,
                                   domain_list = NULL,
                                   iterate = FALSE,
                                   iterations = 50,
                                   normalize_pop_size = TRUE
) {

  # Figure out if we're dealing with a new model or an old one that is
  # being reimplemented. If the proto is not new and usr_funs isn't passed to the new
  # make_ipm call, then restore the old version. If usr_funs are passed,
  # we append them regardless of whether or not the model has been implemented.

  if(isTRUE(attr(proto_ipm, 'implemented')) && rlang::is_empty(usr_funs)) {

    if(!is.na(proto_ipm$usr_funs[[1]][1])) {

      usr_funs  <- proto_ipm$usr_funs[[1]]

    }

  } else if(!rlang::is_empty(usr_funs)) {

    proto_ipm <- .append_usr_funs_to_proto(proto_ipm, usr_funs)

  }

  proto_list <- .initialize_kernels(proto_ipm, iterate)

  others <- proto_list$others
  k_row  <- proto_list$k_row   # k_row is either an NA or a proto_ipm object

  # Initialize the master_environment so these values can all be found at
  # evaluation time

  if(is.null(domain_list)) {
    master_env       <- .make_master_env(others$domain, usr_funs)
  } else {
    master_env       <- .make_master_env(domain_list, usr_funs)
    proto_ipm$domain <- I(list(domain_list))
  }

  temp           <- .prep_di_output(others, k_row,
                                    proto_ipm, iterations,
                                    normalize_pop_size)

  master_env     <- .add_pop_state_to_master_env(temp$pop_state,
                                                 master_env)

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

    # .iterate_model is internal generic - see internal-model_iteration.R

    pop_state <- .iterate_model(proto_ipm,
                                iterators,
                                sub_kern_list,
                                iterations,
                                temp$pop_state,
                                master_env,
                                k_row,
                                normalize_pop_size)
  } else {

    pop_state <- NA_real_

  }


  if(return_all) {
    out_ret <- all_sub_kerns$env_list
  } else {
    out_ret <- NA_character_
  }

  out_seq <- NA_integer_

  iterators     <- set_ipmr_classes(iterators)
  sub_kern_list <- set_ipmr_classes(sub_kern_list)

  out <- list(iterators   = iterators,
              sub_kernels = sub_kern_list,
              env_list    = out_ret,
              env_seq     = out_seq,
              pop_state   = pop_state,
              proto_ipm   = proto_ipm
  )

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)
  class(out)                         <- c('simple_di_det_ipm', 'list')

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
                                          kernel_seq  = NULL,
                                          normalize_pop_size = TRUE,
                                          report_progress = FALSE) {


  if(iterate && all(is.null(kernel_seq) | is.na(kernel_seq))) {

    message("'kernel_seq' not defined. Will generate one internally")

    kernel_seq <- "internal"
  }

  # Work out whether to append usr_funs to proto or to restore them from prior
  # implemenation. Logic is documented in make_ipm.simple_di_det()

  if(isTRUE(attr(proto_ipm, 'implemented')) && rlang::is_empty(usr_funs)) {

    if(!is.na(proto_ipm$usr_funs[[1]][1])) {

      usr_funs  <- proto_ipm$usr_funs[[1]]

    }

  } else if(!rlang::is_empty(usr_funs)) {

    proto_ipm <- .append_usr_funs_to_proto(proto_ipm, usr_funs)

  }

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

  temp       <- .prep_di_output(others, k_row,
                                proto_ipm, iterations,
                                normalize_pop_size)

  master_env <- .add_pop_state_to_master_env(temp$pop_state, master_env)

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

    # .iterate_model is internal generic - see internal-model_iteration.R

    pop_state      <- .iterate_model(proto_ipm,
                                     iterators,
                                     sub_kern_list,
                                     iterations,
                                     kern_seq,
                                     temp$pop_state,
                                     master_env,
                                     k_row,
                                     normalize_pop_size,
                                     report_progress)

    # In order to operate properly, we had to insert an extra entry in
    # lambda to make sure it had the same length as pop_state (it'll always be
    # one fewer in reality though). The first entry is ALWAYS NA, but that doesn't
    # seem user friendly (ipmr-internal friend either, tbh)

    pop_state$lambda <- pop_state$lambda[-1]

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
              proto_ipm   = proto_ipm

  )

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)
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
                                           iterations  = 50,
                                           kernel_seq  = NULL,
                                           normalize_pop_size = TRUE,
                                           report_progress = FALSE) {

  # Work out whether to append usr_funs to proto or to restore them from prior
  # implemenation. Logic is documented in make_ipm.simple_di_det()

  if(isTRUE(attr(proto_ipm, 'implemented')) && rlang::is_empty(usr_funs)) {

    if(!is.na(proto_ipm$usr_funs[[1]][1])) {

      usr_funs  <- proto_ipm$usr_funs[[1]]

    }

  } else if(!rlang::is_empty(usr_funs)) {

    proto_ipm <- .append_usr_funs_to_proto(proto_ipm, usr_funs)

  }

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

  master_env <- .bind_all_constants(env_state   = others$env_state[[1]]$constants,
                                    env_to_bind = master_env)

  temp        <- .prep_di_output(others, k_row, proto_ipm, iterations,
                                 normalize_pop_size)

  # initialize the pop_state vectors in master_env so they can be found
  # at evaluation time

  master_env <- .add_pop_state_to_master_env(temp$pop_state,
                                             master_env)

  # list to hold the possibly returned evaluation environments

  env_list <- list(master_env = master_env)

  kern_seq <- .make_kern_seq(others,
                             c(others$kernel_id,
                               k_row$kernel_id),
                             iterations,
                             kernel_seq)

  for(i in seq_len(iterations)) {

    # Lazy variant makes sure that whatever functions that generate parameter
    # values stochastically are only evaluated 1 time per iteration! This is so
    # multiple parameters meant to come from a joint distribution really come from
    # the joint distribution!

    .bind_iter_var(master_env, i)
    .stoch_progress_message(report_progress, iterations, i)

    sys         <- .make_sub_kernel_simple_lazy(others,
                                                master_env,
                                                return_envs = return_all,
                                                dd = 'n')

    sub_kernels <- sys$ipm_system$sub_kernels

    sys_i       <- .make_k_param_samp(k_row,
                                      sub_kernels,
                                      master_env)

    # .iterate_model is internal generic - see internal-model_iteration.R

    pop_state   <- .iterate_model(proto_ipm,
                                  iterators = sys_i$iterators,
                                  sub_kernels,
                                  current_iteration = i,
                                  kern_seq = kern_seq,
                                  temp$pop_state,
                                  master_env,
                                  k_row,
                                  normalize_pop_size)


    if(return_all) {

      env_list <- c(sys$data_envs, env_list)

    } else {

      env_list <- NA_character_

    }

    temp <- .update_param_simple_output(
      sub_kernels,
      sys_i,
      pop_state,
      env_list,
      master_env,
      temp,
      iterations,
      i
    )

  }

  temp$pop_state$lambda <- temp$pop_state$lambda[-1]

  out <- list(
    iterators   = temp$iterators,
    sub_kernels = temp$sub_kernels,
    env_list    = temp$sub_kernel_envs,
    env_seq     = temp$env_seq,
    pop_state   = temp$pop_state,
    proto_ipm   = proto_ipm
  )

  out$iterators   <- set_ipmr_classes(out$iterators)
  out$sub_kernels <- set_ipmr_classes(out$sub_kernels)

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)
  class(out) <- c('simple_di_stoch_param_ipm', 'list')

  if(return_all) out$data_envs <- env_list

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
                                    iterations  = 50,
                                    normalize_pop_size = TRUE) {

  # Work out whether to append usr_funs to proto or to restore them from prior
  # implemenation. Logic is documented in make_ipm.simple_di_det()

  if(isTRUE(attr(proto_ipm, 'implemented')) && rlang::is_empty(usr_funs)) {

    if(!is.na(proto_ipm$usr_funs[[1]][1])) {

      usr_funs  <- proto_ipm$usr_funs[[1]]

    }

  } else if(!rlang::is_empty(usr_funs)) {

    proto_ipm <- .append_usr_funs_to_proto(proto_ipm, usr_funs)

  }

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

  master_env <- .bind_all_constants(env_state   = others$env_state[[1]],
                                    env_to_bind = master_env)


  temp       <- .prep_di_output(others, k_row, proto_ipm, iterations,
                                normalize_pop_size)

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

    pop_state <- .iterate_model(proto_ipm,
                                k_row,
                                sub_kern_list,
                                iterations,
                                temp$pop_state,
                                master_env,
                                normalize_pop_size)


  }


  # Final bits of housekeeping post-iteration. .prep_other_output deals specifically
  # with items that depend on `return_all` and `iterate` switches, but have
  # length > 1 (ifelse() output always equals length(input))

  sub_kern_list  <- set_ipmr_classes(sub_kern_list)


  temp_other_out <- .prep_other_output(env_list,
                                       kern_seq = NA_integer_,
                                       pop_state,
                                       return_all,
                                       iterate)

  env_ret     <- temp_other_out$env_ret
  env_seq_ret <- temp_other_out$env_seq_ret
  pop_ret     <- temp_other_out$pop_ret

  out <- list(
    iterators   = NA_real_,
    sub_kernels = sub_kern_list,
    env_list    = env_ret,
    env_seq     = env_seq_ret,
    pop_state   = pop_ret,
    proto_ipm   = proto_ipm
  )

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)
  class(out) <- c('general_di_det_ipm', 'list')

  return(out)

}

#' @rdname make_ipm
#'
#' @export

make_ipm.general_di_stoch_kern <- function(proto_ipm,
                                           return_all  = FALSE,
                                           usr_funs    = list(),
                                           ...,
                                           domain_list = NULL,
                                           iterate     = TRUE,
                                           iterations  = 50,
                                           kernel_seq  = NULL,
                                           normalize_pop_size = TRUE,
                                           report_progress = FALSE) {


  if(iterate && all(is.null(kernel_seq) | is.na(kernel_seq))) {

    message("'kernel_seq' not defined. Will generate one internally")

    kernel_seq <- "internal"
  }

  # Work out whether to append usr_funs to proto or to restore them from prior
  # implemenation. Logic is documented in make_ipm.simple_di_det()

  if(isTRUE(attr(proto_ipm, 'implemented')) && rlang::is_empty(usr_funs)) {

    if(!is.na(proto_ipm$usr_funs[[1]][1])) {

      usr_funs  <- proto_ipm$usr_funs[[1]]

    }

  } else if(!rlang::is_empty(usr_funs)) {

    proto_ipm <- .append_usr_funs_to_proto(proto_ipm, usr_funs)

  }

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

  master_env <- .bind_all_constants(env_state   = others$env_state[[1]],
                                    env_to_bind = master_env)


  temp       <- .prep_di_output(others, k_row, proto_ipm, iterations,
                                normalize_pop_size)

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

    kern_seq  <- .make_kern_seq(proto_ipm,
                                sub_kern_list,
                                iterations,
                                kernel_seq)

    pop_state <- .iterate_model(proto_ipm,
                                k_row,
                                sub_kern_list,
                                iterations,
                                kern_seq,
                                temp$pop_state,
                                master_env,
                                normalize_pop_size,
                                report_progress)

    pop_state$lambda <- pop_state$lambda[-1]
  }

  # Final bits of housekeeping post-iteration. .prep_other_output deals specifically
  # with items that depend on `return_all` and `iterate` switches, but have
  # length > 1 (ifelse() output always equals length(input))

  temp_other_out <- .prep_other_output(env_list,
                                       kern_seq,
                                       pop_state,
                                       return_all,
                                       iterate)

  env_ret     <- temp_other_out$env_ret
  env_seq_ret <- temp_other_out$env_seq_ret
  pop_ret     <- temp_other_out$pop_ret

  sub_kern_list  <- set_ipmr_classes(sub_kern_list)

  out <- list(
    iterators   = NA_real_,
    sub_kernels = sub_kern_list,
    env_list    = env_ret,
    env_seq     = env_seq_ret,
    pop_state   = pop_ret,
    proto_ipm   = proto_ipm
  )

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)
  class(out) <- c('general_di_stoch_kern_ipm', 'list')

  return(out)

}


#' @rdname make_ipm
#'
#' @export

make_ipm.general_di_stoch_param <- function(proto_ipm,
                                            return_all = FALSE,
                                            usr_funs = list(),
                                            ...,
                                            domain_list = NULL,
                                            iterate     = TRUE,
                                            iterations  = 50,
                                            kernel_seq  = NULL,
                                            normalize_pop_size = TRUE,
                                            report_progress = FALSE) {

  # Work out whether to append usr_funs to proto or to restore them from prior
  # implemenation. Logic is documented in make_ipm.simple_di_det()

  if(isTRUE(attr(proto_ipm, 'implemented')) && rlang::is_empty(usr_funs)) {

    if(!is.na(proto_ipm$usr_funs[[1]][1])) {

      usr_funs  <- proto_ipm$usr_funs[[1]]

    }

  } else if(!rlang::is_empty(usr_funs)) {

    proto_ipm <- .append_usr_funs_to_proto(proto_ipm, usr_funs)

  }

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

  master_env <- .bind_all_constants(env_state   = others$env_state[[1]]$constants,
                                    env_to_bind = master_env)

  temp       <- .prep_di_output(others, k_row, proto_ipm, iterations,
                                normalize_pop_size)

  # initialize the pop_state vectors in master_env so they can be found
  # at evaluation time

  if(all(is.na(proto_ipm$pop_state[[1]]))) {

    stop("All general_* IPMs must have a 'pop_state' defined.",
         "See ?define_pop_state() for more details." )

  } else {

    master_env  <- .add_pop_state_to_master_env(temp$pop_state, master_env)

  }

  if(!iterate || iterations < 1) {
    stop("All 'general_*_stoch_param' models must be iterated at least once!",
         call. = FALSE)
  }

  kern_seq <- .make_kern_seq(others,
                             c(others$kernel_id,
                               k_row$kernel_id),
                             iterations,
                             kernel_seq)

  for(i in seq_len(iterations)) {

    # add the variable to let the user access the current iteration.

    .bind_iter_var(master_env, i)

    .stoch_progress_message(report_progress, iterations, i)

    # Lazy variant makes sure that whatever functions that generate parameter
    # values stochastically are only evaluated 1 time per iteration! This is so
    # multiple parameters meant to come from a joint distribution really come from
    # the joint distribution!

    sys         <- .make_sub_kernel_general_lazy(others,
                                                 master_env,
                                                 return_envs = return_all)

    sub_kernels <- sys$ipm_system$sub_kernels

    # Generate the pop_state for a single iteration! This is critical to ensuring
    # env_state_funs are only evaluated once per iteration. kern_seq = NULL
    # because the environmental parameters are generated on the fly by the
    # user defined function

    pop_state <- .iterate_model(proto_ipm,
                                k_row,
                                sub_kernels,
                                current_iteration = i,
                                kern_seq = kern_seq,
                                temp$pop_state,
                                master_env,
                                normalize_pop_size)

    if(return_all) {

      env_ret <- sys$ipm_system$env_list

    } else {

      env_ret <- NA_character_

    }

    # Variant that doesn't require an "iterator" slot

    temp         <- .update_param_general_output(sub_kernels,
                                                 pop_state,
                                                 env_ret,
                                                 master_env,
                                                 temp,
                                                 iterations,
                                                 i)

  }

  temp$env_seq <- data.frame(temp$env_seq, stringsAsFactors = FALSE)

  if(!is.null(kern_seq)) {

    temp$env_seq <- cbind(temp$env_seq,
                          kernel_seq = kern_seq)

  }

  temp$pop_state$lambda <- temp$pop_state$lambda[-1]

  out <- list(
    iterators   = NA_real_,
    sub_kernels = temp$sub_kernels,
    env_list    = temp$sub_kernel_envs,
    env_seq     = temp$env_seq,
    pop_state   = temp$pop_state,
    proto_ipm   = proto_ipm
  )

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)

  out$sub_kernels                    <- set_ipmr_classes(out$sub_kernels)

  class(out) <- c('general_di_stoch_param_ipm', 'list')

  return(out)

}

# Density dependent methods----------

make_ipm.simple_dd_det <- function(proto_ipm,
                                   return_all = FALSE,
                                   usr_funs = list(),
                                   ...,
                                   domain_list = NULL,
                                   iterate     = TRUE,
                                   iterations  = 50,
                                   normalize_pop_size = FALSE) {

  # Work out whether to append usr_funs to proto or to restore them from prior
  # implemenation. Logic is documented in make_ipm.simple_di_det()

  if(isTRUE(attr(proto_ipm, 'implemented')) && rlang::is_empty(usr_funs)) {

    if(!is.na(proto_ipm$usr_funs[[1]][1])) {

      usr_funs  <- proto_ipm$usr_funs[[1]]

    }

  } else if(!rlang::is_empty(usr_funs)) {

    proto_ipm <- .append_usr_funs_to_proto(proto_ipm, usr_funs)

  }

  proto_list <- .initialize_kernels_dd(proto_ipm, iterate)

  others <- proto_list$others
  k_row  <- proto_list$k_row

  # Initialize the master_environment so these values can all be found at
  # evaluation time

  if(is.null(domain_list)) {
    master_env <- .make_master_env(others$domain, usr_funs)
  } else {
    master_env <- .make_master_env(domain_list, usr_funs)
  }

  # prepare dd output. I'm not sure how this would be different
  # than di output, so I'm going to keep using the di function for now.
  # will change to something else

  out        <- .prep_di_output(others, k_row, proto_ipm, iterations,
                                normalize_pop_size)
  pop_state  <- out$pop_state

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

    env_list <- list(master_env = master_env)

    sys         <- .make_sub_kernel_simple(others,
                                           env_list,
                                           return_envs = return_all)

    sub_kernels <- sys$sub_kernels

    iterator       <- .make_k_simple(k_row,
                                     proto_ipm,
                                     sub_kernels,
                                     master_env)

    if(return_all) {

      env_ret <- sys$data_envs

    } else {

      env_ret <- NA_character_

    }

    # Use iterate_kerns_simple with iterations = 1. We have to rebuild kernels
    # each time so don't want to pass more than that

    pop_state <- .iterate_model(proto_ipm,
                                iterator,
                                sub_kernels,
                                iterations = 1,
                                current_iteration = i,
                                kern_seq = NULL,
                                pop_state,
                                master_env,
                                k_row,
                                normalize_pop_size)

    names(iterator)    <- paste(names(iterator), i, sep = "_")
    names(sub_kernels) <- paste(names(sub_kernels), i, sep = "_")

    out$iterators   <- purrr::splice(out$iterators, iterator)
    out$sub_kernels <- purrr::splice(out$sub_kernels, sub_kernels)

  }

  out$iterators   <- set_ipmr_classes(out$iterators)
  out$sub_kernels <- set_ipmr_classes(out$sub_kernels)

  out$pop_state   <- pop_state

  attr(out$proto_ipm, 'implemented') <- TRUE
  class(out) <- c('simple_dd_stoch_param_ipm', 'list')

  if(return_all) out$data_envs <- purrr::splice(env_list, out$data_envs)

  return(out)


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
