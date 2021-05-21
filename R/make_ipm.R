#' @title Methods to implement an IPM
#' @rdname make_ipm
#'
#' @description The \code{make_ipm.*} methods convert a \code{proto_ipm} into a
#' set of discretized kernels and population vectors. Methods have different
#' requirements, so be sure to read the parameter documentation. \code{
#' vignette('ipmr-introduction', 'ipmr')}  a more complete introduction.
#'
#' @param proto_ipm A proto_ipm. This should be the
#' output of \code{define_kernel}, or the \code{define_*} functions.
#'
#' @param ... Other arguments passed to methods.
#'
#' @param return_main_env A logical indicating whether to return the main environment
#' for the model. This environment contains the integration mesh, weights, and
#' other potentially useful variables for subsequent analyses. Default is
#' \code{TRUE}.
#'
#' @param return_all_envs A logical indicating whether to return the environments that
#' the kernel expressions are evaluated in. These may be useful for some analyses,
#' such as regression-level sensitivity/elasticity analyses, but can also rapidly
#' increase memory consumption for models with many kernels (e.g. ones with
#' parameter set indices that have many levels, or any \code{*_stoch_param} model).
#' Default is \code{FALSE}.
#'
#' @param domain_list An optional list of new domain information to implement
#' the IPM with.
#'
#' @param usr_funs An optional list of user-specified functions that are passed
#' on to the  model building process. This can help make vital rate expressions
#' more concise and expressive. Names in this list should exactly match the names
#' of the function calls in the \code{...} or \code{formula}.
#'
#' @param iterate A logical indicating whether or not iterate the model before exiting
#' or just return the sub-kernels. Only applies to density-independent, deterministic
#' models and density-independent, stochastic kernel re-sampled models.
#'
#' @param iterations If \code{iterate} is \code{TRUE}, then the number of iterations
#' to run the model for.
#'
#' @param normalize_pop_size A logical indicating whether to re-scale the population
#' vector to sum to 1 before each iteration. Default is \code{TRUE} for
#' \code{*_di_*} methods and \code{FALSE} for \code{*_dd_*} methods.
#'
#' @param kernel_seq For \code{*_stoch_kern} methods, the sequence of kernels
#' to use during the simulation process. It should have the same number of entries
#' as the number of \code{iterations}.
#' This should be a vector containing values of the parameter set indices specified
#' in \code{par_set_indices}, or empty. Support for Markov chains will eventually
#' get implemented. If it is empty, \code{make_ipm} will try to generate a
#' sequence internally using a random selection of the \code{par_set_indices}
#' defined in \code{define_kernel}.
#'
#' @param report_progress A logical indicating whether or not to periodically
#' report progress for a stochastic simulation. Does not apply to deterministic
#' methods. Default is \code{FALSE}.
#'
#' @param iteration_direction Either \code{"right"} (default) or \code{"left"}.
#' This controls the direction of projection. Right iteration will generate
#' the right eigenvector (if it exists), while left iteration generates
#' the left eigenvector. These correspond to the stable trait distributions, and
#' reproductive values, respectively. This parameter is mostly used internally
#' by other functions. Use with care.
#'
#'
#' @return
#'  The \code{make_ipm.*} methods will always return a list of length 5
#' containing the following components:
#'
#' \itemize{
#'   \item{\strong{sub_kernels}}{: a list of arrays specified in \code{define_kernel}.}
#'   \item{\strong{env_list}}{: a list containing the evaluation environments of
#'                            kernel. This will contain the \code{main_env} object
#'                            if \code{return_main_env = TRUE}. It will also contain
#'                            the sub-kernels evaluation environments if
#'                            \code{return_all_envs = TRUE}. }
#'   \item{\strong{env_seq}}{: a character vector with length \code{iterations} of
#'                              kernel indices indicating the order
#'                              in which kernels are to be/were resampled OR
#'                              a matrix with as many columns as stochastic parameters
#'                              and \code{n_iterations} rows.}
#'   \item{\strong{pop_state}}{: population vectors
#'                              stored as a list of arrays. The first dimension
#'                              of each array corresponds to the state variable distribution,
#'                              and the second dimension corresponds to time
#'                              steps.}
#'   \item{\strong{proto_ipm}}{: the \code{proto_ipm} object used to implement
#'                              the model.}
#' }
#'
#' In addition to the list class, each object will have a class comprised of the
#' arguments from  \code{init_ipm} plus \code{'ipm'} pasted together with
#' underscores. This is to facilitate \code{print}, \code{plot}, and
#' \code{lambda} methods. For example, a \code{init_ipm("general", "di", "det")}
#' model will have the class \code{'general_di_det_ipm'} once it has been
#' implemented using \code{make_ipm}.
#'
#' @export

make_ipm <- function(proto_ipm,
                     return_main_env = TRUE,
                     return_all_envs = FALSE,
                     usr_funs = list(),
                     ...) {

  UseMethod('make_ipm')

}


#' @rdname make_ipm
#'
#' @export

make_ipm.simple_di_det <- function(proto_ipm,
                                   return_main_env = TRUE,
                                   return_all_envs = FALSE,
                                   usr_funs = list(),
                                   ...,
                                   domain_list = NULL,
                                   iterate = TRUE,
                                   iterations = 50,
                                   normalize_pop_size = TRUE,
                                   iteration_direction = "right"
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

  proto_list <- .initialize_kernels(proto_ipm, iterate, iteration_direction)

  others <- proto_list$others
  k_row  <- proto_list$k_row   # k_row is either an NA or a proto_ipm object

  # Initialize the main_environment so these values can all be found at
  # evaluation time

  if(is.null(domain_list)) {
    main_env       <- .make_main_env(others$domain, usr_funs,
                                     age_size = .uses_age(others))
  } else {
    main_env       <- .make_main_env(domain_list, usr_funs,
                                     age_size = .uses_age(others))
    proto_ipm$domain <- I(list(domain_list))
  }

  temp           <- .prep_di_output(others, k_row,
                                    proto_ipm, iterations,
                                    normalize_pop_size)

  main_env     <- .add_pop_state_to_main_env(temp$pop_state,
                                                 main_env)

  # construct the kernels from their function definitions

  env_list      <- list(main_env = main_env)

  all_sub_kerns <- .make_sub_kernel_simple(others,
                                           env_list,
                                           return_envs = return_all_envs)

  sub_kern_list <- all_sub_kerns$sub_kernels


  if(iterate) {

    # .iterate_model is internal generic - see internal-model_iteration.R

    pop_state <- .iterate_model(proto_ipm,
                                sub_kern_list,
                                iterations,
                                temp$pop_state,
                                main_env,
                                k_row,
                                normalize_pop_size)

    # Convert pop_state names back to the user-supplied ones

    names(pop_state) <- gsub("^pop_state_", "n_", names(pop_state))

  } else {

    pop_state <- NA_real_

  }


  if(return_all_envs) {

    out_ret <- all_sub_kerns$env_list

  } else if(return_main_env) {

    out_ret <- list(main_env = main_env)

  } else {

    out_ret <- NA_character_

  }

  out_seq <- NA_integer_

  sub_kern_list <- set_ipmr_classes(sub_kern_list)

  out <- list(sub_kernels = sub_kern_list,
              env_list    = out_ret,
              env_seq     = out_seq,
              pop_state   = pop_state,
              proto_ipm   = proto_ipm
  )

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)

  class(out)                         <- c('simple_di_det_ipm',
                                          'list')

  return(out)

}

#' @rdname make_ipm
#'
#' @export
make_ipm.simple_di_stoch_kern <- function(proto_ipm,
                                          return_main_env     = TRUE,
                                          return_all_envs     = FALSE,
                                          usr_funs            = list(),
                                          ...,
                                          domain_list         = NULL,
                                          iterate             = TRUE,
                                          iterations          = 50,
                                          kernel_seq          = NULL,
                                          normalize_pop_size  = TRUE,
                                          report_progress     = FALSE,
                                          iteration_direction = "right") {


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

  proto_list <- .initialize_kernels(proto_ipm, iterate, iteration_direction)

  others <- proto_list$others
  k_row  <- proto_list$k_row

  # Initialize the main_environment so these values can all be found at
  # evaluation time

  if(is.null(domain_list)){
    main_env <- .make_main_env(others$domain, usr_funs,
                               age_size = .uses_age(others))
  } else {
    main_env <- .make_main_env(domain_list, usr_funs,
                               age_size = .uses_age(others))
  }

  temp       <- .prep_di_output(others, k_row,
                                proto_ipm, iterations,
                                normalize_pop_size)

  main_env <- .add_pop_state_to_main_env(temp$pop_state, main_env)

  # construct the kernels from their function defintions

  env_list      <- list(main_env = main_env)

  all_sub_kerns <- .make_sub_kernel_simple(others,
                                           env_list,
                                           return_envs = return_all_envs)

  sub_kern_list <- all_sub_kerns$sub_kernels
  env_list      <- all_sub_kerns$env_list

  kern_seq      <-  .make_kern_seq(proto_ipm,
                                   iterations,
                                   kernel_seq)

  if(iterate) {

    # .iterate_model is internal generic - see internal-model_iteration.R

    pop_state      <- .iterate_model(proto_ipm,
                                     sub_kern_list,
                                     iterations,
                                     kern_seq,
                                     temp$pop_state,
                                     main_env,
                                     k_row,
                                     normalize_pop_size,
                                     report_progress)

    # In order to operate properly, we had to insert an extra entry in
    # lambda to make sure it had the same length as pop_state (it'll always be
    # one fewer in reality though). The first entry is ALWAYS NA, but that doesn't
    # seem user friendly (ipmr-internal friend either, tbh). After that,
    # we need to conver the names of the trait distributions back to the
    # user-supplied ones

    pop_state$lambda <- pop_state$lambda[-1]

    names(pop_state) <- gsub("^pop_state_", "n_", names(pop_state))

  } else {

    pop_state <- NA_real_

  }

  if(return_all_envs) {

    out_ret <- env_list

  } else if(return_main_env) {

    out_ret <- list(main_env = main_env)

  } else {

    out_ret <- NA_character_

  }

  sub_kern_list <- set_ipmr_classes(sub_kern_list)

  out <- list(sub_kernels = sub_kern_list,
              env_list    = out_ret,
              env_seq     = kern_seq,
              pop_state   = pop_state,
              proto_ipm   = proto_ipm

  )

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)

  class(out) <- c('simple_di_stoch_kern_ipm',
                  'list')

  return(out)
}

#' @rdname make_ipm
#'
#' @export

make_ipm.simple_di_stoch_param <- function(proto_ipm,
                                           return_main_env     = TRUE,
                                           return_all_envs     = FALSE,
                                           usr_funs            = list(),
                                           ...,
                                           domain_list         = NULL,
                                           iterate             = TRUE,
                                           iterations          = 50,
                                           kernel_seq          = NULL,
                                           normalize_pop_size  = TRUE,
                                           report_progress     = FALSE,
                                           iteration_direction = "right") {

  # Work out whether to append usr_funs to proto or to restore them from prior
  # implemenation. Logic is documented in make_ipm.simple_di_det()

  if(isTRUE(attr(proto_ipm, 'implemented')) && rlang::is_empty(usr_funs)) {

    if(!is.na(proto_ipm$usr_funs[[1]][1])) {

      usr_funs  <- proto_ipm$usr_funs[[1]]

    }

  } else if(!rlang::is_empty(usr_funs)) {

    proto_ipm <- .append_usr_funs_to_proto(proto_ipm, usr_funs)

  }

  proto_list <- .initialize_kernels(proto_ipm, iterate, iteration_direction)

  others <- proto_list$others
  k_row  <- proto_list$k_row

  # Initialize the main_environment so these values can all be found at
  # evaluation time

  if(is.null(domain_list)) {
    main_env <- .make_main_env(others$domain, usr_funs,
                               age_size = .uses_age(others))
  } else {
    main_env <- .make_main_env(domain_list, usr_funs,
                               age_size = .uses_age(others))
  }

  # Bind env_exprs, constants, and pop_vectors to main_env so that
  # we can always find them and avoid that miserable repitition

  main_env <- .bind_all_constants(env_state   = others$env_state[[1]]$constants,
                                    env_to_bind = main_env)

  temp        <- .prep_di_output(others, k_row, proto_ipm, iterations,
                                 normalize_pop_size)

  # initialize the pop_state vectors in main_env so they can be found
  # at evaluation time

  main_env <- .add_pop_state_to_main_env(temp$pop_state,
                                             main_env)

  # list to hold the possibly returned evaluation environments

  env_list <- list(main_env = main_env)

  kern_seq <- .make_kern_seq(others,
                             iterations,
                             kernel_seq)

  for(i in seq_len(iterations)) {

    # Lazy variant makes sure that whatever functions that generate parameter
    # values stochastically are only evaluated 1 time per iteration! This is so
    # multiple parameters meant to come from a joint distribution really come from
    # the joint distribution!

    .bind_iter_var(main_env, i)
    .stoch_progress_message(report_progress, iterations, i)

    sys         <- .make_sub_kernel_simple_lazy(others,
                                                main_env,
                                                return_envs = return_all_envs,
                                                dd = 'n')

    sub_kernels <- sys$ipm_system$sub_kernels

    # .iterate_model is internal generic - see internal-model_iteration.R

    pop_state   <- .iterate_model(proto_ipm,
                                  sub_kernels,
                                  current_iteration = i,
                                  kern_seq = kern_seq,
                                  temp$pop_state,
                                  main_env,
                                  k_row,
                                  normalize_pop_size)


    if(return_all_envs) {

      sys$env_list$main_env <- NULL

      names(sys$env_list) <- paste(names(sys$env_list), "it", i, sep = "_")

      env_ret             <- c(env_ret, sys$env_list)

    } else if(!return_main_env) {

      env_list <- NA_character_

    }

    names(sub_kernels) <- paste(names(sub_kernels), "it", i, sep = "_")

    temp <- .update_param_output(
      sub_kernels,
      pop_state,
      env_list,
      main_env,
      temp,
      iterations,
      i
    )

  }

  temp$env_seq <- data.frame(temp$env_seq, stringsAsFactors = FALSE)

  if(!is.null(kern_seq)) {

    temp$env_seq <- cbind(temp$env_seq,
                          kernel_seq = kern_seq)

  }

  temp$pop_state$lambda <- temp$pop_state$lambda[-1]
  names(temp$pop_state) <- gsub("^pop_state_", "n_", names(temp$pop_state))

  out <- list(
    sub_kernels = temp$sub_kernels,
    env_list    = env_list,
    env_seq     = temp$env_seq,
    pop_state   = temp$pop_state,
    proto_ipm   = proto_ipm
  )

  out$sub_kernels <- set_ipmr_classes(out$sub_kernels)

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)

  class(out) <- c('simple_di_stoch_param_ipm',
                  'list')

  return(out)


}


#' @rdname make_ipm
#'
#' @export

make_ipm.general_di_det <- function(proto_ipm,
                                    return_main_env     = TRUE,
                                    return_all_envs     = FALSE,
                                    usr_funs            = list(),
                                    ...,
                                    domain_list         = NULL,
                                    iterate             = TRUE,
                                    iterations          = 50,
                                    normalize_pop_size  = TRUE,
                                    iteration_direction = "right") {

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

  proto_list <- .initialize_kernels(proto_ipm, iterate, iteration_direction)

  others <- proto_list$others
  k_row  <- proto_list$k_row

  if(is.null(domain_list)) {
    main_env <- .make_main_env(others$domain,
                               usr_funs,
                               age_size = .uses_age(others))
  } else {
    main_env <- .make_main_env(domain_list, usr_funs,
                               age_size = .uses_age(others))
  }

  # Bind env_exprs, constants, and pop_vectors to main_env so that
  # we can always find them and avoid that miserable repitition

  main_env <- .bind_all_constants(env_state   = others$env_state[[1]],
                                    env_to_bind = main_env)


  temp     <- .prep_di_output(others, k_row, proto_ipm, iterations,
                              normalize_pop_size)

  # Thus far, I think general_* methods will have to use iteration for lambdas
  # as I'm not sure I want to work out the correct cbind(rbind(...)) rules for
  # creating a mega-matrix. So throw an error if pop_state isn't defined
  if(all(is.na(proto_ipm$pop_state[[1]]))) {

    stop("All general_* IPMs must have a 'pop_state' defined.",
         "See ?define_pop_state() for more details." )

  } else {

    main_env  <- .add_pop_state_to_main_env(temp$pop_state, main_env)

  }

  env_list      <- list(main_env = main_env)

  all_sub_kerns <- .make_sub_kernel_general(others,
                                            env_list,
                                            return_envs = return_all_envs)

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
                                main_env,
                                normalize_pop_size)

    names(pop_state) <- gsub("^pop_state_", "n_", names(pop_state))

  }


  # Final bits of housekeeping post-iteration. .prep_other_output deals specifically
  # with items that depend on `return_all_envs` and `iterate` switches, but have
  # length > 1 (ifelse() output always equals length(input))

  sub_kern_list  <- set_ipmr_classes(sub_kern_list)


  temp_other_out <- .prep_other_output(env_list,
                                       kern_seq = NA_integer_,
                                       pop_state,
                                       return_main_env,
                                       return_all_envs,
                                       iterate)

  env_ret     <- temp_other_out$env_ret
  env_seq_ret <- temp_other_out$env_seq_ret
  pop_ret     <- temp_other_out$pop_ret

  out <- list(
    sub_kernels = sub_kern_list,
    env_list    = env_ret,
    env_seq     = env_seq_ret,
    pop_state   = pop_ret,
    proto_ipm   = proto_ipm
  )

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)
  if(inherits(proto_ipm, "age_x_size")) {
    a_s_class <- "age_x_size_ipm"
  } else {
    a_s_class <- NULL
  }
  class(out) <- c('general_di_det_ipm',
                  a_s_class,
                  'list')

  return(out)

}

#' @rdname make_ipm
#'
#' @export

make_ipm.general_di_stoch_kern <- function(proto_ipm,
                                           return_main_env     = TRUE,
                                           return_all_envs     = FALSE,
                                           usr_funs            = list(),
                                           ...,
                                           domain_list         = NULL,
                                           iterate             = TRUE,
                                           iterations          = 50,
                                           kernel_seq          = NULL,
                                           normalize_pop_size  = TRUE,
                                           report_progress     = FALSE,
                                           iteration_direction = "right") {


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

  proto_list <- .initialize_kernels(proto_ipm, iterate, iteration_direction)

  others <- proto_list$others
  k_row  <- proto_list$k_row

  if(is.null(domain_list)) {

    main_env <- .make_main_env(others$domain, usr_funs,
                               age_size = .uses_age(others))

  } else {

    main_env <- .make_main_env(domain_list, usr_funs,
                               age_size = .uses_age(others))

  }

  # Bind env_exprs, constants, and pop_vectors to main_env so that
  # we can always find them and avoid that miserable repitition

  main_env <- .bind_all_constants(env_state   = others$env_state[[1]],
                                    env_to_bind = main_env)


  temp       <- .prep_di_output(others, k_row, proto_ipm, iterations,
                                normalize_pop_size)

  # Thus far, I think general_* methods will have to use iteration for lambdas
  # as I'm not sure I want to work out the correct cbind(rbind(...)) rules for
  # creating a mega-matrix. So throw an error if pop_state isn't defined
  if(all(is.na(proto_ipm$pop_state[[1]]))) {

    stop("All general_* IPMs must have a 'pop_state' defined.",
         "See ?define_pop_state() for more details." )

  } else {

    main_env  <- .add_pop_state_to_main_env(temp$pop_state, main_env)

  }

  env_list      <- list(main_env = main_env)

  all_sub_kerns <- .make_sub_kernel_general(others,
                                            env_list,
                                            return_envs = return_all_envs)

  sub_kern_list <- all_sub_kerns$sub_kernels
  env_list      <- all_sub_kerns$env_list

  if(iterate) {

    kern_seq  <- .make_kern_seq(proto_ipm,
                                iterations,
                                kernel_seq)

    pop_state <- .iterate_model(proto_ipm,
                                k_row,
                                sub_kern_list,
                                iterations,
                                kern_seq,
                                temp$pop_state,
                                main_env,
                                normalize_pop_size,
                                report_progress)

    pop_state$lambda <- pop_state$lambda[-1]
    names(pop_state) <- gsub("^pop_state_", "n_", names(pop_state))

  }

  # Final bits of housekeeping post-iteration. .prep_other_output deals specifically
  # with items that depend on `return_all_envs` and `iterate` switches, but have
  # length > 1 (ifelse() output always equals length(input))

  temp_other_out <- .prep_other_output(env_list,
                                       kern_seq,
                                       pop_state,
                                       return_main_env,
                                       return_all_envs,
                                       iterate)

  env_ret     <- temp_other_out$env_ret
  env_seq_ret <- temp_other_out$env_seq_ret
  pop_ret     <- temp_other_out$pop_ret

  sub_kern_list  <- set_ipmr_classes(sub_kern_list)

  out <- list(
    sub_kernels = sub_kern_list,
    env_list    = env_ret,
    env_seq     = env_seq_ret,
    pop_state   = pop_ret,
    proto_ipm   = proto_ipm
  )

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)
  if(inherits(proto_ipm, "age_x_size")) {
    a_s_class <- "age_x_size_ipm"
  } else {
    a_s_class <- NULL
  }
  class(out) <- c('general_di_stoch_kern_ipm',
                  a_s_class,
                  'list')

  return(out)

}


#' @rdname make_ipm
#'
#' @export

make_ipm.general_di_stoch_param <- function(proto_ipm,
                                            return_main_env     = TRUE,
                                            return_all_envs     = FALSE,
                                            usr_funs            = list(),
                                            ...,
                                            domain_list         = NULL,
                                            iterate             = TRUE,
                                            iterations          = 50,
                                            kernel_seq          = NULL,
                                            normalize_pop_size  = TRUE,
                                            report_progress     = FALSE,
                                            iteration_direction = "right") {

  if(iterate &&
     all(is.null(kernel_seq) | is.na(kernel_seq))  &&
     any(proto_ipm$uses_par_sets)) {

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

  proto_list <- .initialize_kernels(proto_ipm, iterate, iteration_direction)

  others <- proto_list$others
  k_row  <- proto_list$k_row

  # Initialize the main_environment so these values can all be found at
  # evaluation time

  if(is.null(domain_list)) {
    main_env <- .make_main_env(others$domain, usr_funs,
                               age_size = .uses_age(others))
  } else {
    main_env <- .make_main_env(domain_list, usr_funs,
                               age_size = .uses_age(others))
  }

  # Bind env_exprs, constants, and pop_vectors to main_env so that
  # we can always find them and avoid that miserable repitition

  main_env <- .bind_all_constants(env_state   = others$env_state[[1]]$constants,
                                    env_to_bind = main_env)

  temp       <- .prep_di_output(others, k_row, proto_ipm, iterations,
                                normalize_pop_size)

  # initialize the pop_state vectors in main_env so they can be found
  # at evaluation time

  if(all(is.na(proto_ipm$pop_state[[1]]))) {

    stop("All general_* IPMs must have a 'pop_state' defined.",
         "See ?define_pop_state() for more details." )

  } else {

    main_env  <- .add_pop_state_to_main_env(temp$pop_state, main_env)

  }

  if(!iterate || iterations < 1) {
    stop("All 'general_*_stoch_param' models must be iterated at least once!",
         call. = FALSE)
  }

  kern_seq <- .make_kern_seq(others,
                             iterations,
                             kernel_seq)

  env_ret <- list()

  for(i in seq_len(iterations)) {

    # add the variable to let the user access the current iteration.

    .bind_iter_var(main_env, i)

    .stoch_progress_message(report_progress, iterations, i)

    # Lazy variant makes sure that whatever functions that generate parameter
    # values stochastically are only evaluated 1 time per iteration! This is so
    # multiple parameters meant to come from a joint distribution really come from
    # the joint distribution!

    sys         <- .make_sub_kernel_general_lazy(others,
                                                 main_env,
                                                 return_envs = return_all_envs)

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
                                main_env,
                                normalize_pop_size)

    if(return_all_envs) {

      sys$ipm_system$env_list$main_env <- NULL

      names(sys$ipm_system$env_list)   <- paste(names(sys$ipm_system$env_list),
                                                "it", i, sep = "_")

      env_ret                          <- c(env_ret, sys$ipm_system$env_list)

    } else if(return_main_env) {

      env_ret <- list(main_env = main_env)

     } else{

      env_ret <- NA_character_

    }

    names(sub_kernels) <- paste(names(sub_kernels), "it", i, sep = "_")

    temp         <- .update_param_output(sub_kernels,
                                         pop_state,
                                         env_ret,
                                         main_env,
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
  names(temp$pop_state) <- gsub("^pop_state_", "n_", names(temp$pop_state))


  out <- list(
    sub_kernels = temp$sub_kernels,
    env_list    = c(list(main_env = main_env), temp$sub_kernel_envs),
    env_seq     = temp$env_seq,
    pop_state   = temp$pop_state,
    proto_ipm   = proto_ipm
  )

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)

  out$sub_kernels                    <- set_ipmr_classes(out$sub_kernels)
  if(inherits(proto_ipm, "age_x_size")) {
    a_s_class <- "age_x_size_ipm"
  } else {
    a_s_class <- NULL
  }

  class(out) <- c('general_di_stoch_param_ipm',
                  a_s_class,
                  'list')

  return(out)

}

# Density dependent methods----------

#' @rdname make_ipm
#' @export

make_ipm.simple_dd_det <- function(proto_ipm,
                                   return_main_env    = TRUE,
                                   return_all_envs    = FALSE,
                                   usr_funs           = list(),
                                   ...,
                                   domain_list        = NULL,
                                   iterate            = TRUE,
                                   iterations         = 50,
                                   normalize_pop_size = FALSE,
                                   report_progress    = FALSE,
                                   iteration_direction = "right") {


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

  proto_list <- .initialize_kernels(proto_ipm, iterate, iteration_direction)

  others <- proto_list$others
  k_row  <- proto_list$k_row   # k_row is either an NA or a proto_ipm object

  # Initialize the main_environment so these values can all be found at
  # evaluation time

  if(is.null(domain_list)) {
    main_env       <- .make_main_env(others$domain, usr_funs,
                                     age_size = .uses_age(others))
  } else {
    main_env       <- .make_main_env(domain_list, usr_funs,
                                     age_size = .uses_age(others))
    proto_ipm$domain <- I(list(domain_list))
  }

  temp           <- .prep_di_output(others,
                                    k_row,
                                    proto_ipm,
                                    iterations,
                                    normalize_pop_size)

  main_env     <- .add_pop_state_to_main_env(temp$pop_state,
                                             main_env)

  # construct the kernels from their function defintions

  # env_list      <- list(main_env = main_env)

  others        <- .prep_dd_vr_exprs(others)
  k_row         <- .prep_dd_k_exprs(k_row)

  if(iterate) {

    sub_kern_out <- list()
    env_ret      <- list(main_env = main_env)

    for(i in seq_len(iterations)) {

      # add the variable to let the user access the current iteration.

      .bind_iter_var(main_env, i)

      .stoch_progress_message(report_progress, iterations, i)

      # Lazy variant makes sure that whatever functions that generate parameter
      # values stochastically are only evaluated 1 time per iteration! This is so
      # multiple parameters meant to come from a joint distribution really come from
      # the joint distribution!

      sys         <- .make_sub_kernel_simple_lazy(others,
                                                  main_env,
                                                  return_envs = return_all_envs,
                                                  dd = "y")

      sub_kernels <- sys$sub_kernels

      # Generate the pop_state for a single iteration! This is critical to ensuring
      # env_state_funs are only evaluated once per iteration. kern_seq = NULL
      # because the environmental parameters are generated on the fly by the
      # user defined function

      pop_state <- .iterate_model(proto_ipm,
                                  k_row,
                                  sub_kernels,
                                  current_iteration = i,
                                  iterations,
                                  temp$pop_state,
                                  main_env,
                                  normalize_pop_size)

      if(return_all_envs) {

        sys$env_list$main_env <- NULL

        names(sys$env_list) <- paste(names(sys$env_list), "it", i, sep = "_")

        env_ret             <- c(env_ret, sys$env_list)

      } else if(return_main_env) {

        env_ret <- list(main_env = main_env)

      } else{

        env_ret <- NA_character_

      }

      temp         <- .update_param_output(sub_kernels,
                                              pop_state,
                                              env_ret,
                                              main_env,
                                              temp,
                                              iterations,
                                              i)

      names(sub_kernels)           <- paste(names(sub_kernels), "it", i, sep = "_")
      sub_kern_out                 <- c(sub_kern_out, sub_kernels)


    }

    # Remove the NA at 1 - it is a dummy so that purrr::map2 can work
    # properly

    pop_state[grepl('lambda', names(pop_state))] <- lapply(
      pop_state[grepl("lambda", names(pop_state))],
      function(x) {
        out <- x[ , -1, drop = FALSE]
        return(out)
      }
    )
    # Convert pop_state names back to the user-supplied ones

    names(pop_state) <- gsub("^pop_state_", "n_", names(pop_state))

  } else {

    pop_state <- NA_real_

  }

  out_seq      <- NA_integer_

  sub_kern_out <- set_ipmr_classes(sub_kern_out)

  # Not storing iteration kernels for now, though that *should* be
  # fairly easy to change...

  out <- list(sub_kernels = sub_kern_out,
              env_list    = env_ret,
              env_seq     = out_seq,
              pop_state   = pop_state,
              proto_ipm   = proto_ipm
  )

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)

  class(out)                         <- c('simple_dd_det_ipm',
                                          'list')

  return(out)


}

#' @rdname make_ipm
#'
#' @export

make_ipm.simple_dd_stoch_kern <- function(proto_ipm,
                                          return_main_env    = TRUE,
                                          return_all_envs    = FALSE,
                                          usr_funs           = list(),
                                          ...,
                                          domain_list        = NULL,
                                          iterate            = TRUE,
                                          iterations         = 50,
                                          kernel_seq         = NA_character_,
                                          normalize_pop_size = FALSE,
                                          report_progress    = FALSE,
                                          iteration_direction = "right") {


  if(iterate &&
     all(is.null(kernel_seq) | is.na(kernel_seq))) {

    message("'kernel_seq' not defined. Will generate one internally")

    kernel_seq <- "internal"

  }

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

  proto_list <- .initialize_kernels(proto_ipm, iterate, iteration_direction)

  others <- proto_list$others
  k_row  <- proto_list$k_row   # k_row is either an NA or a proto_ipm object

  # Initialize the main_environment so these values can all be found at
  # evaluation time

  if(is.null(domain_list)) {
    main_env       <- .make_main_env(others$domain, usr_funs,
                                     age_size = .uses_age(others))
  } else {
    main_env       <- .make_main_env(domain_list, usr_funs,
                                     age_size = .uses_age(others))
    proto_ipm$domain <- I(list(domain_list))
  }

  temp           <- .prep_di_output(others,
                                    k_row,
                                    proto_ipm,
                                    iterations,
                                    normalize_pop_size)

  main_env     <- .add_pop_state_to_main_env(temp$pop_state,
                                             main_env)

  # construct the kernels from their function defintions

  # env_list      <- list(main_env = main_env)

  others        <- .prep_dd_vr_exprs(others)
  k_row         <- .prep_dd_k_exprs(k_row)

  if(iterate) {

    sub_kern_out <- list()
    env_ret      <- list(main_env = main_env)
    kern_seq     <- .make_kern_seq(others,
                                   iterations,
                                   kernel_seq)

    for(i in seq_len(iterations)) {

      # add the variable to let the user access the current iteration.

      .bind_iter_var(main_env, i)

      .stoch_progress_message(report_progress, iterations, i)

      # Lazy variant makes sure that whatever functions that generate parameter
      # values stochastically are only evaluated 1 time per iteration! This is so
      # multiple parameters meant to come from a joint distribution really come from
      # the joint distribution!

      sys         <- .make_sub_kernel_simple_lazy(others,
                                                  main_env,
                                                  return_envs = return_all_envs,
                                                  dd = "y")

      sub_kernels <- sys$sub_kernels

      # Generate the pop_state for a single iteration! This is critical to ensuring
      # env_state_funs are only evaluated once per iteration. kern_seq = NULL
      # because the environmental parameters are generated on the fly by the
      # user defined function

      pop_state <- .iterate_model(proto_ipm,
                                  k_row,
                                  sub_kernels,
                                  current_iteration = i,
                                  iterations,
                                  kern_seq,
                                  temp$pop_state,
                                  main_env,
                                  normalize_pop_size,
                                  report_progress)

      if(return_all_envs) {

        sys$env_list$main_env <- NULL

        names(sys$env_list) <- paste(names(sys$env_list), "it", i, sep = "_")

        env_ret             <- c(env_ret, sys$env_list)

      } else if(return_main_env) {

        env_ret <- list(main_env = main_env)

      } else{

        env_ret <- NA_character_

      }

      temp         <- .update_param_output(sub_kernels,
                                              pop_state,
                                              env_ret,
                                              main_env,
                                              temp,
                                              iterations,
                                              i)

      names(sub_kernels)           <- paste(names(sub_kernels), "it", i, sep = "_")
      sub_kern_out                 <- c(sub_kern_out, sub_kernels)


    }

    # Remove the NA at 1 - it is a dummy so that purrr::map2 can work
    # properly

    pop_state[grepl('lambda', names(pop_state))] <- lapply(
      pop_state[grepl("lambda", names(pop_state))],
      function(x) {
        out <- x[ , -1, drop = FALSE]
        return(out)
      }
    )
    # Convert pop_state names back to the user-supplied ones

    names(pop_state) <- gsub("^pop_state_", "n_", names(pop_state))

  } else {

    pop_state <- NA_real_

  }

  out_seq      <- kern_seq

  sub_kern_out <- set_ipmr_classes(sub_kern_out)

  # Not storing iteration kernels for now, though that *should* be
  # fairly easy to change...

  out <- list(sub_kernels = sub_kern_out,
              env_list    = env_ret,
              env_seq     = out_seq,
              pop_state   = pop_state,
              proto_ipm   = proto_ipm
  )

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)

  class(out)                         <- c('simple_dd_stoch_kern_ipm',
                                        'list')

  return(out)

}

#' @rdname make_ipm
#'
#' @export

make_ipm.simple_dd_stoch_param <-function(proto_ipm,
                                          return_main_env    = TRUE,
                                          return_all_envs    = FALSE,
                                          usr_funs           = list(),
                                          ...,
                                          domain_list        = NULL,
                                          iterate            = TRUE,
                                          iterations         = 50,
                                          kernel_seq         = NA_character_,
                                          normalize_pop_size = FALSE,
                                          report_progress    = FALSE,
                                          iteration_direction = "right") {


  if(iterate &&
     all(is.null(kernel_seq) | is.na(kernel_seq)) &&
     any(proto_ipm$uses_par_sets)) {

    message("'kernel_seq' not defined. Will generate one internally")

    kernel_seq <- "internal"

  }

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

  proto_list <- .initialize_kernels(proto_ipm, iterate, iteration_direction)

  others <- proto_list$others
  k_row  <- proto_list$k_row   # k_row is either an NA or a proto_ipm object

  # Initialize the main_environment so these values can all be found at
  # evaluation time

  if(is.null(domain_list)) {
    main_env       <- .make_main_env(others$domain, usr_funs,
                                     age_size = .uses_age(others))
  } else {
    main_env       <- .make_main_env(domain_list, usr_funs,
                                     age_size = .uses_age(others))
    proto_ipm$domain <- I(list(domain_list))
  }

  temp           <- .prep_di_output(others,
                                    k_row,
                                    proto_ipm,
                                    iterations,
                                    normalize_pop_size)


  # Bind env_exprs, constants, and pop_vectors to main_env so that
  # we can always find them and avoid that miserable repitition

  main_env     <- .bind_all_constants(env_state   = others$env_state[[1]]$constants,
                                      env_to_bind = main_env)

  main_env     <- .add_pop_state_to_main_env(temp$pop_state,
                                             main_env)

  # construct the kernels from their function defintions

  # env_list      <- list(main_env = main_env)

  others        <- .prep_dd_vr_exprs(others)
  k_row         <- .prep_dd_k_exprs(k_row)

  if(iterate) {

    sub_kern_out <- list()
    env_ret      <- list(main_env = main_env)
    kern_seq     <- .make_kern_seq(others,
                                   iterations,
                                   kernel_seq)

    for(i in seq_len(iterations)) {

      # add the variable to let the user access the current iteration.

      .bind_iter_var(main_env, i)

      .stoch_progress_message(report_progress, iterations, i)

      # Lazy variant makes sure that whatever functions that generate parameter
      # values stochastically are only evaluated 1 time per iteration! This is so
      # multiple parameters meant to come from a joint distribution really come from
      # the joint distribution!

      sys         <- .make_sub_kernel_simple_lazy(others,
                                                  main_env,
                                                  return_envs = return_all_envs,
                                                  dd = "n")

      sub_kernels <- sys$ipm_system$sub_kernels

      # Generate the pop_state for a single iteration! This is critical to ensuring
      # env_state_funs are only evaluated once per iteration. kern_seq = NULL
      # because the environmental parameters are generated on the fly by the
      # user defined function

      pop_state <- .iterate_model(proto_ipm,
                                  k_row,
                                  sub_kernels,
                                  current_iteration = i,
                                  iterations,
                                  kern_seq,
                                  temp$pop_state,
                                  main_env,
                                  normalize_pop_size,
                                  report_progress)

      if(return_all_envs) {

        sys$env_list$main_env <- NULL

        names(sys$env_list) <- paste(names(sys$env_list), "it", i, sep = "_")

        env_ret             <- c(env_ret, sys$env_list)

      } else if(return_main_env) {

        env_ret <- list(main_env = main_env)

      } else{

        env_ret <- NA_character_

      }

      names(sub_kernels)           <- paste(names(sub_kernels),
                                            "it",
                                            i,
                                            sep = "_")

      temp         <- .update_param_output(sub_kernels,
                                                   pop_state,
                                                   env_ret,
                                                   main_env,
                                                   temp,
                                                   iterations,
                                                   i)

    }

    # Remove the NA at 1 - it is a dummy so that purrr::map2 can work
    # properly

    pop_state[grepl('lambda', names(pop_state))] <- lapply(
      pop_state[grepl("lambda", names(pop_state))],
      function(x) {
        out <- x[ , -1, drop = FALSE]
        return(out)
      }
    )
    # Convert pop_state names back to the user-supplied ones

    names(pop_state) <- gsub("^pop_state_", "n_", names(pop_state))

  } else {

    pop_state <- NA_real_

  }

  temp$env_seq <- data.frame(temp$env_seq, stringsAsFactors = FALSE)

  if(!is.null(kern_seq)) {

    temp$env_seq <- cbind(temp$env_seq,
                          kernel_seq = kern_seq)

  }

  out_seq <- temp$env_seq

  sub_kern_out <- set_ipmr_classes(sub_kern_out)

  out <- list(sub_kernels = temp$sub_kernels,
              env_list    = env_ret,
              env_seq     = out_seq,
              pop_state   = pop_state,
              proto_ipm   = proto_ipm
  )

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)

  class(out)                         <- c('simple_dd_stoch_param_ipm',
                                          'list')

  return(out)

}

#' @rdname make_ipm
#'
#' @export

make_ipm.general_dd_det <- function(proto_ipm,
                                    return_main_env = TRUE,
                                    return_all_envs = FALSE,
                                    usr_funs = list(),
                                    ...,
                                    domain_list = NULL,
                                    iterate = TRUE,
                                    iterations = 50,
                                    normalize_pop_size = FALSE,
                                    report_progress = FALSE,
                                    iteration_direction = "right"
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

  proto_list <- .initialize_kernels(proto_ipm, iterate, iteration_direction)

  others <- proto_list$others
  k_row  <- proto_list$k_row   # k_row is either an NA or a proto_ipm object

  # Initialize the main_environment so these values can all be found at
  # evaluation time

  if(is.null(domain_list)) {
    main_env       <- .make_main_env(others$domain, usr_funs,
                                     age_size = .uses_age(others))
  } else {
    main_env       <- .make_main_env(domain_list, usr_funs,
                                     age_size = .uses_age(others))
    proto_ipm$domain <- I(list(domain_list))
  }

  main_env <- .bind_all_constants(env_state   = others$env_state[[1]],
                                  env_to_bind = main_env)

  temp           <- .prep_di_output(others, k_row,
                                    proto_ipm, iterations,
                                    normalize_pop_size)

  if(all(is.na(proto_ipm$pop_state[[1]]))) {

    stop("All general_* IPMs must have a 'pop_state' defined.",
         "See ?define_pop_state() for more details." )

  } else {

    main_env  <- .add_pop_state_to_main_env(temp$pop_state, main_env)

  }

  others        <- .prep_dd_vr_exprs(others)
  k_row         <- .prep_dd_k_exprs(k_row)


  if(iterate) {

    sub_kern_out <- list()
    env_ret      <- list(main_env = main_env)

    for(i in seq_len(iterations)) {


      .bind_iter_var(main_env, i)

      .stoch_progress_message(report_progress, iterations, i)

      # Lazy variant makes sure that whatever functions that generate parameter
      # values stochastically are only evaluated 1 time per iteration! This is so
      # multiple parameters meant to come from a joint distribution really come from
      # the joint distribution!

      sys         <- .make_sub_kernel_simple_lazy(others,
                                                  main_env,
                                                  return_envs = return_all_envs,
                                                  dd = "y")

      sub_kernels <- sys$sub_kernels

      pop_state <- .iterate_model(proto_ipm,
                                  k_row,
                                  sub_kernels,
                                  current_iteration = i,
                                  iterations,
                                  temp$pop_state,
                                  main_env,
                                  normalize_pop_size)

      if(return_all_envs) {

        sys$env_list$main_env <- NULL

        names(sys$env_list) <- paste(names(sys$env_list), "it", i, sep = "_")

        env_ret             <- c(env_ret, sys$env_list)

      } else if(return_main_env) {

        env_ret <- list(main_env = main_env)

      } else{

        env_ret <- NA_character_

      }

      temp         <- .update_param_output(sub_kernels,
                                           pop_state,
                                           env_ret,
                                           main_env,
                                           temp,
                                           iterations,
                                           i)

      names(sub_kernels)           <- paste(names(sub_kernels), "it", i, sep = "_")
      sub_kern_out                 <- c(sub_kern_out, sub_kernels)

    }

    # Remove the NA at 1 - it is a dummy so that purrr::map2 can work
    # properly

    pop_state[grepl('lambda', names(pop_state))] <- lapply(
      pop_state[grepl("lambda", names(pop_state))],
      function(x) {
        out <- x[ , -1, drop = FALSE]
        return(out)
      }
    )
    # Convert pop_state names back to the user-supplied ones

    names(pop_state) <- gsub("^pop_state_", "n_", names(pop_state))

  } else {

    pop_state <- NA_real_

  }


  if(return_all_envs) {

    out_ret <- temp$env_list

  } else if(return_main_env) {

    out_ret <- list(main_env = main_env)

  } else {

    out_ret <- NA_character_

  }

  out_seq      <- NA_integer_

  sub_kern_out <- set_ipmr_classes(sub_kern_out)

  # Not storing iteration kernels for now, though that *should* be
  # fairly easy to change...

  out <- list(sub_kernels = sub_kern_out,
              env_list    = env_ret,
              env_seq     = out_seq,
              pop_state   = pop_state,
              proto_ipm   = proto_ipm
  )

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)

  if(inherits(proto_ipm, "age_x_size")) {
    a_s_class <- "age_x_size_ipm"
  } else {
    a_s_class <- NULL
  }

  class(out)                         <- c('general_dd_det_ipm',
                                          a_s_class,
                                          'list')

  return(out)

}

#' @rdname make_ipm
#'
#' @export

make_ipm.general_dd_stoch_kern <- function(proto_ipm,
                                           return_main_env    = TRUE,
                                           return_all_envs    = FALSE,
                                           usr_funs           = list(),
                                           ...,
                                           domain_list        = NULL,
                                           iterate            = TRUE,
                                           iterations         = 50,
                                           kernel_seq         = NA_character_,
                                           normalize_pop_size = FALSE,
                                           report_progress    = FALSE,
                                           iteration_direction = "right") {


  if(iterate &&
     all(is.null(kernel_seq) | is.na(kernel_seq))) {

    message("'kernel_seq' not defined. Will generate one internally")

    kernel_seq <- "internal"

  }

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

  proto_list <- .initialize_kernels(proto_ipm, iterate, iteration_direction)

  others <- proto_list$others
  k_row  <- proto_list$k_row   # k_row is either an NA or a proto_ipm object

  # Initialize the main_environment so these values can all be found at
  # evaluation time

  if(is.null(domain_list)) {
    main_env       <- .make_main_env(others$domain, usr_funs,
                                     age_size = .uses_age(others))
  } else {
    main_env       <- .make_main_env(domain_list, usr_funs,
                                     age_size = .uses_age(others))
    proto_ipm$domain <- I(list(domain_list))
  }

  temp           <- .prep_di_output(others,
                                    k_row,
                                    proto_ipm,
                                    iterations,
                                    normalize_pop_size)

  main_env     <- .add_pop_state_to_main_env(temp$pop_state,
                                             main_env)

  # construct the kernels from their function defintions

  # env_list      <- list(main_env = main_env)

  others        <- .prep_dd_vr_exprs(others)
  k_row         <- .prep_dd_k_exprs(k_row)

  if(iterate) {

    sub_kern_out <- list()
    env_ret      <- list(main_env = main_env)
    kern_seq     <- .make_kern_seq(others,
                                   iterations,
                                   kernel_seq)

    for(i in seq_len(iterations)) {

      # add the variable to let the user access the current iteration.

      .bind_iter_var(main_env, i)

      .stoch_progress_message(report_progress, iterations, i)

      # Lazy variant makes sure that whatever functions that generate parameter
      # values stochastically are only evaluated 1 time per iteration! This is so
      # multiple parameters meant to come from a joint distribution really come from
      # the joint distribution!

      sys         <- .make_sub_kernel_simple_lazy(others,
                                                  main_env,
                                                  return_envs = return_all_envs,
                                                  dd = "y")

      sub_kernels <- sys$sub_kernels

      # Generate the pop_state for a single iteration! This is critical to ensuring
      # env_state_funs are only evaluated once per iteration. kern_seq = NULL
      # because the environmental parameters are generated on the fly by the
      # user defined function

      pop_state <- .iterate_model(proto_ipm,
                                  k_row,
                                  sub_kernels,
                                  current_iteration = i,
                                  iterations,
                                  kern_seq,
                                  temp$pop_state,
                                  main_env,
                                  normalize_pop_size,
                                  report_progress)

      if(return_all_envs) {

        sys$env_list$main_env <- NULL

        names(sys$env_list) <- paste(names(sys$env_list), "it", i, sep = "_")

        env_ret             <- c(env_ret, sys$env_list)

      } else if(return_main_env) {

        env_ret <- list(main_env = main_env)

      } else{

        env_ret <- NA_character_

      }

      temp         <- .update_param_output(sub_kernels,
                                              pop_state,
                                              env_ret,
                                              main_env,
                                              temp,
                                              iterations,
                                              i)

      names(sub_kernels)           <- paste(names(sub_kernels), "it", i, sep = "_")
      sub_kern_out                 <- c(sub_kern_out, sub_kernels)


    }

    # Remove the NA at 1 - it is a dummy so that purrr::map2 can work
    # properly

    pop_state[grepl('lambda', names(pop_state))] <- lapply(
      pop_state[grepl("lambda", names(pop_state))],
      function(x) {
        out <- x[ , -1, drop = FALSE]
        return(out)
      }
    )
    # Convert pop_state names back to the user-supplied ones

    names(pop_state) <- gsub("^pop_state_", "n_", names(pop_state))

  } else {

    pop_state <- NA_real_

  }

  out_seq      <- kern_seq

  sub_kern_out <- set_ipmr_classes(sub_kern_out)

  # Not storing iteration kernels for now, though that *should* be
  # fairly easy to change...

  out <- list(sub_kernels = sub_kern_out,
              env_list    = env_ret,
              env_seq     = out_seq,
              pop_state   = pop_state,
              proto_ipm   = proto_ipm
  )

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)

  if(inherits(proto_ipm, "age_x_size")) {
    a_s_class <- "age_x_size_ipm"
  } else {
    a_s_class <- NULL
  }

  class(out)                         <- c('general_dd_stoch_kern_ipm',
                                          a_s_class,
                                          'list')

  return(out)

}

#' @rdname make_ipm
#'
#' @export

make_ipm.general_dd_stoch_param <-function(proto_ipm,
                                          return_main_env    = TRUE,
                                          return_all_envs    = FALSE,
                                          usr_funs           = list(),
                                          ...,
                                          domain_list        = NULL,
                                          iterate            = TRUE,
                                          iterations         = 50,
                                          kernel_seq         = NA_character_,
                                          normalize_pop_size = FALSE,
                                          report_progress    = FALSE,
                                          iteration_direction = "right") {


  if(iterate &&
     all(is.null(kernel_seq) | is.na(kernel_seq)) &&
     any(proto_ipm$uses_par_sets)) {

    message("'kernel_seq' not defined. Will generate one internally")

    kernel_seq <- "internal"

  }

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

  proto_list <- .initialize_kernels(proto_ipm, iterate, iteration_direction)

  others <- proto_list$others
  k_row  <- proto_list$k_row   # k_row is either an NA or a proto_ipm object

  # Initialize the main_environment so these values can all be found at
  # evaluation time

  if(is.null(domain_list)) {
    main_env       <- .make_main_env(others$domain, usr_funs,
                                     age_size = .uses_age(others))
  } else {
    main_env       <- .make_main_env(domain_list, usr_funs,
                                     age_size = .uses_age(others))
    proto_ipm$domain <- I(list(domain_list))
  }

  temp           <- .prep_di_output(others,
                                    k_row,
                                    proto_ipm,
                                    iterations,
                                    normalize_pop_size)


  # Bind env_exprs, constants, and pop_vectors to main_env so that
  # we can always find them and avoid that miserable repitition

  main_env     <- .bind_all_constants(env_state   = others$env_state[[1]]$constants,
                                      env_to_bind = main_env)

  main_env     <- .add_pop_state_to_main_env(temp$pop_state,
                                             main_env)

  # construct the kernels from their function defintions

  # env_list      <- list(main_env = main_env)

  others        <- .prep_dd_vr_exprs(others)
  k_row         <- .prep_dd_k_exprs(k_row)

  if(iterate) {

    sub_kern_out <- list()
    env_ret      <- list(main_env = main_env)
    kern_seq     <- .make_kern_seq(others,
                                   iterations,
                                   kernel_seq)

    for(i in seq_len(iterations)) {

      # add the variable to let the user access the current iteration.

      .bind_iter_var(main_env, i)

      .stoch_progress_message(report_progress, iterations, i)

      # Lazy variant makes sure that whatever functions that generate parameter
      # values stochastically are only evaluated 1 time per iteration! This is so
      # multiple parameters meant to come from a joint distribution really come from
      # the joint distribution!

      sys         <- .make_sub_kernel_simple_lazy(others,
                                                  main_env,
                                                  return_envs = return_all_envs,
                                                  dd = "n")

      sub_kernels <- sys$ipm_system$sub_kernels

      # Generate the pop_state for a single iteration! This is critical to ensuring
      # env_state_funs are only evaluated once per iteration. kern_seq = NULL
      # because the environmental parameters are generated on the fly by the
      # user defined function

      pop_state <- .iterate_model(proto_ipm,
                                  k_row,
                                  sub_kernels,
                                  current_iteration = i,
                                  iterations,
                                  kern_seq,
                                  temp$pop_state,
                                  main_env,
                                  normalize_pop_size,
                                  report_progress)

      if(return_all_envs) {

        sys$env_list$main_env <- NULL

        names(sys$env_list) <- paste(names(sys$env_list), "it", i, sep = "_")

        env_ret             <- c(env_ret, sys$env_list)

      } else if(return_main_env) {

        env_ret <- list(main_env = main_env)

      } else{

        env_ret <- NA_character_

      }

      names(sub_kernels) <- paste(names(sub_kernels),
                                  "it",
                                  i,
                                  sep = "_")

      temp               <- .update_param_output(sub_kernels,
                                                 pop_state,
                                                 env_ret,
                                                 main_env,
                                                 temp,
                                                 iterations,
                                                 i)

    }

    # Remove the NA at 1 - it is a dummy so that purrr::map2 can work
    # properly

    pop_state[grepl('lambda', names(pop_state))] <- lapply(
      pop_state[grepl("lambda", names(pop_state))],
      function(x) {
        out <- x[ , -1, drop = FALSE]
        return(out)
      }
    )
    # Convert pop_state names back to the user-supplied ones

    names(pop_state) <- gsub("^pop_state_", "n_", names(pop_state))

  } else {

    pop_state <- NA_real_

  }

  temp$env_seq <- data.frame(temp$env_seq,
                             stringsAsFactors = FALSE)

  if(!is.null(kern_seq)) {

    temp$env_seq <- cbind(temp$env_seq,
                          kernel_seq = kern_seq)

  }

  out_seq <- temp$env_seq

  sub_kern_out <- set_ipmr_classes(sub_kern_out)

  out <- list(sub_kernels = temp$sub_kernels,
              env_list    = env_ret,
              env_seq     = out_seq,
              pop_state   = pop_state,
              proto_ipm   = proto_ipm
  )

  attr(out$proto_ipm, 'implemented') <- TRUE
  attr(out, "iterated")              <- (iterate && iterations >= 1)

  if(inherits(proto_ipm, "age_x_size")) {
    a_s_class <- "age_x_size_ipm"
  } else {
    a_s_class <- NULL
  }

  class(out)                         <- c('general_dd_stoch_param_ipm',
                                          a_s_class,
                                          'list')

  return(out)

}
