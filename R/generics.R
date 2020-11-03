# Print----------------

#' @title Print proto_ipms or *_ipm objects
#' @rdname print_star
#'
#' @param x An object of class \code{proto_ipm} or produced by \code{make_ipm()}.
#' @param ... Ignored
#'
#' @details For printing \code{proto_ipm} objects, suffixes are wrapped in
#' \code{<suffix>} to assist with debugging. These are not carried into the model,
#' just a visual aid.
#'
#' @export

print.proto_ipm <- function(x, ...) {

  n_kerns      <- dim(x)[1]

  cls_switch   <- class(x)[1]

  pretty_class <- .pretty_class(cls_switch)


  msg          <- paste(ifelse(.has_age(x),
                               "An age structured",
                               "A"),
                        pretty_class,
                        "proto_ipm with",
                        n_kerns,
                        "kernels defined:\n")

  msg          <- c(msg, paste(x$kernel_id, collapse = ', '))

  msg          <- c(msg, '\n\nKernel formulae:\n\n')

  cat(msg)

  forms        <- kernel_formulae(x)
  print(forms)

  cat("\n\nVital rates:\n\n")
  vrs          <- vital_rates(x)

  print(vrs)

  cat("\n\nParameter names:\n\n")

  pars <- parameters(x)

  print(names(pars))

  # Check if all parameters are in vital rate expressions, and if all
  # vital rate expression args are found in parameters.

  in_mod <- .all_params_in_mod(x)

  in_mod_tst <- vapply(in_mod, function(x) x$out, logical(1L)) %>%
    all()

  cat("\nAll parameters in vital rate expressions found in 'data_list': ",
      in_mod_tst)

  if(!in_mod_tst) {
    stoch_params <- c(paste("simple_d", c('i', 'd'), "_stoch_param", sep = ""),
                      paste("general_d", c('i', 'd'), "_stoch_param", sep = ""))
    if(inherits(x, stoch_params)) {
      cat("\nDid not test parameters in 'define_env_state()'.\n")
    }

    cat("\nItems in vital rate expressions not found in 'data_list':\n")

    for(i in seq_along(in_mod)) {

      if(!in_mod[[i]]$out) {
        nm <- paste("\n*", names(in_mod)[i], "*\n", sep = "")
        miss_msg <- paste(nm,
                          paste(in_mod[[i]]$missing, collapse = "\n"),
                          sep = "")
        cat(miss_msg, '\n')
      }

    }
    cat("\nDouble check the model definition and 'data_list'!\n")

  }

  cat("\n\nDomains for state variables:\n\n")

  print(domains(x))

  cat("\n\nPopulation states defined:\n\n")

  print(pop_state(x))

  # Add more later -----------

  invisible(x)
}

#' @noRd

.all_params_in_mod <- function(proto) {

  kern_list <- list()

  for(row in seq_len(nrow(proto))) {

    use_proto <- proto[row, ]
    kern_nm   <- proto$kernel_id[row]

    vr_exprs <- vital_rates(use_proto)
    params   <- parameters(proto)

    # Age x size models and/or hier_effs models need expansion before we
    # we can test whether all parameters are in the model

    if(.has_age(proto) | any(proto$has_hier_effs)) {

      all_levs <- .flatten_to_depth(proto$levels_hier_effs[proto$has_hier_effs],
                                    1L) %>%
        .[!duplicated(names(.))]

      hier_nms <- names(all_levs)

      if(length(all_levs) > 1) {
        lev_df <- expand.grid(all_levs,
                              stringsAsFactors = FALSE)
      } else {

        lev_df <- as.data.frame(all_levs)

      }

      # Hold the results

      all_vr_exprs <- list()

      it <- 1

      # Replace suffixes in expressions
      for(i in seq_along(vr_exprs)) {

        base_expr <- rlang::expr_text(vr_exprs[[i]])

        for(j in seq_len(dim(lev_df)[2])) {

          hier_nm <- names(lev_df)[j]

          for(k in seq_len(dim(lev_df)[1])) {

            # Replace suffixes in expression and name
            all_vr_exprs[[it]]      <- gsub(hier_nm,
                                            lev_df[k, j],
                                            base_expr)

            names(all_vr_exprs)[it] <- gsub(hier_nm,
                                            lev_df[k, j],
                                            names(vr_exprs)[i])


            it <- it + 1

          }

        }


      }

      vr_exprs <- lapply(all_vr_exprs, rlang::parse_expr)

      vr_args  <- .vr_args(vr_exprs)

    } else {

      vr_args <- .vr_args(vr_exprs)

    }

    # If a vital rate expression creates another parameter, then remove that
    # from this test, as it shouldn't appear in the data_list anyway.

    to_exclude <- names(vr_exprs)
    domains    <- names(domains(proto))

    vr_args    <- Filter(.can_be_number, vr_args)

    out <- vapply(vr_args[!vr_args %in% to_exclude],
                  function(nm, pars) nm %in% pars,
                  logical(1L),
                  pars = c(names(params), domains))

    if(!all(out)) {
      vr_args <- vr_args[!vr_args %in% to_exclude]
      ind     <- which(!vr_args %in% c(names(params), domains))

      miss_args <- vr_args[ind]
      out       <- FALSE

    } else {

      out       <- TRUE
      miss_args <- NA_character_
    }

    kern_list <- c(kern_list,
                   list(list(out = out,
                        missing = miss_args)))

    names(kern_list)[row] <- kern_nm
  }

  return(kern_list)

}

.can_be_number <- function(x) {
  suppressWarnings(is.na(as.numeric(x)))
}

#' @noRd

.vr_args <- function(vr_exprs) {

  all_args <- lapply(vr_exprs, function(y) {
    temp <- rlang::expr_text(y)
    .args_from_txt(temp)
  }) %>%
    unlist() %>%
    unique()

  return(all_args)

}

#' @noRd

.has_age <- function(x) {

  inherits(x, "age_x_size")

}

#' @noRd

.pretty_class <- function(cls_switch) {
  out <-
    switch(
      cls_switch,

      # IPM classes -----------------------------

      'simple_di_det_ipm'          = "simple, density independent, deterministic",
      'simple_di_stoch_kern_ipm'   = "simple, density independent, stochastic, kernel-resampled",
      'simple_di_stoch_param_ipm'  = "simple, density independent, stochastic, parameter-resampled",

      'simple_dd_det_ipm'          = "simple, density dependent, deterministic",
      'simple_dd_stoch_kern_ipm'   = "simple, density dependent, stochastic, kernel-resampled",
      'simple_dd_stoch_param_ipm'  = "simple, density dependent, stochastic, parameter-resampled",

      'general_di_det_ipm'         = "general, density independent, deterministic",
      'general_di_stoch_kern_ipm'  = "general, density independent, stochastic, kernel-resampled",
      'general_di_stoch_param_ipm' = "general, density independent, stochastic, parameter-resampled",

      'general_dd_det_ipm'         = "general, density dependent, deterministic",
      'general_dd_stoch_kern_ipm'  = "general, density dependent, stochastic, kernel-resampled",
      'general_dd_stoch_param_ipm' = "general, density dependent, stochastic, parameter-resampled",

      # proto_ipm classes -----------------------

      'simple_di_det'          = "simple, density independent, deterministic",
      'simple_di_stoch_kern'   = "simple, density independent, stochastic, kernel-resampled",
      'simple_di_stoch_param'  = "simple, density independent, stochastic, parameter-resampled",

      'simple_dd_det'          = "simple, density dependent, deterministic",
      'simple_dd_stoch_kern'   = "simple, density dependent, stochastic, kernel-resampled",
      'simple_dd_stoch_param'  = "simple, density dependent, stochastic, parameter-resampled",

      'general_di_det'         = "general, density independent, deterministic",
      'general_di_stoch_kern'  = "general, density independent, stochastic, kernel-resampled",
      'general_di_stoch_param' = "general, density independent, stochastic, parameter-resampled",

      'general_dd_det'         = "general, density dependent, deterministic",
      'general_dd_stoch_kern'  = "general, density dependent, stochastic, kernel-resampled",
      'general_dd_stoch_param' = "general, density dependent, stochastic, parameter-resampled"

    )

  return(out)
}

#' @rdname print_star
#' @title Generics for IPM classes
#'
#' @param comp_lambda A logical indicating whether or not to calculate lambdas
#' for the iteration kernels and display them.
#' @param comp_method Either \code{"pop_size"} or \code{"eigen"}. \code{"pop_size"}
#' computes lambda as the ratio of population sizes between successive time steps
#' and then takes the geometric mean. \code{"eigen"} computes the dominant eigenvalue
#' of each iteration kernel. See \code{\link{lambda}} for more details.
#' @param sig_digits The number of significant digits to round to if \code{
#' comp_lambda = TRUE}.
#' @param type_lambda Either \code{'all'} or \code{'stochastic'}. See
#' \code{\link{lambda}} for more details.
#' @param check_conv A logical: for \code{general_*} models, check if population state
#' has converged to asymptotic dynamics? If \code{TRUE} and the model has not
#' converged, a message will be printed. Only applies to \code{*_det}  when
#' \code{comp_method = 'pop_size'}
#' @param ... Ignored
#'
#' @return \code{x} invisibly.
#'
#' @export

print.simple_di_det_ipm <- function(x,
                                    comp_lambda = TRUE,
                                    comp_method = 'eigen',
                                    type_lambda = 'all',
                                    sig_digits = 3,
                                    check_conv = TRUE,
                                    ...) {

  pretty_cls <- .pretty_class(class(x)[1])

  msg <- paste0('A ',
                pretty_cls,
                ' IPM with ',
                length(x$iterators),
                ' iteration kernel(s) and ',
                length(x$sub_kernels),
                ' sub-kernel(s) defined.', sep = "")

  if(comp_lambda) {


    nm_ks  <- names(x$iterators)

    lambdas <- lambda(x, comp_method = comp_method, type_lambda = type_lambda)

    l_msg  <- paste0('\nDeterministic lambda for ', nm_ks,' = ', lambdas, sep = "")

    msg    <- c(msg, l_msg)

    if(comp_method == 'pop_size' &&
       check_conv &&
       ! .is_conv_to_asymptotic(lambdas)) {

      # Captures the name of the model that the user gave rather than
      # just print "x isn't converged"

      mod_nm <- deparse(substitute(x))

      message(
        paste(mod_nm,
              ' has has not converged to asymptotic dynamics!',
              sep = "")
      )

    }

  }

  cat(msg)

  invisible(x)


}

#' @rdname print_star
#' @export

print.simple_di_stoch_kern_ipm <- function(x,
                                           comp_lambda = TRUE,
                                           type_lambda = 'stochastic',
                                           sig_digits = 3,
                                           ...) {

  mod_nm     <- deparse(substitute(x))
  pretty_cls <- .pretty_class(class(x)[1])

  msg <- paste0('A ',
                pretty_cls,
                ' IPM with ',
                length(x$iterators),
                ' iteration kernel(s) and ',
                length(x$sub_kernels),
                ' sub-kernel(s) defined.', sep = "")

  if(comp_lambda) {

    all_lams <- lambda(x, type_lambda = type_lambda)


    l_msg  <- paste0('\nStochastic lambda for ',
                     mod_nm,
                     ' = ',
                     round(all_lams, sig_digits),
                     sep = "" )

    det_lam_msg <- paste('\nCall lambda(',
                         mod_nm,
                         ', type_lambda = "all") for deterministic lambdas\n',
                         'from each iteration.',
                         sep = "")
    msg <- c(msg, l_msg, det_lam_msg)
  }
  cat(msg)

  invisible(x)
}

#' @rdname print_star
#' @export

print.simple_di_stoch_param_ipm <- function(x,
                                            comp_lambda = TRUE,
                                            type_lambda = 'stochastic',
                                            sig_digits  = 3,
                                            ...) {
  pretty_cls <- .pretty_class(class(x)[1])

  mod_nm     <- deparse(substitute(x))

  msg <- paste0('A ',
                pretty_cls,
                ' IPM with ',
                length(x$iterators),
                ' iteration kernel(s) and ',
                length(x$sub_kernels),
                ' sub-kernel(s) defined.', sep = "")

  if(comp_lambda) {

    all_lams <- lambda(x, type_lambda = type_lambda)

    l_msg  <- paste0('\nStochastic lambda for ',
                     mod_nm,
                     ' = ',
                     round(all_lams, sig_digits),
                     sep = "" )

    det_lam_msg <- paste('\nCall lambda(',
                         mod_nm,
                         ', type_lambda = "all") for deterministic lambdas\n',
                         'from each iteration.',
                         sep = "")

    msg <- c(msg, l_msg, det_lam_msg)
  }
  cat(msg)

  invisible(x)

}

#' @rdname print_star
#' @export

print.general_di_det_ipm <- function(x,
                                     comp_lambda = TRUE,
                                     type_lambda = 'last',
                                     sig_digits  = 3,
                                     check_conv  = TRUE,
                                     ...) {

  pretty_cls <- .pretty_class(class(x)[1])

  pops <- x$pop_state[!grepl("lambda", names(x$pop_state))]

  msg <- paste0('A ',
                pretty_cls,
                ' IPM with ',
                length(x$sub_kernels),
                ' sub-kernel(s) and ',
                length(pops),
                ' population vectors defined.', sep = "")

  mod_nm <- deparse(substitute(x))


  if(comp_lambda) {

    ret_lam <- lambda(x, type_lambda = type_lambda)

    l_msg <- paste('\nLambda for the final time step of the model is: ',
                   round(ret_lam, sig_digits),
                   '\nCall lambda(',
                   mod_nm,
                   ', type_lambda = "all") for deterministic lambdas\n',
                   'from each iteration.',
                   sep = "")

  }

  if(check_conv && ! is_conv_to_asymptotic(x)) {

    message(
      paste(mod_nm,
            ' has has not converged to asymptotic dynamics!',
            sep = "")
    )
  }

  msg <- c(msg, l_msg)



  cat(msg)

  invisible(x)

}

#' @rdname print_star
#' @export

print.general_di_stoch_kern_ipm <- function(x,
                                            comp_lambda = TRUE,
                                            type_lambda = 'stochastic',
                                            sig_digits  = 3,
                                            ...) {

  pretty_cls <- .pretty_class(class(x)[1])

  mod_nm     <- deparse(substitute(x))

  msg <- paste0('A ',
                pretty_cls,
                ' IPM with ',
                length(x$sub_kernels),
                ' sub-kernel(s) and ',
                length(x$pop_state),
                ' population vectors defined.',
                sep = "")

  if(comp_lambda) {


    all_lams <- lambda(x, type_lambda = type_lambda)

    l_msg  <- paste0('\nStochastic lambda for ',
                     mod_nm,
                     ' = ',
                     round(all_lams, sig_digits),
                     sep = "" )

    det_lam_msg <- paste('\nCall lambda(',
                         mod_nm,
                         ', type_lambda = "all") for deterministic lambdas\n',
                         'from each iteration.',
                         sep = "")
    msg <- c(msg, l_msg, det_lam_msg)

  }

  # No check for convergence here - there may quite reasonably be no
  # convergence for

  cat(msg)
  invisible(x)
}



#' @rdname print_star
#' @export

print.general_di_stoch_param_ipm <- function(x,
                                             comp_lambda = TRUE,
                                             type_lambda = 'stochastic',
                                             sig_digits  = 3,
                                             ...) {

  pretty_cls <- .pretty_class(class(x)[1])
  mod_nm <- deparse(substitute(x))

  msg <- paste0(
    'A ',
    pretty_cls,
    ' IPM with ',
    length(x$sub_kernels),
    ' sub-kernel(s), ',
    length(x$pop_state),
    ' population vectors, and ',
    dim(x$env_seq)[2],
    ' environmental parameters defined.',
    sep = ""
  )

  env_nms <- vapply(dimnames(x$env_seq)[[2]],
                    FUN = function(x) strsplit(x, '\\.')[[1]][2],
                    FUN.VALUE = character(1L))
  msg <- c(msg,
           paste(
             '\nEnvironmental variables are: ',
             paste(
               paste(env_nms, collapse = ', ')
             )
           ))

  if(comp_lambda) {

    all_lams <- lambda(x, type_lambda = type_lambda)

    l_msg  <- paste0('\nStochastic lambda for ',
                     mod_nm,
                     ' = ',
                     round(all_lams, sig_digits),
                     sep = "" )

    det_lam_msg <- paste('\nCall lambda(',
                         mod_nm,
                         ', type_lambda = "all") for deterministic lambdas\n',
                         'from each iteration.',
                         sep = "")

    msg <- c(msg, l_msg, det_lam_msg)

  }

  # No check for convergence here - there may quite reasonably be no
  # convergence for

  cat(msg)
  invisible(x)
}

#' @export

print.ipmr_vital_rate_exprs <- function(x, ...) {

  proto <- attr(x, "proto")

  if(any(proto$has_hier_effs) | .has_age(proto)) {

    out <- .pretty_print_hier_effs(x, proto)

  } else {

    out <- lapply(x, rlang::expr_text)

  }

  attr(out, "proto") <- NULL

  out <- paste(names(out), out, sep = ": ") %>%
    paste(., collapse = "\n")

  class(out) <- "character"

  cat(out)
  invisible(x)

}

.pretty_print_hier_effs <- function(x, proto) {

  nm_hier_effs <- proto$levels_hier_effs %>%
    .flatten_to_depth(1L) %>%
    names() %>%
    unique()

  if(.has_age(proto)) {

    ages <- proto$levels_ages %>%
      .flatten_to_depth(1L) %>%
      names() %>%
      unique()

    nm_hier_effs <- c(nm_hier_effs, ages)

  }

  nm_hier_effs <- nm_hier_effs[!is.na(nm_hier_effs)]

  for(i in seq_along(nm_hier_effs)) {

    to_print <- paste("<", nm_hier_effs[i], ">", sep = "")

    x <- lapply(x,
                function(y, nm, to_print) {
                  y <- rlang::expr_text(y)
                  gsub(nm, to_print, y)

                },
                nm       = nm_hier_effs[i],
                to_print = to_print)

    names(x) <- gsub(nm_hier_effs[i], to_print, names(x))

  }

  return(x)
}

#' @export

print.ipmr_kernel_exprs <- function(x, ...) {

  proto <- attr(x, "proto")

  if(any(proto$has_hier_effs) | .has_age(proto)) {

    out <- .pretty_print_hier_effs(x, proto)

  } else {

    out <- lapply(x, rlang::expr_text)

  }

  attr(out, "proto") <- NULL

  out <- paste(names(out), out, sep = ": ") %>%
    paste(., collapse = "\n")

  class(out) <- "character"

  cat(out)
  invisible(x)

}

#' @export

print.ipmr_parameters <- function(x, ...) {

  out <- as.numeric(x)
  names(out) <- names(x)

  print(out)
  invisible(x)

}

#' @export
# Will need to update w/ other names for other integration rules as those
# become available

print.ipmr_domains <- function(x, ...) {

  # Strips class and proto attributes

  out <- lapply(x, function(y) {
    y
  })

  out <- paste(names(out), out, sep = ": ") %>%
    paste(., collapse = "\n")

  out <- gsub("c\\(", "", out)
  out <- gsub("\\)", "", out)

  cat(out)
  invisible(x)

}

#' @export
print.ipmr_pop_state <- function(x, ...) {

  out <- lapply(x, function(y) y)

  if(out[[1]] == "No population state defined.") {

    cat(out[[1]])

  } else {

    if(rlang::is_named(x) & !grepl("pop_state_", names(out))) {
      message("Names of population states don't appear to be prefixed with 'n_'.")
    }

    names(out) <- gsub("pop_state_", "n_", names(out))

    out <- paste(names(out), out, sep = ": ") %>%
      paste(., collapse = "\n")

    cat(out)

  }

  invisible(x)

}

# Lambda------------
#' @title Compute the per-capita growth rate for an IPM object
#' @rdname lambda
#'
#' @description Compute the per-capita growth rate for a given model. Can handle
#' stochastic and deterministic models, and has options for thinning.
#'
#' @param ipm An object returned by \code{make_ipm()}.
#' @param comp_method Either \code{"eigen"} or \code{"pop_size"}. \code{"eigen"}
#' is only possible for \code{"simple_*"} methods.
#' @param type_lambda Either \code{'all'}, \code{'last'},
#'  or \code{'stochastic'}. \code{'all'}
#' returns a vector of lambda values for each time step of the simulation (equal
#' in length to the \code{iterations} argument of \code{make_ipm()}).
#' \code{'last'} returns the lambda value for the final timestep.
#' \code{'stochastic'} returns a single value which is the geometric mean
#' of log-per-capita growth rate from each time step.
#' @param ... other arguments passed to methods.
#' @param burn_in The proportion of iterations to discard. Default is 0.1
#' (i.e. first 10\% of iterations in the simulation).
#'
#' @return An array. Rows correspond to time steps, and columns correspond
#' to hierarchical effects (if any).
#'
#' @details There are two possible methods for computing \code{lambda} and these
#' are controlled by the \code{comp_method} argument. Possible values are
#' \code{"eigen"} and \code{"pop_size"}. The first computes the dominant
#' eigenvalue of all entries in \code{ipm$iterators} slot of the \code{*_ipm}
#' object. Since iteration kernels aren't generated for \code{general_*} methods,
#' \code{"pop_size"} is the only possible option for those objects.
#'
#' @export

lambda <- function(ipm, ...) {
  UseMethod('lambda')
}

#' @rdname lambda
#' @export

lambda.simple_di_det_ipm <- function(ipm,
                                     comp_method = "pop_size",
                                     type_lambda = 'last',
                                     ...) {

  # Need to insert arg checking - probably should make it it's own s3 generic.

  .check_lambda_args(ipm, comp_method, type_lambda)

  all_lams <- switch(type_lambda,
                     "all"  = TRUE,
                     'last' = FALSE)

  switch(comp_method,
         'eigen'    = .lambda_eigen(ipm),
         'pop_size' = .lambda_pop_size(ipm, all_lambdas = all_lams))

}

#' @rdname lambda
#' @export

lambda.simple_di_stoch_kern_ipm <- function(ipm,
                                            comp_method = "pop_size",
                                            type_lambda = 'stochastic',
                                            burn_in     = 0.1,
                                            ...) {


  .check_lambda_args(ipm, comp_method, type_lambda)

  all_lams <- switch(type_lambda,
                     "all"        = TRUE,
                     "stochastic" = TRUE,
                     'last'       = FALSE)

  temp <- .lambda_pop_size(ipm, all_lambdas = all_lams)

  if(all_lams) {
    burn_ind <- seq_len(round(length(temp) * burn_in))
  }

  return(
    switch(type_lambda,
           'all'        = temp,
           'last'       = temp,
           'stochastic' = .thin_stoch_lambda(temp, burn_ind))
  )
}

#' @rdname lambda
#' @export

lambda.simple_di_stoch_param_ipm <- function(ipm,
                                             comp_method = "pop_size",
                                             type_lambda = 'stochastic',
                                             burn_in     = 0.1,
                                             ...) {

  .check_lambda_args(ipm, comp_method, type_lambda)

  all_lams <- switch(type_lambda,
                     "all"        = TRUE,
                     "stochastic" = TRUE,
                     'last'       = FALSE)

  temp <- .lambda_pop_size(ipm, all_lambdas = all_lams)

  if(all_lams) {
    burn_ind <- seq_len(round(length(temp) * burn_in))
  }

  return(
    switch(type_lambda,
           'all'        = temp,
           'last'       = temp,
           'stochastic' = .thin_stoch_lambda(temp, burn_ind))
  )

}

#' @rdname lambda
#' @importFrom methods hasArg
#' @export

lambda.general_di_det_ipm <- function(ipm,
                                      comp_method = "pop_size",
                                      type_lambda = 'last',
                                      ...) {

  .check_lambda_args(ipm, comp_method, type_lambda)

  all_lams <- switch(type_lambda,
                     'all'  = TRUE,
                     'last' = FALSE,
                     'stochastic' = stop("Cannot compute stochastic lambda for deterministic IPM",
                                         call. = FALSE))

  out <- .lambda_pop_size(ipm, all_lambdas = all_lams)

  return(out)
}


#' @rdname lambda
#' @export

lambda.general_di_stoch_kern_ipm <- function(ipm,
                                             ...,
                                             comp_method = 'pop_size',
                                             type_lambda = 'stochastic',
                                             burn_in     = 0.1) {

  .check_lambda_args(ipm, comp_method, type_lambda)

  all_lams <- switch(type_lambda,
                     'all'        = TRUE,
                     'stochastic' = TRUE,
                     'last'       = FALSE)

  temp <- .lambda_pop_size(ipm, all_lambdas = all_lams)

  if(all_lams) {
    burn_ind <- seq_len(round(length(temp) * burn_in))
  }

  return(
    switch(type_lambda,
           'all'        = temp,
           'last'       = temp,
           'stochastic' = .thin_stoch_lambda(temp, burn_ind))
  )
}

#' @rdname lambda
#' @export

lambda.general_di_stoch_param_ipm <- function(ipm,
                                              ...,
                                              comp_method = 'pop_size',
                                              type_lambda = 'stochastic',
                                              burn_in     = 0.1) {

  .check_lambda_args(ipm, comp_method, type_lambda)

  all_lams <- switch(type_lambda,
                     'all'        = TRUE,
                     'stochastic' = TRUE,
                     'last'       = FALSE)

  temp <- .lambda_pop_size(ipm, all_lambdas = all_lams)

  if(all_lams) {
    burn_ind <- seq_len(round(length(temp) * burn_in))
  }

  return(
    switch(type_lambda,
           'all'        = temp,
           'last'       = temp,
           'stochastic' = .thin_stoch_lambda(temp, burn_ind))
  )
}

# Plot ---------------
# Authors - whoever wrote the code from the IPM book

#' @title Plot a matrix or an *_ipm object
#' @rdname plot_star
#'
#' @param x,y Either the values of the meshpoints or \code{NULL}. If \code{NULL},
#' then a sequence is generated so that meshpoints are given sequential bin numbers.
#' @param A,ipm A matrix or a result from \code{make_ipm}, or \code{NULL} if \code{x}
#' is specified as the matrix or IPM object.
#' @param col A vector of colors to use for plotting
#' @param bw A logical indicating whether to use a greyscale palette for plotting
#' @param do_contour A logical indicating whether or not draw contour lines
#' on the plot
#' @param do_legend A logical indicating whether to draw a legend for the plot
#' @param ... further arguments passed to legend
#'
#' @return \code{A} or \code{ipm} invisibly
#'
#' @details \code{plot.ipmr_matrix} is intended for internal use only, and it
#' is usually safer to use \code{plot.*_ipm} methods for visualizing kernels.
#' If a K kernel is overwhelmed by information in say, a fecundity sub-kernel,
#' use the \code{exponent} argument in \code{plot.*_ipm} to make it more visually
#' appealing.
#'
#' @importFrom grDevices grey rainbow
#' @importFrom graphics abline axis contour image layout par
#' @export


plot.ipmr_matrix <- function(x = NULL, y = NULL,
                             A,
                             col = grDevices::rainbow(100, start=0.67, end=0),
                             bw = FALSE,
                             do_contour = FALSE,
                             do_legend = FALSE,
                             ...) {

  old_par <- graphics::par('mar')
  on.exit(par(old_par))

  if(do_legend) graphics::layout(mat = cbind(matrix(1, 5, 5), rep(2, 5)))

  graphics::par(mar = c(6, 5, 3, 2))

  if(is.null(x)) x <- seq_len(ncol(A))
  if(is.null(y)) y <- seq_len(nrow(A))

  nx = length(x)
  ny = length(y)
  x1 = c(1.5 * x[1] - 0.5 * x[2], 1.5 * x[nx] - 0.5 * x[nx - 1])
  y1 = c(1.5 * y[1] - 0.5 * y[2], 1.5 * y[ny] - 0.5 * y[ny - 1])

  if(bw) col = grDevices::grey( (200:50) / 200 )

  graphics::image(list(x = x,
                       y = y,
                       z = t(A)),
                  xlim     = x1,
                  ylim     = rev(y1),
                  col      = col,
                  cex.axis = 1.5,
                  cex.lab  = 1.5,
                  bty      = "u",
                  xlab     = 'T',
                  ylab     = 'T + 1',
                  ...)

  graphics::abline(v = range(x1))
  graphics::abline(h = range(y1))

  if(do_contour) graphics::contour(x,
                                   y,
                                   t(A),
                                   nlevels = 5,
                                   labcex  = 1.2,
                                   add     = TRUE)

  if(do_legend) {
    l.y = seq(min(A), max(A),length = 100)
    par(mar = c(6, 2, 3, 1))
    graphics::image(list(x = 1:2,
                         y = l.y,
                         z = rbind(l.y, l.y)),
                    col  = col,
                    bty  = "o",
                    xaxt = "n",
                    yaxt = "n")
    graphics::axis(side     = 2,
                   cex.axis = 1.5,
                   at       = pretty(seq(min(A), max(A), length=10)))
  }

  invisible(A)
}

#' @rdname plot_star
#' @param sub_kernels A logical - also plot the sub-kernels?
#' @param exponent The exponent to raise each kernel to. Setting this to a low
#' number can help visualize kernels that are overwhelmed by a few very large numbers.
#'
#' @export

plot.simple_di_det_ipm <- function(x = NULL, y = NULL,
                                   ipm = NULL,
                                   sub_kernels = FALSE,
                                   col = rainbow(100, start=0.67, end=0),
                                   bw = FALSE,
                                   do_contour = FALSE,
                                   do_legend = FALSE,
                                   exponent = 1,
                                   ...) {

  # This is used so that users can just say plot(my_model) instead of
  # plot(ipm = my_model). ipmr_matrix expects x and y to both be NULL

  if(!is.null(x) && is.null(ipm)){
    ipm <- x
    x   <- NULL
  }

  old_par <- par('mar')
  on.exit(par(old_par))

  dots <- list(...)

  if(sub_kernels) {

    plot_list <- purrr::splice(ipm$iterators, ipm$sub_kernels)

  } else {

    plot_list <- ipm$iterators

  }

  lapply(plot_list, function(ipm) plot.ipmr_matrix(x = x,
                                                   y = y,
                                                   A = ipm ^ exponent,
                                                   col = col,
                                                   bw = bw,
                                                   do_contour = do_contour,
                                                   do_legend = do_legend,
                                                   dots))

  invisible(ipm)
}

#' @rdname plot_star
#' @export

plot.simple_di_stoch_param_ipm <- function(x = NULL, y = NULL,
                                           ipm = NULL,
                                           sub_kernels = FALSE,
                                           col = rainbow(100, start=0.67, end=0),
                                           bw = FALSE,
                                           do_contour = FALSE,
                                           do_legend = FALSE,
                                           exponent = 1,
                                           ...) {

  # This is used so that users can just say plot(my_model) instead of
  # plot(ipm = my_model). ipmr_matrix expects x and y to both be NULL

  if(!is.null(x) && is.null(ipm)){
    ipm <- x
    x   <- NULL
  }

  old_par <- par('mar')
  on.exit(par(old_par))

  dots <- list(...)

  if(sub_kernels) {

    plot_list <- purrr::splice(ipm$iterators, ipm$sub_kernels)

  } else {

    plot_list <- ipm$iterators

  }

  lapply(plot_list, function(ipm) plot.ipmr_matrix(x = x,
                                                   y = y,
                                                   A = ipm ^ exponent,
                                                   col = col,
                                                   bw = bw,
                                                   do_contour = do_contour,
                                                   do_legend = do_legend,
                                                   dots))

  invisible(ipm)
}

#' @rdname plot_star
#' @export

plot.simple_di_stoch_kern_ipm <- function(x = NULL, y = NULL,
                                          ipm = NULL,
                                          sub_kernels = FALSE,
                                          col = rainbow(100, start=0.67, end=0),
                                          bw = FALSE,
                                          do_contour = FALSE,
                                          do_legend = FALSE,
                                          exponent = 1,
                                          ...) {

  # This is used so that users can just say plot(my_model) instead of
  # plot(ipm = my_model). ipmr_matrix expects x and y to both be NULL

  if(!is.null(x) && is.null(ipm)){
    ipm <- x
    x   <- NULL
  }

  old_par <- par('mar')
  on.exit(par(old_par))

  dots <- list(...)

  if(sub_kernels) {

    plot_list <- purrr::splice(ipm$iterators, ipm$sub_kernels)

  } else {

    plot_list <- ipm$iterators

  }

  lapply(plot_list, function(ipm) plot.ipmr_matrix(x = x,
                                                   y = y,
                                                   A = ipm ^ exponent,
                                                   col = col,
                                                   bw = bw,
                                                   do_contour = do_contour,
                                                   do_legend = do_legend,
                                                   dots))

  invisible(ipm)
}

#' @rdname plot_star
#' @inheritParams format_mega_matrix
#'
#' @export

plot.general_di_det_ipm <- function(x = NULL, y = NULL,
                                    ipm = NULL,
                                    mega_mat = NA_character_,
                                    col = NA_character_,
                                    bw = FALSE,
                                    do_contour = FALSE,
                                    do_legend = FALSE,
                                    exponent = 1,
                                    ...) {

  # This is used so that users can just say plot(my_model) instead of
  # plot(ipm = my_model). ipmr_matrix expects x and y to both be NULL

  if(!is.null(x) && is.null(ipm)){
    ipm <- x
    x   <- NULL
  }

  if(is.na(mega_mat)) {

    stop("Plotting general IPMs requires building a 'mega_mat'.\n",
         "Please specify an expression for the 'mega_mat' argument.")
  }

  if(is.na(col)) {
    col <- rainbow(100,
                   start = 0.67,
                   end = 0,
                   alpha = NULL)
  }

  old_par <- par('mar')
  on.exit(par(old_par))

  dots <- list(...)

  mega_mat <- rlang::enquo(mega_mat)

  plot_list <- list(format_mega_matrix(ipm, mega_mat = !! mega_mat))

  lapply(plot_list, function(ipm) plot.ipmr_matrix(x = x,
                                                   y = y,
                                                   A = ipm ^ exponent,
                                                   col = col,
                                                   bw = bw,
                                                   do_contour = do_contour,
                                                   do_legend = do_legend,
                                                   dots))

  invisible(ipm)

}

# right_ev ----------------

#' @rdname eigenvectors
#'
#' @title Compute the standardized left and right eigenvectors via iteration
#'
#' @param ipm Output from \code{make_ipm()}.
#' @param ... other arguments passed to methods
#'
#' @return A list of named numeric vector(s) corresponding to the stable trait distribution
#' function (\code{right_ev}) or the reproductive values for each trait (\code{left_ev}).
#'
#' @export

right_ev <- function(ipm, ...) {

  UseMethod('right_ev')

}

#' @rdname eigenvectors
#' @param n_iterations The number of times to iterate the model to reach
#' convergence. Default is 100.
#'
#' @export

right_ev.simple_di_det_ipm <- function(ipm,
                                       n_iterations = 100,
                                       ...) {

  mod_nm <- deparse(substitute(ipm))

  # Identify state variable name

  pop_nm <- .get_pop_nm_simple(ipm)

  # if it's already been iterated to convergence, we don't have much work to do.

  if(.already_iterated(ipm)) {

    # get index for population vector of final iteration
    final_it <- dim(ipm$pop_state[[1]])[2]

    if(is_conv_to_asymptotic(ipm)) {

      out    <- ipm$pop_state[[1]][ , final_it]
      out_nm <- paste(pop_nm, 'w', sep = "_")

    } else {

      # If it's iterated, but not converged, we can at least use the final
      # iteration which is hopefully closer to convergence than the user-defined
      # initial population vector

      message(
        paste(
          "'",
          mod_nm,
          "'",
          ' did not converge to asymptotic dynamics after ',
          final_it,
          ' iterations.\n',
          'Will re-iterate the model ',
          n_iterations,
          ' times and check for convergence.',
          sep = ""
        )
      )

      init_pop_vec <- ipm$pop_state[[1]][ , final_it]

      pop_nm       <- ipm$proto_ipm$state_var %>%
        unlist() %>%
        unique() %>%
        .[1] %>%
        paste('n_', ., '_t', sep = "")

      pop_nm <- rlang::ensym(pop_nm)

      test_conv <- ipm$proto_ipm %>%
        define_pop_state(!! pop_nm := init_pop_vec) %>%
        make_ipm(iterate    = TRUE,
                 iterations = n_iterations)

      if(is_conv_to_asymptotic(test_conv)) {

        final_it <- dim(test_conv$pop_state[[1]])[2]

        out      <- test_conv$pop_state[[1]][ , final_it]

        out_nm   <- paste(pop_nm, 'w', sep = "_")

        message('model is now converged :)')

      } else {
        warning(
          paste(
            "'",
            mod_nm,
            "'",
            ' did not converge after ',
            final_it + n_iterations,
            ' iterations. Returning NA, please try again with more iterations.',
            sep = ""
          )
        )

        return(NA_real_)

      }
    }

  } else {

    message(
      paste(
        "'",
        mod_nm,
        "'",
        ' has not been iterated yet. ',
        'Generating a population vector using runif() and\niterating the model ',
        n_iterations,
        ' times to check for convergence to asymptotic dynamics',
        sep = ""
      )
    )

    # A model that hasn't been iterated yet. The part below only applies to
    # simple density independent models that might not have an expression
    # for iterating the model defined. For simple models, this always
    # n_sv_t_1 = k %*% n_sv_t, so we can define that internally without
    # any trouble. Other model classes will error during make_ipm() if they
    # don't have a pop_state, so this next part will look quite different for
    # them. Nearly every other model type also requires a minimum of 1 iteration,
    # so we can also be confident that it'll never reach this point for the
    # vast majority of cases.

    # Create variables for internal usage

    pop_states  <- paste('n', pop_nm, c('t', 't_1'), sep = '_')

    k_nm        <- names(ipm$iterators)

    # Construct the call described above to relate pop_state_t_1 to K
    # and pop_state_t

    text_call   <- paste(k_nm, ' %*% ', pop_states[1])

    to_add      <- rlang::list2(!! pop_states[2] := text_call)

    # Insert into proto

    proto_ind   <- which(ipm$proto_ipm$kernel_id == k_nm)

    ipm$proto_ipm$params[[proto_ind]]$formula <- purrr::splice(
      ipm$proto_ipm$params[[proto_ind]]$formula,
      to_add
    )

    # Final step is to generate an initial pop_state. This is always just
    # vector drawn from a random uniform distribution

    len_pop_state <- dim(ipm$iterators[[1]])[1]
    init_pop      <- stats::runif(len_pop_state)

    test_conv     <- ipm$proto_ipm %>%
      define_pop_state(!! pop_states[1] := init_pop) %>%
      make_ipm(iterate = TRUE,
               iterations = n_iterations)

    if(is_conv_to_asymptotic(test_conv)) {

      out    <- test_conv$pop_state[[1]][ , (n_iterations + 1)]
      out_nm <- paste(pop_nm, 'w', sep = "_")

    } else {

      warning(
        paste(
          "'",
          mod_nm,
          "'",
          ' did not converge after ',
          n_iterations,
          ' iterations. Returning NA, please try again with more iterations.',
          sep = ""
        )
      )

      return(NA_real_)
    }

  }

  # Stuff into a list and standardize

  out <- rlang::list2(!! out_nm := (out / sum(out)))

  return(out)

}

#' @rdname eigenvectors
#' @export

right_ev.general_di_det_ipm <- function(ipm,
                                        n_iterations = 100,
                                        mega_mat     = NULL,
                                        mega_vec     = NULL,
                                        keep_mega    = FALSE,
                                        ...) {

  mega_mat  <- rlang::enquo(mega_mat)
  mega_vec  <- rlang::enquo(mega_vec)

  mod_nm    <- deparse(substitute(ipm))

  final_it  <- dim(ipm$pop_state[[1]])[2]

  if(is_conv_to_asymptotic(ipm)) {

    out <- .extract_conv_ev_general(ipm$pop_state)

  } else {

    # Not yet converged, so re-iterate the model and see what we get

    message(
      paste(
        "'",
        mod_nm,
        "'",
        ' did not converge to asymptotic dynamics after ',
        final_it,
        ' iterations.\n',
        'Will re-iterate the model ',
        n_iterations,
        ' times and check for convergence.',
        sep = ""
      )
    )

    use_pop_state       <- ipm$pop_state[names(ipm$pop_state) != 'lambda']

    init_pop_vec        <- lapply(use_pop_state,
                                  function(x, final_it) x[ , final_it],
                                  final_it = final_it)

    names(init_pop_vec) <- gsub('pop_state', 'n', names(init_pop_vec))

    test_conv           <- ipm$proto_ipm %>%
      define_pop_state(
        pop_vectors = init_pop_vec
      ) %>%
      make_ipm(iterate    = TRUE,
               iterations = n_iterations)

    if(is_conv_to_asymptotic(test_conv)) {

      out <- .extract_conv_ev_general(test_conv$pop_state)

    } else {

      warning(
        paste(
          "'",
          mod_nm,
          "'",
          ' did not converge after ',
          n_iterations,
          ' iterations. Returning NA, please try again with more iterations.',
          sep = ""
        )
      )

      return(NA_real_)

    }
  }

  # If the user wants the formatted mega-mat/vector (or we want it internally),
  # Then build those up. Otherwise, just substitute names and return the list
  # of vectors

  if(keep_mega) {

    if(!rlang::quo_is_null(mega_vec) && !rlang::quo_is_null(mega_mat)) {
      vec_out <- .make_mega_vec(mega_vec, out)
      mat_out <- format_mega_matrix(ipm, !! mega_mat)[[1]]
      out     <- list(mega_mat = mat_out, mega_vec = vec_out)

    } else {
      stop("Cannot set 'keep_mega = TRUE' and leave mega_vec and/or mega_mat as NULL.")
    }

  } else {

    names(out) <- gsub('n_', '', names(out))
    names(out) <- paste(names(out), 'w', sep = '_')

  }

  return(out)
}

# left_ev -----------------

#' @rdname eigenvectors
#' @export

left_ev <- function(ipm, ...) {

  UseMethod('left_ev')

}

#' @export
#' @rdname eigenvectors
#' @importFrom stats runif

left_ev.simple_di_det_ipm <- function(ipm, n_iterations = 100, ...) {

  mod_nm  <- deparse(substitute(ipm))

  # Identify state variable name

  pop_nm <- .get_pop_nm_simple(ipm)

  # If it's already iterated, then we need to wrap K with a t()
  if(.already_iterated(ipm)) {

    # Models that are already iterated may have additional things.
    # however, the id of the kernel should match the top level kernel in this
    # list, so we'll use that.

    k_nm  <- names(ipm$iterators)
    k_ind <- ifelse(length(k_nm) > 1,
                    which(k_nm %in% ipm$proto_ipm$kernel_id),
                    1)

    t_k   <- t(ipm$iterators[[k_ind]])

    # next, we pull out the initial population vector and set up
    # something to hold the population state while we iterate. Simple ipms don't
    # have a complicated pop_state structure, so this is straightforward.

    n_row          <- dim(ipm$pop_state[[1]])[1]
    temp_pop_state <- matrix(NA_real_,
                             nrow = n_row,
                             ncol = (n_iterations + 1))

    temp_pop_state[ , 1] <- ipm$pop_state[[1]][ , 1]

    for(i in seq_len(n_iterations)) {

      temp_pop_state[ , (i + 1)] <- t_k %*% temp_pop_state[ , i]

    }

    if(.is_conv_to_asymptotic(temp_pop_state)) {

      out <- temp_pop_state[ , (n_iterations + 1)]
      out_nm <- paste(pop_nm, 'v', sep = "_")

    } else {

      warning(
        paste(
          "'",
          mod_nm,
          "'",
          ' did not converge after ',
          n_iterations,
          ' iterations. Returning NA, please try again with more iterations.',
          sep = ""
        )
      )

      return(NA_real_)

    }

  } else {

    # If not iterated, we need to generate a population vector and t(k)
    message(
      paste(
        "'",
        mod_nm,
        "'",
        ' has not been iterated yet. ',
        'Generating a population vector using runif() and\niterating the model ',
        n_iterations,
        ' times to check for convergence to asymptotic dynamics',
        sep = ""
      )
    )

    t_k    <- t(ipm$iterators[[1]])

    # Create variables for internal usage

    k_nm        <- names(ipm$iterators)

    # Final step is to generate an initial pop_state. This is always just
    # vector drawn from a random uniform distribution

    len_pop_state <- dim(ipm$iterators[[1]])[1]

    temp_pop_state <- matrix(NA_real_,
                             nrow = len_pop_state,
                             ncol = (n_iterations + 1))

    temp_pop_state[ , 1] <- stats::runif(len_pop_state)

    for(i in seq_len(n_iterations)) {

      temp_pop_state[ , (i + 1)] <- t_k %*% temp_pop_state[ , i]

    }

    if(.is_conv_to_asymptotic(temp_pop_state)) {

      out    <- temp_pop_state[ , (n_iterations + 1)]
      out_nm <- paste(pop_nm, 'v', sep = "_")

    } else {

      warning(
        paste(
          "'",
          mod_nm,
          "'",
          ' did not converge after ',
          n_iterations,
          ' iterations. Returning NA, please try again with more iterations.',
          sep = ""
        )
      )

      return(NA_real_)
    }

  }

  out <- rlang::list2(!! out_nm := (out / sum(out)))

  return(out)
}

#' @rdname eigenvectors
#' @param mega_mat A vector of names and/or 0s specifying the relationship
#' between the kernels in the model. The names should correspond to kernel
#' names, with 0s corresponding to sparse areas of the mega-matrix. The names can
#' be either symbols or characters. These functions support suffix expansion as in
#' \code{define_k(ernel)}, so expressions don't need to be re-written for every
#' combination hierarchical effects. DO NOT supply arguments \code{nrow},
#' \code{ncol}, or \code{byrow}. These are set internally based on dimensions
#' of each sub-kernel.
#' @param mega_vec A vector of names specifying the format of the population
#' vector. The names can be either symbols or characters.
#' @param keep_mega A logical. TRUE returns a list with the mega matrix
#' and full eigenvector, FALSE returns a list with each state's contribution
#' to eigenvector. Default is FALSE, mostly for internal usage.
#'
#'
#' @details If the model has already been iterated, then these functions
#' will just extract population state of the final iteration and return
#' that in a named list. Each element of the list is a vector with length
#' \code{>= 1} and corresponds each state variable's portion of the eigenvector.
#'
#' Note that for \code{*_di_stoch_kern_ipm}'s, these generics will create a mean
#' matrix and then compute the left/right eigenvectors for that. For
#' \code{*_di_stoch_param_ipm}'s, it will compute the average environment kernel
#' (e.g. using the means from the \code{env_seq} slot of the IPM).
#'
#' \code{mega_mat} fits the pieces of the model together to create a
#' a transpose of the model to iterate with. Kernel names/0s should be supplied
#' in ROW MAJOR order (as in \code{byrow = TRUE}). \code{ipmr} supplies
#' a helper function for large models \code{format_mega_matrix}
#' to help with age-size models or ones with a lot of hierarchical effects.
#'
#' @examples
#'
#' data(gen_di_det_ex)
#'
#' # mega matrix is specifed as a row-major vector/matrix of symbols. You can also
#' # supply a character vector/matrix. mega_vec is done the same way.
#' # DO NOT supply arguments like nrow or ncol if using \code{matrix(...)}.
#' # The function expects to work these out on its own.
#'
#' ipmr_v <- left_ev(gen_di_det_ex,
#'                   mega_mat     = c(stay_discrete, go_discrete,
#'                                    leave_discrete, P),
#'                   mega_vec = c(b, ht),
#'                   n_iterations = 100)
#'
#' @export

left_ev.general_di_det_ipm <- function(ipm,
                                       mega_mat,
                                       mega_vec,
                                       n_iterations = 100,
                                       keep_mega = FALSE,
                                       ...) {

  # capture expressions and model names

  mega_mat     <- rlang::enquo(mega_mat)
  mega_vec     <- rlang::enquo(mega_vec)

  mod_nm       <- deparse(substitute(ipm))

  # extract kernels and make a mega_k. same for initial population vector
  sub_kernels  <- ipm$sub_kernels


  # Set everything up and iterate the transposed model

  mega_k       <- format_mega_matrix(ipm, !! mega_mat)[[1]]
  mega_pop     <- .make_mega_vec(mega_vec, ipm$pop_state)
  t_k          <- t(mega_k)

  pop_holder   <- matrix(NA_real_,
                         nrow = length(mega_pop),
                         ncol = (n_iterations + 1))

  # Insert t_0 into pop_holder, then iterate!

  pop_holder[ , 1] <- mega_pop

  for(i in seq_len(n_iterations)) {

    pop_holder[ , (i + 1)] <- t_k %*%  pop_holder[ , i]

  }

  # Convert output back to ipmr style list. I *really* hate this mega-matrix
  # method, but genuinely have no idea how else to solve this.

  if(.is_conv_to_asymptotic(pop_holder)) {

    # Convert *back* to standard ipmr pop_state format - so dumb. Get the final
    # iteration and convert it back to a pop_state list

    use_pop    <- pop_holder[ , dim(pop_holder)[2]]

    if(! keep_mega) {

      pop_holder <- .mega_vec_to_list(mega_vec,
                                      use_pop,
                                      ipm$pop_state[names(ipm$pop_state != 'lambda')])

      # I'm foolish and wrote .extract_conv_ev in a way that isn't compatible with
      # the pop_holder format here, so we perform the standardization by hand
      # before returning.

      pop_std <- Reduce('sum', unlist(pop_holder), init = 0)

      out     <- lapply(pop_holder,
                        function(x, pop_std) x / pop_std,
                        pop_std = pop_std)
    } else {

      out_vec <- use_pop / sum(use_pop)
      out_mat <- mega_k

      out <- list(mega_mat = out_mat, mega_vec = out_vec)

    }

  } else {

      warning(
        paste(
          "'",
          mod_nm,
          "'",
        ' did not converge after ',
        n_iterations,
        ' iterations. Returning NA, please try again with more iterations.',
        sep = ""
      )
    )

    return(NA_real_)

  }


  return(out)
}


#' @title Age X Size models
#' @rdname age_x_size
#'
#' @description Helper functions for implementing age x size models
#'
#' @param expr An expression to compute a sum for.
#' @param ... Used internally. Modifying this will likely cause models to break.
#' @param na.rm Not used and ignored.
#'
#' @return A modified expression.
#'
#' @details Note that these functions do not compute a single number, but rather
#' expand an expression to compute the trait distribution. This differs from the
#' usual behavior of the numeric methods for \code{sum, prod} as those
#' return a single number or \code{NA}.
#'
#' @examples
#'
#' \dontrun{
#'
#'define_k(
#'  proto_ipm          = some_age_size_proto_ipm,
#'  name               = "K",
#'  family             = "IPM",
#'  n_wt_0_t_1         = sum(F_age %*% n_wt_age_t),
#'  n_wt_age_t_1       = P_age_minus_1 %*% n_wt_age_minus_1_t,
#'  n_wt_max_age_t_1   = P_max_age %*% n_wt_max_age_t +
#'                       P_max_age_minus_1 %*% n_wt_max_age_minus_1_t,
#'  data_list          = param_list,
#'  states             = list (c("wt")),
#'  has_hier_effs      = FALSE,
#'  levels_ages        = list(age = c(0:20), max_age = 21),
#'  evict_cor          = FALSE
#')
#'
#'
#' }
#'
#' @export

sum.age_size_expr <- function(expr, ..., na.rm = NA){

  expr <- rlang::enexpr(expr)
  return(.all_ages(!! expr, fun = "+", ...))

}

#' @rdname age_x_size
#' @export

prod.age_size_expr <- function(expr, ..., na.rm = NA){

  expr <- rlang::enexpr(expr)
  return(.all_ages(!! expr, fun = "*", ...))

}

#' @noRd

.all_ages <- function(expr, fun, ...) {

  expr <- rlang::enexpr(expr)

  ages <- list(...) %>%
    unlist()

  new_text <- character(length(ages))

  for(i in seq_along(ages)) {

    new_text[i] <- gsub("age", ages[i], expr)

  }

  # add spaces around "fun" so output is more human readable. mostly for
  # my own debugging purposes.
  # NB: cosndier updating to use call_modify instead of literal
  # subbing + expansion on the expression

  out_text <- paste(new_text, collapse = paste(" ",
                                               as.character(fun),
                                               " ",
                                               sep = ""))

  return(out_text)

}



# diagnose ---------------

diagnose <- function(ipm, ...) {

  UseMethod('diagnose')

}


# update -------------------

#' @rdname update
#' @title Update expressions in an IPM
#'
#' @description Provides methods to readily update vital rate functions and kernel
#' formulae, as well as pass new parameter values and/or functions to the model.
#'
#' @param object A \code{proto_ipm} or output from \code{make_ipm}.
#' @param kernels The names of the kernels to modify.
#' @param formulae The new kernel formulae (if any). This should be a named expression
#' or list of expressions.
#' @param dots The new vital rate expressions (e.g. the \code{...} of \code{define_(k)ernel}).
#' This should be a named expression or list of expressions.
#' @param parameters A list of new parameters to include in the model. Names that
#' overlap with the current parameter list will overwrite the current values, while
#' names that are not in the current parameter list will be appended to it.
#' @param functions New user defined functions to include in the model.
#' @param ... not used

update.proto_ipm <- function(object,
                             kernels,
                             formulae,
                             dots,
                             parameters,
                             functions,
                             ...) {

  # Define ME!!!!

}
