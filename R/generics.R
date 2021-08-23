# Print----------------

#' @title Print proto_ipms or *_ipm objects
#' @rdname print_star
#'
#' @param x An object of class \code{proto_ipm} or produced by \code{make_ipm()}.
#' @param ... Ignored
#'
#' @details For printing \code{proto_ipm} objects, indices are wrapped in
#' \code{<index>} to assist with debugging. These are not carried into the model,
#' just a visual aid.
#'
#' @export

print.proto_ipm <- function(x, ...) {

  n_kerns      <- dim(x)[1]

  cls_switch   <- class(x)[1]

  pretty_class <- .pretty_class(cls_switch)


  msg          <- paste(ifelse(.uses_age(x),
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
  vrs          <- vital_rate_exprs(x)

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
    cat("\nThis may be because 'define_domains()' has not yet been called.\n",
        "If it has been called, then double check the model definition and",
        "'data_list'.\n")

  }

  cat("\n\nDomains for state variables:\n\n")

  print(domains(x))

  cat("\n\nPopulation states defined:\n\n")

  print(pop_state(x))

  cat("\n\nInternally generated model iteration procedure:\n\n")

  if(!all(is.na(x$domain))){

    k_row <- .init_iteration(x, TRUE, direction = "right")

    print(kernel_formulae(k_row))

  } else {

    cat("Not enough information to generate an iteration procedure.\n",
        "define_impl() and define_pop_state() must be called first to generate one.")

  }

  # Add more later -----------

  invisible(x)
}

#' @noRd

.all_params_in_mod <- function(proto) {

  kern_list <- list()

  for(row in seq_len(nrow(proto))) {

    use_proto <- proto[row, ]
    kern_nm   <- proto$kernel_id[row]

    vr_exprs <- vital_rate_exprs(use_proto)
    params   <- parameters(proto)

    # Age x size models and/or par_sets models need expansion before we
    # we can test whether all parameters are in the model

    if(.uses_age(proto) | any(proto$uses_par_sets)) {

      all_levs <- .flatten_to_depth(proto$par_set_indices[proto$uses_par_sets],
                                    1L) %>%
        .[!duplicated(names(.))]

      par_set_nms <- names(all_levs)

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

          par_set_nm <- names(lev_df)[j]

          for(k in seq_len(dim(lev_df)[1])) {

            # Replace suffixes in expression and name
            all_vr_exprs[[it]]      <- gsub(par_set_nm,
                                            lev_df[k, j],
                                            base_expr)

            names(all_vr_exprs)[it] <- gsub(par_set_nm,
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
    domains    <- names(domains(proto)) %>%
      vapply(function(x) paste(x, c("_1", "_2"), sep = ""),
             character(2L)) %>%
      as.vector()

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

.uses_age <- function(x) {

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
#' @param sig_digits The number of significant digits to round to if \code{
#' comp_lambda = TRUE}.
#' @param type_lambda Either \code{'all'} or \code{'stochastic'}. See
#' \code{\link{lambda}} for more details.
#' @param check_conv A logical: for deterministic models, check if population state
#' has converged to asymptotic dynamics? If \code{TRUE} and the model has not
#' converged, a message will be printed.
#' @param ... Ignored
#'
#' @return \code{x} invisibly.
#'
#' @export

print.simple_di_det_ipm <- function(x,
                                    comp_lambda = TRUE,
                                    type_lambda = 'last',
                                    sig_digits = 3,
                                    check_conv = TRUE,
                                    ...) {

  pretty_cls <- .pretty_class(class(x)[1])

  msg <- paste0('A ',
                pretty_cls,
                ' IPM with ',
                length(x$sub_kernels),
                ' sub-kernel(s) defined.', sep = "")

  if(comp_lambda) {

    lambdas <- lambda(x, type_lambda = type_lambda)

    l_msg  <- paste0('\nDeterministic lambda = ',
                     round(lambdas, sig_digits),
                     sep = "")

    msg    <- c(msg, l_msg)

    if(check_conv &&
       !all(is_conv_to_asymptotic(x))) {

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

print.simple_dd_det_ipm <- function(x,
                                    comp_lambda = TRUE,
                                    type_lambda = 'last',
                                    sig_digits = 3,
                                    ...) {

  mod_nm <- deparse(substitute(x))

  pretty_cls <- .pretty_class(class(x)[1])

  msg <- paste0('A ',
                pretty_cls,
                ' IPM with ',
                length(x$sub_kernels),
                ' sub-kernel(s) defined.', sep = "")

  if(comp_lambda) {

    lambdas <- lambda(x, type_lambda = type_lambda)

    l_msg  <- paste('\nLambda for the final time step of the model is: ',
                    round(lambdas, sig_digits),
                    sep = "")

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

print.simple_dd_stoch_kern_ipm <- print.simple_di_stoch_kern_ipm

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

print.simple_dd_stoch_param_ipm <- print.simple_di_stoch_param_ipm

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

  if(check_conv && !all(is_conv_to_asymptotic(x))) {

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

print.general_dd_det_ipm <- function(x,
                                     comp_lambda = TRUE,
                                     type_lambda = 'last',
                                     sig_digits  = 3,
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

    msg <- c(msg, l_msg)

  }




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
                length(x$pop_state[!grepl("lambda", names(x$pop_state))]),
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
print.general_dd_stoch_kern_ipm <- print.general_di_stoch_kern_ipm

#' @rdname print_star
#' @export

print.general_di_stoch_param_ipm <-  function(x,
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
    length(x$pop_state[!grepl("lambda", names(x$pop_state))]),
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

#' @rdname print_star
#' @export
print.general_dd_stoch_param_ipm <- print.general_di_stoch_param_ipm

#' @export

print.ipmr_vital_rate_exprs <- function(x, ...) {

  proto <- attr(x, "proto")

  if(any(proto$uses_par_sets) | .uses_age(proto)) {

    out <- .pretty_print_par_sets(x, proto)

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

.pretty_print_par_sets <- function(x, proto) {

  nm_par_sets <- proto$par_set_indices %>%
    .flatten_to_depth(1L) %>%
    names() %>%
    unique()

  if(.uses_age(proto)) {

    ages <- proto$age_indices %>%
      .flatten_to_depth(1L) %>%
      names() %>%
      unique()

    nm_par_sets <- c(nm_par_sets, ages)

  }

  nm_par_sets <- nm_par_sets[!is.na(nm_par_sets)]

  x <- lapply(x, rlang::expr_text)

  for(i in seq_along(nm_par_sets)) {

    to_print <- paste("<", nm_par_sets[i], ">", sep = "")

    x <- lapply(x,
                function(y, nm, to_print) {
                  gsub(nm, to_print, y)
                },
                nm       = nm_par_sets[i],
                to_print = to_print)

    names(x) <- gsub(nm_par_sets[i], to_print, names(x))

  }

  return(x)
}

#' @export

print.ipmr_kernel_exprs <- function(x, ...) {

  proto <- attr(x, "proto")

  if(any(proto$uses_par_sets) | .uses_age(proto)) {

    out <- .pretty_print_par_sets(x, proto)

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

    if(rlang::is_named(x) && !grepl("pop_state_", names(out))) {
      message("Names of population states don't appear to be prefixed with 'n_'.")
    }

    names(out) <- gsub("pop_state_", "n_", names(out))

    out <- paste(names(out), out, sep = ": ") %>%
      paste(., collapse = "\n")

    cat(out)


  }

  invisible(x)

}

#' @export
print.ipmr_mp_mesh <- function(x, ...) {

  uses <- Filter(Negate(is.null), x)

  msg <- vapply(uses, function(z) {

    if(length(z) > 1) {

      paste("Lower bound: ", range(z)[1], ", Upper bound: ", range(z)[2],
            ", # of Meshpoints: ", sqrt(length(z)),
            sep = "")

    } else {
      as.character(z)
    }

  },
  character(1L)) %>%
    paste(names(.), " - ", ., sep = "") %>%
    paste(., collapse = "\n")

  cat(msg)

  invisible(x)

}

#' @export

print.ipmr_matrix <- function(x, ...) {

  rng <- range(x)

  msg <- paste("Minimum value: ",
               round(rng[1], 5),
               ", maximum value: ",
               round(rng[2], 5), sep = "")

  pp_tst <- all(x >= 0)

  msg <- c(msg,
           paste("\nAll entries greater than or equal to 0: ",
                 pp_tst ,
                 sep = ""))

  cat("\n", msg, "\n")

  invisible(x)
}

#' @export
print.ipmr_vital_rate_funs <- function(x, ...) {

  if(length(x) == 1 && x == "No vital rate functions specified") {
    cat("No vital rates specified\n")
    invisible(x)
    return()
  }
  vr_nms <- names(x)

  cls    <- vapply(x,
                   function(x) .mat_cls(x),
                   character(1L))

  rngs <- vapply(x,
                 function(y) round(range(y), 4),
                 numeric(2L))

  out <- vapply(seq_along(cls),
                function(it, vr_nm, cls, rngs) {
                  paste(
                    paste(vr_nms[it],
                        " (not yet discretized)",
                        sep = ""),
                    cls[it],
                    sep = ": "
                  ) %>%
                    paste(" with minimum value: ", rngs[1, it],
                          " and maximum value: ", rngs[2, it],
                          sep = "")
                },
                vr_nm = vr_nms, cls = cls, rngs = rngs,
                FUN.VALUE = character(1L)) %>%
    paste(., collapse = "\n")

  out <- paste(out, '\n', sep = "")

    cat(out)
    invisible(x)
}

.pretty_dim <- function(x) {

  col <- ncol(x)
  row <- nrow(x)

  cls <- switch(class(x)[1],
                "CC" = "kernel",
                "DC" = "column vector",
                "CD" = "row vector",
                "DD" = "matrix")

  paste("A", row, "x", col, cls, sep = " ")

}

.mat_cls <- function(x) {

  .pretty_dim(x)

}

# Lambda------------
#' @title Compute the per-capita growth rate for an IPM object
#' @rdname lambda
#'
#' @description Compute the per-capita growth rate for a given model. Can handle
#' stochastic and deterministic models, and has the option to discard burn in for
#' stochastic models.
#'
#' @param ipm An object returned by \code{make_ipm()}.
#' @param type_lambda Either \code{'all'}, \code{'last'},
#'  or \code{'stochastic'}. \code{'all'}
#' returns a vector of lambda values for each time step of the simulation (equal
#' in length to the \code{iterations} argument of \code{make_ipm()}).
#' \code{'last'} returns the lambda value for the final timestep.
#' \code{'stochastic'} returns a single value, which by default is
#' \code{mean(log(lambda(ipm, type_lambda = "all")))}, with the proportion of
#' \code{burn_in} iterations removed from the beginning of the simulation. Set
#' \code{log} to \code{FASLE} to get the arithmetic mean for stochastic models.
#' @param ... other arguments passed to methods.
#' @param burn_in The proportion of iterations to discard. Default is 0.1
#' (i.e. first 10\% of iterations in the simulation).
#' @param log Return lambda on the log scale? This is \code{TRUE} by default for
#' stochastic models, and \code{FALSE} for deterministic models.
#'
#' @return When \code{type_lambda = "all"}, an array. Rows correspond to time
#' steps, and columns correspond to parameter sets (if any). For other types,
#' a numeric vector.
#'
#' @export

lambda <- function(ipm, ...) {
  UseMethod('lambda')
}

#' @rdname lambda
#' @export

lambda.simple_di_det_ipm <- function(ipm,
                                     type_lambda = 'last',
                                     log = FALSE,
                                     ...) {

  .check_lambda_args(ipm, type_lambda)

  all_lams <- switch(type_lambda,
                     "all"  = TRUE,
                     'last' = FALSE)

  out <- .lambda_pop_size(ipm, all_lambdas = all_lams)

  if(log) out <- log(out)

  return(out)

}

#' @rdname lambda
#' @export

lambda.simple_di_stoch_kern_ipm <- function(ipm,
                                            type_lambda = 'stochastic',
                                            burn_in     = 0.1,
                                            log         = TRUE,
                                            ...) {

  .check_lambda_args(ipm, type_lambda)

  all_lams <- switch(type_lambda,
                     "all"        = TRUE,
                     "stochastic" = TRUE,
                     'last'       = FALSE)

  temp <- .lambda_pop_size(ipm, all_lambdas = all_lams)

  if(all_lams) {
    burn_ind <- seq_len(round(length(temp) * burn_in))
  }

  out <- switch(type_lambda,
                'all'        = temp,
                'last'       = temp,
                'stochastic' = .thin_stoch_lambda(temp, burn_ind, log))

  if(type_lambda == "stochastic") {
    if(log) {
      message("log(lambda) is returned by default for stochastic models. Set ",
              "'log = FALSE' for lambda on linear scale.")
    }
  }

  return(out)

}

#' @rdname lambda
#' @export

lambda.simple_di_stoch_param_ipm <- function(ipm,
                                             type_lambda = 'stochastic',
                                             burn_in     = 0.1,
                                             log         = TRUE,
                                             ...) {

  .check_lambda_args(ipm, type_lambda)

  all_lams <- switch(type_lambda,
                     "all"        = TRUE,
                     "stochastic" = TRUE,
                     'last'       = FALSE)

  temp <- .lambda_pop_size(ipm, all_lambdas = all_lams)

  if(all_lams) {
    burn_ind <- seq_len(round(length(temp) * burn_in))
  }

  out <- switch(type_lambda,
                'all'        = temp,
                'last'       = temp,
                'stochastic' = .thin_stoch_lambda(temp, burn_ind, log))

  if(type_lambda == "stochastic") {
    if(log) {
      message("log(lambda) is returned by default for stochastic models. Set ",
              "'log = FALSE' for lambda on linear scale.")
    }
  }

  return(out)


}

#' @rdname lambda
#' @importFrom methods hasArg
#' @export

lambda.general_di_det_ipm <- function(ipm,
                                      type_lambda = 'last',
                                      log = FALSE,
                                      ...) {

  .check_lambda_args(ipm, type_lambda)

  all_lams <- switch(type_lambda,
                     'all'  = TRUE,
                     'last' = FALSE,
                     'stochastic' = stop("Cannot compute stochastic lambda for deterministic IPM",
                                         call. = FALSE))

  out <- .lambda_pop_size(ipm, all_lambdas = all_lams)

  if(log) out <- log(out)

  return(out)
}


#' @rdname lambda
#' @export

lambda.general_di_stoch_kern_ipm <- function(ipm,
                                             ...,
                                             type_lambda = 'stochastic',
                                             burn_in     = 0.1,
                                             log = TRUE) {

  .check_lambda_args(ipm, type_lambda)

  all_lams <- switch(type_lambda,
                     'all'        = TRUE,
                     'stochastic' = TRUE,
                     'last'       = FALSE)

  temp <- .lambda_pop_size(ipm, all_lambdas = all_lams)

  if(all_lams) {
    burn_ind <- seq_len(round(length(temp) * burn_in))
  }

  out <- switch(type_lambda,
                'all'        = temp,
                'last'       = temp,
                'stochastic' = .thin_stoch_lambda(temp, burn_ind, log))

  if(type_lambda == "stochastic") {
    if(log) {
      message("log(lambda) is returned by default for stochastic models. Set ",
              "'log = FALSE' for lambda on linear scale.")
    }
  }


  return(out)
}

#' @rdname lambda
#' @export

lambda.general_di_stoch_param_ipm <- function(ipm,
                                              ...,
                                              type_lambda = 'stochastic',
                                              burn_in     = 0.1,
                                              log = TRUE) {

  .check_lambda_args(ipm, type_lambda)

  all_lams <- switch(type_lambda,
                     'all'        = TRUE,
                     'stochastic' = TRUE,
                     'last'       = FALSE)

  temp <- .lambda_pop_size(ipm, all_lambdas = all_lams)

  if(all_lams) {
    burn_ind <- seq_len(round(length(temp) * burn_in))
  }

  out <- switch(type_lambda,
                'all'        = temp,
                'last'       = temp,
                'stochastic' = .thin_stoch_lambda(temp, burn_ind, log))

  if(type_lambda == "stochastic") {
    if(log) {
      message("log(lambda) is returned by default for stochastic models. Set ",
              "'log = FALSE' for lambda on linear scale.")
    }
  }

  return(out)

}

#' @rdname lambda
#' @export

lambda.simple_dd_det_ipm <- function(ipm,
                                     type_lambda = "all",
                                     ...,
                                     log = FALSE) {

  .check_lambda_args(ipm, type_lambda)

  all_lams <- switch(
    type_lambda,
    'all'  = TRUE,
    'last' = FALSE,
    'stochastic' = stop("Cannot compute stochastic lambda for deterministic IPM",
                        call. = FALSE))

  out <- .lambda_pop_size(ipm, all_lambdas = all_lams)

  if(log) out <- log(out)

  return(out)
}

#' @rdname lambda
#' @export

lambda.simple_dd_stoch_kern_ipm <- function(ipm,
                                            ...,
                                            type_lambda = 'stochastic',
                                            burn_in     = 0.1,
                                            log = TRUE) {

  .check_lambda_args(ipm, type_lambda)

  all_lams <- switch(type_lambda,
                     'all'        = TRUE,
                     'stochastic' = TRUE,
                     'last'       = FALSE)

  temp <- .lambda_pop_size(ipm, all_lambdas = all_lams)

  if(all_lams) {
    burn_ind <- seq_len(round(length(temp) * burn_in))
  }

  out <- switch(type_lambda,
                'all'        = temp,
                'last'       = temp,
                'stochastic' = .thin_stoch_lambda(temp, burn_ind, log))

  if(type_lambda == "stochastic") {

    if(log) {
      message("log(lambda) is returned by default for stochastic models. Set ",
              "'log = FALSE' for lambda on linear scale.")
    }
  }

  return(out)

}

#' @rdname lambda
#' @export

lambda.simple_dd_stoch_param_ipm <- function(ipm,
                                             ...,
                                             type_lambda = 'stochastic',
                                             burn_in     = 0.1,
                                             log = TRUE) {

  .check_lambda_args(ipm, type_lambda)

  all_lams <- switch(type_lambda,
                     'all'        = TRUE,
                     'stochastic' = TRUE,
                     'last'       = FALSE)

  temp <- .lambda_pop_size(ipm, all_lambdas = all_lams)

  if(all_lams) {
    burn_ind <- seq_len(round(length(temp) * burn_in))
  }

  out <- switch(type_lambda,
                'all'        = temp,
                'last'       = temp,
                'stochastic' = .thin_stoch_lambda(temp, burn_ind, log))

  if(type_lambda == "stochastic") {
    if(log) {
      message("log(lambda) is returned by default for stochastic models. Set ",
              "'log = FALSE' for lambda on linear scale.")
    }
  }

  return(out)
}

#' @rdname lambda
#' @export

lambda.general_dd_det_ipm <- function(ipm,
                                      type_lambda = 'last',
                                      ...,
                                      log = FALSE) {

  .check_lambda_args(ipm, type_lambda)

  all_lams <- switch(type_lambda,
                     'all'  = TRUE,
                     'last' = FALSE,
                     'stochastic' = stop("Cannot compute stochastic lambda for deterministic IPM",
                                         call. = FALSE))

  out <- .lambda_pop_size(ipm, all_lambdas = all_lams)

  if(log) out <- log(out)

  return(out)
}

#' @rdname lambda
#' @export

lambda.general_dd_stoch_kern_ipm <- function(ipm,
                                            ...,
                                            type_lambda = 'stochastic',
                                            burn_in     = 0.1,
                                            log = TRUE) {

  .check_lambda_args(ipm, type_lambda)

  all_lams <- switch(type_lambda,
                     'all'        = TRUE,
                     'stochastic' = TRUE,
                     'last'       = FALSE)

  temp <- .lambda_pop_size(ipm, all_lambdas = all_lams)

  if(all_lams) {
    burn_ind <- seq_len(round(length(temp) * burn_in))
  }

  out <- switch(type_lambda,
                'all'        = temp,
                'last'       = temp,
                'stochastic' = .thin_stoch_lambda(temp, burn_ind, log))

  if(type_lambda == "stochastic") {
    if(log) {
      message("log(lambda) is returned by default for stochastic models. Set ",
              "'log = FALSE' for lambda on linear scale.")
    }
  }

  return(out)

}

#' @rdname lambda
#' @export

lambda.general_dd_stoch_param_ipm <- function(ipm,
                                             ...,
                                             type_lambda = 'stochastic',
                                             burn_in     = 0.1,
                                             log = TRUE) {

  .check_lambda_args(ipm, type_lambda)

  all_lams <- switch(type_lambda,
                     'all'        = TRUE,
                     'stochastic' = TRUE,
                     'last'       = FALSE)

  temp <- .lambda_pop_size(ipm, all_lambdas = all_lams)

  if(all_lams) {
    burn_ind <- seq_len(round(length(temp) * burn_in))
  }

  out <- switch(type_lambda,
                'all'        = temp,
                'last'       = temp,
                'stochastic' = .thin_stoch_lambda(temp, burn_ind, log))

  if(type_lambda == "stochastic") {
    if(log) {
      message("log(lambda) is returned by default for stochastic models. Set ",
              "'log = FALSE' for lambda on linear scale.")
    }
  }

  return(out)
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
#' @param contour_cex A numeric specifying how large to make labels for the
#' contour lines.
#' @param do_legend A logical indicating whether to draw a legend for the plot
#' @param ... further arguments passed to legend
#'
#' @return \code{A} or \code{ipm} invisibly
#'
#' @details
#' If an IPM kernel is overwhelmed by information in say, a fecundity sub-kernel,
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
                             contour_cex = 1,
                             ...) {

  if(missing(A)) A <- x
  x <- seq_len(ncol(A))
  y <- seq_len(nrow(A))

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
                  bty      = "u",
                  ...)

  graphics::abline(v = range(x1))
  graphics::abline(h = range(y1))

  if(do_contour) graphics::contour(x,
                                   y,
                                   t(A),
                                   nlevels = 5,
                                   labcex  = contour_cex,
                                   add     = TRUE)

  if(do_legend) {
    l.y = seq(min(A), max(A),length = 100)
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
#' @param exponent The exponent to raise each kernel to. Setting this to a low
#' number can help visualize kernels that are overwhelmed by a few very large numbers.
#'
#' @export

plot.simple_di_det_ipm <- function(x = NULL, y = NULL,
                                   ipm = NULL,
                                   col = rainbow(100, start=0.67, end=0),
                                   bw = FALSE,
                                   do_contour = FALSE,
                                   do_legend = FALSE,
                                   exponent = 1,
                                   ...) {

  dots <- list(...)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(par(old_par))

  # This is used so that users can just say plot(my_model) instead of
  # plot(ipm = my_model). ipmr_matrix expects x and y to both be NULL

  if(!is.null(x) && is.null(ipm)){
    ipm <- x
    x   <- NULL
  }


  plot_list <- ipm$sub_kernels

  plt_seq   <- seq_along(plot_list)

  canvas_dims <- .ncol_nrow(plt_seq)

  if(do_legend) {

    graphics::layout(mat = cbind(matrix(1, 5, 5), rep(2, 5)))
    par(mar = c(6, 2, 3, 1))

  } else {


    graphics::par(mar = c(6, 5, 3, 2),
                  mfrow = c(canvas_dims$nrow,
                            canvas_dims$ncol))
  }


  for(i in seq_along(plot_list)){

    use_kern <- plot_list[[i]]
    nm       <- names(plot_list)[i]

    plot.ipmr_matrix(x = x,
                     y = y,
                     A = use_kern ^ exponent,
                     col = col,
                     bw = bw,
                     do_contour = do_contour,
                     do_legend = do_legend,
                     dots,
                     main = nm)
  }

  invisible(ipm)
}

#' @rdname plot_star
#' @export

plot.simple_di_stoch_param_ipm <- function(x = NULL, y = NULL,
                                           ipm = NULL,
                                           col = rainbow(100, start=0.67, end=0),
                                           bw = FALSE,
                                           do_contour = FALSE,
                                           do_legend = FALSE,
                                           exponent = 1,
                                           ...) {

  dots <- list(...)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(par(old_par))

  # This is used so that users can just say plot(my_model) instead of
  # plot(ipm = my_model). ipmr_matrix expects x and y to both be NULL

  if(!is.null(x) && is.null(ipm)){
    ipm <- x
    x   <- NULL
  }


  plot_list <- ipm$sub_kernels

  plt_seq   <- seq_along(plot_list)

  canvas_dims <- .ncol_nrow(plt_seq)

  if(do_legend) {

    graphics::layout(mat = cbind(matrix(1, 5, 5), rep(2, 5)))
    par(mar = c(6, 2, 3, 1))

  } else {


    graphics::par(mar = c(6, 5, 3, 2),
                  mfrow = c(canvas_dims$nrow,
                            canvas_dims$ncol))
  }


  for(i in seq_along(plot_list)){

    use_kern <- plot_list[[i]]
    nm       <- names(plot_list)[i]

    plot.ipmr_matrix(x = x,
                     y = y,
                     A = use_kern ^ exponent,
                     col = col,
                     bw = bw,
                     do_contour = do_contour,
                     do_legend = do_legend,
                     dots,
                     main = nm)
  }


  invisible(ipm)
}

#' @rdname plot_star
#' @export

plot.simple_di_stoch_kern_ipm <- function(x = NULL, y = NULL,
                                          ipm = NULL,
                                          col = rainbow(100, start=0.67, end=0),
                                          bw = FALSE,
                                          do_contour = FALSE,
                                          do_legend = FALSE,
                                          exponent = 1,
                                          ...) {

  dots <- list(...)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(par(old_par))

  # This is used so that users can just say plot(my_model) instead of
  # plot(ipm = my_model). ipmr_matrix expects x and y to both be NULL

  if(!is.null(x) && is.null(ipm)){
    ipm <- x
    x   <- NULL
  }


  plot_list   <- ipm$sub_kernels

  plt_seq     <- seq_along(plot_list)

  canvas_dims <- .ncol_nrow(plt_seq)

  if(do_legend) {

    graphics::layout(mat = cbind(matrix(1, 5, 5), rep(2, 5)))
    par(mar = c(6, 2, 3, 1))

  } else {


    graphics::par(mar = c(6, 5, 3, 2),
                  mfrow = c(canvas_dims$nrow,
                            canvas_dims$ncol))
  }


  for(i in seq_along(plot_list)){

    use_kern <- plot_list[[i]]
    nm       <- names(plot_list)[i]

    plot.ipmr_matrix(x = x,
                     y = y,
                     A = use_kern ^ exponent,
                     col = col,
                     bw = bw,
                     do_contour = do_contour,
                     do_legend = do_legend,
                     dots,
                     main = nm)
  }


  invisible(ipm)
}

#' @rdname plot_star
#' @inheritParams format_mega_kernel
#'
#' @export

plot.general_di_det_ipm <- function(x = NULL, y = NULL,
                                    ipm = NULL,
                                    mega_mat = NA_character_,
                                    col = rainbow(100,
                                                  start = 0.67,
                                                  end = 0),
                                    bw = FALSE,
                                    do_contour = FALSE,
                                    do_legend = FALSE,
                                    exponent = 1,
                                    ...) {
  dots <- list(...)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(par(old_par))

  # This is used so that users can just say plot(my_model) instead of
  # plot(ipm = my_model). ipmr_matrix expects x and y to both be NULL

  if(!is.null(x) && is.null(ipm)){
    ipm <- x
    x   <- NULL
  }

  mega_mat <- rlang::enquo(mega_mat)

  if(rlang::quo_text(mega_mat) == "NA") {

    stop("Plotting general IPMs requires building a 'mega_mat'.\n",
         "Please specify an expression for the 'mega_mat' argument.")
  }

  plot_list <- format_mega_kernel(ipm, mega_mat = !! mega_mat)
  plt_seq   <- seq_along(plot_list)
  canvas_dims <- .ncol_nrow(plt_seq)

  if(do_legend) {

    graphics::layout(mat = cbind(matrix(1, 5, 5), rep(2, 5)))
    par(mar = c(6, 2, 3, 1))

  } else {


    graphics::par(mar = c(6, 5, 3, 2),
                  mfrow = c(canvas_dims$nrow,
                            canvas_dims$ncol))
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

.ncol_nrow <- function(plt_seq) {

  out <- list(ncol = NA,
              nrow = NA)

  if(length(plt_seq) > 25) {

    out$ncol <- 5
    out$nrow <- 5

  } else {

    out$ncol <- out$nrow <- round(sqrt(length(plt_seq)))

  }

  return(out)

}

# right_ev ----------------

#' @rdname eigenvectors
#'
#' @title Compute the standardized left and right eigenvectors via iteration
#'
#' @param ipm Output from \code{make_ipm()}.
#' @param ... Other arguments passed to methods
#'
#' @return A list of named numeric vector(s) corresponding to the stable trait distribution
#' function (\code{right_ev}) or the reproductive values for each trait (\code{left_ev}).
#'
#' @export

right_ev <- function(ipm, ... ) {

  UseMethod('right_ev')

}

#' @rdname eigenvectors
#' @param iterations The number of times to iterate the model to reach
#' convergence. Default is 100.
#' @param tolerance Tolerance to evaluate convergence to asymptotic dynamics.
#'
#' @export

right_ev.simple_di_det_ipm <- function(ipm,
                                       iterations = 100,
                                       tolerance = 1e-10,
                                       ...) {

  mod_nm <- deparse(substitute(ipm))

  # if it's already been iterated to convergence, we don't have much work to do.

  if(.already_iterated(ipm)) {

    # get index for population vector of final iteration
    final_it <- dim(ipm$pop_state[[1]])[2]

    if(all(is_conv_to_asymptotic(ipm,
                                 tolerance = tolerance))) {

      out    <- .extract_conv_ev_simple(ipm$pop_state)

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
          iterations,
          ' times and check for convergence.',
          sep = ""
        )
      )

      init_pop_vec <- lapply(ipm$pop_state[!grepl("lambda", names(ipm$pop_state))],
                             function(x, final_it) {
                               x[ , final_it]
                             },
                             final_it = ncol(ipm$pop_state[[1]]))
#
#       pop_nm       <- ipm$proto_ipm$state_var %>%
#         unlist() %>%
#         unique() %>%
#         .[1] %>%
#         paste('n_', ., sep = "")
#
#       pop_nm <- rlang::ensym(pop_nm)

      test_conv <- ipm$proto_ipm %>%
        define_pop_state(!!! init_pop_vec) %>%
        make_ipm(iterate    = TRUE,
                 iterations = iterations)

      if(all(is_conv_to_asymptotic(test_conv,
                                   tolerance = tolerance))) {

        out    <- .extract_conv_ev_simple(test_conv$pop_state)

        message('model is now converged :)')

      } else {
        warning(
          paste(
            "'",
            mod_nm,
            "'",
            ' did not converge after ',
            final_it + iterations,
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
        iterations,
        ' times to check for convergence to asymptotic dynamics',
        sep = ""
      )
    )

    # A model that hasn't been iterated yet

    # Create variables for internal usage
    pop_nm <- .get_pop_nm_simple(ipm)

    pop_states  <- paste('n', pop_nm, sep = '_')

    # Final step is to generate an initial pop_state. This is always just
    # vector drawn from a random uniform distribution

    len_pop_state <- dim(ipm$sub_kernels[[1]])[1]
    init_pop      <- stats::runif(len_pop_state)

    # Drop _t for defining the initial population vector so define_pop_state
    # doesn't complain

    test_conv     <- ipm$proto_ipm %>%
      define_pop_state(!! pop_states[1] := init_pop) %>%
      make_ipm(iterate = TRUE,
               iterations = iterations,
               normalize_pop_size = TRUE)

    if(all(is_conv_to_asymptotic(test_conv,
                                 tolerance = tolerance))) {

      out    <- .extract_conv_ev_simple(test_conv$pop_state)

    } else {

      warning(
        paste(
          "'",
          mod_nm,
          "'",
          ' did not converge after ',
          iterations,
          ' iterations. Returning NA, please try again with more iterations.',
          sep = ""
        )
      )

      return(NA_real_)
    }

  }

  # Stuff into a list and standardize

  # out <- rlang::list2(!! out_nm := (out / sum(out)))

  out <- Filter(function(x) !any(is.na(x)), out)

  # Replaces the leading n_ with a trailing _v
  names(out) <- substr(names(out), 3, nchar(names(out)))
  names(out) <- paste(names(out), 'w', sep = '_')
  class(out) <- "ipmr_w"

  return(out)

}

#' @rdname eigenvectors
#' @param burn_in The proportion of early iterations to discard from the
#' stochastic simulation
#' @export

right_ev.simple_di_stoch_kern_ipm <- function(ipm,
                                              burn_in = 0.25,
                                              ...) {

  mod_nm <- deparse(substitute(ipm))

  # Identify state variable name

  pop_nm <- .get_pop_nm_simple(ipm)

  out <- ipm$pop_state[!grepl("lambda", names(ipm$pop_state))]

  burn_in_seq <- seq_len(ncol(out[[1]])) %>%
    .[seq.int(from = ceiling(ncol(out[[1]]) * burn_in),
              to   = ncol(out[[1]]),
              by   = 1)]

  out <- lapply(out,
                function(x, burn_seq) x[ , burn_seq],
                burn_seq = burn_in_seq)

  # Normalize just in case. easy for simple IPMs where we know there's only
  # one trait contributing to total population size. Tricker in general IPMs.

  out[[1]] <- out[[1]] / colSums(out[[1]])

  # Replaces the leading n_ with a trailing _v
  names(out) <- substr(names(out), 3, nchar(names(out)))
  names(out) <- paste(names(out), 'w', sep = '_')
  class(out) <- "ipmr_w"

  return(out)

}

#' @rdname eigenvectors
#' @export

right_ev.simple_di_stoch_param_ipm <- right_ev.simple_di_stoch_kern_ipm


#' @rdname eigenvectors
#' @export

right_ev.general_di_det_ipm <- function(ipm,
                                        iterations = 100,
                                        tolerance = 1e-10,
                                        ...) {

  mod_nm    <- deparse(substitute(ipm))

  final_it  <- dim(ipm$pop_state[[1]])[2]

  if(all(is_conv_to_asymptotic(ipm,
                           tolerance = tolerance))) {

    out <- .extract_conv_ev_general(ipm$pop_state, ipm$proto_ipm)

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
        iterations,
        ' times and check for convergence.',
        sep = ""
      )
    )

    use_pop_state       <- ipm$pop_state[!grepl("lambda", names(ipm$pop_state))]

    init_pop_vec        <- lapply(use_pop_state,
                                  function(x, final_it) x[ , final_it],
                                  final_it = final_it)

    # Remove NAs generated by the unsubstituted pop_state place holders
    init_pop_vec <- Filter(f = function(y) !any(is.na(y)), x = init_pop_vec)

    names(init_pop_vec) <- gsub('pop_state', 'n', names(init_pop_vec))

    test_conv           <- ipm$proto_ipm %>%
      define_pop_state(
        pop_vectors = init_pop_vec
      ) %>%
      make_ipm(iterate    = TRUE,
               iterations = iterations)

    if(all(is_conv_to_asymptotic(test_conv,
                                 tolerance = tolerance))) {

      out <- .extract_conv_ev_general(test_conv$pop_state, test_conv$proto_ipm)

    } else {

      warning(
        paste(
          "'",
          mod_nm,
          "'",
          ' did not converge after ',
          iterations,
          ' iterations. Returning NA, please try again with more iterations.',
          sep = ""
        )
      )

      return(NA_real_)

    }
  }

  out <- Filter(function(x) !any(is.na(x)), out)

  # Replaces the leading n_ with a trailing _v
  names(out) <- substr(names(out), 3, nchar(names(out)))
  names(out) <- paste(names(out), 'w', sep = '_')
  class(out) <- "ipmr_w"


  return(out)
}

#' @rdname eigenvectors
#' @export

right_ev.general_di_stoch_kern_ipm <- function(ipm,
                                               burn_in = 0.25,
                                               ...) {

  mod_nm   <- deparse(substitute(ipm))

  out      <- ipm$pop_state[!grepl("lambda", names(ipm$pop_state))]

  burn_in_seq <- seq_len(ncol(out[[1]])) %>%
    .[seq.int(from = ceiling(ncol(out[[1]]) * burn_in),
              to   = ncol(out[[1]]),
              by   = 1)]

  out <- lapply(out,
                function(x, burn_seq) x[ , burn_seq, drop = FALSE],
                burn_seq = burn_in_seq)

  # Normalize everything just in case

  pop_sizes <- lapply(out, colSums) %>%
    do.call(what = ".add", args = .)

  for(i in seq_along(out)) {

    for(j in seq_len(ncol(out[[1]]))) {

      out[[i]][ , j] <- out[[i]][ ,j] / pop_sizes[j]

    }
  }

  out <- Filter(function(x) !any(is.na(x)), out)

  # Replaces the leading n_ with a trailing _v
  names(out) <- substr(names(out), 3, nchar(names(out)))
  names(out) <- paste(names(out), 'w', sep = '_')
  class(out) <- "ipmr_w"


  return(out)
}

#' @rdname eigenvectors
#' @export

right_ev.general_di_stoch_param_ipm <- right_ev.general_di_stoch_kern_ipm

# left_ev -----------------

#' @rdname eigenvectors
#' @export

left_ev <- function(ipm, ...) {

  UseMethod('left_ev')

}

#' @export
#' @rdname eigenvectors
#' @importFrom stats runif

left_ev.simple_di_det_ipm <- function(ipm,
                                      iterations = 100,
                                      tolerance = 1e-10,
                                      ...) {

  mod_nm <- deparse(substitute(ipm))

  # Identify state variable name

  pop_nm <- .get_pop_nm_simple(ipm)

  if(any(ipm$proto_ipm$uses_par_sets)) {

    ps_nms <- .flatten_to_depth(ipm$proto_ipm$par_set_indices, 1L) %>%
      names() %>%
      unique()
    ps_nms <- ps_nms[ps_nms != "levels"]

    ps_suff <- paste(ps_nms, collapse = "_")

    pop_nm <- paste(pop_nm, ps_suff, sep = "_")
  }

  # Create variables for internal usage

  pop_states  <- paste('n', pop_nm, sep = '_')

  # Final step is to generate an initial pop_state. This is always just
  # vector drawn from a random uniform distribution

  len_pop_state <- dim(ipm$sub_kernels[[1]])[1]
  init_pop      <- stats::runif(len_pop_state)

  test_conv     <- ipm$proto_ipm %>%
    define_pop_state(!! pop_states[1] := init_pop) %>%
    make_ipm(iterate = TRUE,
             iterations = iterations,
             iteration_direction = "left",
             normalize_pop_size = TRUE)

  if(all(is_conv_to_asymptotic(test_conv,
                               tolerance = tolerance))) {

    out    <- .extract_conv_ev_simple(test_conv$pop_state)
    # out_nm <- paste(pop_nm, 'v', sep = "_")

  } else {

    warning(
      paste(
        "'",
        mod_nm,
        "'",
        ' did not converge after ',
        iterations,
        ' iterations. Returning NA, please try again with more iterations.',
        sep = ""
      )
    )

    return(NA_real_)
  }

  # Stuff into a list and standardize

  out        <- Filter(function(x) !any(is.na(x)), out)
  names(out) <- substr(names(out), 3, nchar(names(out)))
  names(out) <- paste(names(out), 'v', sep = '_')
  class(out) <- "ipmr_v"

  return(out)

}

#' @rdname eigenvectors
#' @param kernel_seq The sequece of parameter set indices used to select kernels
#' during the iteration procedure. If \code{NULL}, will use the sequence stored
#' in the \code{ipm} object. Should usually be left as \code{NULL}.
#' @export

left_ev.simple_di_stoch_kern_ipm <- function(ipm,
                                             iterations = 10000,
                                             burn_in    = 0.25,
                                             kernel_seq = NULL,
                                             ...) {

  mod_nm <- deparse(substitute(ipm))

  # Identify state variable name

  pop_nm <- .get_pop_nm_simple(ipm)

  # Create variables for internal usage

  pop_states  <- paste('n', pop_nm, c('t', 't_1'), sep = '_')

  # Final step is to generate an initial pop_state. This is always just
  # vector drawn from a random uniform distribution

  len_pop_state <- nrow(ipm$sub_kernels[[1]])
  init_pop      <- stats::runif(len_pop_state)
  if(is.null(kernel_seq)) kernel_seq <- ipm$env_seq

  # Drop _t for defining the initial population vector so define_pop_state
  # doesn't complain

  pop_states[1] <- gsub("_t$", "", pop_states[1])

  test_conv     <- ipm$proto_ipm %>%
    define_pop_state(!! pop_states[1] := init_pop) %>%
    make_ipm(iterate = TRUE,
             iterations = iterations,
             iteration_direction = "left",
             normalize_pop_size = TRUE,
             kernel_seq = kernel_seq)

    out    <- test_conv$pop_state
    out$lambda <- NULL

    burn_in_seq <- seq_len(ncol(out[[1]])) %>%
      .[seq.int(from = ceiling(ncol(out[[1]]) * burn_in),
                to   = ncol(out[[1]]),
                by   = 1)]

    out <- lapply(out,
                  function(x, burn_seq) x[ , burn_seq],
                  burn_seq = burn_in_seq)

    out_nm <- paste(pop_nm, 'v', sep = "_")

    names(out) <- out_nm
    class(out) <- "ipmr_v"

    return(out)
}

#' @rdname eigenvectors
#'
#' @section \strong{Deterministic eigenvectors}:
#'  For \code{right_ev}, if the model has already been iterated and has
#' converged to asymptotic dynamics, then it will just extract the final
#' population state and return that in a named list. Each element of the list
#' is a vector with length \code{>= 1} and corresponds each state variable's
#' portion of the eigenvector.
#' If the model has been iterated, but has not yet converged to asymptotic dynamics,
#' \code{right_ev} will try to iterate it further using the final population state
#' as the starting point. The default number of iterations is 100, and can be
#' adjusted using the \code{iterations} argument.
#' If the model hasn't been iterated, then \code{right_ev} will try iterating it
#' for \code{iterations} number of time steps and check for convergence. In the
#' latter two cases, if the model still has not converged to asymptotic dynamics,
#' it will return \code{NA} with a warning.
#'
#' For \code{left_ev}, the transpose iteration (\emph{sensu} Ellner & Rees 2006,
#' Appendix A) is worked out based on the \code{state_start} and \code{state_end}
#' in the model's \code{proto_ipm} object. The model is then iterated for
#' \code{iterations} times to produce a standardized left eigenvector.
#'
#' @section \strong{Stochastic eigenvectors}:
#' \code{left_ev} and \code{right_ev} return different things for stochastic models.
#' \code{right_ev} returns the trait distribution through time from the stochastic
#' simulation (i.e. \code{ipm$pop_state}), and normalizes it such that the
#' distribution at each time step integrates to 1 (if it is not already).
#' It then discards the first \code{burn_in * iterations} time steps of the
#' simulation to eliminate transient dynamics. See Ellner, Childs, & Rees 2016,
#' Chapter 7.5 for more details.
#'
#' \code{left_ev} returns a similar result as \code{right_ev}, except the trait
#' distributions are the result of left multiplying the kernel and trait
#'  distribution. See Ellner, Childs, & Rees 2016, Chapter 7.5 for more
#' details.
#'
#' @export

left_ev.general_di_det_ipm <- function(ipm,
                                       iterations = 100,
                                       tolerance = 1e-10,
                                       ...) {

  mod_nm    <- deparse(substitute(ipm))

  use_pop_state       <- ipm$pop_state[!grepl("lambda", names(ipm$pop_state))]

  init_pop_vec       <- lapply(use_pop_state, function(x) x[ , 1])

  # Remove NAs generated by the unsubstituted pop_state place holders
  init_pop_vec <- Filter(f = function(y) !any(is.na(y)), x = init_pop_vec)

  names(init_pop_vec) <- gsub('pop_state', 'n', names(init_pop_vec))

  test_conv           <- ipm$proto_ipm %>%
    define_pop_state(
      pop_vectors = init_pop_vec
    ) %>%
    make_ipm(iterate    = TRUE,
             iterations = iterations,
             iteration_direction = "left",
             normalize_pop_size = TRUE)

  if(all(is_conv_to_asymptotic(test_conv, tolerance = tolerance))) {

    out <- .extract_conv_ev_general(test_conv$pop_state, test_conv$proto_ipm)

  } else {

    warning(
      paste(
        "'",
        mod_nm,
        "'",
        ' did not converge after ',
        iterations,
        ' iterations. Returning NA, please try again with more iterations.',
        sep = ""
      )
    )

    return(NA_real_)

  }

  out <- Filter(function(x) !any(is.na(x)), out)
  # Replaces the leading n_ with a trailing _v
  names(out) <- substr(names(out), 3, nchar(names(out)))
  names(out) <- paste(names(out), 'v', sep = '_')

  class(out) <- "ipmr_v"

  return(out)
}

#' @rdname eigenvectors
#' @export

left_ev.general_di_stoch_kern_ipm <- function(ipm,
                                              iterations = 10000,
                                              burn_in    = 0.25,
                                              kernel_seq = NULL,
                                              ...) {

  mod_nm    <- deparse(substitute(ipm))

  use_pop_state      <- ipm$pop_state[!grepl("lambda", names(ipm$pop_state))]

  init_pop_vec       <- lapply(use_pop_state, function(x) x[ , 1])

  names(init_pop_vec) <- gsub('pop_state', 'n', names(init_pop_vec))

  if(is.null(kernel_seq)) kernel_seq <- ipm$env_seq

  test_conv           <- ipm$proto_ipm %>%
    define_pop_state(
      pop_vectors = init_pop_vec
    ) %>%
    make_ipm(iterate    = TRUE,
             iterations = iterations,
             iteration_direction = "left",
             normalize_pop_size = TRUE,
             kernel_seq = kernel_seq)

  out <- test_conv$pop_state[!grepl('lambda', names(test_conv$pop_state))]

  burn_in_seq <- seq_len(ncol(out[[1]])) %>%
    .[seq.int(from = ceiling(ncol(out[[1]]) * burn_in),
              to   = ncol(out[[1]]),
              by   = 1)]

  out <- lapply(out,
                function(x, burn_seq) x[ , burn_seq, drop = FALSE],
                burn_seq = burn_in_seq)

  # Replaces the leading n_ with a trailing _v
  names(out) <- substr(names(out), 3, nchar(names(out)))
  names(out) <- paste(names(out), 'v', sep = '_')

  class(out) <- "ipmr_v"

  return(out)
}

#' @rdname eigenvectors
#' @export

# Needs to re-iterate the model, but use the same kernel sequence that is used
# in right_ev/make_ipm. In effect, we want a deterministic model using the
# sub-kernels slot in lieu of parameters.

left_ev.general_di_stoch_param_ipm <- function(ipm,
                                               iterations = 10000,
                                               burn_in = 0.25,
                                               kernel_seq = NULL,
                                               ...) {

  mod_nm      <- deparse(substitute(ipm))

  sub_kernels <- ipm$sub_kernels
  proto_ipm   <- ipm$proto_ipm

  usr_its <- ncol(ipm$pop_state[[1]]) - 1

  if(usr_its < iterations) {

    stop("'iterations' should be less than or equal to the number of iterations in 'make_ipm()'!")

  }

  # Set up the model iteration expressions and initial population state.
  # This is a leaner version of make_ipm.general_di_stoch_param, because
  # we don't need all the user functions or other stuff - just the basic
  # infrastructure to house model iteration evaluation.

  proto_list <- .initialize_kernels(proto_ipm,
                                    iterate = TRUE,
                                    iter_dir = "left")

  others <- proto_list$others
  k_row  <- proto_list$k_row

  main_env <- .make_main_env(others$domain,
                             proto_ipm$usr_funs[[1]],
                             age_size = .uses_age(others))


  temp      <- .prep_di_output(others, k_row, proto_ipm,
                               iterations, normal= TRUE)
  main_env  <- .add_pop_state_to_main_env(temp$pop_state, main_env)

  # Throw this in just in case we have parameter set indices in addition
  # to continuous environmental variation.

  kern_seq <- .make_kern_seq(others, iterations, kernel_seq)

  # Create an index to select sub-kernels based on which iteration they were
  # generated in.

  kern_it_ind <- vapply(names(sub_kernels), function(x) {
    as.integer(strsplit(x, "(_it_)")[[1]][2])
  }, integer(1L))

  for(i in seq_len(iterations)) {

    # This selects sub-kernels by matching the _it_XXX suffix to the current
    # iteration.
    use_kerns <- sub_kernels[kern_it_ind == i]

    # Now, drop the _it_XXX suffix and pretend they're just regular kernels
    # generated by .make_sub_kernel_general_lazy.

    rm_regex <- paste("_it_", i, '$', sep = "")
    names(use_kerns) <- gsub(rm_regex, "", names(use_kerns))

    # We don't care about any of the other output from model iteration this,
    # so skip all the naming/stashing other outputs.

    pop_state <- .iterate_model(proto_ipm = proto_ipm,
                                k_row = k_row,
                                sub_kern_list = use_kerns,
                                current_iteration = i,
                                kern_seq = kern_seq,
                                pop_state = temp$pop_state,
                                main_env = main_env,
                                normal = TRUE)

    temp$pop_state <- pop_state
  }

  out <- temp$pop_state
  out[grepl("lambda", names(out))] <- NULL

  burn_in_seq <- seq_len(ncol(out[[1]])) %>%
    .[seq.int(from = ceiling(ncol(out[[1]]) * burn_in),
              to   = ncol(out[[1]]),
              by   = 1)]

  out <- lapply(out,
                function(x, burn_seq) x[ , burn_seq, drop = FALSE],
                burn_seq = burn_in_seq)

  names(out) <- gsub("pop_state_", "", names(out))
  names(out) <- paste(names(out), 'v', sep = '_')

  class(out) <- "ipmr_v"

  return(out)
}

#' @rdname eigenvectors
#' @export

left_ev.simple_di_stoch_param_ipm <- left_ev.general_di_stoch_param_ipm

#' @noRd
# Helper to be used in do.call(".add", X) where X is a list with potentially more
# than two elements. "+" won't work in this case, and sum() will return a single
# number.

.add <- function(...) {
  dots <- list(...)
  Reduce("+", x = dots, init = 0)
}



