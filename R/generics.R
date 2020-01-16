# Print----------------

#' @title Print proto_ipms or *_ipm objects
#' @rdname print_star
#'
#' @param x An object of class \code{proto_ipm}
#' @param ... Ignored
#'
#' @export

print.proto_ipm <- function(x, ...) {

  n_kerns      <- dim(x)[1]

  cls_switch   <- class(x)[1]

  pretty_class <- .pretty_class(cls_switch)

  msg          <- paste("A",
                        pretty_class,
                        "proto_ipm with",
                        n_kerns,
                        "kernels defined:\n")

  msg          <- c(msg, paste(x$kernel_id, collapse = ', '))

  msg          <- c(msg, '\n')

  # Add more later -----------
  cat(msg)

  invisible(x)
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

      # Perhaps add in the age x state_var as well, but work out what that
      # even looks like first


    )

  return(out)
}

#' @rdname print_star
#' @title Generics for IPM classes
#'
#' @param x An object produced by \code{make_ipm}.
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

    lambdas <- lambda(x, comp_method = comp_method, type_lambda = 'all')

    l_msg  <- paste0('\nDeterministic lambda for ', nm_ks,' = ', lambdas, sep = "" )

    msg    <- c(msg, l_msg)

    if(comp_method == 'pop_size' &&
       check_conv &&
       !is_conv_to_asymptotic(lambdas)) {

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

  print(msg)

  invisible(x)


}

#' @rdname print_star
#' @export

print.simple_di_stoch_kern_ipm <- function(x,
                                           comp_lambda = TRUE,
                                           comp_method = c('pop_size'),
                                           type_lambda = 'all',
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
    all_lams <- lambda(x, comp_method = comp_method, type_lambda = type_lambda)

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
                                            comp_method = c('pop_size'),
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

    all_lams <- lambda(x, comp_method = comp_method, type_lambda = type_lambda)

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
#' @param all_lambdas Return lambdas for every time step or just the final one?
#' @export

print.general_di_det_ipm <- function(x,
                                     comp_lambda = TRUE,
                                     comp_method = 'pop_size',
                                     type_lambda = 'all',
                                     sig_digits  = 3,
                                     all_lambdas = FALSE,
                                     check_conv  = TRUE,
                                     ...) {

  pretty_cls <- .pretty_class(class(x)[1])

  msg <- paste0('A ',
                pretty_cls,
                ' IPM with ',
                length(x$sub_kernels),
                ' sub-kernel(s) and ',
                length(x$pop_state),
                ' population vectors defined.', sep = "")

  if(comp_lambda) {

    all_lams <- lambda(x, comp_method = 'pop_size', type_lambda = type_lambda)

    if(!all_lambdas) {

      ret_lam  <- all_lams[length(all_lams)]

      l_msg <- paste('\nLambda for the final time step of the model is: ',
                     round(ret_lam, sig_digits),
                     sep = "")

    } else {

      ret_lam <- all_lams

      l_msg <- paste('\nLambda for time step ',
                     seq_len(length(ret_lam)),
                     ' is: ',
                     round(ret_lam, sig_digits),
                     sep = "")

    }

    if(check_conv && !is_conv_to_asymptotic(all_lams)) {

      # Captures the name of the model that the user gave rather than
      # just print "x isn't converged"

      mod_nm <- deparse(substitute(x))


      message(
        paste(mod_nm,
              ' has has not converged to asymptotic dynamics!',
              sep = "")
      )
    }

    msg <- c(msg, l_msg)

  }

  cat(msg)

  invisible(x)

}

#' @rdname print_star
#' @export

print.general_di_stoch_kern_ipm <- function(x,
                                            comp_lambda = TRUE,
                                            comp_method = 'pop_size',
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

    all_lams <- lambda(x, comp_method = 'pop_size', type_lambda = type_lambda)

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
                                             comp_method = 'pop_size',
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

    all_lams <- lambda(x, comp_method = 'pop_size', type_lambda = type_lambda)

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

# Lambda------------
#' @title Compute the per-capita growth rate for an IPM object
#' @rdname lambda
#'
#' @param ipm An object returned by \code{make_ipm()}.
#' @param comp_method Either \code{"eigen"} or \code{"pop_size"}. \code{"eigen"}
#' is only possible for \code{"simple_*"} methods.
#' @param type_lambda Either \code{'all'} or \code{'stochastic'}. \code{'all'}
#' returns a vector of lambda values for each time step of the simulation (equal
#' in length to the \code{iterations} argument of \code{make_ipm()}).
#' \code{'stochastic'} returns a single value which is the geometric mean
#' of log-per-capita growth rate from each time step.
#' @param ... other arguments passed to methods.
#'
#' @return A single numeric vector or single value.
#'
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
#'
#' @export

lambda.simple_di_det_ipm <- function(ipm,
                                     comp_method = c("eigen", "pop_size"),
                                     type_lambda = 'all',
                                     ...) {

  switch(comp_method,
         'eigen'    = .lambda_eigen(ipm),
         'pop_size' = .lambda_pop_size(ipm, all_lambdas = TRUE))

}

#' @rdname lambda
#'
#' @export

lambda.simple_di_stoch_kern_ipm <- function(ipm,
                                            comp_method = c("eigen",
                                                            "pop_size"),
                                            type_lambda = c('all',
                                                            'stochastic'),
                                            ...) {

  temp <- switch(comp_method,
                 'eigen'    = .lambda_eigen(ipm),
                 'pop_size' = .lambda_pop_size(ipm, all_lambdas = TRUE))

  return(
    switch(type_lambda,
           'all'        = temp,
           'stochastic' = .geom_mean(log(temp)))
  )
}

#' @rdname lambda
#'
#' @export

lambda.simple_di_stoch_param_ipm <- function(ipm,
                                             comp_method = c("eigen", "pop_size"),
                                             type_lambda = c('all',
                                                             'stochastic'),
                                             ...) {

  temp <- switch(comp_method,
                 'eigen'    = .lambda_eigen(ipm),
                 'pop_size' = .lambda_pop_size(ipm, all_lambdas = TRUE))

  return(
    switch(type_lambda,
           'all'        = temp,
           'stochastic' = .geom_mean(log(temp)))
  )

}

#' @rdname lambda
#' @export
#'
lambda.general_di_det_ipm <- function(ipm, type_lambda = 'all', ...) {

  temp <- .lambda_pop_size(ipm, all_lambdas = TRUE)

  return(
    switch(type_lambda,
           'all'        = temp,
           'stochastic' = .geom_mean(log(temp)))
  )
}


#' @rdname lambda
#'
#' @export
#'
lambda.general_di_stoch_kern_ipm <- function(ipm,
                                             ...,
                                             type_lambda = c('all', 'stochastic')) {

  temp <- .lambda_pop_size(ipm, all_lambdas = TRUE)

  return(
    switch(type_lambda,
           'all'        = temp,
           'stochastic' = .geom_mean(log(temp)))
  )
}

#' @rdname lambda
#'
#' @export
#'
lambda.general_di_stoch_param_ipm <- function(ipm,
                                              ...,
                                              type_lambda = c('all', 'stochastic')) {

  temp <- .lambda_pop_size(ipm, all_lambdas = TRUE)

  return(
    switch(type_lambda,
           'all'        = temp,
           'stochastic' = .geom_mean(log(temp)))
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

  if(is.null(x)) x = seq_len(ncol(A))
  if(is.null(y)) y = seq_len(nrow(A))

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


# right_ev ----------------

#' @rdname eigenvectors
#'
#' @title Compute the standardized left and right eigenvectors via iteration
#'
#' @param ipm Output from \code{make_ipm()}.
#' @param ... other arguments passed to methods
#'
#' @value A list of named numeric vector(s) corresponding to the stable trait distribution
#' function (\code{right_ev}) or the reproductive values for each trait (\code{left_ev}).
#'
#' @details If the model has already been iterated, then these functions
#' will just extract population state of the final iteration and return
#' that in a named list. Each element of the list is a vector with length
#' \code{>= 1} and corresponds each state variable's portion of the eigenvector.
#'
#' Note that \code{right/left_ev} methods only exist for deterministic IPMs. To
#'  see if there is a quasi-stable distribution for stochastic models, use
#' \code{\link{qsd_converge()}}.
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

  pop_nm      <- ipm$proto_ipm$state_var %>%
    unlist() %>%
    unique() %>%
    .[1]

  # if it's already been iterated to convergence, we don't have much work to do.

  if(.already_iterated(ipm)) {

    # get index for population vector of final iteration
    final_it <- dim(ipm$pop_state[[1]])[2]

    if(is_conv_to_asymptotic(ipm$pop_state[[1]])) {

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

      if(is_conv_to_asymptotic(test_conv$pop_state[[1]])) {

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
    init_pop      <- runif(len_pop_state)

    test_conv     <- ipm$proto_ipm %>%
      define_pop_state(!! pop_states[1] := init_pop) %>%
      make_ipm(iterate = TRUE,
               iterations = n_iterations)

    if(is_conv_to_asymptotic(test_conv$pop_state[[1]])) {

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
                                        ...) {

  mod_nm    <- deparse(substitute(ipm))

  final_it <- dim(ipm$pop_state[[1]])[2]

  if(is_conv_to_asymptotic(ipm$pop_state)) {

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

    init_pop_vec        <- lapply(ipm$pop_state,
                                  function(x, final_it) x[ , final_it],
                                  final_it = final_it)

    names(init_pop_vec) <- gsub('pop_state', 'n', names(init_pop_vec))

    test_conv           <- ipm$proto_ipm %>%
      define_pop_state(
        pop_vectors = init_pop_vec
      ) %>%
      make_ipm(iterate    = TRUE,
               iterations = n_iterations)

    if(is_conv_to_asymptotic(test_conv$pop_state)) {

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

  return(out)
}

# left_ev -----------------

#' @rdname eigenvectors
#' @export

left_ev <- function(ipm, ...) {

  UseMethod('left_ev')

}

left_ev.simple_di_det_ipm <- function(ipm, n_iterations = 100, ...) {

  mod_nm  <- deparse(substitute(ipm))

  # Identify state variable name

  pop_nm <- ipm$proto_ipm$state_var %>%
    unlist() %>%
    unique() %>%
    .[1]

  # If it's already iterated, then we need to wrap K with a t()
  if(.already_iterated(ipm)) {

    # Models that are already iterated may have additional things.
    # however, the id of the kernel should match the top level kernel in this
    # list, so we'll use that.

    k_nm  <- names(ipm$iterators)
    k_ind <- ifelse(length(k_nm) > 1,
                    which(k_nm %in% ipm$proto_ipm$kernel_id),
                    1)

    t_k   <- t(ipm$iterators[[1]])

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

    if(is_conv_to_asymptotic(temp_pop_state)) {

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

    temp_pop_state[ , 1] <- runif(len_pop_state)

    for(i in seq_len(n_iterations)) {

      temp_pop_state[ , (i + 1)] <- t_k %*% temp_pop_state[ , i]

    }

    if(is_conv_to_asymptotic(temp_pop_state)) {

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


left_ev.general_di_det_ipm <- function(ipm,
                                       n_iterations = 100,
                                       ...) {




  param_ind <- vapply(ipm$proto_ipm$params[[proto_ind]]$formula,
                      function(x) grepl(k_nm,
                                        x),
                      logical(1L))

  text_call <- ipm$proto_ipm$params[[proto_ind]]$formula[[param_ind]]
  mod_call  <- gsub(k_nm,
                    paste('t(', k_nm, ')', sep = ""),
                    text_call)

  ipm$proto_ipm$params[[proto_ind]]$formula[[param_ind]] <- mod_call

  iterated <- ipm$proto_ipm %>%
    make_ipm(iterate    = TRUE,
             iterations = n_iterations)



}

# ev helper functions -------------

# Checks if model has already been iterated, saving us from re-iterating a model
# when we can just use the pop_state slot

.already_iterated <- function(ipm) {

  return(
    ! all(
      is.na(
        ipm$pop_state
      )
    )
  )


}

#' @noRd

# We need the total pop size to standardize our vectors, but we don't really
# know the proper order in which to combine them to create a single output.
# thus, we keep them in a list with entries corresponding to states, and
# use Reduce to compute the total size. The final output will be the list
# of states standardized by this total population size.

.extract_conv_ev_general <- function(pop_state) {

  final_it <- dim(pop_state[[1]])[2]

  temp <- lapply(pop_state,
                 function(x, final_it) {

                   x[ , final_it]

                 },
                 final_it = final_it)

  pop_std <- Reduce('sum', unlist(temp), init = 0)

  out <- lapply(temp,
                function(x, pop_std) x / pop_std,
                pop_std = pop_std)

  return(out)
}


#' @rdname qsd
#' @title Compute the quasi-stable trait distribution for a stochastic model
#'
#' @param ipm The output from \code{make_ipm}
#' @param p_nms The name of the survival/growth kernel(s)
#' @param tol The distance between each kernels true eigenvector and the
#' quasi-stable distribution that is acceptable. Computed as
#' \code{0.5 * sum(abs(pop_size - right_ev))}. Default is \code{1e-7}
#' @param start_life The index of the population vector to start on. Should
#' be a named list where names correspond to population states and entries
#' are length 1 integer vectors. Default is 1.
#' @param n_steps The number of times to iterate the model before checking
#' for convergence. Default is 1000.
#' @param ... other arguments passed to methods.
#'

qsd_converge <- function(ipm, p_nms, tol, start_life, n_steps, ...) {

  UseMethod('qsd_converge')

}

qsd_converge.simple_di_stoch_kern_ipm <- function(ipm,
                                                  p_nms,
                                                  tol = 1e-7,
                                                  start_life = rep(1,
                                                                   length(p_nms)),
                                                  n_steps = 1000,
                                                  ...) {

  mat_ps <- ipm$sub_kernels[p_nms]

  r_evs  <- lapply(mat_ps,
                   function(x) Re(eigen(x)$vectors[ , 1])
  )

  ns     <- lapply(mat_ps,
                   function(x) rep(0, dim(x)[1]))

  ns     <- lapply(seq_along(ns),
                   function(ind, evs, pop_states, start_lifes) {

                     pop_states[[ind]][start_lifes[ind]] <- sum(evs[[ind]])
                   },
                   evs         = r_evs,
                   pop_states  = ns,
                   start_lifes = start_life)

  dist <- numeric(n_steps)

  # Resume here------------------
  for(i in seq_len(n_steps)) {
    # p <- n / sum(n)
    dist[i] <- 0.5 * (sum(abs(n - r_ev)))
    n <- mat_p %*% n
  }

  if(min(dist, na.rm = TRUE) < tol) {
    out <- which.min(dist < tol)
  } else {
    out <- NA_integer_
  }

  return(out)
}
# Sensitivitiy ---------------

#' @rdname sensitivity
#' @title Compute sensitivity
#'
#' @param ipm Output from \code{make_ipm()}.
#' @param level The level to compute sensitivity at. \code{"kernel"} computes
#' the model wide sensitivity surface and returns that. \code{"vital_rate"}
#' computes lambda's sensitivity to specific vital rates, \code{"parameter"} lambda's
#' sensitivity to each parameter. Use \code{subset} to specify a subset of
#' vital rates or parameters to compute values for.
#' @param subset A character vector corresponding to the \code{"vital_rates"} or
#' \code{"parameters"} you want to compute  sensitivities for. Can save time if
#' only a few parameters or vital rates are of interest.
#' @param ... other arguments passed to methods
#'
#' @return A list. Contains either kernel response surfaces, vital rate level
#' sensitivity values, or parameter level sensitivity values. The class is always
#' \code{c("ipmr_sensitivity", "list")} (for internal usage and plot methods).
#'
#' @export

sensitivity <- function(ipm,
                        level  = c('kernel', 'vital_rate', 'parameter'),
                        subset = NA_character_,
                        ...) {

  UseMethod('sensitivity')

}

# Elastictiy   ---------------

#' @rdname elasticity
#' @title Compute elasticity
#'
#' @param ipm Output from \code{make_ipm()}.
#' @param level The level to compute elasticity at. \code{"kernel"} computes
#' the model wide elasticity surface and returns that. \code{"vital_rate"}
#' computes lambda's elasticity to specific vital rates, \code{"parameter"} lambda's
#' elasticity to each parameter. Use \code{subset} to specify a subset of
#' vital rates or parameters to compute values for.
#' @param subset A character vector corresponding to the \code{"vital_rates"} or
#' \code{"parameters"} you want to compute  sensitivities for. Can save time if
#' only a few parameters or vital rates are of interest.
#' @param ... other arguments passed to methods
#'
#' @return A list. Contains either kernel response surfaces, vital rate level
#' elasticity values, or parameter level elasticity values. The class is always
#' \code{c("ipmr_elasticity", "list")} (for internal usage and plot methods).
#'
elasticity <- function(ipm,
                       level = c('kernel', 'vital_rate', 'parameter'),
                       subset = NA_character_,
                       ...) {

  UseMethod('elasticity')

}

# diagnose ---------------

diagnose <- function(ipm, ...) {

  UseMethod('diagnose')

}


