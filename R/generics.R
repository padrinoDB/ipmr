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
#' @param all_lambdas A logical: print lambdas for each individual kernel/time step,
#' or just the final one? Only applies to \code{comp_method = "pop_size"}. When
#' \code{TRUE}, returns the ratio of population sizes for every iteration of the
#' model. When \code{FALSE}, only returns the ratio of population sizes for the
#' final iteration. \code{comp_method = 'eigen'} always returns the lambda values
#' for every single kernel in the \code{x$iterators} slot.
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
                                    all_lambdas = FALSE,
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

    lambda <- lambda(x, comp_method = comp_method, all_lambdas = all_lambdas)

    l_msg  <- paste0('\nDeterministic lambda for ', nm_ks,' = ', lambda, sep = "" )

    msg    <- c(msg, l_msg)

    if(comp_method == 'pop_size' &&
       check_conv &&
       !is_conv_to_asymptotic(lambda)) {

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
                                           comp_method = c('pop_size',
                                                           'eigen'),
                                           all_lambdas = FALSE,
                                           sig_digits = 3,
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
    lambda <- lambda(x, comp_method = comp_method, all_lambdas = all_lambdas)

    l_msg  <- paste0('\nDeterministic lambda for ',
                     nm_ks,
                     ' = ',
                     round(lambda, sig_digits),
                     sep = "" )

    msg <- c(msg, l_msg)
  }
  cat(msg)

  invisible(x)
}

#' @rdname print_star
#' @export

print.simple_di_stoch_param_ipm <- function(x,
                                            comp_lambda = TRUE,
                                            comp_method = c('pop_size',
                                                            'eigen'),
                                            all_lambdas = FALSE,
                                            sig_digits  = 3,
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
    lambda <- lambda(x, comp_method = comp_method, all_lambdas = all_lambdas)

    l_msg  <- paste0('\nDeterministic lambda for ',
                     nm_ks,
                     ' = ',
                     round(lambda, sig_digits),
                     sep = "" )

    msg <- c(msg, l_msg)
  }
  cat(msg)

  invisible(x)

}

#' @rdname print_star
#' @export

print.general_di_det_ipm <- function(x,
                                     comp_lambda = TRUE,
                                     comp_method = 'pop_size',
                                     sig_digits  = 3,
                                     all_lambdas = TRUE,
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

    all_lams <- lambda(x, comp_method = 'pop_size', all_lambdas = TRUE)

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
                                            sig_digits  = 3,
                                            all_lambdas = TRUE,
                                            ...) {

  pretty_cls <- .pretty_class(class(x)[1])

  msg <- paste0('A ',
                pretty_cls,
                ' IPM with ',
                length(x$sub_kernels),
                ' sub-kernel(s) and ',
                length(x$pop_state),
                ' population vectors defined.',
                sep = "")

  if(comp_lambda) {

    all_lams <- lambda(x, comp_method = 'pop_size', all_lambdas = TRUE)

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

    msg <- c(msg, l_msg)

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
                                             sig_digits  = 3,
                                             all_lambdas = TRUE,
                                             ...) {

  pretty_cls <- .pretty_class(class(x)[1])

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

    all_lams <- lambda(x, comp_method = 'pop_size', all_lambdas = TRUE)

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

    msg <- c(msg, l_msg)

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
#' @param all_lambdas A logical: return lambdas for each individual kernel, or
#' a single number? Only applies to \code{comp_method = "pop_size"}. When
#' \code{TRUE}, returns the ratio of population sizes for every iteration of the
#' model. When \code{FALSE}, only returns the ratio of population sizes for the
#' final iteration.
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
                                     all_lambdas = TRUE,
                                     ...) {

  switch(comp_method,
         'eigen'    = .lambda_eigen(ipm),
         'pop_size' = .lambda_pop_size(ipm, all_lambdas = all_lambdas))

}

#' @rdname lambda
#'
#' @export

lambda.simple_di_stoch_kern_ipm <- function(ipm,
                                            comp_method = c("eigen", "pop_size"),
                                            all_lambdas = TRUE,
                                            ...) {

  switch(comp_method,
         'eigen'    = .lambda_eigen(ipm),
         'pop_size' = .lambda_pop_size(ipm, all_lambdas = all_lambdas))
}

#' @rdname lambda
#'
#' @export

lambda.simple_di_stoch_param_ipm <- function(ipm,
                                             comp_method = c("eigen", "pop_size"),
                                             all_lambdas = TRUE,
                                             ...) {

  switch(comp_method,
         'eigen'    = .lambda_eigen(ipm),
         'pop_size' = .lambda_pop_size(ipm, all_lambdas = all_lambdas))
}

#' @rdname lambda
#' @export
#'
lambda.general_di_det_ipm <- function(ipm, all_lambdas = TRUE, ...) {
  .lambda_pop_size(ipm, all_lambdas = all_lambdas)
}


#' @rdname lambda
#'
#' @export
#'
lambda.general_di_stoch_kern_ipm <- function(ipm, ..., all_lambdas = TRUE) {

  .lambda_pop_size(ipm, all_lambdas = all_lambdas)

}

#' @rdname lambda
#'
#' @export
#'
lambda.general_di_stoch_param_ipm <- function(ipm, ..., all_lambdas = TRUE) {

  .lambda_pop_size(ipm, all_lambdas = all_lambdas)

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


# `[` ----------------
# ^^ I think these will exist... but hold on

# Sensitivitiy ---------------


# Elastictiy   ---------------


# diagnose_ipm ---------------



