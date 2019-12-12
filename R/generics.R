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
  out <- switch(cls_switch,
                'simple_di_det'          = "simple, density independent, deterministic",
                'simple_di_stoch_kern'   = "simple, density independent, stochastic, kernel-resampled",
                'simple_di_stoch_param'  = "simple, density independent, stochastic, parmaeter-resampled",

                'simple_dd_det'          = "simple, density dependent, deterministic",
                'simple_dd_stoch_kern'   = "simple, density dependent, stochastic, kernel-resampled",
                'simple_dd_stoch_param'  = "simple, density dependent, stochastic, parmaeter-resampled",

                'general_di_det'         = "general, density independent, deterministic",
                'general_di_stoch_kern'  = "general, density independent, stochastic, kernel-resampled",
                'general_di_stoch_param' = "general, density independent, stochastic, parmaeter-resampled",

                'general_dd_det'         = "general, density dependent, deterministic",
                'general_dd_stoch_kern'  = "general, density dependent, stochastic, kernel-resampled",
                'general_dd_stoch_param' = "general, density dependent, stochastic, parmaeter-resampled"

                # Perhaps add in the age x state_var as well, but work out what that
                # even looks like first
  )

  return(out)
}

#' @rdname print_star
#' @title Generics for IPM classes
#'
#' @param x An object produced by \code{make_ipm}.
#' @param compute_lambda A logical indicating whether or not to calculate lambdas
#' for the iteration kernels and display them
#' @param sig_digits The number of significant digits to round to if \code{
#' compute_lambda = TRUE}.
#' @param ... Ignored
#'
#' @return \code{x} invisibly.
#'
#' @export

print.simple_di_det_ipm <- function(x, compute_lambda = TRUE,
                                    sig_digits = 3,
                                    ...) {
  msg <- paste0('A simple, density independent, deterministic IPM with ',
                length(x$iterators),
                ' iteration kernel(s) and ',
                length(x$sub_kernels),
                ' sub-kernel(s) defined.', sep = "")

  if(compute_lambda) {

    nm_ks  <- names(x$iterators)

    lambda <- .det_lambda(x)

    l_msg  <- paste0('\nDeterministic lambda for ', nm_ks,' = ', lambda, sep = "" )

    msg    <- c(msg, l_msg)

  }

  cat(msg)

  invisible(x)


}

#' @rdname print_star
#' @param lambda_type If \code{compute_lambda} is \code{TRUE}, then either
#' \code{"stochastic"} or \code{"deterministic"}. \code{"deterministic"} will return
#' the dominant eigenvalue of the iteration kernel from each iteration. \code{
#' "stochastic"} depends on \code{compute_type}.
#' @param compute_type Either \code{"pop_size"} or \code{"eigen"}. \code{"pop_size"}
#' computes lambda as the ratio of population sizes between successive time steps
#' and then takes the geometric mean. \code{"eigen"} computes the dominant eigenvalue
#' of each iteration kernel and then takes the geometric mean. For large population
#' vectors, \code{"pop_size"} will likely be substantially faster. Note that
#' option \code{"pop_size"} is only possible if an initial population vector was
#' supplied when constructing the IPM.
#'  @export
print.simple_di_stoch_kern_ipm <- function(x,
                                           compute_lambda = TRUE,
                                           lambda_type = c("stochastic",
                                                           "deterministic"),
                                           compute_type = c('pop_size',
                                                            'eigen'),
                                           sig_digits = 3, ...) {

  msg <- paste0('A simple, density independent, deterministic IPM with ',
                length(x$iterators),
                ' iteration kernel(s) and ',
                length(x$sub_kernels),
                ' sub-kernel(s) defined.', sep = "")

  if(compute_lambda){
    nm_ks  <- names(x$iterators)
    lambda <- switch(lambda_type,
                     'stochastic'    = switch(compute_type,
                                              'pop_size' = .stoch_lambda_pop_size(x),
                                              'eigen'    = .stoch_lambda_eigen(x)),
                     'deterministic' = .det_lambda(x))


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
                                            compute_lambda = TRUE,
                                            lambda_type = c("stochastic",
                                                            "deterministic"),
                                            compute_type = c('pop_size',
                                                             'eigen'),
                                            sig_digits = 3,
                                            ...) {
  cat('not yet implemented')
  invisible(x)
}


# Lambda------------
#' @title Compute the per-capita growth rate for an IPM object
#' @rdname lambda
#'
#' @param ipm An object returned by \code{make_ipm()}.
#' @param ... other arguments passed to methods.
#'
#' @return A single numeric vector or single value.
#'
#'
#' @details Determinstic lambda is computed as the dominant eigenvalue of the
#' iteration kernel(s) where possible. For \code{simple_*_stoch_kern_ipm} models,
#' a vector containing the dominant eigenvalue of each entry in
#' \code{ipm$iterators}. For \code{simple_*_det_ipm} models, a single value.
#'
#'  NOTE ^^^^ Needs clarification for simple_dd_det models - that can either be
#'  pop_size or eigen. implement the make_ipm methods before returning!
#'
#' Because \code{ipmr} doesn't construct the iteration
#' kernel for \code{general_*} IPMs the way that they are typically implemented in
#' interactive sessions (e.g. by r/cbind()ing discrete and continuous stages together),
#' lambda must be computed by dividing successive population sizes
#' by their prior sizes. Thus, only
#' \code{type = 'stochastic'} and \code{comp_method = "pop_size"} are available
#' for these models.
#'
#' @export

lambda <- function(ipm, ...) {
  UseMethod('lambda')
}

#' @rdname lambda
#' @param type Either \code{"stochastic"} or \code{"deterministic"}.
#' \code{"stochastic"} also has two types - \code{"eigen"} and \code{"pop_size"}.
#' See details for more information.
#'
#' @export

lambda.simple_di_det_ipm <- function(ipm, type = "deterministic", ...) {

  .det_lambda(ipm)

}

#' @rdname lambda
#' @param comp_method Either \code{"eigen"} or \code{"pop_size"}. \code{"eigen"}
#' is not possible except for \code{"simple_*_stoch_kern"} and \code{"simple_*_det"}.
#' @param all_lambdas A logical to return lamdas for each individual kernel, or a single
#' number. For \code{comp_method = 'pop_size'} and \code{all_lambdas = FALSE},
#' this will return the final value of computed lambdas (e.g. presumably when
#' convergence has been reached). For \code{comp_method = 'eigen'} and
#' \code{all_lambdas = FALSE}, it will return the geometric mean of the
#' dominant eigenvalues for each projection matrix. For \code{all_lambdas = TRUE},
#' it will always return a numeric vector of lambdas computed either from the ratio
#' of population sizes or dominant eigenvalues of each kernel.
#'
#'
#' @export

lambda.simple_di_stoch_kern_ipm <- function(ipm,
                                            type = "stochastic",
                                            comp_method = c("eigen", "pop_size"),
                                            all_lambdas = TRUE,
                                            ...) {

  switch(comp_method,
         'eigen'    = .stoch_lambda_eigen(ipm, all_lambdas = all_lambdas),
         'pop_size' = .stoch_lambda_pop_size(ipm, all_lambdas = all_lambdas))
}

#' @rdname lambda
#'
#' @export

lambda.simple_di_stoch_param_ipm <- function(ipm, ..., all_lambdas = TRUE) {

  .stoch_lambda_pop_size(ipm, all_lambdas = all_lambdas)

}

#' @rdname lambda
#' @export
#'
lambda.general_di_det_ipm <- function(ipm, type, ...) {
  .stoch_lambda_pop_size(ipm)
}


#' @rdname lambda
#'
#' @export
#'
lambda.general_di_stoch_kern_ipm <- function(ipm, ..., all_lambdas = TRUE) {

  .stoch_lambda_pop_size(ipm, all_lambdas = all_lambdas)

}

#' @rdname lambda
#'
#' @export
#'
lambda.general_di_stoch_param_ipm <- function(ipm, ...) {

  .stoch_lambda_pop_size(ipm)

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



