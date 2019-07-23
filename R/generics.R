# Print----------------

#' Methods for proto_ipms
#' @rdname proto_generics
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

#' @rdname ipm-generics
#' @title Generics for IPM classes
#'
#' @param x An object produced by \code{make_ipm}.
#' @param compute_lambda A logical indicating whether or not to calculate lambdas
#' for the iteration kernels and display them
#' @param sig_digits The number of significant digits to round to if \code{
#' comput_lambda = TRUE}.
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

#' @rdname ipm-generics
#' @inheritParams print.simple_di_det
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
                     lambda,
                     sep = "" )

    msg <- c(msg, l_msg)
  }
  cat(msg)

  invisible(x)
}

#' @rdname ipm-generics
#' @inheritParams print.simple_di_stoch_kern
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
#' @return A single numeric value.
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
#' by their prior sizes and taking the geometric mean of those ratios. Thus, only
#' \code{type = 'stochastic'} and \code{comp_method = "pop_size"} are available
#' for these models.
#'
#' @export

lambda <- function(ipm, ...) {
  UseMethod('lambda')
}

#' @rdname lambda
#' @inheritParams lambda
#' @param type Either \code{"stochastic"} or \code{"deterministic"}.
#' \code{"stochastic"} also has two types - \code{"eigen"} and \code{"pop_size"}.
#' See details for more information.
#'
#' @export

lambda.simple_di_det_ipm <- function(ipm, type = "deterministic", ...) {

  .det_lambda(ipm)

}

#' @rdname lambda
#' @inheritParams lambda
#' @param comp_method Either \code{"eigen"} or \code{"pop_size"}. \code{"eigen"}
#' is not possible except for \code{"simple_*_stoch_kern"} and \code{"simple_*_det"}
#'
#' @export

lambda.simple_di_stoch_kern_ipm <- function(ipm,
                                            type = "stochastic",
                                            comp_method = c("eigen", "pop_size"),
                                            ...) {

  switch(comp_method,
         'eigen'    = .stoch_lambda_eigen(ipm),
         'pop_size' = .stoch_lambda_pop_size(ipm))
}

#' @rdname lambda
#' @inheritParams lambda
#'
#' @export

lambda.simple_di_stoch_param_ipm <- function(ipm, ...) {

  .stoch_lambda_pop_size(ipm)

}

#' @rdname lambda
#' @inheritParams lambda
#' @export
#'
lambda.general_di_det_ipm <- function(ipm, type, ...) {
  .stoch_lambda_pop_size(ipm)
}

#' @rdname lambda
#' @inheritParams lambda.simple_di_stoch_kern_ipm
#'
#' @export
#'
lambda.general_di_stoch_param_ipm <- function(ipm, ...) {

  .stoch_lambda_pop_size(ipm)

}




# Plot ---------------


# `[` ----------------
# ^^ I think these will exist... but hold on

# Sensitivitiy ---------------


# Elastictiy   ---------------


# diagnose_ipm ---------------



