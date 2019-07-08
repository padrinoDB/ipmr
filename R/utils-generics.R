#' Generic methods for \code{proto_ipm}s
#' @rdname proto_generics
#'
#' @param x An object of class \code{proto_ipm}
#' @param ... Ignored
#'
#' @export

print.proto_ipm <- function(x, ...) {

  n_kerns <- dim(x)[1]

  cls_switch <- class(x)[1]

  pretty_class <- .pretty_class(cls_switch)

  msg <- paste("A", pretty_class, "proto_ipm with", n_kerns, "kernels defined:\n")

  msg <- c(msg, paste(x$kernel_id, collapse = ', '))

  msg <- c(msg, '\n')

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

