#' @rdname define_star
#'
#' @inheritParams define_kernel
#' @param kernel_names A character vector with the names of the kernels
#' that parameters are being defined for.
#'
#' @param kernel_impl_list A named list. Names correspond to kernel names. Each
#' kernel should have 3 slots defined - the \code{int_rule} (integration rule),
#' the \code{state_start} (the domain the kernel begins on), and the \code{state_end}
#' (the domain the kernel ends on). For more complicated models, it is usually
#' safest to use \code{make_impl_args_list} to generate this.
#'
#' @param state_start The name of the state variable for the kernel that the
#' kernel acts on at time \emph{t}.
#'
#' @param state_end The name of the state variable that the kernel produces
#' at time \emph{t+1}.
#'
#' @param int_rule The integration rule to be used for the kernel. The default is
#' "midpoint". "trapezoid" and "g-l" (Gauss-Legendre) will be implemented as well.
#'
#' @export

define_impl <- function(proto_ipm,
                        kernel_impl_list) {

  cls     <- class(proto_ipm)
  kernels <- names(kernel_impl_list)

  for(i in seq_along(kernel_impl_list)) {

    proto_ind                     <- which(kernels[i] == proto_ipm$kernel_id)

    proto_ipm$int_rule[proto_ind] <- kernel_impl_list[[i]]$int_rule

    domain_info                   <- .state_to_domain_info(kernel_impl_list[[i]]$state_start,
                                                           kernel_impl_list[[i]]$state_end,
                                                           proto_ipm)


    proto_ipm$domain[proto_ind]   <- list(domain_info)

  }

  class(proto_ipm) <- cls

  return(proto_ipm)
}

#' @rdname define_star
#' @export

make_impl_args_list <- function(kernel_names,
                                int_rule,
                                state_start,
                                state_end) {

  ln <- length(kernel_names)

  if(ln != length(int_rule)) {

    warning("Assuming that all kernels are implemented with the same",
            " 'int_rule'.",
            call. = FALSE)

    i_rule   <- replicate(ln, int_rule, simplify = FALSE) %>%
      unlist()

    int_rule <- i_rule[seq_len(ln)]
  }
  if(ln != length(state_start)) {

    warning("Assuming that all kernels are implemented with the same",
            " 'state_start'.",
            call. = FALSE)

    ds_rule   <- replicate(ln, state_start, simplify = FALSE) %>%
      unlist()

    state_start <- ds_rule[seq_len(ln)]
  }

  if(ln != length(state_end)) {

    warning("Assuming that all kernels are implemented with the same",
            " 'state_end'.",
            call. = FALSE)

    de_rule <- replicate(ln, state_end, simplify = FALSE) %>%
      unlist()

    state_end <- de_rule[seq_len(ln)]
  }

  out <- vector('list', length = ln)

  for(i in seq_along(kernel_names)) {
    out[[i]]$int_rule  <- int_rule[i]
    out[[i]]$state_start <- state_start[i]
    out[[i]]$state_end   <- state_end[i]
  }

  names(out) <- kernel_names

  return(out)

}


#' @noRd
.state_to_domain_info <- function(state_start, state_end, proto_ipm) {

  # Generate dummy vectors to hold domain information.

  start_state_info <- rep(NA_real_, 3)

  end_state_info <- rep(NA_real_, 3)

  out <- rlang::list2(!!state_start := start_state_info,
                      !!state_end   := end_state_info)

  return(out)
}
