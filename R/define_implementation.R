#' @rdname define_star
#'
#' @inheritParams define_kernel
#' @param kernel_names A character vector with the names of the kernels
#' that parameters are being defined for.
#' @param kernel_impl_list A named list. Names correspond to kernel names. Each
#' kernel should have 3 slots defined - the \code{int_rule} (integration rule),
#' the \code{dom_start} (the domain the kernel begins on), and the \code{dom_end}
#' (the domain the kernel ends on). It is safest to use \code{make_impl_args_list}
#' to generate this.
#' @param dom_start The name of the state variable for the kernel at time \emph{t}.
#' @param dom_end The name of the state variable for the kernel at time \emph{t+1}.
#' This is usually the same as \code{dom_start}, but general IPMs
#' with discrete states or IPMs that move from one state to another (e.g. tree
#' seedling going from a height domain at \emph{t} to a DBH domain at \emph{t+1})
#' may have another value here. For cases with a discrete stage, kernels moving
#' individuals from discrete to continuous should have a state variable entered
#' here and an \code{NA} for \code{dom_start}. For kernels moving from continuous
#' to discrete, vice versa. For discrete to discrete, both are \code{NA}.
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

    domain_info                   <- .state_to_domain_info(kernel_impl_list[[i]]$dom_start,
                                                           kernel_impl_list[[i]]$dom_end)

    proto_ipm$domain[proto_ind]   <- list(domain_info)

  }

  class(proto_ipm) <- cls

  return(proto_ipm)
}

#' @rdname define_star
#' @export

make_impl_args_list <- function(kernel_names,
                                int_rule,
                                dom_start,
                                dom_end,
                                ...) {

  ln <- length(kernel_names)

  if(ln != length(int_rule)) {

    warning("Assuming that all kernels are implemented with the same",
            " 'int_rule'.",
            call. = FALSE)

    i_rule   <- replicate(ln, int_rule, simplify = FALSE) %>%
      unlist()

    int_rule <- i_rule[seq_len(ln)]
  }
  if(ln != length(dom_start)) {

    warning("Assuming that all kernels are implemented with the same",
            " 'dom_start'.",
            call. = FALSE)

    ds_rule   <- replicate(ln, dom_start, simplify = FALSE) %>%
      unlist()

    dom_start <- ds_rule[seq_len(ln)]
  }

  if(ln != length(dom_end)) {

    warning("Assuming that all kernels are implemented with the same",
            " 'dom_end'.",
            call. = FALSE)

    de_rule <- replicate(ln, dom_end, simplify = FALSE) %>%
      unlist()

    dom_end <- de_rule[seq_len(ln)]
  }

  out <- vector('list', length = ln)

  for(i in seq_along(kernel_names)) {
    out[[i]]$int_rule  <- int_rule[i]
    out[[i]]$dom_start <- dom_start[i]
    out[[i]]$dom_end   <- dom_end[i]
  }

  names(out) <- kernel_names

  return(out)

}


#' @noRd
.state_to_domain_info <- function(dom_start, dom_end) {

  # match names, then get info. Otherwise, generate an NA. the domain name
  # will always be first entry.

  if(!is.na(dom_start)) {

    start_state_info <- rep(NA_real_, 3)
    dom_start        <- paste(dom_start, "_1", sep = "")

  } else {

    start_state_info <- NA_real_
    dom_start        <- 'start_not_applicable'

  }

  if(!is.na(dom_end)) {

    end_state_info <- rep(NA_real_, 3)
    dom_end        <- paste(dom_end, "_2", sep = "")

  } else {

    end_state_info <- NA_real_
    dom_end        <- 'end_not_applicable'

  }

  out <- rlang::list2(!!dom_start := start_state_info,
                      !!dom_end   := end_state_info)
  return(out)
}
