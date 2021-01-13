#' @rdname stoch_helpers
#'
#' @title Helpers for stochastic models
#'
#' @description Various functions to compute useful quantities for stochastic
#' models
#'
#' @param ipm A stochastic model created by \code{make_ipm()}.
#'
#' @return \code{mean_kernel} returns a list of mean sub-kernels
#' for the model.
#'
#' @export

mean_kernel <- function(ipm) {

  cls <- class(ipm)[1]

  if(grepl("stoch_kern", cls)) {

    class(ipm) <- c(class(ipm), "stoch_kern")

  } else if(grepl("stoch_param", cls)) {

    class(ipm) <- c(class(ipm), "stoch_param")

  } else {

    stop("Cannot compute mean kernels for deterministic models.",
         call. = FALSE)
  }

  .mean_kernel(ipm)


}

#' @noRd

.mean_kernel <- function(ipm) {

  UseMethod(".mean_kernel")
}

#' @noRd

.mean_kernel.stoch_kern <- function(ipm) {

  proto     <- ipm$proto_ipm

  base_nms  <- proto$kernel_id

  sub_kerns <- ipm$sub_kernels

  out       <- list()

  for(i in seq_along(unique(base_nms))) {

    kern_nm <- unique(base_nms)[i]
    p_row   <- proto[proto$kernel_id == kern_nm, ]

    # If there aren't hierarachical effects, then we don't really
    # need to compute anything - there's just 1 version of the sub-kernel.

    if(!p_row$has_hier_effs) {

      use_nm    <- paste("mean_", kern_nm, sep = "")

      mean_kern <- rlang::list2(!!use_nm := sub_kerns[kern_nm])

      out       <- c(out, mean_kern)

      next
    }

    # otherwise - we generate exact names for each set of sub-kernels, extract them
    # from the sub_kernel list, and then compute the point-wise mean.

    levs <- .make_hier_levels(p_row$levels_hier_effs)
    kern_nms <- character(length(levs))
    to_sub   <- names(p_row$levels_hier_effs[[1]]) %>%
      .[!"to_drop" %in% .] %>%
      paste(collapse = "_")

    for(j in seq_along(levs)) {

      kern_nms[j] <- gsub(to_sub, levs[j], kern_nm)

    }

    use_kerns <- sub_kerns[kern_nms]

    mean_kern <- mean_kernel_impl(use_kerns)

    names(mean_kern) <- paste("mean_", kern_nm, sep = "")

    out <- c(out, mean_kern)
  }

  return(out)
}

#' @noRd

.mean_kernel.stoch_param <- function(ipm) {

  proto     <- ipm$proto_ipm

  base_nms  <- proto$kernel_id

  sub_kerns <- ipm$sub_kernels

  out       <- list()

  n_its     <- ncol(ipm$pop_state[[1]]) - 1

  for(i in seq_along(unique(base_nms))) {

    kern_nm <- unique(base_nms)[i]
    p_row   <- proto[proto$kernel_id == kern_nm, ]

    # If there aren't hierarachical effects, then we still need to compute
    # the mean of all iterations. The kernels will always have _it_x appended
    # to them to distinguish them from each other, so we create those labels,
    # then subset the kernel list w exact name matching

    if(!p_row$has_hier_effs) {

      use_nm    <- paste("mean_", kern_nm, sep = "")

      kern_nms  <- paste(kern_nm, "it", seq(1, n_its, by = 1), sep = "_")

      use_kerns <- sub_kerns[kern_nms]

      mean_kern <- mean_kernel_impl(use_kerns)

      names(mean_kern) <- use_nm

      out       <- c(out, mean_kern)

      next
    }

    # otherwise - we generate exact names for each set of sub-kernels, extract them
    # from the sub_kernel list, and then compute the point-wise mean.

    levs <- paste(ipm$env_seq$kernel_seq,
                  "it",
                  seq(1, n_its, by = 1),
                  sep = "_")

    to_sub   <- names(p_row$levels_hier_effs[[1]]) %>%
      .[!"to_drop" %in% .] %>%
      paste(collapse = "_")

    kern_nms <- character(n_its)

    for(j in seq_along(levs)) {
      kern_nms[j] <- gsub(to_sub, levs[j], kern_nm)
    }
    use_kerns <- list()

    for(j in seq_along(kern_nms)) {
      use_kerns <- c(use_kerns,
                     sub_kerns[grepl(kern_nms[j],
                               names(sub_kerns))])
    }
    mean_kern <- mean_kernel_impl(use_kerns)

    names(mean_kern) <- paste("mean_", kern_nm, sep = "")

    out <- c(out, mean_kern)
  }

  return(out)
}
