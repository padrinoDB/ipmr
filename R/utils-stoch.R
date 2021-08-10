#' @title Mean kernels for stochastic models
#'
#' @description This function computes mean sub-kernels for stochastic
#' parameter re-sampled and stochastic kernel re-sampled models.
#'
#' @param ipm A stochastic model created by \code{make_ipm()}.
#'
#' @details For \code{*_stoch_kern} models, this computes the element-wise
#' mean for each sub-kernel across all the different levels of
#' \code{par_set_indices}. For models where not all sub-kernels contain
#' parameter set indices, sub-kernels that do not have varying
#' parameters are included in the output and are identical to their input.
#'
#' For \code{*_stoch_param} models, this computes the element-wise mean for each
#' sub-kernel created by the iteration procedure.
#'
#' @return A list of mean sub-kernels for the model.
#'
#' @export

mean_kernel <- function(ipm) {
  UseMethod("mean_kernel")
}

#' @export

mean_kernel.ipmr_ipm <- function(ipm) {

  cls <- class(ipm)[1]

  if(grepl("stoch_kern", cls)) {

    class(ipm) <- c(class(ipm), "stoch_kern")

  } else if(grepl("stoch_param", cls)) {

    class(ipm) <- c(class(ipm), "stoch_param")

  } else {

    stop("Cannot compute mean kernels for deterministic models.",
         call. = FALSE)
  }

  out <- .mean_kernel(ipm) %>%
    set_ipmr_classes()

  return(out)

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

    # If there aren't par_setarachical effects, then we don't really
    # need to compute anything - there's just 1 version of the sub-kernel.

    if(!p_row$uses_par_sets) {

      use_nm    <- paste("mean_", kern_nm, sep = "")

      mean_kern <- rlang::list2(!!use_nm := sub_kerns[[kern_nm]])

      out       <- c(out, mean_kern)

      next
    }

    # otherwise - we generate exact names for each set of sub-kernels, extract them
    # from the sub_kernel list, and then compute the point-wise mean.

    levs     <- .make_par_set_indices(p_row$par_set_indices)

    kern_nms <- character(length(levs))

    to_sub   <- names(p_row$par_set_indices[[1]]) %>%
      .[!. %in% "to_drop"] %>%
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

    # If there aren't par_setarachical effects, then we still need to compute
    # the mean of all iterations. The kernels will always have _it_x appended
    # to them to distinguish them from each other, so we create those labels,
    # then subset the kernel list w exact name matching

    if(!p_row$uses_par_sets) {

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

    to_sub   <- names(p_row$par_set_indices[[1]]) %>%
      .[!"to_drop" %in% .] %>%
      paste(collapse = "_")

    kern_nms <- character(n_its)

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
