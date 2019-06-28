#' @title Eviction correction
#' @rdname eviction
#'
#' @description Various helpers to correct for unintentional eviction (Williams
#' et al. 2012).
#'
#' @param discretized_kernel The kernel or function that needs correcting.
#' @param n_mesh_p The number of meshpoints for the kernel being corrected.
#'
#' @return A matrix of the same dimension as the input.
#'
#' @references
#' Williams JL, Miller TEX & Ellner SP, (2012). Avoiding unintentional eviction
#' from integral projection models.Ecology 93(9): 2008-2014.
#'
#' @export

truncated_distributions <- function(discretized_kernel,
                                    n_mesh_p) {

  if(is.matrix(discretized_kernel)){
    dim_in <- dim(discretized_kernel)
  } else {

    discretized_kernel <- matrix(discretized_kernel,
                                 nrow = n_mesh_p,
                                 ncol = n_mesh_p)
    dim_in <- c(n_mesh_p, n_mesh_p)

  }

  out <- discretized_kernel /
    matrix(
      as.vector(
        apply(
          discretized_kernel,
          2,
          sum
        )
      ),
      nrow = dim_in[1],
      ncol = dim_in[2],
      byrow = TRUE
    )
  return(out)

}

rescale_kernel <- function(discretized_kernel, n_mesh_p) {
  return(
    truncated_distributions(discretized_kernel = discretized_kernel,
                            n_mesh_p = n_mesh_p)
  )
}


# Internal helpers for make_ipm() methods

#' @noRd
.correct_eviction <- function(evict_fun, kernel_env) {

  # Set the quosure environment for eval_tidy and get the name of the symbol
  # being modified
  evict_correction <- rlang::quo_set_env(evict_fun,
                                         kernel_env)
  nm <- strsplit(rlang::quo_text(evict_correction), '\\(|,|\\)')[[1]][2]

  # Need to evaluate before bindidng as the symbols are no longer uniquely named
  assign(nm,
         rlang::eval_tidy(evict_correction),
         envir = kernel_env)

  invisible(kernel_env)
}
