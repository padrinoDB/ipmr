#' @title Eviction correction
#' @rdname eviction
#'
#' @description Various helpers to correct for unintentional eviction (Williams
#' et al. 2012).
#'
#' @param discretized_kernel The kernel or function that needs correcting.
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
