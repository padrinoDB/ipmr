#' @title Eviction correction
#' @rdname eviction
#'
#' @description Various helpers to correct for unintentional eviction (Williams
#' et al. 2012).
#'
#' @param discretized_mat The kernel or function that needs correcting.
#'
#' @return A matrix of the same dimension as the input.
#'
#' @references
#' Williams JL, Miller TEX & Ellner SP, (2012). Avoiding unintentional eviction
#' from integral projection models.Ecology 93(9): 2008-2014.
#'
#' @export

truncated_distributions <- function(discretized_kernel) {

  dim_in <- dim(discretized_kernel)

  out <- discretized_mat / matrix(
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

rescale_kernel <- function(discretized_kernel) {
  return(
    truncated_distributions(discretized_kernel = discretized_kernel)
  )
}
