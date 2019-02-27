# General exported utilities

#' @rdname vec-operations
#'
#' @title Perform the \emph{vec} operation and it's inverse
#'
#' @description Unstacks a matrix into a single column vector. Each successive
#' column is appended to the bottom of the vector
#'
#' @param mat The matrix
#'
#' @return \code{vec}: A matrix with dimensions nrow(mat) * ncol(mat) X 1
#'
#' \code{inv_vec}: A matrix with dimensions \code{sqrt(length(mat))} X \code{sqrt(length(mat))} or
#' dimensions specified by \code{nrow} and \code{ncol}
#'
#' @examples
#'
#' mat <- matrix(1:20, nrow = 4, ncol = 5)
#'
#' v <- vec(mat)
#'
#' v
#'
#' inv_vec(v)
#'
#'
#' @export

vec <- function(mat) {
  matrix(mat, ncol = 1)
}

#' @rdname vec-operations
#'
#' @inheritParams vec
#' @param square A logical indicating whether to return a square matrix. For IPMs,
#' this should usually be \code{TRUE}. If set to \code{FALSE}, \code{nrow} and \code{ncol}
#' must be specified
#' @param nrow The number of rows for the new matrix
#' @param ncol The number of columns for the new matrix
#'
#' @export

inv_vec <- function(mat, square = TRUE, nrow = NULL, ncol = NULL) {

  if(!square & (is.null(nrow) | is.null(ncol))) {
    stop("Must specify 'nrow' and 'ncol' when square = FALSE")
  }

  if(square) {
    mat_dim <- sqrt(length(mat))

    out <- matrix(mat, nrow = mat_dim, ncol = mat_dim)
  } else {
    out <- matrix(mat, nrow = nrow, ncol = ncol)
  }

  return(out)

}
