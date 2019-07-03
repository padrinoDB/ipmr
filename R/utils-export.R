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

#' @title Helpers for IPM construction
#' @inheritParams define_kernel
#' @param kernels The kernel names you that correspond to given population state vector.
#' Suffix expansion is supported.
#' @param ... Named expressions. See details for more information on their usage in
#' each \code{define_*} function.
#'
#' @details
#' These are helper functions to define certain types of IPM classes. It is recommended
#' to use them within the IPM construction pipeline somewhere between \code{add_kernels}
#' and \code{make_ipm}. The order in which they are called does not matter provided
#' all values on the left hand side of the \code{...} appear in the right hand side
#' of another expression used elsewhere. Expressions where the left hand side does
#' not appear anywhere else will be ignored during construction, and so will not
#' be present in the returned \code{ipm} object.
#'
#' Each \code{define_*} function takes dots in the same way: named expressions
#' are captured and evaluated later in the appropriate contexts with their respective
#' kernels. Suffix expansion is supported so that hierarchical models with
#' year/plot/what-have-you effects do not need to  be rewritten 10s or 100s of times.
#'
#' @return All \code{define_*} functions return a proto_ipm.
#'
#' @rdname define_star
#' @export

define_pop_state <- function(proto_ipm, ...) {


}

#' @rdname define_star
#' @inheritParams define_pop_state
#' @export

define_state_vars <- function(proto_ipm, ...) {

  # DEFINE ME
}

#' @inheritParams define_pop_state
#' @param data_list A list of named values that correspond to constants in the formula.
#' You do not need to specify vectors corresponding to the domains here.
#' @rdname define_star
#' @export

define_env_state <- function(proto_ipm, ..., data_list) {

  # DEFINE ME
}


#' @rdname define_star
#' @inheritParams define_pop_state
#' @export

define_hier_effs <- function(proto_ipm, ...) {

  # DEFINE ME
}
