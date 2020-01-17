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

    out     <- matrix(mat, nrow = mat_dim, ncol = mat_dim)

  } else {

    out     <- matrix(mat, nrow = nrow, ncol = ncol)

  }

  return(out)

}

#' @title Helpers for IPM construction
#' @inheritParams define_kernel
#' @param ... Named expressions. See details for more information on their usage in
#' each \code{define_*} function.
#'
#' @param pop_vectors If the population vectors are already pre-defined (i.e. are
#' not defined by a function passed to \code{...}), then they can
#' be passed as a named list here.
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
#' Most \code{define_*} functions takes dots in the same way: named expressions
#' are captured and evaluated later in the appropriate contexts with their respective
#' kernels. Suffix expansion is supported so that hierarchical models with
#' year/plot/what-have-you effects do not need to  be rewritten 10s or 100s of times.
#'
#' The one exception to this is \code{define_domains}. It takes named numeric
#' vectors of length 3 where the name corresponds to the
#' state variable, the first entry is the lower bound of the domain, the second
#' is the upper bound of the domain, and the third entry is the number of
#' meshpoints (applies to \code{int_rule = "midpoint"} only! Other rules are not
#' implemented yet).
#'
#' \code{define_impl} is meant to help distinguish the process of
#' generating the kernels' mathematical form from their implementation details.
#' It takes a \code{proto_ipm} object and returns a modified one containing the
#' implementation information.
#'
#' \code{make_impl_args_list} helps generate that
#' information in the correct format. It is usually easiest to call \code{make_impl_args_list}
#' before calling \code{init_ipm} and then substituting that variable into
#' the call to \code{define_impl}. Alternatively, one can do something like
#' \code{define_impl(kernel_impl_list = make_impl_arg_list(...))} within the course
#' of the model definitition pipeline.
#'
#' @return All \code{define_*} functions return a proto_ipm. \code{make_impl_args_list}
#' returns a list, and so must be used within a call to \code{define_impl} or
#' before initiating the piped model creation procedure.
#'
#' @rdname define_star
#' @importFrom rlang is_empty
#' @export

define_pop_state <- function(proto_ipm, ..., pop_vectors = list()) {

  pop_quos            <- rlang::enquos(...)

  temp                <- rlang::list2(!!! pop_quos, !!! pop_vectors)

  out                 <- Filter(Negate(rlang::is_empty), temp)

  names(out)          <- gsub('^n_', 'pop_state_', names(out))

  proto_ipm$pop_state <- list(out)

  return(proto_ipm)

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

  env_quos            <- rlang::enquos(...)

  out                 <- list(env_quos = unlist(env_quos),
                              constants = data_list)

  proto_ipm$env_state <- list(out)

  return(proto_ipm)
}


#' @rdname define_star
#' @inheritParams define_pop_state
#' @export

define_hier_effs <- function(proto_ipm, ...) {

  # DEFINE ME
}


#' @rdname fun-mult-helpers
#'
#' @param ... Symbols representing functions - usually things that appear on the
#' left hand side of \code{formula} and/or \code{...} in \code{define_kernel} and
#' \code{define_k}.
#'
#' @export
right_mult <- function(...) {

  to_mult <- list(...)
  dims <- vapply(to_mult,
                 function(x) {
                   if(is.matrix(x)) {
                     dim(x)
                   } else {
                     c(NA_integer_, NA_integer_)
                   }
                 },
                 integer(2))

  id_dim <- max(dims, na.rm = TRUE)

  init   <- diag(id_dim) # Identity matrix as initial starting point

  Reduce('%*%', to_mult, init = init)
}

# #' @rdname fun-mult-helpers
# #' @inheritParams left_mult
# #' @export
#
# left_mult <- function(...) {
#   to_mult <- rev(list(...))
#   init <- diag(dim(to_mult[[length(to_mult)]])[1]) # Identity matrix as initial starting point
#
#   Reduce('%*%', to_mult, init = init)
#
# }

#' @title Function multiplication helpers
#' @rdname fun-mult-helpers
#'
#' @description Helpers for multiplying functions - follows the terminology from
#' Ellner, Rees & Childs (2016), Table 3.1. \code{s_g_mult} is a special function
#' that correctly multiplies survival and growth
#'
#' @param s A vector of survival probabilities
#' @param g A discretized growth kernel
#'
#' @export

s_g_mult <- function(s, g) {

  return(t(s * t(g)))

}

#' @title Raise a matrix to a power
#' @rdname matrix-power
#'
#' @description Raises a matrix \code{x} to the \code{y}-th power. \code{x ^ y} computes
#' element wise powers, whereas this computes \emph{y - 1} matrix multiplications.
#' \code{mat_power(x, y)} is identical to \code{x \%^\% y}.
#'
#' @param x A numeric or integer matrix.
#' @param y An integer.
#'
#' @return A matrix.
#'
#' @export
#'

`%^%` <- function(x, y) {

  if(!is_square(x)) {
    stop('not implemented for non-square matrices')
  }

  if(!is.integer(y)) {

    warning("`%^%` is coercing second argument to an integer",
            call. = FALSE)

    y <- as.integer(y)

  }

  init_dim <- dim(x)[1]

  use_list <- lapply(seq_len(y), function(a, b) b, b = x)

  init_i   <- diag(init_dim)

  out <- Reduce('%*%', use_list, init = init_i)

  return(out)

}

#' @rdname matrix-power
#'
#' @export

mat_power <- function(x, y) {

  return(x %^% y)

}

#' @title Format a mega-matrix
#'
#' @param ipm Output from \code{make_ipm}.
#' @param mega_mat A vector with symbols and/or 0s representing the matrix blocks.
#' They should be specified in ROW MAJOR order! See examples.
#' @param presets Either empty or one of 'age-size' or 'hier_effs'. Currently
#' not implemented
#'
#' @return A large matrix
#'
#' @examples

#' data(gen_di_det_ex)
#'
#' big_k <- format_mega_matrix(gen_di_det_ex,
#'                             mega_mat = c(0, go_discrete,
#'                                          leave_discrete, P))
#'
#'
#'
#' @export

format_mega_matrix <- function(ipm, mega_mat, presets = NULL) {

  mega_mat    <- rlang::enquo(mega_mat)
  sub_kernels <- ipm$sub_kernels

  out         <- .make_mega_mat(mega_mat, sub_kernels)

  return(out)

}
