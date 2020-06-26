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
#' @param ... Named expressions. See Details for more information on their usage in
#' each \code{define_*} function.
#'
#' @param pop_vectors If the population vectors are already pre-defined (i.e. are
#' not defined by a function passed to \code{...}), then they can
#' be passed as a named list here.
#'
#' @details
#' These are helper functions to define IPMs. They are used after defining the kernels,
#' but before calling \code{make_ipm()} They are meant to be called in the
#' following order:
#'
#' \enumerate{
#'
#'   \item \code{define_impl()}
#'
#'   \item \code{define_domains()}
#'
#'   \item \code{define_pop_state()}
#'
#'   \item \code{define_env_state()}
#'
#' }
#'
#' The order requirement is so that information is correctly matched to each kernel.
#' Below are specific details on the way each one takes \code{...}.
#'
#' \strong{\code{define_impl}}
#'
#' This has two arguments - \code{proto_ipm} (the model object you wish to work with),
#' and the \code{kernel_impl_list}. The \code{kernel_impl_list} has a very specific format.
#' Names of the list should be kernel names, and each kernel should have 3 entries:
#' \code{int_rule}, \code{dom_start}, and \code{dom_end}. See examples.
#'
#' \strong{\code{define_domains}}
#'
#' If the \code{int_rule = "midpoint"}, these are vectors of length 3 where the name corresponds to the
#' state variable, the first entry is the lower bound of the domain, the second
#' is the upper bound of the domain, and the third entry is the number of
#' meshpoints. Other \code{int_rule}s are not yet implemented, so for now this is the
#' only format they can take. See examples.
#'
#' \strong{\code{define_pop_state}}
#'
#' This takes either calls to functions in the \code{...}, or a pre-generated
#' list of vectors in the \code{pop_vectors}. The names used
#' for each entry in \code{...} and/or for the \code{pop_vectors} should be
#' \code{n_<state_variable>}. See examples.
#'
#' \strong{\code{define_env_state}}
#'
#' Takes expressions that generate values for environmental covariates at each
#' iteration of the model in \code{...}. The \code{data_list} should contain any
#' parameters that the function uses, as well as the function itself. The
#' functions should return named lists. Names in that list can be referenced in
#' vital rate expressions and/or kernel formulas.
#'
#' @return All \code{define_*} functions return a proto_ipm. \code{make_impl_args_list}
#' returns a list, and so must be used within a call to \code{define_impl} or
#' before initiating the model creation procedure.
#'
#' @examples
#'
#' # Example with kernels named "P" and "F", and a domain "z"
#'
#' kernel_impl_list <- list(P = list(int_rule = "midpoint", dom_start = "z", dom_end = "z"),
#'                          F = list(int_rule = "midpoint", dom_start = "z", dom_end = "z"))
#'
#' # an equivalent version using make_impl_args_list
#'
#' kernel_impl_list <- make_impl_args_list(
#'      kernel_names = c("P", "F"),
#'      int_rule     = c("midpoint", "midpoint"),
#'      dom_start    = c("z", "z"),
#'      dom_end      = c("z", "z")
#' )
#'
#' # define_domains
#'
#' lower_bound <- 1
#' upper_bound <- 100
#' n_meshpoints <- 50
#'
#' \dontrun{
#'
#' define_domains(proto_ipm, c(lower_bound, upper_bound, n_meshpoints))
#'
#' # define_pop_state with a state variable named "z". Note that "n_" is prefixed
#' # to denote that it is a population state function!
#'
#' define_pop_state(proto_ipm, n_z = runif(100))
#'
#' # alternative, we can make a list before starting to make the IPM
#'
#' pop_vecs <- list(n_z = runif(100))
#'
#' define_pop_state(proto_ipm, pop_vectors = pop_vecs)
#'
#' # define_env_state. Generates a random draw from a known distribution
#' # of temperatures.
#'
#' env_sampler <- function(temp_mean, temp_sd) {
#'
#'   temp <- rnorm(1, temp_mean, temp_sd)
#'
#'   return(list(temp = temp))
#'
#' }
#'
#' env_pars <- list(temp_mean = 12, temp_sd = 2)
#'
#' define_env_state(
#'  proto_ipm,
#'  env_values = env_sampler(env_pars$temp_mean, env_pars$temp_sd),
#'  data_list  = list(env_sampler = env_sampler,
#'                    env_pars    = env_pars)
#'
#' )
#'
#' }
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

#' @inheritParams define_pop_state
#' @param data_list A list of named values that correspond to constants in the formula.
#' This behaves in the same way as the \code{data_list} argument to \code{define_kernel}.
#' @rdname define_star
#' @export

define_env_state <- function(proto_ipm, ..., data_list = list()) {

  env_quos            <- rlang::enquos(...)

  out                 <- list(env_quos = unlist(env_quos),
                              constants = data_list)

  proto_ipm$env_state <- list(out)

  return(proto_ipm)
}

#' @title Predict methods in ipmr
#' @rdname predict_methods
#'
#' @description This function is used when a \code{predict} method is incorporated
#' into the vital rate expressions of a kernel. Generally, ipmr can handle this
#' without any additional user effort, but some model classes will fail (often
#' with an obscure error message).
#' When this happens, \code{use_vr_model} can ensure that model object is
#' correctly represented in the \code{data_list}.
#'
#' @param model A fitted model object representing a vital rate. Primarily used to avoid
#' writing the mathematical expression for a vital rate, and using a \code{predict()}
#' method instead.
#'
#' @return A model object with a \code{"flat_protect"} attribute.
#'
#' @details ipmr usually recognizes model objects passed into the \code{data_list} argument
#' automatically. Unfortunately, sometimes it'll miss one, and the user will need
#' to manually protect it from the standard build process. This function
#' provides a wrapper around that process.
#'
#' Wrap a model object in \code{use_vr_model} when building the \code{data_list}
#' to pass to \code{define_kernel}.
#'
#'
#' @examples
#'
#' \dontrun{
#'
#' grow_mod <- lm(size_next ~ size, data = grow_data)
#' surv_mod <- glm(surv ~ size, data = surv_data, family = binomial())
#'
#' data_list <- list(
#'   grow_mod = use_vr_model(grow_mod),
#'   surv_mod = use_vr_model(surv_mod),
#'   recruit_mean = 20,
#'   recruit_sd   = 5
#' )
#'}
#' @export

use_vr_model <- function(model) {


  attr(model, "flat_protect") <- TRUE

  return(model)

}


#' @title Right multiplication
#'
#' @description Performs right multiplication.
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

  init   <- diag(id_dim) # Identity matrix as starting point

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
#' They should be specified in ROW MAJOR order! Can also be a character
#' string specifying the call. See examples.
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
#' char_call <- "c(0, go_discrete, leave_discrete, P)"
#'
#' big_k <- format_mega_matrix(gen_di_det_ex, mega_mat = char_call)
#'
#'
#'
#' @export

format_mega_matrix <- function(ipm, mega_mat, presets = NULL) {

  mega_mat    <- rlang::enquo(mega_mat)

  if(!rlang::quo_is_call(mega_mat)) {

    text     <- rlang::eval_tidy(mega_mat)
    exprr    <- rlang::parse_expr(text)
    mega_mat <- rlang::enquo(exprr)

  }

  sub_kernels <- ipm$sub_kernels

  out         <- .make_mega_mat(mega_mat, sub_kernels)

  return(out)

}

#' @title Age X Size models
#' @rdname age_x_size
#'
#' @description Helper functions for implementing age x size models
#'
#' @param expr The expression to modify. See details.
#' @param fun  The function to apply over \code{expr}.
#' @param ... Used internally. Modifying this will likely cause models to break.
#'
#' @return A modified expression.
#'
#' @details Currently only \code{all_ages} is exported. This is a helper
#' to indicate that all of the levels of age should be included in a given expression,
#' rather than split out into separate levels. The \code{fun} is used to combine
#' the results of each individual evaluation of \code{expr}.
#'
#' The most common use case for \code{all_ages} will be when defining the size
#' distribution of new recruits in an age structured model, when many \code{F} kernels
#' are multiplied by their respective age X size distributions, and then summed.
#' Note that using \code{sum} is almost never correct, as this will always return
#' a single number whereas the size distribution should be a vector! Generally,
#' use \code{+} instead.
#'
#' @export

all_ages <- function(expr, fun, ...) {

  ages <- list(...) %>%
    unlist()

  new_text <- character(length(ages))

  cur_text  <- deparse(substitute(expr), width.cutoff = 500L)

  for(i in seq_along(ages)) {

    new_text[i] <- gsub("age", ages[i], cur_text)

  }

  # add spaces around "fun" so output is more human readable. mostly for
  # my own debugging purposes.
  out_text <- paste(new_text, collapse = paste(" ",
                                               as.character(fun),
                                               " ",
                                               sep = ""))

  # out      <- rlang::parse_expr(out_text)

  return(out_text)

}
