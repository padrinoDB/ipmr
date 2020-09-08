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
#' @param mega_mat A vector with symbols, I's, and/or 0s representing the matrix blocks.
#' They should be specified in ROW MAJOR order! Can also be a character
#' string specifying the call. Hierarchical syntax is supported. When used,
#' \code{format_mega_matrix} will produce as many mega-matrices as there are
#' combinations of \code{levels_hier_effs} in the \code{proto_ipm}.
#'
#' @return A list containing a large matrix or many large matrices (when used with
#' hierarchical syntax). The names in the former case will be \code{"mega_matrix"}
#' and in the latter case, the level of the hierarchical effect.
#'
#' @details \code{I} and \code{0} represent identity matrices and 0 matrices,
#' respectively. They can be used to fill in blocks that represent either, without
#' having to create those separately and append them to the model object. The function
#' will work out the correct dimensions for both internally, and there is no
#' restriction on the number that may be used in a given call.
#'
#' @examples
#' data(gen_di_det_ex)
#'
#' big_k <- format_mega_matrix(gen_di_det_ex,
#'                             mega_mat = c(0, go_discrete,
#'                                          leave_discrete, P))
#'
#' char_call <- c(0, "go_discrete", "leave_discrete", "P")
#'
#' big_k <- format_mega_matrix(gen_di_det_ex, mega_mat = char_call)
#'
#' # Now, with an Identity matrix instead of a 0
#'
#' big_k <- format_mega_matrix(gen_di_det_ex,
#'                             mega_mat = c(I, go_discrete,
#'                                          leave_discrete, P))
#'
#'
#'
#' @export

format_mega_matrix <- function(ipm, mega_mat) {

  mega_mat    <- rlang::enquo(mega_mat)

  if(!rlang::quo_is_call(mega_mat)) {

    text     <- rlang::eval_tidy(mega_mat)

    if(length(text) > 1) {

      exprr <- syms(text)
      exprr <- rlang::call2("c", !!! exprr)

    } else{

      exprr    <- rlang::parse_expr(text)

    }

    mega_mat <- rlang::enquo(exprr)

  }

  if(any(ipm$proto_ipm$has_hier_effs)) {

    levs <- ipm$proto_ipm$levels_hier_effs[ipm$proto_ipm$has_hier_effs] %>%
      .flatten_to_depth(1L) %>%
      .[!duplicated(names(.))]

    if("drop_levels" %in% names(levs)) {

      ind <- which(names(levs) != "drop_levels")

      use_levs <- levs[ind]

    } else {

      use_levs <- levs

    }

    use_levs <- expand.grid(use_levs, stringsAsFactors = FALSE)

    out_nms <- temp <- character(dim(use_levs)[1])

    base_expr <- rlang::quo_text(mega_mat)
    base_name <- names(use_levs) %>% paste(collapse = "_")

    it <- 1

    for(i in seq_len(dim(use_levs)[2])) {

      nm <- names(use_levs)[i]

      for(j in seq_len(dim(use_levs)[1])) {

        if(i > 1) {
          base_expr <- temp[it]
          base_name <- out_nms[it]
        }

        temp[it] <- gsub(nm, use_levs[j, i], base_expr)
        out_nms[it] <- gsub(nm, use_levs[j, i], base_name)

        it <- it + 1

      } # end single var substitution - reset counter to modify exprs w/ adtl hier_effs if present

      it <- 1

    }

    if("drop_levels" %in% names(levs)) {

      # Need to use fuzzy matching for temp because the level is already
      # appended to each kernel name. Don't have the same problem for
      # out_nms, as that is just vector with the exact levels

      for(i in seq_along(levs$drop_levels)) {

        temp    <- temp[!grepl(levs$drop_levels[i], temp)]
      }

      out_nms <- out_nms[!out_nms %in% levs$drop_levels]

    }


    mega_mat <- as.list(temp) %>%
      lapply(function(x) {
        y <- rlang::parse_expr(x)
        rlang::enquo(y)
      })



  } else {

    mega_mat <- list(mega_mat)
    out_nms  <- "mega_matrix"

  }

  out         <- lapply(mega_mat,
                        function(x, ipm) {
                          .make_mega_mat(ipm, x)
                        },
                        ipm = ipm)

  names(out) <- out_nms

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
#' a single number whereas we typically want the trait distribution function
#' (i.e. a vector). It is usually best to use \code{+} instead.
#'
#' @examples
#'
#' \dontrun{
#'
#'define_k(
#'  proto_ipm          = some_age_size_proto_ipm,
#'  name               = "K",
#'  family             = "IPM",
#'  n_wt_0_t_1         = all_ages(expr = F_age %*% n_wt_age_t, fun = "+"),
#'  n_wt_age_t_1       = P_age_minus_1 %*% n_wt_age_minus_1_t,
#'  n_wt_max_age_t_1   = P_max_age %*% n_wt_max_age_t +
#'                       P_max_age_minus_1 %*% n_wt_max_age_minus_1_t,
#'  data_list          = param_list,
#'  states             = list (c("wt")),
#'  has_hier_effs      = FALSE,
#'  levels_ages        = list(age = c(0:20), max_age = 21),
#'  evict_cor          = FALSE
#')
#'
#'
#' }
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


#' @title Accessor functions for proto_ipm objects
#' @rdname proto_accessors
#'
#' @description Functions that access specific slots of a \code{proto_ipm} for
#' user examination
#'
#' @param proto_ipm A \code{proto_ipm} object.
#'
#' @return A named list of either numbers (for domains), or bare expressions (for
#' vital rates and kernel formulae).
#'
#' @export

domains <- function(proto_ipm) {

  out <- lapply(proto_ipm$domain, function(x) x)

  out <- lapply(out,
                function(x) {

                  temp <- lapply(x, function(y) {
                    if(!all(is.na(y))) {
                      names(y) <- c("lower_bound",
                                    "upper_bound",
                                    "n_meshpoints")

                    }

                    return(y)
                  }
                  )

                  return(temp)
                }
  ) %>%
    .flatten_to_depth(1L) %>%
    Filter(f = Negate(is.na), x = .) %>%
    .[!duplicated(names(.)) & !is.na(names(.))]

  return(out)

}

#' @rdname proto_accessors
#' @importFrom stats setNames
#' @export

vital_rate_functions <- function(proto_ipm) {

  out <- lapply(proto_ipm$params, function(x) x$vr_text) %>%
    stats::setNames(c("")) %>%
    lapply(function(x)
           if(any(is.na(x) | is.null(x) | rlang::is_empty(x))) {
             return(NULL)
           }  else {
             return(x)
           }
    ) %>%
    Filter(Negate(is.null), x = .) %>%
    Filter(Negate(rlang::is_empty), x = .) %>%
    .flatten_to_depth(1L) %>%
    lapply(rlang::parse_expr)

  out <- out[!duplicated(names(out))]

  return(out)
}

#' @rdname proto_accessors
#' @export

kernel_formulae <- function(proto_ipm) {

  out <- lapply(proto_ipm$params, function(x) x$formula) %>%
    .flatten_to_depth(1L) %>%
    Filter(f = Negate(is.na), x = .) %>%
    lapply(rlang::parse_expr)

  return(out)
}
