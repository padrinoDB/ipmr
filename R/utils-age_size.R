
#' @title Age X Size models
#' @rdname age_x_size
#'
#' @description Helper functions for implementing age x size models. These are
#' only intended for use within \code{define_kernel()}.
#'
#' @param expr An expression to compute a sum for.
#' @param ... Used internally. Modifying this will likely cause models to break.
#' @param na.rm Not used and ignored.
#'
#' @return A modified expression.
#'
#' @details Note that these functions do not compute a single number, but rather
#' expand an expression to compute the trait distribution. This differs from the
#' usual behavior of the numeric methods for \code{sum, prod} as those
#' return a single number or \code{NA}.
#'

sum.age_size_expr <- function(expr, ..., na.rm = NA){

  expr <- rlang::enexpr(expr)
  return(.all_ages(!! expr, fun = "+", ...))

}

#' @rdname age_x_size

prod.age_size_expr <- function(expr, ..., na.rm = NA){

  expr <- rlang::enexpr(expr)
  return(.all_ages(!! expr, fun = "*", ...))

}

#' @noRd

.all_ages <- function(expr, fun, ...) {

  expr <- rlang::enexpr(expr)

  ages <- list(...) %>%
    unlist()

  new_text <- character(length(ages))

  for(i in seq_along(ages)) {

    new_text[i] <- gsub("age", ages[i], expr)

  }

  # add spaces around "fun" so output is more human readable. mostly for
  # my own debugging purposes.
  # NB: cosndier updating to use call_modify instead of literal
  # subbing + expansion on the expression

  out_text <- paste(new_text, collapse = paste(" ",
                                               as.character(fun),
                                               " ",
                                               sep = ""))

  return(out_text)

}


