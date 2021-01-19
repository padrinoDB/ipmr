#' @noRd

sum.age_size_expr <- function(expr, ..., na.rm = NA){

  expr <- rlang::enexpr(expr)
  return(.all_ages(!! expr, fun = "+", ...))

}

#' @noRd

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


