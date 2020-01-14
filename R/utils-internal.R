# Reduces the depth/nestedness of a list to "depth".

.flatten_to_depth <- function(to_flatten, depth) {
  if(rlang::is_empty(to_flatten) |
     !rlang::is_list(to_flatten)) {
    return(to_flatten)
  }

  if(.depth(to_flatten) == depth) {
    return(to_flatten)
  } else {
    to_flatten <- purrr::flatten(to_flatten)
    .flatten_to_depth(to_flatten, depth)
  }
}

# Finds depth of current list

.depth <- function(l, start = 0) {
  if(!is.list(l)) {
    return(start)
  } else {
    max(
      unlist(
        lapply(l, .depth, start = start + 1)
      )
    )
  }
}

.drop_duplicated_names_and_splice <- function(to_drop) {

  if(!any(duplicated(names(to_drop)))) {

    test <- vapply(to_drop, function(x) all(is.na(x)), logical(1))

    if(!any(test)) {

      return(to_drop)

    } else {
      to_drop <- to_drop[!test]
      return(to_drop)
    }
  }

  temp   <- to_drop[rlang::have_name(temp)]
  temp_2 <- temp[!duplicated(names(temp))]

  ind    <- vapply(temp_2, is.list) %>%
    which()

  if(length(ind) > 0) {

    out         <- temp_2[-ind]

    for(i in ind) {

      to_splice <- temp_2[[i]]
      out       <- purrr::splice(out, to_splice)

    }

  } else {

    out         <- temp_2
  }

  return(out)

}

#' @noRd
# x = an pbject from make_ipm()

.det_lambda <- function(x) {

  return(
    vapply(x$iterators,
           function(y) Re(eigen(y)$values[1]),
           numeric(1L))
  )

}


# Lambda helpers----------
#' @noRd
# vec = numeric vector

.geom_mean <- function(vec) {
  n   <- length(vec)

  out <- prod(vec)

  return(out ^ (1/n))
}

#' @noRd

.lambda_pop_size <- function(x, all_lambdas = TRUE) {

  pops  <- x$pop_state
  n_its <- dim(pops[[1]])[2]
  temp  <- numeric(n_its - 1)

  for(i in seq(2, n_its, 1)) {

    tot_pop_size_t <- lapply(pops, function(x, it) {
      sum(x[ ,it])
    },
    it = (i - 1)) %>%
      unlist() %>%
      sum()

    tot_pop_size_t_1 <- lapply(pops, function(x, it) {
      sum(x[ ,it])
    },
    it = i) %>%
      unlist() %>%
      sum()

    temp[(i - 1)] <- tot_pop_size_t_1 / tot_pop_size_t

  }

  if(!all_lambdas) {
    return(temp[length(temp)])
  } else {
    return(temp)
  }

}

.lambda_eigen <- function(x) {

  eigs <- .det_lambda(x)

  return(eigs)


}

#' @noRd

is_square <- function(x) {

  dim(x)[1] == dim(x)[2]

}


#' @noRd
# Checks for convergence to asymptotic dynamics. Supports either
# eigenvalues or population vectors. For population vectors, assumes each column
# of the  matrix represents a single population vector

is_conv_to_asymptotic <- function(x, tol = 1e-7) {

  if(is.matrix(x)) {

    # Standardize columns first
    x <- apply(x, 2, FUN = function(y) y / sum(y))

    n_col     <- end_ind <- dim(x)[2]
    start_ind <- n_col - 1

    start_val <- x[ , start_ind]
    end_val   <- x[ , end_ind]

  } else if(is.vector(x)) {

    len       <- end_ind <- length(x)
    start_ind <- len - 1

    start_val <- x[start_ind]
    end_val   <- x[end_ind]

  }

  return(
    isTRUE(
      all.equal(
        start_val, end_val, tolerance = tol
      )
    )
  )

}
