
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

  vapply(x$iterators,
         function(y) Re(eigen(y)$values[1]),
         numeric(1))

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

.stoch_lambda_pop_size <- function(x, all_lambdas = FALSE) {

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

.stoch_lambda_eigen <- function(x, all_lambdas = FALSE) {

  eigs <- .det_lambda(x)

  if(!all_lambdas) {
    return(.geom_mean(eigs))
  } else {
    return(eigs)
  }

}
