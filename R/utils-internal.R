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


#' @noRd
# Checks if model has already been iterated, saving us from re-iterating a model
# when we can just use the pop_state slot

.already_iterated <- function(ipm) {

  return(
    ! all(
      is.na(
        ipm$pop_state
      )
    )
  )


}
