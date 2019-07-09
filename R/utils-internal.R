
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
    return(to_drop)
  }

  temp   <- to_drop[rlang::have_name(temp)]
  temp_2 <- temp[!duplicated(names(temp))]
  ind <- vapply(temp_2, is.list) %>%
    which()

  if(length(ind) > 0) {
    out <- temp_2[-ind]

    for(i in ind) {
      to_splice <- temp_2[[i]]
      out <- purrr::splice(out, to_splice)
    }
  } else {
    out <- temp_2
  }

  return(out)


}
