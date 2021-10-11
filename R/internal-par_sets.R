# Internal par_setarchical effect functions

.split_par_sets <- function(proto_ipm) {

  kerns <- which(proto_ipm$uses_par_sets)

  # Create a place to hold the output - either kernels with no parameter sets
  # or a new proto

  if(length(kerns) != dim(proto_ipm)[1]) {

    out <- proto_ipm[-kerns, ]

  } else {

    out <- .init_ipm(class(proto_ipm)[1])

  }

  par_set_rows <- proto_ipm[kerns, ]

  all_par_sets <- .flatten_to_depth(par_set_rows$par_set_indices, 1L)


  for(i in seq_len(dim(par_set_rows)[1])) {

    levs   <- par_set_rows$par_set_indices[[i]]

    if("drop_levels" %in% names(levs)) {

      ind <- which(names(levs) != "drop_levels")

      levs <- levs[ind]

    }

    levels <- lapply(levs, eval) %>%
      expand.grid(stringsAsFactors = FALSE)

    temp   <- .expand_par_sets(par_set_rows[i, ], levels)

    out    <- rbind(out, temp)

  }


  if("drop_levels" %in% names(all_par_sets)) {

    out <- .drop_par_set_indices(out, all_par_sets)

  }

  return(out)
}

#'@noRd

.drop_par_set_indices <- function(proto_ipm, all_par_sets) {

  to_drop <- all_par_sets[!duplicated(names(all_par_sets))]
  to_drop <- to_drop$drop_levels

  for(i in seq_along(to_drop)) {

    proto_ipm <- proto_ipm[!grepl(to_drop[i], proto_ipm$kernel_id), ]

  }

  return(proto_ipm)

}

.expand_par_sets <- function(rows, levels) {

  # Place holder, this will ultimately get rbind'ed
  new_proto <- rows[0, ]

  for(i in seq_len(dim(rows)[1])) {

    # Pure garbage, but it works. Creates a row for every level/comibnation
    # of levels, next we substitute in everything

    temp      <- .expand_proto(rows, levels, i)

    for_bind  <- .sub_levels(temp, levels)

    new_proto <- rbind(new_proto, for_bind)

  }

  return(new_proto)
}

.sub_levels <- function(proto, levels) {

  it <- 1

  for(j in seq_len(dim(levels)[2])) { # variable names loop

    nm <- names(levels)[j]

    for(k in seq_len(dim(levels)[1])) { # variable values loop

      proto[it, 'kernel_id']            <- gsub(nm,
                                                levels[k, j],
                                                proto[it, 'kernel_id'])

      # Kernels with multiple population vectors require a named list of
      # `formula`, whereas single population vectors will just be a scalar
      # that is then turned back into a later on

      if(is.list(proto$params[[it]]$formula)) {
        proto$params[[it]]$formula <- purrr::map(
          proto$params[[it]]$formula,
          .f = function(x, level, nm) {

            gsub(nm, level, x)

          },
          level = levels[k,j],
          nm    = nm
        )

        names(proto$params[[it]]$formula) <- purrr::map_chr(
          names(proto$params[[it]]$formula),
                .f = function(x, level, nm) {

                  gsub(nm, level, x)

                },
                level = levels[k,j],
                nm = nm

        )

      } else {

        proto$params[[it]]$formula        <- gsub(
          nm,
          levels[k, j],
          proto$params[[it]]$formula
        )

      }

      proto$params[[it]]$vr_text        <- purrr::map(
        proto$params[[it]]$vr_text,
        .f = function(x, level, nm) {

          gsub(nm, level, x)

        },
        level = levels[k,j],
        nm = nm
      )


      names(proto$pop_state[[it]])      <- purrr::map_chr(
        names(proto$pop_state[[it]]),
        .f = function(x, level, nm) {
          gsub(nm, level, x)
        },
        level = levels[k, j],
        nm = nm
      )


      # We want a character VECTOR of names, not a list. Hence the map_chr instead
      # of standard purrr::map here

      names(proto$params[[it]]$vr_text) <- purrr::map_chr(
        names(proto$params[[it]]$vr_text),
        .f = function(x, level, nm) {
          gsub(nm, level, x)
        },
        level = levels[k, j],
        nm = nm
      )


      # Case when multiple parameter sets exist. The first iteration peels
      # off a layer of list from the quosure, so we need to make sure we don't
      # grab that at the second time of asking.

      if(j > 1) {
        use_ev_fun <- proto$evict_fun[[it]]
      } else {
        use_ev_fun <- proto$evict_fun[[it]][[1]]
      }

      temp <- rlang::quo_text(use_ev_fun) %>%
        gsub(pattern = nm, replacement = levels[k, j], x = .) %>%
        rlang::parse_expr()


      proto$evict_fun[[it]]   <-  rlang::enquo(temp)

      names(proto$params)[it] <- gsub(nm,
                                      levels[k, j],
                                      names(proto$params)[it])


      if(it == dim(levels)[1]){
        it <- 1
      } else {
        it <- it + 1
      }

    } # end var names loop
  } # end proto loop

  return(proto)
}


.expand_proto <- function(rows, levels, i) {

  do.call("rbind",
          replicate(dim(levels)[1],
                    rows[i, ],
                    simplify = FALSE))

}

.make_par_set_indices <- function(par_sets) {

  par_sets <- .flatten_to_depth(par_sets, 1L) %>%
    .[names(.) != "levels"]

  par_sets <- par_sets[!duplicated(names(par_sets))]

  if(length(par_sets) == 1) {

    levs <- as.data.frame(unlist(par_sets),
                          stringsAsFactors = FALSE)
    drop <- FALSE

  } else {

    if("drop_levels" %in% names(par_sets)) {

      to_drop <- par_sets$drop_levels
      par_sets <- par_sets[!names(par_sets) %in% "drop_levels"]
      drop <- TRUE

    } else {

      drop <- FALSE

    }

    levs <- lapply(par_sets, eval) %>%
      expand.grid(stringsAsFactors = FALSE)

  }

  out <- character(length = dim(levs)[1])

  for(i in seq_along(out)) {

    out[i] <- paste(levs[i, ], collapse = "_")

  }

  if(drop){
    out <- out[!out %in% to_drop]
  }

  return(out)

}
