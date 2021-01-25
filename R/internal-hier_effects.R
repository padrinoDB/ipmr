# Internal hierarchical effect functions

.split_hier_effs <- function(proto_ipm) {

  kerns <- which(proto_ipm$has_hier_effs)

  # Create a place to hold the output - either kernels with no hierarchical effects
  # or a new proto

  if(length(kerns) != dim(proto_ipm)[1]) {

    out <- proto_ipm[-kerns, ]

  } else {

    out <- .init_ipm(class(proto_ipm)[1])

  }

  hier_rows <- proto_ipm[kerns, ]

  all_hier_effs <- .flatten_to_depth(hier_rows$levels_hier_effs, 1L)


  for(i in seq_len(dim(hier_rows)[1])) {

    levs   <- hier_rows$levels_hier_effs[[i]]

    if("drop_levels" %in% names(levs)) {

      ind <- which(names(levs) != "drop_levels")

      levs <- levs[ind]

    }

    levels <- lapply(levs, eval) %>%
      expand.grid(stringsAsFactors = FALSE)

    temp   <- .expand_hier_effs(hier_rows[i, ], levels)

    out    <- rbind(out, temp)

  }


  if("drop_levels" %in% names(all_hier_effs)) {

    out <- .drop_levels_hier_effs(out, all_hier_effs)

  }

  return(out)
}

#'@noRd

.drop_levels_hier_effs <- function(proto_ipm, all_hier_effs) {

  to_drop <- all_hier_effs[!duplicated(names(all_hier_effs))]
  to_drop <- to_drop$drop_levels

  for(i in seq_along(to_drop)) {

    proto_ipm <- proto_ipm[!grepl(to_drop[i], proto_ipm$kernel_id), ]

  }

  return(proto_ipm)

}

.expand_hier_effs <- function(rows, levels) {

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


      # Case when multiple hierarchical effects exist. The first iteration peels
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

.make_hier_levels <- function(hier_effs) {

  hier_effs <- .flatten_to_depth(hier_effs, 1L) %>%
    .[!names(.) == "levels"]

  if(length(hier_effs) == 1) {

    levs <- as.data.frame(unlist(hier_effs),
                          stringsAsFactors = FALSE)
    drop <- FALSE

  } else {

    if("drop_levels" %in% names(hier_effs)) {

      to_drop <- hier_effs$drop_levels
      hier_effs <- hier_effs[-c("drop_levels")]
      drop <- TRUE

    } else {

      drop <- FALSE

    }

    levs <- lapply(hier_effs, eval) %>%
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
