# Internal hierarchical effect functions



.split_hier_effs <- function(proto_ipm) {

  kerns <- which(proto_ipm$has_hier_effs)

  # Create a place to hold the output - either kernels with no hierarchical effects
  # or a new proto
  if(length(kerns) != dim(proto_ipm)[1]) {

    out <- proto_ipm[-kerns, ]

  } else {

    out <- init_ipm(class(proto_ipm)[1])

  }

  hier_rows <- proto_ipm[kerns, ]

  for(i in seq_len(dim(hier_rows)[1])) {

    levs   <- hier_rows$levels_hier_effs[[i]]

    levels <- lapply(levs, eval) %>%
      expand.grid(stringsAsFactors = FALSE)

    temp   <- .expand_hier_effs(hier_rows[i, ], levels)

    out    <- rbind(out, temp)

  }

  return(out)
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

      # So far, the only things I can identify that require subbing is
      # the kernel_id, kernel formula, and vital rate exprs. I fear I will
      # come to regret this comment

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

  if(length(hier_effs) == 1) {

    levs <- as.data.frame(unlist(hier_effs),
                          stringsAsFactors = FALSE)

  } else {

    levs <- lapply(hier_effs, eval) %>%
      expand.grid(stringsAsFactors = FALSE)

  }

  out <- character(length = dim(levs)[1])

  for(i in seq_along(out)) {

    out[i] <- paste(levs[i, ], collapse = "_")
  }

  return(out)

}
