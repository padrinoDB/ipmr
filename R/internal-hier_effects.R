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

  # Place holder, this will ultimate get rbind'ed
  new_proto <- rows[0, ]

  for(i in seq_len(dim(rows)[1])) {

    # Unadulterated garbage, but it works. Creates a row for every level/comibnation
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

      proto$params[[it]]$formula        <- gsub(nm,
                                                levels[k, j],
                                                proto$params[[it]]$formula)

      proto$params[[it]]$vr_text        <- purrr::map(proto$params[[it]]$vr_text,
                                                      .f = function(x, level, nm) {
                                                        gsub(nm, level, x)
                                                      },
                                                      level = levels[k,j],
                                                      nm = nm)

      names(proto$params[[it]]$vr_text) <- purrr::map_chr(names(proto$params[[it]]$vr_text),
                                                          .f = function(x, level, nm) {
                                                            gsub(nm, level, x)
                                                          },
                                                          level = levels[k, j],
                                                          nm = nm)


      temp <- rlang::quo_text(proto$evict_fun[[it]][[1]]) %>%
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
