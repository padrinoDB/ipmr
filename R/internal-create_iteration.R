#' @noRd

.init_iteration <- function(others, iterate, direction) {

  UseMethod(".init_iteration")

}

#' @noRd

.init_iteration.default <- function(others, iterate, direction) {

  if(!iterate) return(NA_character_)

  fun        <- .iter_fun(direction)

  all_states <- unique(unlist(others$state_var))

  start_end  <- .find_start_end(others)

  temp       <- .create_iter_exprs(others, start_end, all_states, fun)

  # if every state only has 1 expression operating on it, then we're done!
  # otherwise, we're going to have collapse the ones with a + to generate
  # a complete expression for it

  out <- .combine_iter_exprs(temp$out, temp$all_states)

  kernel_id <- paste(c("K", temp$all_names), collapse = "_")

  temp_proto <- define_k(
    others,
    name = kernel_id,
    family = "IPM",
    !!! out,
    data_list = list(),
    states = list(unique(unlist(others$state_var))),
    has_hier_effs = !is.null(temp$all_names),
    levels_hier_effs = temp$all_levels,
    evict_cor = FALSE,
    evict_fun = NULL,
    integrate = FALSE
  )

  ret <- temp_proto[nrow(temp_proto), ]

  return(ret)


}

#' @noRd

.find_start_end <- function(others) {

  start_end        <- others$domain
  names(start_end) <- others$kernel_id

  start_end <- lapply(start_end,
                      function(x) names(x))

  return(start_end)

}

#' @noRd

.create_iter_exprs <- function(others, start_end, all_states, fun) {

  UseMethod(".create_iter_exprs")

}

#' @noRd

.create_iter_exprs.default <- function(others, start_end, all_states, fun) {

  # Need to generate a deterministic simulation across all levels of hier_effs
  # if needed.

  out <- list()

  if(any(others$has_hier_effs)) {

    all_levels <- .flatten_to_depth(others$levels_hier_effs, 1L) %>%
      .[!duplicated(names(.))] %>%
      .[!is.na(.)]

    all_names <- paste(names(all_levels), collapse = "_")

    if("to_drop" %in% names(all_names)) {
      all_names <- all_names[all_names != "to_drop"]
    }


    if(grepl("_det$", class(others)[1])) {

      start_end <- lapply(start_end,
                          function(x, nms) {
                            paste(x, nms, sep = '_')
                          }, nms = all_names)

      all_states <- paste(all_states, all_names, sep = "_")

    }

  } else {

    all_levels <- list()
    all_names <- NULL

  }

  for(i in seq_along(others$kernel_id)) {

    end_state   <- paste("n_", start_end[[i]][2], "_t_1", sep = "")
    start_state <- paste("n_", start_end[[i]][1], "_t", sep = "")

    args <- list(kernel = rlang::sym(others$kernel_id[i]),
                 vectr  = rlang::ensym(start_state))

    iter_expr <- rlang::call2(fun, !!! args)

    out[[i]]      <- iter_expr
    names(out)[i] <- end_state

  }

  return(list(out = out,
              all_names = all_names,
              all_levels = all_levels,
              all_states = all_states))

}


.combine_iter_exprs <- function(temp, all_states) {

  UseMethod(".combine_iter_exprs")
}

#' @noRd

.combine_iter_exprs.default <- function(temp, all_states) {

  out <- list()

  if(length(temp) == length(all_states)) {

    out <- temp

  } else {

    for(i in seq_along(all_states)) {

      target_nm <- paste("n_", all_states[i], "_t_1", sep = "")

      to_collapse <- temp[names(temp) %in% target_nm] %>%
        lapply(rlang::expr_text)

      out[[i]] <- do.call("paste", list(to_collapse, collapse = " + ")) %>%
        rlang::parse_expr()

      names(out)[i] <- target_nm

    }
  }

  return(out)

}

#' @noRd

.combine_iter_exprs.age_x_size_exprs <- function(temp) {

  out <- list()

  temp       <- temp$out

  if(!anyDuplicated(names(temp))) {

    out <- temp

  } else {

    all_nms <- unique(names(temp))

    for(i in seq_along(all_nms)) {

      target_nm <- all_nms[i]

      to_collapse <- temp[names(temp) %in% target_nm] %>%
        lapply(rlang::expr_text)

      out[[i]] <- do.call("paste", list(to_collapse, collapse = " + ")) %>%
        rlang::parse_expr()

      names(out)[i] <- target_nm

    }
  }

  return(out)

}

#' @noRd

.iter_fun <- function(direction) {

  fun <- switch(direction,
                "left"  = "left_mult",
                "right" = "right_mult",
                stop("Can only perform left or right multiplication for model iteration!",
                     call. = FALSE)
                )

  return(fun)

}

.init_iteration.age_x_size <- function(others, iterate, direction) {

  if(!iterate) return(NA_character_)

  fun        <- .iter_fun(direction)

  all_states <- unique(unlist(others$state_var))

  start_end  <- .find_start_end(others)

  temp        <- .create_iter_exprs(others, start_end, all_states, fun)

  levels_ages <- .flatten_to_depth(others$levels_ages, 1L) %>%
    .[!duplicated(names(.))]

  temp        <- .handle_max_age(temp, levels_ages, start_end, fun)
  class(temp) <- "age_x_size_exprs"

  # if every state only has 1 expression operating on it, then we're done!
  # otherwise, we're going to have collapse the ones with a + to generate
  # a complete expression for it

  out <- .combine_iter_exprs(temp)

  kernel_id <- paste(c("K", temp$all_names), collapse = "_")

  temp_proto <- define_k(
    others,
    name = kernel_id,
    family = "IPM",
    !!! out,
    data_list = list(),
    states = list(unique(unlist(others$state_var))),
    has_hier_effs = !is.null(temp$all_names),
    levels_hier_effs = temp$all_levels,
    levels_ages = levels_ages,
    evict_cor = FALSE,
    evict_fun = NULL,
    integrate = FALSE
  )

  ret <- temp_proto[nrow(temp_proto), ]

  return(ret)

}

.handle_max_age <- function(expr_list, levels_ages, start_end, fun) {

  if(!"max_age" %in% names(levels_ages)) return(expr_list)

  # Generate a kernel name for the max age slot. This is the survival/growth
  # kernel (arbitrarily named) with "max_age_minus_1" tacked on. The trait distribution
  # it operates on is the

  kerns <- lapply(names(expr_list$out),
                  function(x) {
                    grepl("_age_", x)
                    }) %>%
    unlist() %>%
    expr_list$out[.] %>%
    lapply(function(y, start_end){
      args <- rlang::call_args(y)
      kern <- rlang::expr_text(args$kernel)
      sv   <- rlang::expr_text(args$vectr)

      kern <- c(gsub("age", "max_age", kern), # P_max_age_minus_1 %*% n_max_age_minus_1_t
               gsub("age_minus_1", "max_age", kern)) # P_max_age %*% n_max_age_t
      sv   <- c(gsub("age", "max_age", sv),
               gsub("age_minus_1", "max_age", sv))

      out <- list(kern = kern,
                  sv   = sv)
      return(out)
    },
    start_end = start_end)

  end_states <- gsub("_age_", "_max_age_", names(kerns))

  out <- list()

  it <- 1
  for(i in seq_along(kerns)) {

    if(length(kerns[[i]]$kern) != length(kerns[[i]]$sv)) {
      stop("Error generating iteration format for max_age!")
    }

    for(j in seq_along(kerns[[i]]$kern)) {

      args <- list(kernel = rlang::sym(kerns[[i]]$kern[j]),
                   vectr  = rlang::sym(kerns[[i]]$sv[j]))

      out[[it]] <- rlang::call2(fun, !!!args)

      names(out)[it] <- end_states[i]

      it <- it + 1

    }

  }

  expr_list$out <- c(expr_list$out, out)

  return(expr_list)

}

.create_iter_exprs.age_x_size <- function(others, start_end, all_states, fun) {

  # Need to generate a deterministic simulation across all levels of hier_effs
  # if needed.

  out <- list()

  if(any(others$has_hier_effs)) {

    all_levels <- .flatten_to_depth(others$levels_hier_effs, 1L) %>%
      .[!duplicated(names(.))] %>%
      .[!is.na(.)]

    all_names <- paste(names(all_levels), collapse = "_")

    if("to_drop" %in% names(all_names)) {
      all_names <- all_names[all_names != "to_drop"]
    }

    if(grepl("_det$", class(others)[1])) {

      start_end <- lapply(start_end,
                          function(x, nms) {
                            paste(x, nms, sep = '_')
                          }, nms = all_names)

      all_states <- paste(all_states, all_names, sep = "_")

    }

  } else {

    all_levels <- list()
    all_names <- NULL

  }

  # Things get a little hairy here. For any non-fecundity expression,
  # we want end_state_age_t_1 = P_age_minus_1 %*% start_state_minus_1_t.
  # The fecundity expressions will always have end_state_0_t_1, so we do not lag
  # those. The ones that do not have _0_ in the name require an age_minus_1 for
  # start state and the kernel name.

  for(i in seq_along(others$kernel_id)) {

    end_state   <- paste("n_", start_end[[i]][2], "_t_1", sep = "")

    if(!grepl("_0_", end_state)) {

      start_state <- paste("n_", start_end[[i]][1], "_t", sep = "") %>%
        gsub(pattern = "age", replacement = "age_minus_1", x = .)

      use_kernel  <- gsub("age", "age_minus_1", others$kernel_id[i])

      args <- list(kernel = rlang::sym(use_kernel),
                   vectr  = rlang::sym(start_state))

      iter_expr <- rlang::call2(fun, !!! args)

    } else {

      start_state <- paste("n_", start_end[[i]][1], "_t", sep = "")
      use_kernel  <- others$kernel_id[i]

      args <- list(kernel = rlang::sym(use_kernel),
                   vectr  = rlang::sym(start_state))

      iter_expr <- rlang::call2(fun, !!! args)

      iter_expr <- rlang::call2("sum", iter_expr)

    }

    out[[i]]      <- iter_expr
    names(out)[i] <- end_state

  }

  return(list(out = out,
              all_names = all_names,
              all_levels = all_levels,
              all_states = all_states))

}
