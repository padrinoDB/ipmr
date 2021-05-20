# Age x size internal functions
# These work by evaluating the _minus_ and _plus_ keywords in each expression's
# name/content

# These should only appear in .make_age_k_row(...). "others" object is treated
# just like any other par_setarchical model, for now. This will require modification
# for more complicated models like density dependent and/or time lagged ones
# where vital rates in the sub-kernels themselves definitely will require the
# t_minus_1 or n_minus_1 notation.

#' @noRd

`%minus%` <- function(a, b) {

  return(a - b)

}

#' @noRd

`%plus%` <- function(a, b) {

  return(a + b)

}

#' @noRd

.expand_age_kern_forms <- function(x, ages) {

  if(is.character(x)) x <- rlang::parse_expr(x)

  temp     <- rlang::call_standardise(x)
  new_arg  <- rlang::call_args(temp)[[1]]
  new_arg  <- rlang::expr_text(new_arg)
  class(new_arg) <- c("age_size_expr", "character")
  temp[[2]] <- NULL

  new_call <- rlang::call_modify(temp,
                                 expr = new_arg,
                                 ages = ages)

  out <- rlang::eval_tidy(new_call)

  return(out)
}

.make_age_k_row <- function(k_row) {

  forms <- k_row$params[[1]]$formula

  ages  <- k_row$age_indices %>%
    .flatten_to_depth(1L)

  calls <- vapply(forms,
                  function(x) {

                    if(is.character(x)) x <- rlang::parse_expr(x)
                    rlang::call_name(x)

                  }, character(1L))

  # In the case where we have an n_0 = f_1 %*% n_1 + f_2 %*% n_2 etc,
  # We need to treat this differently from other par_sets. Normally, they
  # get one call per level of the par_set_eff. in this case, we need one call
  # expanded to include all levels of the par_set_eff. This is done w/ methods
  # for generics sum and prod (and maybe others later).
  # On the other hand, we also need to generate one expression
  # for each n_x_t_1 = p_x_minus_1 %*% n_x_minus_1_t

  a_s_methods <- c("sum", "prod")

  if(any(calls %in% a_s_methods)) {

    n_0_ind  <- which(calls %in% a_s_methods)
    n_0_call <- forms[n_0_ind]

    n_0_form <- lapply(
      n_0_call,
      function(x, ages) {
        .expand_age_kern_forms(x, ages)

      },
      ages = ages
    )

    reg_calls <- forms[-c(n_0_ind)]

    # We do not want to create an expression for 0 if we already have one

    use_ages <- unlist(ages) %>% .[. != 0]

  } else {

    reg_calls <- forms
    n_0_form <- NULL

    use_ages <- unlist(ages)

  }

  reg_holder <- vector("list",
                       length = length(use_ages))

  # Split ages into max_age if present. We need to modify this first because
  # the fuzzy matching w/ gsub will mess up "max_age" before it reaches that
  # variable name.

  # Perhaps this is a case for more precise matching of names/values.
  # Get implemented first, then worry about this.

  if("max_age" %in% names(ages)) {

    max_age  <- ages$max_age
    use_ages <- use_ages[use_ages != max_age]

    max_age_call <- reg_calls[grepl("max_age", names(reg_calls))]

    temp         <- lapply(max_age_call,
                           function(x, max_age) {
                             place_holder <- gsub("max_age", max_age, x)
                             out          <- .update_age(place_holder)
                             return(out)
                           },
                           max_age = max_age)

    reg_holder[length(reg_holder)]        <- unlist(temp, use.names = FALSE)
    names(reg_holder)[length(reg_holder)] <- gsub("max_age", max_age, names(temp))

    reg_calls <- reg_calls[!grepl("max_age", names(reg_calls))]
  }

  for(i in seq_along(use_ages)) {

    temp <- lapply(reg_calls,
                   function(x, use_age) {

                     place_holder <- gsub("age", use_age, x)
                     out          <- .update_age(place_holder)
                     return(out)

                   },
                   use_age = use_ages[i])

    reg_holder[[i]]   <- unlist(temp, use.names = FALSE)
    names(reg_holder)[i] <- gsub("age", use_ages[i], names(temp))

  }

  k_row$params[[1]]$formula <- c(reg_holder, n_0_form)

  return(k_row)

}

.update_age <- function(txt) {

  # Max width in deparse is apparently 500L. I don't think inserting \n's will
  # make a difference come eval time, but using this to be safe.

  if(!is.character(txt)) txt <- deparse(substitute(txt), width.cutoff = 500L)

  # For now, only supporting age + x, age - x. We can add more if there's demand,
  # but I have no idea why anyone would want them.

  supported_funs <- c("minus", "plus")

  reg_expr       <- paste("[0-9]+", supported_funs, "[0-9]+",
                          sep = "_")

  inds <- vapply(reg_expr,
                 function(regs, txt) {
                   temp <- regexpr(regs, txt)
                   start <- temp[1]
                   end   <- start + attr(temp, "match.length")

                   out <- c(start, (end - 1)) %>% as.integer()

                   return(out)
                 },
                 integer(2L),
                 txt = txt)

  inds <- apply(inds,
                2,
                function(x) {
                  if(x[1] == -1 && x[2] == -3) {
                    x <- NULL
                  }
                  return(x)
                })

  ret_lgl <- vapply(inds, function(x) any(is.null(x)), logical(1L))

  if(all(ret_lgl)) return(txt)

  names(inds) <- supported_funs

  inds <- Filter(Negate(is.null), inds)

  # In this case, there can only be two substitutions (e.g. plus and/or minus)
  # step 1 pulls out the age_minus_n part. Part 2 replaces underscores with
  # infix notation and evaluates it. part 3 replaces the age_minus_n part
  # with the evaluated form (e.g. 4 instead of 5_minus_1)

  for(i in seq_along(inds)){

    rep_text   <- substr(txt, start = inds[[i]][1], stop = inds[[i]][2])
    to_replace <- gsub("_", "%", rep_text)
    val        <- eval(parse(text = to_replace))

    txt      <- gsub(rep_text, val, txt)

  }

  return(txt)
}


#' @noRd
# Basically .split_par_sets, but uses the age column instead of par_sets column

.split_sub_kern_ages <- function(proto_ipm) {

  kerns <- which(proto_ipm$uses_age)

  # Create a place to hold the output - either kernels with no parameter sets
  # or a new proto
  if(length(kerns) != dim(proto_ipm)[1]) {

    out <- proto_ipm[-kerns, ]

  } else {

    out <- .init_ipm(class(proto_ipm)[1])

  }

  par_set_rows <- proto_ipm[kerns, ]

  dup_test  <- gsub("(_age)|(_max_age)", "", par_set_rows$kernel_id)

  is_duped  <- duplicated(dup_test) | duplicated(dup_test, fromLast = TRUE)

  for(i in seq_len(dim(par_set_rows)[1])) {

    levs   <- par_set_rows$age_indices[[i]]

    levs   <- list(age = unlist(levs, use.names = FALSE))

    levels <- lapply(levs, eval)

    temp   <- .expand_age_sub_kerns(par_set_rows[i, ], levels, is_duped[i])

    out    <- rbind(out, temp)

  }

  return(out)

}


.expand_age_sub_kerns <- function(rows, ages, is_duped) {

  # Place holder, this will ultimately get rbind'ed
  new_proto <- rows[0, ]

  for(i in seq_len(dim(rows)[1])) {

    targ_val  <- .get_age_targets_values(rows[i, ], ages, is_duped)

    temp      <- .expand_proto(rows, data.frame(age = targ_val$use_ages), i)

    for_bind  <- .sub_age(temp, targ_val$use_ages, targ_val$target)

    new_proto <- rbind(new_proto, for_bind)

  }

  return(new_proto)
}

#' @noRd
# Function to partition age values correctly across kernels. In the case
# where P_a is defined by P_a AND P_max_a, we want to make sure P_a is defined
# for a < max_a, and P_max_a for a == max_a. On the other hand, if F_a is defined
# identically for all a, then that should use all a for substitution.

.get_age_targets_values <- function(rows, ages, is_duped) {


  if(grepl("max_age", rows$kernel_id)){

    use_ages <- max(ages$age)

    target <- "max_age"

  } else if(!grepl("max_age", rows$kernel_id) & is_duped) {

    use_ages <- ages$age[-which.max(ages$age)]

    target   <- "age"

  } else {

    use_ages <- ages$age

    target   <- "age"

  }

  return(list(use_ages = use_ages,
              target   = target))
}

.sub_age <- function(proto, use_ages, target) {

  it <- 1

  # the target in the vital rate expressions should be either "age" or "max_age".
  # however, pop state should  really always be "age".

  nm    <- target
  ps_nm <- "age"

  for(i in seq_along(use_ages)) { # variable names loop

    proto[it, 'kernel_id']  <- gsub(nm,
                                    use_ages[i],
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
          level = use_ages[i],
          nm    = nm
        )

        names(proto$params[[it]]$formula) <- purrr::map_chr(
          names(proto$params[[it]]$formula),
          .f = function(x, level, nm) {

            gsub(nm, level, x)

          },
          level = use_ages[i],
          nm = nm

        )

      } else {

        proto$params[[it]]$formula <- gsub(
          nm,
          use_ages[i],
          proto$params[[it]]$formula
        )

      }

      proto$params[[it]]$vr_text <- purrr::map(
        proto$params[[it]]$vr_text,
        .f = function(x, level, nm) {

          gsub(nm, level, x)

        },
        level = use_ages[i],
        nm = nm
      )


      names(proto$pop_state[[it]])      <- purrr::map_chr(
        names(proto$pop_state[[it]]),
        .f = function(x, level, nm) {
          gsub(nm, level, x)
        },
        level = use_ages[i],
        nm = ps_nm
      )


      # We want a character VECTOR of names, not a list. Hence the map_chr instead
      # of standard purrr::map here

      names(proto$params[[it]]$vr_text) <- purrr::map_chr(
        names(proto$params[[it]]$vr_text),
        .f = function(x, level, nm) {
          gsub(nm, level, x)
        },
        level = use_ages[i],
        nm = nm
      )


      # Case when multiple parameter sets exist. The first iteration peels
      # off a layer of list from the quosure, so we need to make sure we don't
      # grab that at the second time of asking.

      if(i > 1) {
        use_ev_fun <- proto$evict_fun[[it]]
      } else {
        use_ev_fun <- proto$evict_fun[[it]][[1]]
      }

      temp <- rlang::quo_text(use_ev_fun) %>%
        gsub(pattern = nm, replacement = use_ages[i], x = .) %>%
        rlang::parse_expr()


      proto$evict_fun[[it]]   <-  rlang::enquo(temp)

      names(proto$params)[it] <- gsub(nm,
                                      use_ages[i],
                                      names(proto$params)[it])


      if(it == length(use_ages)) {
        it <- 1
      } else {
        it <- it + 1
      }

    } # end ages names loop

  return(proto)
}


