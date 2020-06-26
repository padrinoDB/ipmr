# Age x size internal functions
# These work by evaluating the _minus_ and _plus_ keywords in each expression's
# name/content

# These should only appear in .make_age_k_row(...). "others" object is treated
# just like any other hierarchical model, for now. This will require modification
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

.make_age_k_row <- function(k_row) {

  forms <- k_row$params[[1]]$formula

  ages  <- k_row$levels_ages %>%
    .flatten_to_depth(1L)

  calls <- vapply(forms,
                  function(x) {

                    if(is.character(x)) x <- rlang::parse_expr(x)
                    rlang::call_name(x)

                  }, character(1L))

  # In the case where we have an n_0 = f_1 %*% n_1 + f_2 %*% n_2 etc,
  # We need to treat this differently from other hier_effs. Normally, they
  # get one call per level of the hier_eff. in this case, we need one call
  # expanded to include all levels of the hier_eff, and collapsed by the "fun"
  # passed to all_ages. On the other hand, we also need to generate one expression
  # for each n_x_t_1 = p_x_minus_1 %*% n_x_minus_1_t

  if(any(calls == "all_ages")) {

    n_0_ind  <- which(calls == "all_ages")
    n_0_call <- forms[n_0_ind]

    n_0_form <- lapply(
      n_0_call,
      function(x, ages) {
        if(is.character(x)) x <- rlang::parse_expr(x)
        temp     <- rlang::call_standardise(x)
        new_call <- rlang::call_modify(temp,
                                       ages = ages)

        out <- rlang::eval_tidy(new_call)
        return(out)

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
# Basically .split_hier_effs, but uses the age column instead of hier_effs column

.split_sub_kern_ages <- function(proto_ipm) {

  kerns <- which(proto_ipm$has_age)

  # Create a place to hold the output - either kernels with no hierarchical effects
  # or a new proto
  if(length(kerns) != dim(proto_ipm)[1]) {

    out <- proto_ipm[-kerns, ]

  } else {

    out <- init_ipm(class(proto_ipm)[1])

  }

  hier_rows <- proto_ipm[kerns, ]

  for(i in seq_len(dim(hier_rows)[1])) {

    levs   <- hier_rows$levels_ages[[i]]

    # Combine the max age w/ the others if it is supplied. This only
    # applies to iteration kernel expressions (I think).

    levs   <- list(age = unlist(levs, use.names = FALSE))

    levels <- lapply(levs, eval) %>%
      expand.grid(stringsAsFactors = FALSE)

    temp   <- .expand_hier_effs(hier_rows[i, ], levels)

    out    <- rbind(out, temp)

  }

  return(out)


}
