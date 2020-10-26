# Reduces the depth/nestedness of a list to "depth".

#' @noRd

.flatten_to_depth <- function(to_flatten, depth) {

  protected <- lapply(to_flatten,
                      function(x) isTRUE(attr(x, "flat_protect"))) %>%
    unlist()

  if(any(protected)) {

    keep <- to_flatten[protected]

    to_flatten <- to_flatten[! protected]

  } else {

    # Make sure "keep" exists, but append nothing

    keep <- NULL

  }

  if(rlang::is_empty(to_flatten) |
     ! rlang::is_list(to_flatten)) {

    return(c(to_flatten, keep))

  }

  if(.depth(to_flatten) == depth) {

    return(c(to_flatten, keep))

  } else {

    to_flatten <- purrr::flatten(to_flatten)

    to_flatten <- c(to_flatten, keep)

    .flatten_to_depth(to_flatten, depth)

  }
}

# Finds depth of current list
#' @noRd

.depth <- function(l, start = 0) {

  if(!is.list(l) |
     isTRUE(attr(l, "flat_protect"))) {

    return(start)

  } else {
    max(
      unlist(
        lapply(l, .depth, start = start + 1)
      )
    )
  }
}


#' @noRd

# Protects model objects from getting flattened. This is useful for user-specified
# models where they use a `predict()` method and pass the model object into
# the data list slot.

.protect_model <- function(obj) {

  # If the user has already defined this via use_vr_model, then
  # we don't want to do anything to it

  if(!is.null(attr(obj, "flat_protect"))) {
    return(obj)
  }

  supported_models <- .supported_models()

  if(inherits(obj, supported_models)) {

    attr(obj, "flat_protect") <- TRUE

  } else {

    attr(obj, "flat_protect") <- FALSE

  }

  return(obj)

}

#' @noRd
# vector of model classes that ipmr will allow to be used in predict() expressions.
# I actually don't think some of these even have predict methods, but just being
# safe.

.supported_models <- function() {

  c("lm",
    "glm",
    "gam",
    "gls",
    "lme",
    "betareg",
    "biglm",
    "glmnet",
    "gamlss",
    "nls",
    "MCMCglmm",
    "rjags",
    "merMod",
    "brmsfit",
    "stanfit")

}

#' @noRd

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


#' @noRd

.stoch_progress_message <- function(report_progress, iterations, iteration) {

  if(report_progress) {

    perc_prog <- iteration / iterations * 100

    if(perc_prog %% 10L == 0){

      message("\nModel progress: ", iteration, "/", iterations)
    }

  }

  invisible(TRUE)

}


#' @noRd
# Function to recursively extract argument names from text representations
# of an expression. This returns a character vector of all arguments a call,
# and does not track what the function call actually is.

.args_from_txt <- function(txt) {

  expr <- rlang::parse_expr(txt)

  if(rlang::is_atomic(expr) || rlang::is_symbol(expr)) {

    return(unlist(txt))

  } else {

    temp <- rlang::call_args(expr)
    tst  <- vapply(temp,
                   function(x) rlang::is_atomic(x) || rlang::is_symbol(x),
                   logical(1L))

    if(all(tst)) {

      out <- vapply(temp,
                    function(x) rlang::expr_text(x),
                    character(1L),
                    USE.NAMES = FALSE)

      return(unlist(out))

    } else {

      temp <- lapply(temp, rlang::expr_text)
      lapply(temp, .args_from_txt) %>%
        unlist()
    }
  }
}


# Utilities for sum() in vr_exprs

#' @noRd

.check_sum_calls <- function(text) {

  grepl("sum\\(", text)

}

# Adds a call to divide by n_meshpoints if needed. nearly every vital rate expression
# is evaluated on the domain [z', z], rather than just, say, z. This means there
# are duplicated values for univariate functions, which for most cases doesn't
# matter and helps us vectorize all calculations. However, if the user calls
# sum(n_z_t * b_z), this will result in either length mismatches, or worse, a silent
# failure whereby the sum of b_z, a function clearly intended to be univariate,
# is internally considered bivariate.
#' @noRd

.prep_sum_calls <- function(text, proto_ipm, main_env) {

  # Get vital rate names and domain names.

  vr_nms <- names(vital_rates(proto_ipm))
  vrs    <- vital_rates(proto_ipm)
  do_nms <- names(domains(proto_ipm))

  # Now get argument names from target expression. We'll check to see if
  # any of those match names of other vital rate expressions in this kernel.

  all_args <- .args_from_txt(text)

  vr_tst <- vapply(vr_nms,
                   function(x, args) x %in% args,
                   logical(1L),
                   args = all_args)

  if(any(vr_tst)) {

    sub_vrs    <- vrs[vr_tst] %>%
      lapply(rlang::expr_text)

    sub_vrs    <- .find_base_vrs(sub_vrs, all_vrs = vrs, do_nms)

    vr_nms <- vr_nms[vr_tst]

    dos    <- vapply(sub_vrs, function(x) .vr_domain(x, do_nms), character(1L)) %>%
      unique()

    for(i in seq_along(dos)) {

      to_add <- paste(vr_nms[i], " / n_", dos[i], sep = "")

      text <- gsub(vr_nms[i], to_add, text)

    }
  }

  return(text)
}

.txt_has_domain <- function(text, domains) {

  temp <- vapply(domains,
                 function(x, text) any(grepl(x, text)),
                 logical(1L),
                 text = text)

  return(any(temp))
}

.find_base_vrs <- function(targets, all_vrs, domains) {

  force(domains)

  all_nms  <- names(all_vrs)
  all_vrs  <- lapply(all_vrs, function(x) {
    if(is.character(x)) {
      return(x)
    } else {
      return(rlang::expr_text(x))
    }
  })

  all_args <- .args_from_txt(unlist(targets))

  if(any(.txt_has_domain(all_args, domains))) {

    return(targets)

  }

  if(any(all_nms %in% all_args)) {

    temp_vrs <- all_vrs[which(all_nms %in% all_args)]

    lapply(temp_vrs,
           function(x, all_vrs, domains) .find_base_vrs(x,
                                                        all_vrs = all_vrs,
                                                        domains = domains),
           all_vrs = all_vrs,
           domains = domains)

  } else {

    return(targets)

  }

}


#' @noRd
# Infer the domain(s) of a vital rate expression

.vr_domain <- function(text, domains) {

  domains <- unique(domains)

  ind <- logical(length(domains))

  for(i in seq_along(domains)) {

    ind[i] <- grepl(domains[i], text)

  }

  out <- domains[ind]

  return(out)

}
