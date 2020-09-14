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


