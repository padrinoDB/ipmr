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

    # Wrapped with suppressWarnings() because the only one that ever comes out
    # is:
    # Warning message:
    # In max(unlist(lapply(l, .depth, start = start + 1))) :
    #   no non-missing arguments to max; returning -Inf
    # This is innocuous and doesn't seem to affect performance in any way, only
    # disconcerts those unaccustomed to its usage.

    suppressWarnings(
      max(
        unlist(
          lapply(l, .depth, start = start + 1)
        )
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
    attr(obj, "na_ok")        <- TRUE

  } else {

    attr(obj, "flat_protect") <- FALSE

    if(is.na(obj)) {
      attr(obj, "na_ok") <- FALSE
    } else {
      attr(obj, "na_ok") <- TRUE
    }

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
# Checks if model has already been iterated, saving us from re-iterating a model
# when we can just use the pop_state slot

.already_iterated <- function(ipm) {

  return(
    attr(ipm, "iterated")
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

  vr_nms <- names(vital_rate_exprs(proto_ipm))
  vrs    <- vital_rate_exprs(proto_ipm)
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

    sub_vrs  <- c(sub_vrs, .find_base_vrs(sub_vrs, vrs, do_nms))


    text <- .update_sum(text, sub_vrs, proto_ipm)

  }

  return(text)

}

#' @noRd

.txt_has_domain <- function(text, domains) {

  temp <- vapply(domains,
                 function(x, text) any(grepl(x, text)),
                 logical(1L),
                 text = text)

  return(any(temp))
}

#' @noRd

.find_base_vrs <- function(targets, all_vrs, domains) {

  all_nms  <- names(all_vrs)
  all_vrs  <- lapply(all_vrs, function(x) {
    if(is.character(x)) {
      return(x)
    } else {
      return(rlang::expr_text(x))
    }
  })

  all_args <- unlist(lapply(targets, .args_from_txt))

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

  if(length(out) == 0) return("")

  return(out)

}

#' @noRd

.expr_type <- function(x) {
  if (rlang::is_syntactic_literal(x)) {
    "constant"
  } else if (rlang::is_symbol(x)) {
    "symbol"
  } else if (rlang::is_call(x)) {
    "call"
  } else if (rlang::is_pairlist(x)) {
    "pairlist"
  } else {
    rlang::type_of(x)
  }
}

#' @noRd

.switch_expr <- function(x, ...) {
  switch(.expr_type(x),
         ...,
         stop("Don't know how to handle type ", typeof(x), call. = FALSE)
  )
}

#' @noRd

.calls_from_txt <- function(txt) {

  if(is.character(txt)){
    expr <- rlang::parse_expr(txt)
  } else {
    expr <- txt
  }

  .switch_expr(expr,
               constant = character(),
               symbol = character(),
               call = {

                 nm <- as.character(expr[[1]])
                 children <- lapply(as.list(expr[-1]), .all_calls) %>%
                   unlist()
                 c(nm, children)

               })
}

#' @noRd

.all_calls <- function(txt) {

  unique(.calls_from_txt(txt))

}

#' @noRd
#' @importFrom rlang child_env as_environment

.update_sum_txt <- function(txt, vr_exprs, nms, calls, proto_ipm) {

  sym_list <- as.list(nms) %>%
    setNames(nms) %>%
    .[!duplicated(names(.))]

  sym_env  <- rlang::as_environment(sym_list, parent = f_env)
  sum_env  <- rlang::child_env(.parent   = sym_env,
                               sum       = .sum2,
                               proto_ipm = proto_ipm,
                               vr_exprs  = vr_exprs)

  expr <- rlang::parse_expr(txt)
  expr <- rlang::quo(!!expr)

  eval_env <- rlang::child_env(sum_env)
  expr <- rlang::quo_set_env(expr, eval_env)

  out <- rlang::eval_tidy(expr)

  return(out)

}

#' @noRd
#' @importFrom rlang env_parents env_names

# This is the function that is evaluated instead of sum() to produce
# the correct result. It is bound as a function
.sum2 <- function(expr) {

  # Track down sum_env object and make it available here. We need this because
  # it contains the proto_ipm and vr_exprs objects that we were unable to pass
  # to sum(). Luckily, it will always be the parent of the caller environment,
  # which we can track down.

  ev_env  <- rlang::env_get(env = caller_env(),
                            nm = rlang::env_names(rlang::caller_env())[2])

  p_list  <- rlang::env_parents(ev_env)
  sum_env <- p_list[[2]]

  vr_exprs  <- rlang::env_get(env = sum_env,
                              nm = "vr_exprs",
                              inherit = FALSE,
                              default = stop("couldn't find 'vr_exprs'"))
  proto_ipm <- rlang::env_get(env = sum_env,
                              nm = "proto_ipm",
                              inherit = FALSE,
                              default = stop("couldn't find 'proto_ipm'"))

  # Now, we need to extract domain names and vr_expr names so that we can
  # search for the correct domain. The goal here is to find out whether or not
  # any of the terms in sum(...) are a function of z_1/z_2, so that we know
  # to divide them by n_z_1/ n_z_2.
  # NB: If we find no domains, then we don't actually do anything except re-wrap
  # the expression in sum(...) and return that instead.

  do_nms    <- names(domains(proto_ipm))

  vr_nms    <- names(vr_exprs)

  dos <- vapply(vr_exprs,
                function(x, do_nms) .vr_domain(x, do_nms),
                character(1L),
                do_nms = do_nms) %>%
    unique() %>%
    Filter(f = function(y) y != "", x = .)

  if(length(dos) > 0) {

    for(i in seq_along(dos)) {

      to_add <- paste(vr_nms[i], " / n_", dos[i], sep = "")

      text <- gsub(vr_nms[i], to_add, expr)
    }

  }

  expr <- paste("sum(", text, ")")

  return(expr)

}

#' @noRd

.update_sum <- function(txt, all_vrs,  proto_ipm) {

  all_nms   <- .args_from_txt(txt)
  all_calls <- .calls_from_txt(txt)
  all_calls <- all_calls[all_calls != "sum"]

  txt <- .update_sum_txt(txt,
                         all_vrs,
                         all_nms,
                         all_calls,
                         proto_ipm)

  return(txt)

}

#' @noRd
#' @importFrom rlang new_function caller_env exprs empty_env

.unary_op <- function(left, right) {
  rlang::new_function(
    rlang::exprs(e1 = ),
    rlang::expr(
      paste0(!!left, e1, !!right)
    ),
    rlang::caller_env()
  )
}

#' @noRd

.binary_op <- function(sep) {
  rlang::new_function(
    rlang::exprs(e1 = , e2 = ),
    rlang::expr(
      paste0(e1, !!sep, e2)
    ),
    rlang::caller_env()
  )
}

f_env <- rlang::child_env(
  .parent = rlang::empty_env(),
  `+` = .binary_op(" + "),
  `-` = .binary_op(" - "),
  `*` = .binary_op(" * "),
  `/` = .binary_op(" / "),
  `^` = .binary_op("^"),

  paste = paste,
  `[` = .unary_op("[", "]"),
  `(` = .unary_op("(", ")"),
  as.vector = .unary_op("as.vector(", ")"),

  # Other math functions
  sqrt = .unary_op("sqrt(", ")"),
  sin  = .unary_op("sin(", ")"),
  cos  = .unary_op("cos(", ")"),
  tan  = .unary_op("tan(", ")"),
  log  = .unary_op("log(", ")"),
  abs  = .unary_op("abs(", ")")
)

