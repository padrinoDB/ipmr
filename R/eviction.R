#' @title Eviction correction
#' @rdname eviction
#'
#' @description Various helpers to correct for unintentional eviction (Williams
#' et al. 2012). \code{rescale_kernel} is an alias for
#' \code{truncated_distributions}.
#'
#' @param fun The cumulative density function to use. For example, could be
#' \code{"norm"} to correct a Gaussian density function, or \code{"lnorm"} to
#' correct a log-normal density function.
#' @param param The parameter/vital rate being modified
#' @param L Optionally, the name of the lower bound of the domain. Otherwise,
#' the function will try to find the \code{param} inside of the \code{proto_ipm}
#' object under construction and use that to infer the domain of the function.
#' @param U Optionally, the name of the upper bound of the domain. Otherwise,
#' the function will try to find the \code{param} inside of the \code{proto_ipm}
#' object under construction and use that to infer the domain of the function.
#' @param ... Only used for internal modification - do not use!
#'
#' @return A matrix of the same dimension as the input.
#'
#' @references
#' Williams JL, Miller TEX & Ellner SP, (2012). Avoiding unintentional eviction
#' from integral projection models.Ecology 93(9): 2008-2014.
#'
#' @importFrom rlang caller_env
#' @export

truncated_distributions <- function(fun,
                                    param,
                                    L = NA,
                                    U = NA,
                                    ...) {

  proto <- ..1

  if(is.na(L) || is.na(U)) {

    LU <- .get_bounds_from_proto(param, proto)
    L <- LU[1]
    U <- LU[2]

  }

  proto <- .sub_new_param_call(fun, param, L, U, proto)

  return(proto)
}

# Internal helpers for eviction correction functions

#' @noRd
#' @importFrom rlang call_args

.sub_new_param_call <- function(fun, param, L, U, proto) {

  fun <- paste('p', fun, sep = "")
  param_form   <- .get_param_form(param, proto)
  fixed_params <- rlang::call_args(rlang::parse_expr(param_form))[-1] %>%
    unlist() %>%
    as.character() %>%
    paste(collapse = ', ')

  denom_1 <- paste(fun, '(', U, ', ', fixed_params, ')', sep = "")
  denom_2 <- paste(fun, '(', L, ', ', fixed_params, ')', sep = "")

  final_form <- paste(param_form,
                      ' / ',
                      '(',
                      denom_1,
                      ' - ',
                      denom_2,
                      ')',
                      sep = "")

  out <- .insert_final_form(param, final_form, proto)

  return(out)

}

.insert_final_form <- function(param, final_form, proto) {

  ind <- which(names(proto$params[[1]]$vr_text) == param)
  proto$params[[1]]$vr_text[ind] <- final_form

  return(proto)
}

#' @noRd

.get_bounds_from_proto <- function(param, proto) {

  # browser()
  # Get parameter function form
  param_form <- .get_param_form(param, proto)

  # Next, infer the state variable and construct the L, U

  svs <- unique(unlist(proto$state_var))

  sv_ind <- vapply(svs,
                   FUN = function(x) grepl(x, param_form),
                   logical(1))

  sv <- svs[sv_ind]

  out <- paste(c("L", "U"), sv, '2', sep = "_")

  return(out)

}

#' @noRd

.get_param_form <- function(param, proto) {

  all_params <- .flatten_to_depth(proto$params, 1)

  ind <- which(names(all_params) %in% param)

  if(length(ind) == 1) {

    param_form <- all_params[[ind]]

  } else{

    param_form <- all_params[ind]

  }
  return(param_form)
}

#' @noRd

.correct_eviction <- function(proto) {

  # If the quosure contains a call of ~list(...), that means it contains multiple
  # eviction functions. we need to loop over those individually. If it's a single
  # one, then just carry on as planned

  # First, get the actual function name in the outermost layer
  ev_call_nm <- rlang::call_name(
    rlang::quo_squash(
      unlist(proto$evict_fun)[[1]]
    )
  )

  # If it's a list of calls, then evaluate them one at a time.b
  if(ev_call_nm == 'list') {

    expr_list <- rlang::call_args(
      rlang::quo_squash(unlist(proto$evict_fun)[[1]])
    )

    temp <- lapply(expr_list, function(x){
      rlang::as_quosure(x,
                        env = rlang::empty_env())
    })

    for(i in seq_along(temp)) {

      ev_fun <- temp[[i]]

      if(grepl('truncated_distributions', rlang::quo_text(ev_fun))) {

        text <- gsub(')$', ', proto = proto[i, ])', rlang::quo_text(ev_fun))
        rep_expr <- rlang::parse_expr(text)
        ev_fun <- rlang::quo_set_expr(ev_fun, rep_expr)

      }

      evict_fun <- rlang::quo_set_env(ev_fun, rlang::caller_env())

      out <- rlang::eval_tidy(evict_fun)

    }

  # The single call gets evaluated in the standard (or rather, haphazard) way
  # I came up with when writing the simple_* methods.
  # This whole thing is RIPE for refactoring.

  } else if(ev_call_nm == 'truncated_distributions' |
            ev_call_nm == 'rescale_kernel') {

    evict_fun <- unlist(proto$evict_fun)[[1]]

    text <- gsub(')$', ', proto = proto[i, ])', rlang::quo_text(evict_fun))
    rep_expr <- rlang::parse_expr(text)
    evict_fun <- rlang::quo_set_expr(evict_fun, rep_expr)



    evict_fun <- rlang::quo_set_env(evict_fun, rlang::caller_env())

    out <- rlang::eval_tidy(evict_fun)

  }

  return(out)

}

.quo_has_list <- function(quo_)  {
  grepl('list\\(', quo_)
}
