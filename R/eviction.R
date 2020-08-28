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
#' @param param The parameter/vital rate being modified. If this is a vector, the
#' distribution specified in \code{fun} will be recycled.
#' @param state The state variable used in the kernel that is being discretized.
#' @param ... Only used for internal modification - do not use!
#'
#' @return For \code{truncated_distributions} and \code{rescale_kernel}, a
#' modified function call with that re-scales the probability density
#' function based on the cumulative density function.
#'
#' For \code{discrete_extrema}, a numeric vector with modified entries based
#' on the discretization process.
#'
#' @note Neither of these functions are intended for use outside of
#' \code{define_kernel}, as they rely on internally generated variables to
#' work inside of \code{make_ipm}.
#'
#' @references
#' Williams JL, Miller TEX & Ellner SP, (2012). Avoiding unintentional eviction
#' from integral projection models.Ecology 93(9): 2008-2014.
#'
#' @importFrom rlang caller_env
#' @export

truncated_distributions <- function(fun,
                                    param,
                                    ...) {

  proto <- ..1

  for(i in seq_along(param)) {

    if(length(fun) != length(param)) {
      warning("length of 'fun' in 'truncated_distributions()' is not equal to ",
              "length of 'param'. Recycling 'fun'.",
              call. = FALSE)

      use_fun <- fun
    } else {
      use_fun <- fun[i]
    }

    LU    <- .get_bounds_from_proto(param[i], proto)
    L     <- LU[1]
    U     <- LU[2]

    proto <- .sub_new_param_call(use_fun, param[i], L, U, proto)
  }

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

  # Get state variable and the names corresponding to its bounds

  svs <- lapply(proto$domain, function(x) names(x)) %>%
    unlist() %>%
    unique()

  sv_ind <- vapply(svs,
                   function(x) grepl("_2", x),
                   logical(1L))

  sv     <- svs[sv_ind]

  out    <- paste(c("L", "U"), sv, sep = "_")

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
#' @importFrom rlang expr_text

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

 if(ev_call_nm == 'truncated_distributions' |
    ev_call_nm == 'rescale_kernel') {

    evict_fun <- unlist(proto$evict_fun)[[1]]

    # This doesn't get evaluated at any point - just fixes a note from R CMD check
    i <- NULL

    evict_fun <- rlang::call_modify(evict_fun, proto = rlang::expr(proto[i, ]))

    evict_fun <- rlang::quo_set_env(evict_fun, rlang::caller_env())

    out <- rlang::eval_tidy(evict_fun)

 } else if(ev_call_nm == "discrete_extrema") {

   evict_fun <- unlist(proto$evict_fun)[[1]]

   target <- rlang::call_args(evict_fun)[[1]][1]
   d_z    <- rlang::call_args(evict_fun)[[2]][1]

   d_z    <- paste("d_", d_z, sep = "")

   vr_nms <- names(proto$params[[1]]$vr_text)

   if(!target %in% vr_nms) {
     stop("Cannot find '", target, "' in vital rate expressions.")
   } else {
     ind <- which(vr_nms == target)
   }

   original_form <- proto$params[[1]]$vr_text[ind] %>%
     lapply(rlang::parse_expr)

   modified_form <- lapply(original_form,
                           function(x, d_z) {

                             temp <- rlang::call2(.fn = "discrete_extrema",
                                                  x, d_z)

                             return(rlang::expr_text(temp))

                           },
                           d_z = rlang::ensym(d_z))

   proto$params[[1]]$vr_text[ind] <- modified_form

   return(proto)

 } else {

   stop("This type of 'evict_fun' is not yet implemented", call. = FALSE)

 }

  return(out)

}

#' @rdname eviction
#' @param ncol,nrow The number of rows or column that the final form of the iteration
#' matrix should have. This is only necessary for rectangular kernels with
#' class \code{"CC"}. \code{make_ipm} works out the correct dimensions for
#' \code{"DC"} and \code{"CD"} kernels internally.
#' @param state The state variable used in the kernel that is being discretized.

#' @export

discrete_extrema <- function(fun, state, ncol = NULL, nrow = NULL) {

  if(is.null(ncol) && is.null(nrow)) {
    ncol <- nrow <- sqrt(length(fun))
  }

  temp <- matrix(fun, ncol = ncol, nrow = nrow, byrow = TRUE) * state

  top_seq <- seq(1, ncol / 2, by = 1)
  bot_seq <- seq(max(top_seq + 1), ncol, by = 1)


  for(i in top_seq) {

    temp[1 ,i] <- temp[1, i] + (1 - sum(temp[ , i]))

  }

  for(i in bot_seq) {

    temp[nrow ,i] <- temp[nrow, i] + (1 - sum(temp[ , i]))

  }

  # We need to transpose before stripping the dimensions, because otherwise
  # it will create a vector from columns first, whereas our domains generate
  # rows first.

  out <- t(temp)

  dim(out) <- NULL

  return(out)

}
