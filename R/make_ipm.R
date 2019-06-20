#' @title Methods to implement an IPM
#'
#' @description The \code{make_ipm.*} methods convert a \code{proto_ipm} into a
#' set of discretized kernels and population vectors. Methods have different
#' requirements, so carefully read the parameter documentation and the
#' \code{vignette("implementation", package = "ipmr")}.
#'
#' @param proto_ipm The proto_ipm object you wish to implement. This should be the
#' output of \code{add_kernel}, \code{add_K}, or the \code{define_*} functions.
#' @param ... Other arguments passed to methods
#'
#' @return The \code{make_ipm.*det} methods will always return a list of length 4
#' containing the following components:
#'
#' \itemize{
#'   \item{\strong{iterators}}{: iteration kernel(s) (if specified by \code{add_K}),
#'                             otherwise contains \code{NA}.}
#'   \item{\strong{sub_kernels}}{: the sub_kernels specified in \code{add_kernel}.}
#'   \item{\strong{pop_state}}{: population vectors stored as an instance of the
#'                              \code{pop_state} class.}
#'   \item{\strong{proto_ipm}}{: the \code{proto_ipm} object used to implement
#'                              the model.}
#' }
#'  The \code{make_ipm.*stoch} methods will always return a list of length 5
#' containing the following components:
#'
#' \itemize{
#'   \item{\strong{iterators}}{: iteration kernel(s) (if specified by \code{add_K}),
#'                             otherwise contains \code{NA}.}
#'   \item{\strong{sub_kernels}}{: the sub_kernels specified in \code{add_kernel}.}
#'   \item{\strong{pop_state}}{: population vectors stored as an instance of the
#'                              \code{pop_state} class.}
#'   \item{\strong{env_state}}{: a matrix with dimension \code{n_iterations} X 1 of
#'                              kernel indices indicating the order
#'                              in which kernels are to be/were resampled OR
#'                              a matrix with as many columns as stochastic parameters
#'                              \code{n_iterations} rows.}
#'   \item{\strong{proto_ipm}}{: the \code{proto_ipm} object used to implement
#'                              the model.}
#' }
#'
#'
#'
#'
#'
#' @author Sam Levin
#'
#' @export

make_ipm <- function(proto_ipm, ...) {
  UseMethod('make_ipm')
}


#' @inheritParams make_ipm
#'
#' @export

make_ipm.simple_di_det <- function(proto_ipm, domain_list = NULL) {

  # Make sure we clean up the shimmed environment from the search path
  on.exit(detach("domain_env", character.only = TRUE))

  # Split out K from others so it isn't evaluated until we're ready. If it
  # isn't there, then proceed as usual

  K_row <- which(grepl("K", proto_ipm$kernel_id))

  if(length(k_row) > 0) {
    k_row <- proto_ipm[K_row, ]
    sub_kernels <- proto_ipm[-c(K_row), ]
  } else {
    sub_kernels <- proto_ipm
  }

  # If vital rates are fit with a hierarchical model of any kind,
  # then split those out into their respective years/plots/what-have-you
  # BE SURE TO WRITE VIGNETTE ON THIS SYNTAX

  if(.has_hier_effs(others) | .has_hier_effs(k_row)) {
    others <- .split_hier_effs(others)
    k_row <- .split_hier_effs(k_row)
  }

  # Initialize the domain_environment so these values can all be found at
  # evaluation time
  if(is.null(domain_list)){
    domain_env <- .generate_domain_env(others$domain)
  } else {
    domain_env <- .generat_domain_env(domain_list)
  }
  # Loop over the kernels for evaluation
  sub_kern_list <- list()

  for(i in seq_len(dim(others)[1])) {

    param_tree <- others$params[[i]]

    # kern_env inherits from domain_env so that those variables are
    # findable at evaluation time

    kern_env <- .generate_kernel_env(param_tree$params, domain_env)

    kern_quos <- .parse_vr_formulae(param_tree$vr_text,
                                    kern_env)
    kern_form <- .parse_kern_formula(param_tree$formula,
                                     kern_env)

    rlang::env_bind_lazy(kern_env,
                         !!! kern_quos,
                         .eval_env = kern_env)

    if(others$evict[i]) {
      rlang::env_bind_lazy(kern_env,
                           !!! evict_type)
    }


    sub_kern_list[[i]] <- rlang::env_get(kern_env, "formula")
    names(sub_kern_list)[i] <- others$kernel_id[i]
    class(sub_kern_list[[i]]) <- others$params[[i]]$family
  }

  if(length(K_row) > 0) {
    K <- make_k(k_row, proto_ipm, sub_kern_list)
  } else {
    K <- NA_character_
  }

  out <- .generate_ipm_output(K,
                              sub_kern_list,
                              pop_state,
                              proto_ipm)

  return(out)

}

make_ipm.simple_di_stoch_kern <- function(proto_ipm, ...) {

  # DEFINE ME
}

make_ipm.simple_di_stoch_param <- function(proto_ipm, ...) {

  # DEFINE ME
}

make_ipm.general_di_det <- function(proto_ipm, ...) {

  # DEFINE ME
}

make_ipm.general_di_stoch_kern <- function(proto_ipm, ...) {

  # DEFINE ME
}

make_ipm.general_di_stoch_param <- function(proto_ipm, ...) {

  # DEFINE ME
}

# Density dependent methods----------

make_ipm.simple_dd_det <- function(proto_ipm, ...) {

  # DEFINE ME
}

make_ipm.simple_dd_stoch_kern <- function(proto_ipm, ...) {

  # DEFINE ME
}

make_ipm.simple_dd_stoch_param <- function(proto_ipm, ...) {

  # DEFINE ME
}

make_ipm.general_dd_det <- function(proto_ipm, ...) {

  # DEFINE ME
}

make_ipm.general_dd_stoch_kern <- function(proto_ipm, ...) {

  # DEFINE ME
}

make_ipm.general_dd_stoch_param <- function(proto_ipm, ...) {

  # DEFINE ME
}
