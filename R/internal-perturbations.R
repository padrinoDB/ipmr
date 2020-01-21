# Internal generics for sensitivity and elasticity

# Sensitivity internal generics --------------
#' @noRd

.sens_lam_kern <- function(ipm, n_iterations, ...) {
  UseMethod(".sens_lam_kern")
}

#' @noRd

.sens_lam_kern.simple_di_det_ipm <- function(ipm,
                                             n_iterations,
                                             ...) {

  r_ev <- right_ev(ipm, n_iterations) %>%
    unlist()
  l_ev <- left_ev(ipm, n_iterations) %>%
    unlist()

  d_zs <- .get_d_z_vars(ipm)

  vw        <- sum(r_ev * l_ev * d_zs)
  pert_kern <- outer(l_ev, r_ev)

  out <- pert_kern / vw

  return(out)

}

#' @noRd

.sens_lam_kern.general_di_det <- function(ipm,
                                          n_iterations,
                                          mega_mat,
                                          mega_vec,
                                          ...) {

  r_ev <- right_ev(ipm, n_iterations, keep_mega = FALSE)
  l_ev <- left_ev(ipm,
                  mega_mat,
                  mega_vec,
                  n_iterations,
                  keep_mega = FALSE)

  d_zs <- .match_d_zs

}



# Elasticty internal generics ------------------
#' @noRd

.elas_lam_kern <- function(ipm, n_iterations, ...) {
  UseMethod('.elas_lam_kern')
}

#' @noRd

.elas_lam_kern.simple_di_det_ipm <- function(ipm,
                                             n_iterations,
                                             ...) {

  sens_kern <- .sens_lam_kern(ipm, n_iterations)

  lam       <- lambda(ipm, comp_method = 'eigen')

  d_zs      <- .get_d_z_vars(ipm)

  # simple_di_det's only have 1 iterator!

  K         <- ipm$iterators[[1]]

  out       <- sens_kern * (K / d_zs) / lam

  return(out)
}

#' @noRd

.elas_lam_kern.general_di_det_ipm <- function(ipm,
                                              n_iterations,
                                              mega_mat,
                                              mega_vec) {

  sens_kern <- .sens_lam_kern(ipm, n_iterations, mega_mat, mega_vec)

  all_lams  <- lambda(ipm, type_lambda = 'all')
  use_lam   <- all_lams[length(all_lams)]

  # Won't work at the moment because we need to make sure each portion of the
  # kernel is divided by it's respective d_z. I think that'll have to come
  # using a combination the mega_mat formula and the kernel/K expressions in
  # the IPM's proto
  d_zs      <- .match_d_zs(ipm)

  # Above-mentioned step should go here, then make the mega_k using the altered
  # sub-kernels

  sub_kern_funs <- .iteration_mat_to_fun(ipm$sub_kernels,
                                         d_zs)

  use_k     <- .make_mega_mat(mega_mat, sub_kern_funs)



  out       <- sens_kern * (use_k / d_zs) / use_lam

  return(out)
}

# Helpers us ed by both sensitivity and elasticity -----------

#' @noRd

# returns a list where names correspond to kernels and entries correspond
# to the appropriate d_z. This is probably a bit shaky - it uses fuzzy matching
# of d_Sv in every kernel/sub-kernel expression. It's hard to imagine anyone would
# have a pattern matching it that *wasn't* referencing a d_z, but it *is* possible
# I guess. Will want to come up with a more robust implementation later.

.match_d_zs <- function(ipm) {

  # Extract the d_z values. This is a named vector. We'll search for their
  # names in the list that contains all of the model formulae in it.

  d_z_vec               <- .get_d_z_vars(ipm)
  d_z_nms               <- names(d_z_vec)

  form_expr_list        <- lapply(ipm$proto_ipm$params,
                                  function(x) .parse_recursive(x$formula))

  # This will return a list of the same structure as form_expr_list,
  # but it will have symbols representing the d_z_nms inserted. We can
  # then eval_tidy this list using the d_z_vec as a data mask to sub in the
  # the actual values that need to be divided/multiplied by.

  d_z_in_kerns          <- lapply(form_expr_list,
                                  function(x, d_z_nms) .find_d_zs(x, d_z_nms),
                                  d_z_nms = d_z_nms)


  # Continue here. Currently working out how to retain list structure/depth
  # to know how deeply, say, a K-expr is buried in with a d_z
  out <- rlang::eval_tidy()

}

.find_d_zs <- function(form_expr, d_z_nms) {

  out <- list()
  for(i in seq_along(d_z_nms)) {

    matches <- .match_dz_to_kern(form_expr,
                                 d_z_nms[i])

    out <- purrr::splice(out, matches)
  }
  return(out)
}

#' @noRd
#' @importFrom rlang is_callable

.match_dz_to_kern <- function(form_expr, nm) {

  if(is.integer(form_expr) ||
     is.numeric(form_expr)) {
    return(1L)
  }

  if(rlang::is_callable(form_expr)) {

    cl_args <- rlang::call_args(form_expr)

    if(any(as.character(cl_args) == nm)) {
      return(rlang::sym(nm))
    } else {
      return(1L)
    }

  } else if(is.list(form_expr)) {

    lapply(form_expr,
           function(x, nm) .match_dz_to_kern(x, nm),
           nm = nm)

  } else {

    stop("Sam fucked up, email him with a reproducible example.")
  }

}

#' @noRd

.get_d_z_vars <- function(ipm) {

  # Makes the master_env and retrieves the d_z variables from it

  if(all(is.na(ipm$env_list))){
    master_env <- .make_master_env(ipm$proto_ipm$domain,
                                   ipm$proto_ipm$usr_funs[[1]])
  } else {
    master_env <- ipm$env_list$master_env
  }

  d_z_vars   <- ipm$proto_ipm$state_var %>%
    unlist() %>%
    unique()

  d_z_var_nms  <- paste('d_', d_z_vars, sep = "")

  d_z_vars <- rlang::env_get_list(master_env, nms = d_z_var_nms)
  d_zs     <- unlist(d_z_vars)


  names(d_zs) <- d_z_var_nms

  return(d_zs)
}

#' @noRd

.parse_recursive <- function(x) {

  if(!is.list(x) && is.character(x)) {
    # Success case
    rlang::parse_expr(x)

  } else if(!is.character(unlist(x))) {
    # Error - Unsure how to proceed
    stop('Error in parsing kernel expressions')

  } else {
    # Recursive case - if x itself is a list, then further down the rabbit hole
    lapply(x, .parse_recursive)

  }
}

#' @noRd

.iteration_mat_to_fun <- function(sub_kernels,
                                  d_zs) {
  stop("not implemented")
}
