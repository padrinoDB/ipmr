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

.sens_lam_kern.general_di_det_ipm <- function(ipm,
                                              n_iterations,
                                              mega_mat,
                                              mega_vec,
                                              ...) {


  r_ev <- right_ev(ipm, n_iterations, keep_mega = FALSE)
  l_ev <- left_ev(ipm,
                  !! mega_mat,
                  !! mega_vec,
                  n_iterations,
                  keep_mega = FALSE)

  # Don't think this works for multiple continuous state variables unless
  # z is rectangular (e.g. won't work for size_type_1 and size_type_2 sort
  # of domains, but will for size x quality).

  d_zs        <- .get_d_z_vars(ipm)
  big_dz      <- Reduce("*", d_zs, init = 1)

  # .sens_kern_lam_general matches d_zs to their respective portions
  # of r/l_ev, and then generates a mega_vec out of each and computes
  # the sensitivity kernel using the standard formula

  out         <- .sens_kern_lam_general(r_ev, l_ev, big_dz, mega_vec)

  return(out)

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

  sens_kern <- .sens_lam_kern(ipm,
                              n_iterations,
                              mega_mat,
                              mega_vec)

  all_lams  <- lambda(ipm, type_lambda = 'all')
  use_lam   <- all_lams[length(all_lams)]

  # Won't work for Ks that comprise multiple state variables that aren't crossed
  # with each other. For example, size X quality or size X age is fine, but
  # something like K_size_type_1 and K_size_type_2 will fail I think. An older form
  # of this segregated K/h by matching the K to the appropriate h, but that produce
  # inaccurate results for models with discrete states. See
  # experimental_funs/possibly_retired_funs.txt for that implementation.

  d_zs      <- .get_d_z_vars(ipm)

  d_zs      <- Reduce("*", d_zs, init = 1)

  k_fun         <- .make_mega_mat(mega_mat, ipm$sub_kernels) / d_zs

  out           <- sens_kern * k_fun / use_lam

  return(out)
}

# Helpers us ed by both sensitivity and elasticity -----------

#' @noRd

.sens_kern_lam_general <- function(r_ev, l_ev, d_zs, mega_vec) {

  names(r_ev) <- gsub('_w', "", names(r_ev))
  names(l_ev) <- gsub('_v', "", names(l_ev))

  r_ev_vec   <- .make_mega_vec(mega_vec, r_ev)
  l_ev_vec   <- .make_mega_vec(mega_vec, l_ev)

  # Already multiplied by d_zs

  vw        <- sum(r_ev_vec * l_ev_vec * d_zs)
  pert_kern <- outer(l_ev_vec, r_ev_vec)

  out <- pert_kern / vw

  return(out)
}

# I think almost everything below is getting retired, but need to reinspect
# this implementation if/when I do get it working to see who's currently
# surplus to requirements. It could be that we *do* need it for IPMs
# with multiple, uncrossed continuous states...

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


  # The unlist %>% as.list pattern seems dumb, but it basically generates
  # a named list of depth 1 that preserves names of each entry regardless of
  # depth. Basically, just want to gsub() out the periods for underscores

  out <- lapply(d_z_in_kerns,
                function(x, data_list) .eval_list(x, data_list),
                data_list = as.list(d_z_vec)) %>%
    unlist() %>%
    as.list()

  names(out) <- gsub('\\.', '_', names(out))

  return(out)
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

.eval_list <- function(x, data_list) {

  if(rlang::is_callable(x) ||
     rlang::is_bare_numeric(x)) {
    out <- rlang::eval_tidy(x, data = data_list)
    return(out)
  } else if(is.character(x)) {
    stop('all entries in form_expr_list should be symbols or numbers by now')
  } else {
    lapply(x,
           function(x, data_list) .eval_list(x, data_list),
           data_list = data_list)
  }

}

#' @noRd
#' @importFrom rlang is_callable is_bare_numeric

.match_dz_to_kern <- function(form_expr, nm) {

  if(rlang::is_bare_numeric(form_expr)) {
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

  # Gnarly pipe. first extract names of domains used in the model. Some of these
  # will be start/end_not_applicable because the kernels they're associated with
  # are DC, CD, or DD (and have no continuous state). We just want to find unique
  # names, but those will have _1/_2 appended, so find those variables only.
  # then filter out the unneeded ones, split open the continuous states to just
  # get the variable they describe, and then unique them to make sure we don't
  # duplicate anything.

  rm_disc_state_nms <- function(x) {
    grepl('_1|_2', x)
  }

  d_z_vars   <- lapply(ipm$proto_ipm$domain,
                       names) %>%
    unlist() %>%
    unique() %>%
    Filter(f = rm_disc_state_nms, x = .) %>%
    vapply(FUN = function(x) strsplit(x, '_')[[1]][1],
           FUN.VALUE = character(1L)) %>%
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

  # just want to work with kernels we need. split out the others so they can
  # be appended before exit

  is_not_1 <- function(x) x != 1
  use_dzs <- Filter(is_not_1, d_zs)

  nm_in_dz <- names(sub_kernels) %in% names(use_dzs)

  use_sub_kerns <- sub_kernels[nm_in_dz]
  append_to_out <- sub_kernels[!nm_in_dz]

  use_sub_kerns <- use_sub_kerns[names(use_dzs)]

  # Perform the conversion and restore the names so .make_mega_mat can find the
  # right objects.

  kern_funs <- lapply(seq_along(use_sub_kerns),
                      function(i, sub_kern_list, dz_list) {
                        out <- sub_kern_list[[i]] / dz_list[[i]]
                        return(out)
                      },
                      sub_kern_list = use_sub_kerns,
                      dz_list       = use_dzs)

  names(kern_funs) <- names(use_sub_kerns)

  out <- purrr::splice(kern_funs, append_to_out)

  return(out)
}







