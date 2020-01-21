# Internal generics for sensitivity and elasticity

#' @noRd

.sens_lam_kern <- function(ipm, n_iterations, ...) {
  UseMethod(".sens_lam_kern")
}

.sens_lam_kern.simple_di_det_ipm <- function(ipm, n_iterations,) {


  r_ev <- right_ev(ipm, n_iterations)
  l_ev <- left_ev(ipm, mega_mat)


}
