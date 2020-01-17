#' @rdname sim_di_det
#'
#' @title Simple deterministic IPM example
#'
#' @format A simple deterministic IPM with the following slots:
#' \describe{
#'
#'   \item{iterators}{The computed iteration kernel, named \code{K}.}
#'   \item{sub_kernels}{The computed sub-kernels, named \code{P} and \code{F}.}
#'   \item{env_list}{Empty.}
#'   \item{env_seq}{Empty.}
#'   \item{pop_state}{Empty.}
#'   \item{proto_ipm}{The \code{proto_ipm} object used to implement the model.}
#'
#' }
'sim_di_det_ex'

#' @rdname gen_di_det
#'
#' @title A general deterministic IPM example
#'
#' @format A general deterministic IPM with the following slots:
#' \describe{
#'
#'   \item{iterators}{Empty.}
#'   \item{sub_kernels}{The computed sub-kernels for the model, named \code{P},
#'   \code{go_discrete}, \code{stay_discrete}, and \code{leave_discrete}.}
#'   \item{env_list}{Empty.}
#'   \item{env_seq}{A vector of 1s. Not particularly useful for deterministic IPMs,
#'   but critical for reproducing stochastic ones.}
#'   \item{pop_state}{A list of length 2, with names \code{pop_state_b} and
#'   \code{pop_state_ht}.}
#'   \item{proto_ipm}{The \code{proto_ipm} used to implement the model.}
#'
#' }
'gen_di_det_ex'
