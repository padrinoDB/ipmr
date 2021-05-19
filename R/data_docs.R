#' @rdname sim_di_det
#'
#' @title Simple deterministic IPM example
#'
#' @format A simple deterministic IPM with the following slots:
#' \describe{
#'
#'   \item{sub_kernels}{The computed sub-kernels, named \code{P} and \code{F}.}
#'   \item{env_list}{Empty}
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
#'   \item{sub_kernels}{The computed sub-kernels for the model, named \code{P},
#'   \code{go_discrete}, \code{stay_discrete}, and \code{leave_discrete}.}
#'   \item{env_list}{Empty}
#'   \item{env_seq}{Contains \code{NA}. Not particularly useful for deterministic IPMs,
#'   but critical for reproducing stochastic ones.}
#'   \item{pop_state}{A list of length 2, with names \code{n_b} and
#'   \code{n_ht}.}
#'   \item{proto_ipm}{The \code{proto_ipm} used to implement the model.}
#'
#' }
'gen_di_det_ex'

#' @rdname raw_data_ex
#'
#' @title Raw demographic data to construct an example IPM
#'
#' @format 288 observations of 10 variables
#' \describe{
#'
#'   \item{id}{Individual identification number}
#'   \item{size}{Surface area in square meters of each individual at time \emph{t}.}
#'   \item{flower_n}{If the plant is reproductive, the number of flowers it made.}
#'   \item{log_size}{Log transformed \code{size}.}
#'   \item{repro}{Either 0 or 1 to indicate whether the plant is reproductive.}
#'   \item{size_next}{Surface area in square meters of each individual at time \emph{t + 1}.}
#'   \item{flower_n_next}{If the plant is reproductive at \emph{t + 1}, the number of
#'   flowers it made.}
#'   \item{survival}{Either 0 or 1 to indicate whether a plant at \emph{t} survives to \emph{t + 1}.}
#'   \item{log_size_next}{Log transformed \code{size_next}.}
#'   \item{repro_next}{Either 0 or 1 to indicate whether a plant is reproductive at \emph{t + 1}.}
#'}
'iceplant_ex'

#' @rdname proto_ex
#' @title A \code{proto_ipm} for a monocarpic perennial
#'
#' @format A \code{proto_ipm} for a simple IPM of \emph{Oenothera glazioviana}.
#' The parameters are from Ellner, Childs, & Rees (2016), Chapter 2, and the data
#' are from Kachi & Hirose (1985). Parameter values can be accessed with
#' \code{parameters(monocarp_proto)}, vital rate expressions can be accessed with
#' \code{vital_rate_exprs(monocarp_proto)}, etc.
#'
#'
#' @references
#' Kachi, H., & Hirose, T. (1985). Population dynamics of _Oenothera glazioviana_
#' in a sand-dune system with special reference to the adaptive significance of
#' size-dependent reproduction. Journal of Ecology 73: 887-901.
#' https://doi.org/10.2307/2260155
#'
#' Ellner, S.P., Childs, D.Z., Rees, M. (2016) Data-driven modelling of
#' structured populations: a practical guide to the integral projection model.
#' Basel, Switzerland: Springer International Publishing AG
"monocarp_proto"
