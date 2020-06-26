#' @title Initialize an IPM
#'
#' @description This is always the first step in constructing an IPM with \code{ipmr}.
#' All you need for this is to know what type of IPM you want to construct - the rest comes later
#' with \code{define_kernel}, \code{define_k}, \code{make_ipm}, and associated helper functions.
#'
#' @param model_class The type of IPM. See \code{Details} for more information on classes
#' @param has_age A logical indicating whether the model has age structure. Default
#' is \code{FALSE}
#'
#' @return An object with classes \code{"proto_ipm"} and \code{model_class}. If
#' \code{has_age = TRUE}, then an \code{"age_x_size"} class is also added.
#'
#' @details Combinations of \code{simple} or \code{general}, \code{dd} or \code{di},
#' and \code{det} or \code{stoch} are separated by an underscore to tell \code{ipmr}
#' the type of model you will be constructing (e.g. \code{"simple_di_det"} for a
#' deterministic IPM with a single continuous state variable and no discrete stages
#' or density dependence).
#'
#' Within \code{stoch} model types, there are two additional options to append to
#' the class: \code{"kern"} or \code{"param"}. These distinguish between models that
#' use kernel resampling vs those that use parameter resampling (\emph{sensu} Metcalf et al.
#' 2015). Below are quick definitions, more detailed explanations can be found
#' in the \code{vignettes("ipmr-introduction", package = 'ipmr')}.
#'
#' \itemize{
#'
#'   \item{One of}
#'   \itemize{
#'     \item{\code{simple}}{: an IPM with a single continuous state variable that does not include
#'     any discrete stages. Simple IPMs can still be stochastic and/or density dependent.}
#'     \item{\code{general}}{: an IPM with more than one continuous state variable
#'     and/or a model that includes discrete stages.}
#'  }
#'  \item{and one of}
#'  \itemize{
#'     \item{\code{dd}}{: used to denote a density dependent IPM.}
#'     \item{\code{di}}{: used to denote a density independent IPM.}
#'  }
#'  \item{and one of}
#'  \itemize{
#'     \item{\code{det}}{: used to denote a deterministic IPM. Use this when you are
#'     building a single iteration kernel.}
#'     \item{\code{stoch}}{: used to denote a stochastic IPM. Stochasticity can be implemented
#'     in two ways in \code{ipmr}: \code{"kern"} resampling, and \code{"param"} resampling.}
#'  }
#'  \item{If using \code{det}, this should be omitted. If using \code{stoch}, then one of the following: }
#'  \itemize{
#'     \item{\code{kern}}{: used to denote an IPM that uses kernel resampling. Briefly,
#'     these models build all of the iteration kernels ahead of time and then choose one
#'     at random or in a user-specified order as they move from iteration to iteration. The
#'     user-specified population vector is multiplied by the chosen kernel and the result
#'     is multiplied by the next kernel for the desired number of iterations.}
#'     \item{\code{param}}{: used to denote parameter resampling. This samples distributions
#'     for each parameter based on user-specified functions supplied to \code{define_env_state()}.
#'     This will be a bit slower than \code{"kern"} resampling because iteration
#'     matrices need to be reconstructed from new parameters at every time step.}
#'}
#'}
#'
#' @references Metcalf et al. (2015). Statistical modelling of annual variation for inference on stochastic
#' population dynamics using Integral Projection Models. Methods in Ecology and Evolution, 6: 1007-1017
#'
#' @export

init_ipm <- function(model_class, has_age = FALSE) {

  out <- data.frame(
    id               = character(0L),
    kernel_id        = character(0L),
    domain           = character(0L),
    state_var        = character(0L),
    int_rule         = character(0L),
    evict            = logical(0L),
    evict_type       = character(0L),
    pop_state        = character(0L),
    env_state        = character(0L),
    hier_effs        = logical(0L),
    levels_hier_effs = character(0L),
    has_age          = logical(0L),
    levels_ages      = character(0L),
    params           = character(0L)
  )

  if(has_age) model_class <- c(model_class, "age_x_size")

  class(out) <- c(model_class, 'proto_ipm', class(out))

  return(out)
}

