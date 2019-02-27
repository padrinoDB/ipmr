#' @title Initialize an IPM
#'
#' @description This is always the first step in constructing an IPM with \code{ipmr}.
#' All you need for this is to know what type of IPM you want to construct - the rest comes later
#' with \code{add_kernel()}, \code{add_K}, \code{make_ipm}, and associated helper functions.
#'
#' @param class The type of IPM. See \code{Details} for more information on classes
#' @param ... Additional arguments passed to methods (???)
#'
#' @return An object of classes \code{proto_ipm} and \code{class}
#'
#' @details Combinations of \code{simple} or \code{general}, \code{dd} or \code{di},
#' and \code{det} or \code{stoch} are separated by an underscore to tell \code{ipmr}
#' the type of model you will be constructing (e.g. \code{"simple_di_det"} for a
#' deterministic IPM with a single continuous state variable and no discrete stages
#' or density dependence).
#'
#' Within \code{stoch} model types, there are two additional options to append to
#' the class: \code{"kern"} or code{"param"}. These distinguish between models that
#' use kernel resampling vs those that use parameter resampling (\emph{sensu} Metcalf et al.
#' 2015). Below are quick definitions, more detailed explanations can be found
#' in the \code{vignettes("classes", package = 'ipmr')}.
#'
#' \itemize{
#'   \item{\code{simple}}{: an IPM with a single continuous state variable that does not include
#'   any discrete stages. Simple IPMs can still be stochastic and/or density dependent.}
#'   \item{\code{general}}{: an IPM with more than one continuous state variable
#'   and/or a model that includes discrete stages.}
#'
#'   \item{\code{dd}}{: used to denote a density dependent IPM.}
#'   \item{\code{di}}{: used to denote a density independent IPM.}
#'
#'   \item{\code{det}}{: used to denote a deterministic IPM. Use this when you are
#'   building a single iteration kernel.}
#'   \item{\code{stoch}}{: used to denote a stochastic IPM. Stochasticity can be implemented
#'   in two ways in \code{ipmr}: \code{"kern"} resampling, and \code{"param"} resampling.}
#'
#'   \item{\code{kern}}{: used to denote an IPM that uses kernel resampling. Briefly,
#'   these models build all of the iteration kernels ahead of time and then choose one
#'   at random or in a user-specified order as they move from iteration to iteration. The
#'   user-specified population vector is multiplied by the chosen kernel and the result
#'   is multiplied by the next kernel for the desired number of iterations or until
#'   convergence is reached (depending on what you request following \code{make_ipm()}.}
#'   \item{\code{param}}{: used to denote parameter resampling. This generates distributions
#'   for each parameter based on user-specified functions supplied to \code{define_stoch_params()},
#'   and then resamples those distributions at each iteration. This will be a bit slower than
#'   \code{"kern"} resampling because iteration matrices need to be reconstructed
#'   from new parameters at every time step. The \code{vignette("classes", package = "ipmr")}
#'   contains more details and a number of suggestions for further reading on how to avoid this
#'   slow down.}
#'
#'}
#'
#' @references Metcalf et al. (2015). Statistical modelling of annual variation for inference on stochastic
#' population dynamics using Integral Projection Models. Methods in Ecology and Evolution, 6: 1007-1017
#'
#' \strong{UPDATE WITH MORE REFERENCES/MAKE COMPLETE LATER}
#'
#' @export

init_ipm <- function(class, ...) {

  out <- data.frame(
    id = character(0L),
    kernel_id = character(0L),
    domain = character(0L),
    state_var =  character(0L),
    int_rule = character(0L),
    evict = logical(0L),
    evict_type = character(0L),
    pop_state =  character(0L),
    env_state =  character(0L),
    params =  character(0L)
  )

  class(out) <- c(class, 'proto_ipm', class(out))

  return(out)
}

