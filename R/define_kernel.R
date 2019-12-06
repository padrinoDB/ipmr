#' @title Functions to initialize and define IPM kernels
#' @rdname kernel-definitions
#'
#' @description Adds a new kernel to the \code{proto_ipm} structure. Different
#' classes of IPMs may have many or only a few kernels. Each one requires its
#' own call to this function, though there are some exceptions, namely for kernels
#' derived from hierarchical models (e.g. vital rate models fit across plots and years).
#' See Details for more info
#'
#' @param proto_ipm The name of object you wish to append the new kernel to.
#' @param name The name of the new kernel - can be anything as long as it is
#' unique to the model.
#' @param formula A bare expression specifying the form of the kernel. See below
#' for examples
#' @param family The type of kernel. Options are \code{"CC"} for continuous to continuous
#' transitions, \code{"DC"} for discrete to continuous (e.g. emergence from a seedbank),
#' \code{"CD"} for continuous to discrete (e.g. entering a seedbank), and \code{"DD"} for
#' discrete to discrete (e.g. stasis in a seedbank).
#' @param ... For \code{define_kernel}, set of named expressions that correspond
#' to vital rates in \code{formula}. For \code{define_K}, a set of named expressions
#' that relate the population state at T + 1 to the population state at T. Alternatively,
#' can be an expression where the left hand side is the name of the kernel specified
#' in \code{name} and the right hand side only describes the structure of the
#' iteration kernel. In all cases, suffix expansion of hierarchical models is supported.
#' Specific details for each case are in their respecitive function's sections
#' below.
#' @param data_list A list of named values that correspond to constants in the formula.
#' You do not need to specify vectors corresponding to the domains here.
#' @param states A character vector containing the names of each state variable used in
#' the kernel.
#' @param has_hier_effs A logical indicating whether or not the kernel and/or its
#' underlying vital rates are structured with hierarchical effects. If so and you
#' specify either the functional forms or the exact parameter values you want to use,
#' the kernels will be automatically split according to the notation in the vital
#' rates and kernel formulae and multiple kernels will be built. See the vignette
#' on the syntax for this feature for more details (\code{vignettes(
#' 'hierarchical-notation', package = 'ipmr')}).
#' @param levels_hier_effs A named list with vectors corresponding the various levels
#' the hierarchical variable can take. Entries in this list should be a single
#' vector and should be character or integer typed.
#' @param evict A logical indicating whether an eviction correction should be applied
#' to the kernel. It is generally recommended to use a function in the vital rate
#' definition as opposed to using this option Default is \code{FALSE}.
#' @param evict_fun If \code{evict == TRUE}, then a function that corrects for it.
#' Currently, the only implemented function is \code{truncated_distributions}.
#' This works by modifying the functional form of the vital rate expressions
#' \code{...}, and so supplying your own function in this slot will be difficult.
#' One can also specify \code{usr_funs} function that performs the correction
#' during the numerical implementation of the model itself. In that case,
#' set \code{evict} to \code{FALSE}.
#'
#' @section \code{define_k}:
#'
#' The preferred method of defining a \code{K} kernel is to use the left
#' hand side of the \code{...} to reference the population vectors that the right
#' hand side creates (e.g. \code{n_T_1 = (P + F) * n_T}). This enables powerful
#' iteration-based methods to work properly. On the other hand, these iteration
#' based methods can be quite time consuming, and many applications only require
#' an iteration matrix while not necessarily requiring the population vectors (e.g.
#' calculations of deterministic population growth rate). In those cases, the
#' \code{...} can contain something like \code{K = P + F}. In this case, the left
#' hand side of the expression should match the \code{name} argument to \code{define_k}.
#'
#' @section \code{define_kernel}:
#'
#' Describe vital rate expression syntax - include suffix expansion of hier_effs
#'
#' @return All functions described here return a \code{proto_ipm}.
#'
#' @importFrom purrr map_chr map
#' @importFrom rlang := list2 enquo enquos parse_expr parse_exprs quo_text
#' quo_is_null is_quosure
#'
#' @export


define_kernel <- function(proto_ipm,
                          name,
                          formula,
                          family,
                          ...,
                          data_list = list(),
                          states,
                          has_hier_effs = FALSE,
                          levels_hier_effs = list(),
                          evict = FALSE,
                          evict_fun = NULL) {

  cls <- class(proto_ipm)
  # Capture formulas and convert to text
  formula <- rlang::enquo(formula)
  vr_quos <- rlang::enquos(...,
                           .named = TRUE,
                           .homonyms = "error",
                           .check_assign = TRUE)

  # make sure eviction function is correctly specified

  evict_fun <- rlang::enquo(evict_fun)
  evict_fun <- .check_evict_fun(evict, evict_fun)

  # protos store text rather than quos. They get converted back into quos later
  form_text <- rlang::quo_text(formula)
  vr_text <- lapply(vr_quos, rlang::quo_text)

  # retain names
  names(vr_text) <- names(vr_quos)

  # Param_tree should always contain these four entries, regardless of class.
  # pop_states and env_states get defined separately.

  param_tree <- list(formula = form_text,
                     family = family,
                     vr_text = vr_text,
                     params = data_list)

  if(!methods::hasArg(levels_hier_effs)) levels_hier_effs <- list(levels = NA)

  temp <- data.frame(
    id = 'A1',
    kernel_id = name,
    domain = I(list(NA_character_)),
    state_var = I(rlang::list2(!!name := states)),
    int_rule = NA_character_,
    evict = evict,
    evict_fun = I(list(evict_fun)),
    pop_state = I(list(NA_character_)),
    env_state = I(list(NA_character_)),
    has_hier_effs = has_hier_effs,
    levels_hier_effs = I(rlang::list2(levels_hier_effs)),
    params = I(rlang::list2(!! name := param_tree)),
    stringsAsFactors = FALSE
  )

  out <- rbind(proto_ipm,
               temp,
               stringsAsFactors = FALSE)

  class(out) <- cls
  return(out)
}


#' @rdname kernel-definitions
#'
#'
#' @export

define_k <- function(proto_ipm,
                     name,
                     family,
                     ...,
                     data_list = list(),
                     states,
                     has_hier_effs = FALSE,
                     levels_hier_effs = list(),
                     evict = FALSE,
                     evict_fun = NULL) {

  cls <- class(proto_ipm)
  forms <- rlang::enquos(...,
                         .named = TRUE,
                         .homonyms = "error",
                         .check_assign = TRUE)

  .check_k_def(proto_ipm,
               name = name,
               family = family)

  # make sure eviction function is correctly specified
  evict_fun <- rlang::enquo(evict_fun)
  evict_fun <- .check_evict_fun(evict, evict_fun)

  # protos store text rather than quos. They get converted back into quos later
  forms_text <- lapply(forms, function(x) {
    temp <- rlang::quo_text(x)
    out <- gsub('n_t$', 'pop_state_t', temp)
    return(out)
  })

  names(forms_text) <- gsub('^n_', 'pop_state_', names(forms_text))

  # retain names
  names(forms_text) <- names(forms)

  # Param_tree should always contain these four entries, regardless of class.
  # pop_states and env_states get defined separately. In the case of define_k,
  # vr_text is no longer relevant - expressions all get tossed into formula and
  # will be evaluated in the kernel environment anyway.

  param_tree <- list(formula = forms_text,
                     family = family,
                     vr_text = NA_character_,
                     params = data_list)

  if(!methods::hasArg(levels_hier_effs)) levels_hier_effs <- list(levels = NA)

  temp <- data.frame(
    id = 'A1',
    kernel_id = name,
    domain = I(list(NA_character_)),
    state_var = I(rlang::list2(!!name := states)),
    int_rule = NA_character_,
    evict = evict,
    evict_fun = I(list(evict_fun)),
    pop_state = I(list(NA_character_)),
    env_state = I(list(NA_character_)),
    has_hier_effs = has_hier_effs,
    levels_hier_effs = I(rlang::list2(levels_hier_effs)),
    params = I(rlang::list2(!! name := param_tree)),
    stringsAsFactors = FALSE
  )

  out <- rbind(proto_ipm,
               temp,
               stringsAsFactors = FALSE)

  class(out) <- cls
  return(out)


  return(out)
}


.check_k_def <- function(proto_ipm, name, family) {

  if(!grepl('K|k', name)) {
    stop("'name' passed to define_k must be of the form K, k, K_effects, or k_effects.")
  }

  # Possible family types. These need to be updated eventually, this is mostly
  # a placeholder
  families <- c('IPM')

  if(!family %in% families) {
    stop("'family' should be one of the following options: ", families)
  }

  invisible(proto_ipm)

}

.check_evict_fun <- function(evict, fun) {


  # need to supply a function if you want to correct for eviction!

  if(evict & rlang::quo_is_null(fun)) {

    stop('"evict" is TRUE but no fun supplied!')

    # if we have one, we need to get the name of the corrected object

  } else if(!rlang::quo_is_null(fun)) {

    text <- rlang::quo_text(fun)
    nm <- strsplit(text, '\\(|,|\\)')[[1]][2]

    fun <- list(fun)
    names(fun) <- nm

  } else if(rlang::quo_is_null(fun)) {
    fun <- NA_character_
  }



  return(fun)
}
