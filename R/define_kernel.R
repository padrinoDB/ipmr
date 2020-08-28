#' @title Functions to initialize and define IPM kernels
#' @rdname kernel-definitions
#'
#' @description Adds a new kernel to the \code{proto_ipm} structure.
#'
#' @param proto_ipm The name of the model.
#' @param name The name of the new kernel. There are only two rules: 1. For
#' \code{define_k}, the name must start with \code{"K"} or \code{"k"}.
#' 2. For \code{define_kernel}, the name cannot start with \code{"K"} or
#' \code{"k"}.
#' @param formula A bare expression specifying the form of the kernel. See Details.
#' @param family The type of kernel. Options are \code{"CC"} for continuous to continuous
#' transitions, \code{"DC"} for discrete to continuous (e.g. emergence from a seedbank),
#' \code{"CD"} for continuous to discrete (e.g. entering a seedbank), and \code{"DD"} for
#' discrete to discrete (e.g. stasis in a seedbank).
#' @param ... For \code{define_kernel}, a set of named expressions that correspond
#' to vital rates in \code{formula}. For \code{define_k}, a set of named expressions
#' that relate the population state at T + 1 to the population state at T. Alternatively,
#' can be an expression where the left hand side is the name of the kernel specified
#' in \code{name} and the right hand side only describes the structure of the
#' iteration kernel. In all cases, suffix expansion of hierarchical models is supported.
#' Specific details for each case are in their respecitive function's sections
#' below.
#' @param data_list A list of named values that correspond to constants in the formula.
#' You do not need to specify vectors corresponding to the domains here.
#' @param states A list with character vector containing the names of each state
#' variable used in the kernel.
#' @param has_hier_effs A logical indicating whether or not the kernel and/or its
#' underlying vital rates are structured with hierarchical effects. If so and you
#' specify either the functional forms or the exact parameter values you want to use,
#' the kernels will be automatically split according to the notation in the vital
#' rates and kernel formulae and multiple kernels will be built. See the vignette
#' on the syntax for this feature for more details (\code{vignettes(
#' 'hierarchical-notation', package = 'ipmr')}).
#' @param levels_hier_effs A named list with vectors corresponding to the various levels
#' the hierarchical variable can take. The names should match the suffixes used
#' in the vital rate expressions.
#' @param evict_cor A logical indicating whether an eviction correction should be applied
#' to the kernel.
#' @param evict_fun If \code{evict_cor== TRUE}, then a function that corrects for it.
#' Currently, the only implemented function is \code{truncated_distributions}.
#' This works by modifying the functional form of the vital rate expressions
#' \code{...}, and so supplying your own function in this slot will be difficult.
#' You can also either specify a \code{usr_funs} function that performs the correction
#' during the implementation of the model itself, or incorporate the correction
#' into individual vital rate expressions. In either of those cases,
#' set \code{evict_cor} to \code{FALSE}.
#'
#'
#' @details
#' Different classes of IPMs may have many or only a few kernels. Each
#' one requires its own call to this function, though there are some exceptions,
#' namely for kernels derived from hierarchical models (e.g. vital rate models
#' fit across plots and years). \code{remove_k} is a helper for re-building models
#' but with a different \code{K} kernel.
#'
#' Internally, \code{ipmr} generates the \code{levels_hier_effs} using \code{expand.grid}.
#' This  means that every combination of the values in \code{levels_hier_effs}
#' must exist in the \code{data_list}, or the model will fail with an error along
#' the lines of \code{Error in eval_tidy: object 'x_yz' not found}. In order to
#' exclude levels that do not exist in your data, you can add a vector to the list
#' in \code{levels_hier_effs} called \code{drop_levels}. This should contain the
#' values you wish to exclude as a character vector. For example, if data are
#' collected for \code{sites = c("a", "b", "c")}, and \code{years = c(2005:2008)},
#' but there is no data from site \code{"a"} in \code{2007}, we can use
#' \code{levels_hier_effs = list(site = c("a", "b", "c"), year = c(2005:2008),
#' drop_levels = c("a_2007"))}.
#'
#' \strong{\code{define_k}}
#'
#' The preferred method of defining a \code{K} kernel is to use the left
#' hand side of the \code{...} to reference the population vectors that the right
#' hand side creates (e.g. \code{n_T_1 = (P + F) \%*\% n_T}). This enables
#' iteration-based methods to work properly. On the other hand, these iterations
#' can be quite time consuming, and some applications only require
#' an iteration matrix while not necessarily requiring the population vectors (e.g.
#' calculations of deterministic population growth rate). In those cases, the
#' \code{...} can contain something like \code{K = P + F}. In this case, the left
#' hand side of the expression should match the \code{name} argument to \code{define_k}.
#'
#' \strong{\code{define_kernel}}
#'
#' \code{define_kernel} generates most of the information needed to create an IPM
#' kernel. There are a few requirements - \code{name}, \code{family},
#' \code{formula}, and \code{states} must not be empty. The \code{formula} should be an expression
#' for how vital rates produce a kernel (e.g. \code{formula = S * G}). The
#' \code{...} should be a set of named expressions that correspond to vital rate
#' expressions (e.g. \code{G = dnorm(size_2, mean_size, sd_size)}). See the vignettes
#' for more information on how to get started.
#'
#' @return All functions described here return a  \code{proto_ipm}.
#'
#' @importFrom purrr map_chr map
#' @importFrom rlang := list2 enquo enquos parse_expr parse_exprs quo_text
#' quo_is_null is_quosure
#' @importFrom methods hasArg
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
                          levels_ages      = list(),
                          evict_cor= FALSE,
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
  evict_fun <- .check_evict_fun(evict_cor, evict_fun)

  # protos store text rather than quos. They get converted back into quos later
  form_text <- rlang::quo_text(formula)
  vr_text <- lapply(vr_quos, rlang::quo_text)

  # retain names
  names(vr_text) <- names(vr_quos)

  # Param_tree should always contain these four entries, regardless of class.
  # pop_states and env_states get defined separately. .protect_model detects
  # model objects and keeps them from getting flattened beyond

  data_list <- lapply(data_list, .protect_model)

  param_tree <- list(formula = form_text,
                     family = family,
                     vr_text = vr_text,
                     params = data_list)

  if(!methods::hasArg(levels_hier_effs)) levels_hier_effs <- list(levels = NA)

  temp <- data.frame(
    id               = 'A1',
    kernel_id        = name,
    domain           = I(list(NA_character_)),
    state_var        = I(rlang::list2(!!name := states)),
    int_rule         = NA_character_,
    evict            = evict_cor,
    evict_fun        = I(list(evict_fun)),
    pop_state        = I(list(NA_character_)),
    env_state        = I(list(NA_character_)),
    has_hier_effs    = has_hier_effs,
    levels_hier_effs = I(rlang::list2(levels_hier_effs)),
    has_age          = ifelse(.has_age(proto_ipm), TRUE, FALSE),
    levels_ages      = I(rlang::list2(levels_ages)),
    params           = I(rlang::list2(!! name := param_tree)),
    usr_funs         = I(list(NA_character_)),
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
#' @export

define_k <- function(proto_ipm,
                     name,
                     family,
                     ...,
                     data_list = list(),
                     states,
                     has_hier_effs = FALSE,
                     levels_hier_effs = list(),
                     levels_ages      = list(),
                     evict_cor = FALSE,
                     evict_fun = NULL) {

  UseMethod("define_k")

}

#' @rdname kernel-definitions
#'
#' @export

define_k.default <- function(proto_ipm,
                             name,
                             family,
                             ...,
                             data_list = list(),
                             states,
                             has_hier_effs = FALSE,
                             levels_hier_effs = list(),
                             evict_cor = FALSE,
                             evict_fun = NULL) {

  cls <- class(proto_ipm)
  forms <- rlang::enquos(...,
                         .named = TRUE,
                         .homonyms = "error",
                         .check_assign = TRUE)

  name_suff <- substr(name, 2, nchar(name))
  name      <- toupper(substr(name, 1, 1)) %>%
    paste(., name_suff, sep = "")

  .check_k_def(proto_ipm,
               name = name,
               family = family)

  # make sure eviction function is correctly specified
  evict_fun <- rlang::enquo(evict_fun)
  evict_fun <- .check_evict_fun(evict_cor, evict_fun)

  # protos store text rather than quos. They get converted back into quos later
  forms_text <- lapply(forms, function(x) {
    temp <- rlang::quo_text(x)
    out <- gsub('n_t$', 'pop_state_t', temp)
    return(out)
  })


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
    kernel_id        = name,
    domain           = I(list(NA_character_)),
    state_var        = I(rlang::list2(!! name := states)),
    int_rule         = NA_character_,
    evict            = evict_cor,
    evict_fun        = I(list(evict_fun)),
    pop_state        = I(list(NA_character_)),
    env_state        = I(list(NA_character_)),
    has_hier_effs    = has_hier_effs,
    levels_hier_effs = I(rlang::list2(levels_hier_effs)),
    has_age          = FALSE,
    levels_ages      = I(list(NA_character_)),
    params           = I(rlang::list2(!! name := param_tree)),
    usr_funs         = I(list(NA_character_)),
    stringsAsFactors = FALSE
  )

  out <- rbind(proto_ipm,
               temp,
               stringsAsFactors = FALSE)

  class(out) <- cls
  return(out)

}

#' @rdname kernel-definitions
#' @param levels_ages A list with possible two entries: 1. \code{"age"}: the range
#' of ages, and, optionally, 2. \code{"max_age"}: The maximum age for the model.
#'
#' @export

define_k.age_x_size <- function(proto_ipm,
                                name,
                                family,
                                ...,
                                data_list = list(),
                                states,
                                has_hier_effs = FALSE,
                                levels_hier_effs = list(),
                                levels_ages      = list(),
                                evict_cor = FALSE,
                                evict_fun = NULL) {
  cls <- class(proto_ipm)
  forms <- rlang::enquos(...,
                         .named = TRUE,
                         .homonyms = "error",
                         .check_assign = TRUE)

  name_suff <- substr(name, 2, nchar(name))
  name      <- toupper(substr(name, 1, 1)) %>%
    paste(., name_suff, sep = "")

  .check_k_def(proto_ipm,
               name = name,
               family = family)

  # make sure eviction function is correctly specified
  evict_fun <- rlang::enquo(evict_fun)
  evict_fun <- .check_evict_fun(evict_cor, evict_fun)

  # protos store text rather than quos. They get converted back into quos later
  forms_text <- lapply(forms, function(x) {
    temp <- rlang::quo_text(x)
    out <- gsub('n_t$', 'pop_state_t', temp)
    return(out)
  })


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
    kernel_id        = name,
    domain           = I(list(NA_character_)),
    state_var        = I(rlang::list2(!! name := states)),
    int_rule         = NA_character_,
    evict            = evict_cor,
    evict_fun        = I(list(evict_fun)),
    pop_state        = I(list(NA_character_)),
    env_state        = I(list(NA_character_)),
    has_hier_effs    = has_hier_effs,
    levels_hier_effs = I(rlang::list2(levels_hier_effs)),
    has_age          = TRUE,
    levels_ages      = I(rlang::list2(levels_ages)),
    params           = I(rlang::list2(!! name := param_tree)),
    usr_funs         = I(list(NA_character_)),
    stringsAsFactors = FALSE
  )

  out <- rbind(proto_ipm,
               temp,
               stringsAsFactors = FALSE)

  class(out) <- cls
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

#' @rdname kernel-definitions
#' @export

remove_k <- function(proto_ipm) {

  # Using substr() %>% grepl() to avoid a case where P_sink or something
  # flags a P kernel for removal. kernel_ids for Ks must start with K | k,
  # so just searching the first characters *should* be safer. Will fail if
  # for some reason, the beginning of another kernel's name is also a K.

  K_test   <- vapply(proto_ipm$kernel_id,
                     function(x) substr(x, 1, 1),
                     character(1L))

  keep_ind <- !grepl('K|k', K_test)

  out      <- proto_ipm[keep_ind, ]

  return(out)

}

.check_evict_fun <- function(evict_cor, fun) {


  # need to supply a function if you want to correct for eviction!

  if(evict_cor && rlang::quo_is_null(fun)) {

    stop('"evict_cor" is TRUE but no fun supplied!')

    # if we have one, we need to get the name of the corrected object

  } else if(!rlang::quo_is_null(fun)) {

    temp <- rlang::call_args(fun)

    nm   <- temp[[length(temp)]][1]

    fun <- list(fun)
    names(fun) <- nm

  } else if(rlang::quo_is_null(fun)) {
    fun <- NA_character_
  }

  return(fun)
}
