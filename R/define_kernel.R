#' @title Functions to initialize and define IPM kernels
#' @rdname kernel-definitions
#'
#' @description Adds a new kernel to the \code{proto_ipm} structure.
#'
#' @param proto_ipm The name of the model.
#' @param name The name of the new kernel.
#' @param formula A bare expression specifying the form of the kernel.
#' @param family The type of kernel. Options are \code{"CC"} for continuous to continuous
#' transitions, \code{"DC"} for discrete to continuous (e.g. emergence from a seedbank),
#' \code{"CD"} for continuous to discrete (e.g. entering a seedbank), and \code{"DD"} for
#' discrete to discrete (e.g. stasis in a seedbank).
#' @param ... A set of named expressions that correspond
#' to vital rates in \code{formula}. Parameter set index syntax is supported.
#' @param data_list A list of named values that correspond to constants in the formula
#' and vital rate expressions in \code{...}.
#' @param states A list with character vector containing the names of each state
#' variable used in the kernel.
#' @param uses_par_sets A logical indicating whether or not the parameters in the kernel and/or its
#' underlying vital rates are derived from sets. See the
#' introduction vignette for this feature for more details
#' (\code{vignettes(ipmr-introduction', package = 'ipmr')}, and
#' \code{vignettes( index-notation', package = 'ipmr')}).
#' @param par_set_indices A named list with vectors corresponding to the values
#' the index variable can take. The names should match the suffixes used
#' in the vital rate expressions.
#' @param evict_cor A logical indicating whether an eviction correction should be applied
#' to the kernel.
#' @param evict_fun If \code{evict_cor = TRUE}, then a function that corrects for it.
#' Currently, only \code{truncated_distributions} and \code{discrete_extrema} are
#' possible.
#' @param age_indices If \code{init_ipm(uses_age = TRUE)}, a list with possibly
#' 2 entries: 1. \code{"age"}: the range
#' of possible ages in the model and, optionally, 2. \code{"max_age"}: the maximum
#' age individuals in the model can attain. Otherwise, not used.
#' @param integrate For \code{simple_*} models, this controls whether a \code{"d_z"}
#' is automatically appended to the \code{formula} argument. When \code{TRUE},
#' this automatically generates \code{formula * d_z}. There may be some cases where
#' this behavior is not desirable. Set this to \code{FALSE} and specify the correct
#' form if needed. The default is \code{TRUE}. This argument is ignored for
#' all \code{general_*} models.
#'
#'
#' @details
#' Different classes of IPMs may have many or only a few kernels. Each
#' one requires its own call to \code{define_kernel}, though there are some exceptions,
#' namely for kernels derived for models derived from parameter sets (e.g. vital
#' rate models fit across plots and years).
#'
#' A much more complete overview of how to generate kernels is provided in
#' \code{vignette("ipmr-introduction", "ipmr")}.
#'
#' @return A \code{proto_ipm}.
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
                          uses_par_sets = FALSE,
                          par_set_indices = list(),
                          age_indices      = list(),
                          evict_cor= FALSE,
                          evict_fun = NULL,
                          integrate = TRUE) {

  cls <- class(proto_ipm)

  simple <- any(grepl("simple", cls))

  if(missing(family) && simple) family <- "CC"

  integrate <- integrate && simple

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

  # Param_tree should always contain these five entries, regardless of class.
  # pop_states and env_states get defined separately. .protect_model detects
  # model objects and keeps them from getting flattened beyond

  data_list <- lapply(
    data_list,
    function(x) {
      x <- .protect_model(x)

      na_test <- attr(x, "na_ok")

      if(! na_test) {
        warning("'data_list' in 'define_kernel()' contains NAs. Is this correct?",
                call. = FALSE)
      }

      return(x)
    })

  param_tree <- list(formula = form_text,
                     family = family,
                     vr_text = vr_text,
                     params = data_list,
                     integrate = integrate)

  if(!methods::hasArg(par_set_indices)) par_set_indices <- list(levels = NA)

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
    uses_par_sets    = uses_par_sets,
    par_set_indices = I(rlang::list2(par_set_indices)),
    uses_age          = ifelse(.uses_age(proto_ipm), TRUE, FALSE),
    age_indices      = I(rlang::list2(age_indices)),
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

#' @noRd

.define_k <- function(proto_ipm,
                     name,
                     family,
                     ...,
                     data_list = list(),
                     states,
                     uses_par_sets = FALSE,
                     par_set_indices = list(),
                     age_indices      = list(),
                     evict_cor = FALSE,
                     evict_fun = NULL,
                     integrate = FALSE) {

  UseMethod(".define_k")

}

#' @noRd

.define_k.default <- function(proto_ipm,
                             name,
                             family,
                             ...,
                             data_list = list(),
                             states,
                             uses_par_sets = FALSE,
                             par_set_indices = list(),
                             evict_cor = FALSE,
                             evict_fun = NULL,
                             integrate = FALSE) {

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
                     params = data_list,
                     integrate = integrate)

  if(!methods::hasArg(par_set_indices)) par_set_indices <- list(levels = NA)

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
    uses_par_sets    = uses_par_sets,
    par_set_indices = I(rlang::list2(par_set_indices)),
    uses_age          = FALSE,
    age_indices      = I(list(NA_character_)),
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

#' @noRd

.define_k.age_x_size <- function(proto_ipm,
                                name,
                                family,
                                ...,
                                data_list = list(),
                                states,
                                uses_par_sets = FALSE,
                                par_set_indices = list(),
                                age_indices      = list(),
                                evict_cor = FALSE,
                                evict_fun = NULL,
                                integrate = FALSE) {
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
                     params = data_list,
                     integrate = integrate)

  if(!methods::hasArg(par_set_indices)) par_set_indices <- list(levels = NA)

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
    uses_par_sets    = uses_par_sets,
    par_set_indices = I(rlang::list2(par_set_indices)),
    uses_age          = TRUE,
    age_indices      = I(rlang::list2(age_indices)),
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

#' @noRd

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
