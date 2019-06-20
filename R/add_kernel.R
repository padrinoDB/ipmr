#' @title Add a new kernel to the IPM
#'
#' @description Adds a new kernel to the \code{proto_ipm} structure. Different
#' classes of IPMs may have many or only a few kernels. Each one requires its
#' own call to this function, though there are some exceptions, namely for kernels
#' derived from hierarchical models (e.g. vital rate models fit across plots and years).
#' See Details for more info
#'
#' @param proto_ipm The name of object you wish to append the new kernel to.
#' @param name The name of the new kernel. Survival and growth kernels should be
#' named \code{"P"} for now. Fecundity should be named \code{"F"}, and clonal
#' kernels should be named \code{"C"} (again, just for now. \strong{ need to work on
#' generalizing that).}
#' @param formula An bare expression specifying the form of the kernel. See below
#' for examples
#' @param family The type of kernel. Options are \code{"CC"} for continuous to continuous
#' transitions, \code{"DC"} for discrete to continuous (e.g. emergence from a seedbank),
#' \code{"CD"} for continuous to discrete (e.g. entering a seedbank), and \code{"DD"} for
#' discrete to discrete (e.g. stasis in a seedbank).
#' @param ... A set of named expressions that correspond to vital rates in \code{formula}.
#' @param data_list A list of named values that correspond to constants in the formula.
#' You do not need to specify vectors corresponding to the domains here.
#' @param state_list A list containing the names of each state variable used in
#' the kernel.
#' @param dom_start The name of the state variable for the kernel at time \emph{t}.
#' @param dom_end The name of the state variable for the kernel at time \emph{t+1}.
#' This is usually the same as \code{dom_start}, but general IPMs
#' with discrete classes or IPMs that move from one state to another (e.g. tree
#' seedling going from a height domain to a DBH domain at T+1) may have another
#' value here. For cases with a discrete stage, kernels moving individuals from
#' discrete to continuous should have a state variable entered here and an \code{NA}
#' for \code{dom_start}. For kernels moving from continuous to discrete, vice versa.
#' For discrete to discrete, both are \code{NA}.
#' @param int_rule The integration rule to be used for the kernel. The default is
#' "midpoint". "trapezoid" and "g-l" (Gauss-Legendre) will be implemented as well.
#' If "g-l", additional arguments need to be supplied (\strong{Work on this later!!}).
#' @param evict A logical indicating whether an eviction correction should be applied
#' to the kernel. Default is \code{TRUE}.
#' @param evict_fun If \code{evict == TRUE}, then a function that corrects for it.
#' Currently, \code{truncated_distributions} and \code{rescale_kernel} are the
#' only built-in functions, but user specified ones should work as well.
#'
#' @details \strong{BLAH BLAH BLAH}
#'
#' @importFrom purrr map_chr map
#' @importFrom rlang := list2 enquo enquos parse_expr parse_exprs quo_text
#' quo_is_null
#'
#' @export


add_kernel <- function(proto_ipm, name, formula, family,
                       ..., data_list, state_list,
                       dom_start, dom_end,
                       int_rule = "midpoint",
                       evict = TRUE,
                       evict_fun = NULL) {

  # Capture formulas and convert to text
  formula <- rlang::enquo(formula)
  vr_quos <- rlang::enquos(...,
                           .named = TRUE,
                           .homonyms = "error",
                           .check_assign = TRUE)

  evict_fun <- rlang::enquo(evict_fun)
  form_text <- rlang::quo_text(formula)
  vr_text <- lapply(vr_quos, rlang::quo_text)

  # need to supply a function if you want to correct for eviction!
  if(evict & rlang::quo_is_null(evict_fun)) {
    stop('"evict" is TRUE but no evict_fun supplied!')

  # if we have one, we need to get the name of the corrected object
  } else if(!rlang::quo_is_null(evict_fun)) {
    text <- rlang::quo_text(evict_fun)
    nm <- strsplit(text, '\\(|\\)')[[1]][2]

    evict_fun <- list(evict_fun)
    names(evict_fun) <- nm
  }


  # retain names
  names(vr_text) <- names(vr_quos)

  # convert state list into usable domain information
  domain_info <- .get_state_info(state_list, dom_start, dom_end)

  # Param_tree should always contain these four entries, regardless of class.
  # pop_states and env_states get defined separately.

  param_tree <- list(formula = form_text,
                     family = family,
                     vr_text = vr_text,
                     params = data_list)

  out <- data.frame(
    id = 'A1',
    kernel_id = name,
    domain = I(rlang::list2(!! name := domain_info)),
    state_var = I(rlang::list2(!! name := state_list)),
    int_rule = int_rule,
    evict = evict,
    evict_fun = I(list(evict_fun)),
    pop_state = I(list(NA_character_)),
    env_state = I(list(NA_character_)),
    params = I(rlang::list2(!! name := param_tree))
  )

  rbind(proto_ipm, out)
}


#' @noRd
.get_state_info <- function(state_list, dom_start, dom_end) {

  # match names, then get info. Otherwise, generate an NA. the domain name
  # will always be first entry.

  if(!is.na(dom_start)) {
    state_ind <- which(grepl(dom_start, purrr::map_chr(state_list, ~.x[1])))
    start_state_info <- rep(NA_real_, 3)

    dom_start <- paste(dom_start, "_1", sep = "")
  } else {
    start_state_info <- NA_real_
    dom_start <- 'start_not_applicable'
  }

  if(!is.na(dom_end)) {
    state_ind <- which(grepl(dom_end, purrr::map_chr(state_list, ~.x[1])))
    end_state_info <- rep(NA_real_, 3)

    dom_end <- paste(dom_end, "_2", sep = "")
  } else {
    end_state_info <- NA_real_
    dom_end <- 'end_not_applicable'
  }

  out <- rlang::list2(!!dom_start := start_state_info,
                      !!dom_end := end_state_info)
  return(out)
}

