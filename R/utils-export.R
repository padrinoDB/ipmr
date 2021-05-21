# General exported utilities

#' @title Helpers for IPM construction
#' @inheritParams define_kernel
#' @param ... Named expressions. See Details for more information on their usage in
#' each \code{define_*} function.
#'
#' @param pop_vectors If the population vectors are already pre-defined (i.e. are
#' not defined by a function passed to \code{...}), then they can
#' be passed as a named list here.
#'
#' @details
#' These are helper functions to define IPMs. They are used after defining the kernels,
#' but before calling \code{make_ipm()} They are meant to be called in the
#' following order:
#'
#' \enumerate{
#'
#'   \item \code{define_impl()}
#'
#'   \item \code{define_domains()}
#'
#'   \item \code{define_pop_state()}
#'
#'   \item \code{define_env_state()}
#'
#' }
#'
#' The order requirement is so that information is correctly matched to each kernel.
#' Below are specific details on the way each works.
#'
#' \strong{\code{define_impl}}
#'
#' This has two arguments - \code{proto_ipm} (the model object you wish to work with),
#' and the \code{kernel_impl_list}. The format of the \code{kernel_impl_list} is:
#' names of the list should be kernel names, and each kernel should have 3 entries:
#' \code{int_rule}, \code{state_start}, and \code{state_end}. See examples.
#'
#' \strong{\code{define_domains}}
#'
#' If the \code{int_rule = "midpoint"}, the \code{...} entries are vectors of
#' length 3 where the name corresponds to the
#' state variable, the first entry is the lower bound of the domain, the second
#' is the upper bound of the domain, and the third entry is the number of
#' meshpoints. Other \code{int_rule}s are not yet implemented, so for now this is the
#' only format they can take. See examples.
#'
#' \strong{\code{define_pop_state}}
#'
#' This takes either calls to functions in the \code{...}, or a pre-generated
#' list of vectors in the \code{pop_vectors}. The names used
#' for each entry in \code{...} and/or for the \code{pop_vectors} should be
#' \code{n_<state_variable>}. See examples.
#'
#' \strong{\code{define_env_state}}
#'
#' Takes expressions that generate values for environmental covariates at each
#' iteration of the model in \code{...}. The \code{data_list} should contain any
#' parameters that the function uses, as well as the function itself. The
#' functions should return named lists. Names in that list can be referenced in
#' vital rate expressions and/or kernel formulas.
#'
#' @return All \code{define_*} functions return a proto_ipm. \code{make_impl_args_list}
#' returns a list, and so must be used within a call to \code{define_impl} or
#' before initiating the model creation procedure.
#'
#' @examples
#'
#' # Example with kernels named "P" and "F", and a domain "z"
#'
#' kernel_impl_list <- list(P = list(int_rule = "midpoint",
#'                                   state_start = "z",
#'                                   state_end = "z"),
#'                          F = list(int_rule = "midpoint",
#'                                   state_start = "z",
#'                                   state_end = "z"))
#'
#' # an equivalent version using make_impl_args_list
#'
#' kernel_impl_list <- make_impl_args_list(
#'      kernel_names = c("P", "F"),
#'      int_rule     = c("midpoint", "midpoint"),
#'      state_start  = c("z", "z"),
#'      state_end    = c("z", "z")
#' )
#'
#' data(sim_di_det_ex)
#'
#' proto_ipm <- sim_di_det_ex$proto_ipm
#'
#' # define_domains
#'
#' lower_bound <- 1
#' upper_bound <- 100
#' n_meshpoints <- 50
#'
#'
#' define_domains(proto_ipm, c(lower_bound, upper_bound, n_meshpoints))
#'
#' # define_pop_state with a state variable named "z". Note that "n_" is prefixed
#' # to denote that it is a population state function!
#'
#' define_pop_state(proto_ipm, n_z = runif(100))
#'
#' # alternative, we can make a list before starting to make the IPM
#'
#' pop_vecs <- list(n_z = runif(100))
#'
#' define_pop_state(proto_ipm, pop_vectors = pop_vecs)
#'
#' # define_env_state. Generates a random draw from a known distribution
#' # of temperatures.
#'
#' env_sampler <- function(env_pars) {
#'
#'   temp <- rnorm(1, env_pars$temp_mean, env_pars$temp_sd)
#'
#'   return(list(temp = temp))
#'
#' }
#'
#' env_pars <- list(temp_mean = 12, temp_sd = 2)
#'
#' define_env_state(
#'  proto_ipm,
#'  env_values = env_sampler(env_pars),
#'  data_list  = list(env_sampler = env_sampler,
#'                    env_pars    = env_pars)
#'
#' )
#'
#'
#' @rdname define_star
#' @importFrom rlang is_empty
#' @export

define_pop_state <- function(proto_ipm, ..., pop_vectors = list()) {

  pop_quos            <- rlang::enquos(...)

  temp                <- rlang::list2(!!! pop_quos, !!! pop_vectors)

  out                 <- Filter(Negate(rlang::is_empty), temp)

  # Catch a few cases where names me be defined incorrectly:

  nm_test             <- vapply(names(out), function(x) substr(x, 1, 1), character(1L))

  if(any(nm_test != "n")) {

    stop("All population state names must start with 'n_'",
         " (e.g. 'n_<stateVariable>)'.",
         call. = FALSE)

  }

  nm_test             <- vapply(names(out),
                                function(x) grepl("*_t$", x) | grepl("*_t_1$", x),
                                logical(1L))

  if(any(nm_test)) {

    stop("Detected '_t' attached to end of name supplied in 'define_pop_state()'.",
         "\nVariables in define_pop_state() automatically have '_t' and '_t_1'",
         " appended to them.\nPlease remove these suffixes!.")

  }

  names(out)          <- gsub('^n_', 'pop_state_', names(out))

  proto_ipm$pop_state <- list(out)

  return(proto_ipm)

}

#' @inheritParams define_pop_state
#' @param data_list A list of named values that contain data used in the expressions
#' in \code{...} in \code{define_env_state()}.
#'
#' @rdname define_star
#' @export

define_env_state <- function(proto_ipm, ..., data_list = list()) {

  env_quos            <- rlang::enquos(...)

  out                 <- list(env_quos = unlist(env_quos),
                              constants = data_list)

  proto_ipm$env_state <- list(out)

  return(proto_ipm)
}

#' @title Predict methods in ipmr
#' @rdname predict_methods
#'
#' @description This function is used when a \code{predict} method is incorporated
#' into the vital rate expressions of a kernel. Generally, ipmr can handle this
#' without any additional user effort, but some model classes will fail (often
#' with an obscure error message).
#' When this happens, \code{use_vr_model} can ensure that model object is
#' correctly represented in the \code{data_list}.
#'
#' @param model A fitted model object representing a vital rate. Primarily used to avoid
#' writing the mathematical expression for a vital rate, and using a \code{predict()}
#' method instead.
#'
#' @return A model object with a \code{"flat_protect"} attribute.
#'
#' @details ipmr usually recognizes model objects passed into the \code{data_list} argument
#' automatically. Unfortunately, sometimes it'll miss one, and the user will need
#' to manually protect it from the standard build process. This function
#' provides a wrapper around that process. Additionally, please file a bug
#' report here: \url{https://github.com/levisc8/ipmr/issues} describing what type
#' of model you are trying to use so it can be added to later versions of the
#' package.
#'
#' Wrap a model object in \code{use_vr_model} when building the \code{data_list}
#' to pass to \code{define_kernel}.
#'
#'
#' @examples
#'
#' data(iceplant_ex)
#'
#' grow_mod <- lm(log_size_next ~ log_size, data = iceplant_ex)
#' surv_mod <- glm(survival ~ log_size, data = iceplant_ex, family = binomial())
#'
#' data_list <- list(
#'   grow_mod = use_vr_model(grow_mod),
#'   surv_mod = use_vr_model(surv_mod),
#'   recruit_mean = 20,
#'   recruit_sd   = 5
#' )
#'
#' @export

use_vr_model <- function(model) {


  attr(model, "flat_protect") <- TRUE

  return(model)

}


#' @title Right/left multiplication
#'
#' @description Performs right and left multiplication.
#'
#' @param kernel,vectr \code{kernel} should be a bivariate kernel, \code{vectr}
#' should be a univariate trait distribution.
#'
#' @return \code{left_mult} returns \code{t(kernel) \%*\% vectr}. \code{right_mult}
#' returns \code{kernel \%*\% vectr}.
#'
#'
#' @export

right_mult <- function(kernel, vectr) {

  kernel %*% vectr

}

#' @rdname right_mult
#' @export

left_mult <- function(kernel, vectr) {

  t(kernel) %*% vectr

}

#' @title Raise a matrix to a power
#' @rdname matrix-power
#'
#' @description Raises a matrix \code{x} to the \code{y}-th power. \code{x ^ y} computes
#' element wise powers, whereas this computes \emph{y - 1} matrix multiplications.
#' \code{mat_power(x, y)} is identical to \code{x \%^\% y}.
#'
#' @param x A numeric or integer matrix.
#' @param y An integer.
#'
#' @return A matrix.
#'
#' @export
#'

`%^%` <- function(x, y) {

  if(!is_square(x)) {
    stop('not implemented for non-square matrices')
  }

  if(!is.integer(y)) {

    warning("`%^%` is coercing second argument to an integer",
            call. = FALSE)

    y <- as.integer(y)

  }

  init_dim <- dim(x)[1]

  use_list <- lapply(seq_len(y), function(a, b) b, b = x)

  init_i   <- diag(init_dim)

  out <- Reduce('%*%', use_list, init = init_i)

  return(out)

}

#' @rdname matrix-power
#'
#' @export

mat_power <- function(x, y) {

  return(x %^% y)

}

#' @rdname make_iter_kernel
#' @export

format_mega_kernel <- function(ipm, ...) {

  UseMethod("format_mega_kernel")

}

#' @rdname make_iter_kernel
#' @export

format_mega_kernel.default <- function(ipm, mega_mat, ...) {

  mega_mat    <- rlang::enquo(mega_mat)

  if(!rlang::quo_is_call(mega_mat)) {

    text     <- rlang::eval_tidy(mega_mat)

    if(length(text) > 1) {

      exprr <- syms(text)
      exprr <- rlang::call2("c", !!! exprr)

    } else{

      exprr    <- rlang::parse_expr(text)

    }

    mega_mat <- rlang::enquo(exprr)

  }

  if(any(ipm$proto_ipm$uses_par_sets)) {

    levs <- ipm$proto_ipm$par_set_indices[ipm$proto_ipm$uses_par_sets] %>%
      .flatten_to_depth(1L) %>%
      .[!duplicated(names(.))]

    if("drop_levels" %in% names(levs)) {

      ind <- which(names(levs) != "drop_levels")

      use_levs <- levs[ind]

    } else {

      use_levs <- levs

    }

    use_levs <- expand.grid(use_levs, stringsAsFactors = FALSE)

    out_nms <- temp <- character(dim(use_levs)[1])

    base_expr <- rlang::quo_text(mega_mat)
    base_name <- names(use_levs) %>% paste(collapse = "_")

    it <- 1

    for(i in seq_len(dim(use_levs)[2])) {

      nm <- names(use_levs)[i]

      for(j in seq_len(dim(use_levs)[1])) {

        if(i > 1) {
          base_expr <- temp[it]
          base_name <- out_nms[it]
        }

        temp[it] <- gsub(nm, use_levs[j, i], base_expr)
        out_nms[it] <- gsub(nm, use_levs[j, i], base_name)

        it <- it + 1

      } # end single var substitution - reset counter to modify exprs w/ adtl par_sets if present

      it <- 1

    }

    if("drop_levels" %in% names(levs)) {

      # Need to use fuzzy matching for temp because the level is already
      # appended to each kernel name. Don't have the same problem for
      # out_nms, as that is just vector with the exact levels

      for(i in seq_along(levs$drop_levels)) {

        temp    <- temp[!grepl(levs$drop_levels[i], temp)]
      }

      out_nms <- out_nms[!out_nms %in% levs$drop_levels]

    }

    mega_mat <- as.list(temp) %>%
      lapply(function(x) {
        y <- rlang::parse_expr(x)
        rlang::enquo(y)
      })


    out_nms <- paste("mega_matrix", out_nms, sep = "_")

  } else {

    mega_mat <- list(mega_mat)
    out_nms  <- "mega_matrix"

  }

  out         <- lapply(mega_mat,
                        function(x, ipm) {
                          .make_mega_mat(ipm, x)
                        },
                        ipm = ipm)

  names(out) <- out_nms

  return(out)

}

#' @rdname make_iter_kernel
#' @export

format_mega_kernel.age_x_size_ipm <- function(ipm,
                                              name_ps,
                                              f_forms,
                                              ...) {

  kern_ids <- ipm$proto_ipm$kernel_id

  # Remove the iteration procedure from the id list. We shouldn't need that
  # to generate a block matrix from the sub kernels anyway.

  fams <- vapply(ipm$proto_ipm$params,
                 function(x) x$family,
                 character(1L))

  if(any(fams == "IPM")) {
    k_ind <- which(fams == "IPM")
    kern_ids <- kern_ids[-c(k_ind)]
  }

  clean_ids <- vapply(kern_ids,
                      function(x) strsplit(x, "_")[[1]][1],
                      character(1L))

  clean_p   <- strsplit(name_ps, "_")[[1]][1]
  p_ind     <- which(clean_ids == clean_p)

  # Get all ages. If there's  max age, we'll create an object with that
  # info called add_p. if there isn't then add_p is just a column of 0s

  ages     <- ipm$proto_ipm$age_indices %>%
    .flatten_to_depth(1L) %>%
    .[!duplicated(names(.))]

  all_ages <- unlist(ages)

  tot_len  <- length(all_ages)

  add_p <- matrix("0", nrow = (tot_len - 1), ncol = 1)

  Ps <- matrix(data = NA_character_,
                nrow = (tot_len - 1),
                ncol = (tot_len - 1))

  Fs <- matrix(NA_character_,
               nrow = 1,
               ncol = tot_len)

  # Next, get the base kernel name. We need this because we only want to modify
  # age here, not any additional par_sets (if present). The ps are a bit easier
  # because there really only should be 1 of those. However, the Fs may have a
  # some additional terms (e.g. F + C), so we have to get both of those, then
  # iterate over each one using the base expression and substituting as needed.

  base_p  <- kern_ids[p_ind]

  f_terms <- .args_from_txt(f_forms)

  clean_fs <- vapply(f_terms,
                     function(x) strsplit(x, "_")[[1]][1],
                     character(1L))

  base_fs <- kern_ids[clean_ids %in% clean_fs]

  # Now, modify the f_forms expression to use the actual kernel names
  # as they appear in the kernel_ids column

  for(i in seq_along(base_fs)) {
    f_forms <- gsub(f_terms[i], base_fs[i], f_forms)
  }

  diag_ps <- character(tot_len - 1)

  for(i in seq_along(all_ages)) {

    diag_ps[i] <- gsub("age", all_ages[i], base_p)
    Fs[ , i]   <- gsub("age", all_ages[i], f_forms)

  }

  # now, handle the max_age case.
  if("max_age" %in% names(ages)) {

    max_p <- diag_ps[length(diag_ps)]
    add_p[(tot_len - 1), ] <- max_p
    diag_ps <- diag_ps[-c(tot_len)]
  }

  diag(Ps) <- diag_ps

  mega_mat <- rbind(
    Fs,
    cbind(Ps, add_p)
  )

  # Insert 0s for impossible transitions. If only we could grow younger... ::sigh::
  mega_mat[is.na(mega_mat)] <- "0"

  # finally, convert back to vector that format_mega.default expects

  mega_mat <- as.vector(t(mega_mat))

  out <- format_mega_kernel.default(ipm, mega_mat = mega_mat)

  return(out)
}

# Accessors functions ----------------------

#' @rdname accessors
#' @export

domains <- function(object) {
  UseMethod("domains")
}

#' @title Accessor functions for (proto_)ipm objects
#' @rdname accessors
#'
#' @description Functions that access slots of a \code{*_ipm} (including
#'  \code{proto_ipm}). \code{default} methods correspond to \code{*_ipm} objects.
#'
#' @param object A \code{proto_ipm} or object created by \code{make_ipm()}.
#'
#' @return Depending on the class of \code{object}, a list
#' with types numeric or character.
#'
#' @details The \code{*.default} method corresponds to output from \code{make_ipm()},
#' and the \code{*.proto_ipm} methods correspond to outputs from \code{define_*}.
#'
#' When using \code{kernel_formulae<-} and \code{vital_rates_exprs<-}, the right
#' hand side of the expression must be wrapped in \code{new_fun_form}. See
#' examples.
#'
#' Note that when using \code{vital_rate_funs}, unless the vital rate expression
#' explicitly contains an expression for integration, these functions
#' \strong{are not yet integrated!} This is useful for things like sensitivity
#' and elasticity analysis, but care must be taken to not use these values
#' incorrectly.
#'
#' @examples
#'
#' data(gen_di_det_ex)
#'
#' proto <- gen_di_det_ex$proto_ipm
#'
#' # Create a new, iterated IPM
#' new_ipm <- make_ipm(proto, iterate = TRUE,
#'                     iterations = 100, return_all_envs = TRUE)
#'
#' vital_rate_exprs(new_ipm)
#' kernel_formulae(new_ipm)
#' vital_rate_funs(new_ipm)
#'
#' domains(new_ipm)
#' parameters(new_ipm)
#'
#' # Usage is the same for proto_ipm's as *_ipm's
#'
#' vital_rate_exprs(proto)
#' kernel_formulae(proto)
#'
#' domains(proto)
#' parameters(proto)
#'
#' int_mesh(new_ipm)
#'
#' # Setting new parameters, vital rate expressions, and kernel formulae
#' # only works on proto_ipm's.
#'
#' # This replaces the "g_int" parameter and leaves the rest untouched
#'
#' parameters(proto) <- list(g_int = 1.5)
#'
#' # This creates a new g_z parameter and leaves the rest of parameters untouched
#' parameters(proto) <- list(g_z = 2.2)
#'
#' # setting a new vital rate or kernel expression requires wrapping the
#' # right-hand side in a call to new_fun_form(). new_fun_form uses expressions
#' # with the same format as ... in define_kernel()
#'
#' vital_rate_exprs(proto,
#'                  kernel = "P",
#'                  vital_rate = "g_mu") <- new_fun_form(g_int + g_z + g_slope * ht_1)
#'
#' kernel_formulae(proto, kernel = "stay_discrete") <- new_fun_form(g_z * d_ht)
#'
#' @export

domains.proto_ipm <- function(object) {

  out <- lapply(object$domain, function(x) x)

  out <- lapply(out,
                function(x) {

                  temp <- lapply(x, function(y) {
                    if(!all(is.na(y))) {
                      names(y) <- c("lower_bound",
                                    "upper_bound",
                                    "n_meshpoints")

                    }

                    return(y)
                  }
                  )

                  return(temp)
                }
  ) %>%
    .flatten_to_depth(1L) %>%
    lapply(function(x) {
      if(any(is.na(x))) {
        return(NULL)
      } else {
        return(x)
      }
    }) %>%
    Filter(f = Negate(is.null), x = .) %>%
    .[!duplicated(names(.)) & !is.na(names(.))]

  class(out) <- c("ipmr_domains", "list")
  attr(out, "proto") <- object

  return(out)

}

#' @rdname accessors
#' @export

domains.default <- function(object) {

  domains(object$proto_ipm)

}

#' @rdname accessors
#' @export

vital_rate_exprs <- function(object) {

  UseMethod("vital_rate_exprs")

}

#' @rdname accessors
#' @importFrom stats setNames
#' @export

vital_rate_exprs.proto_ipm <- function(object) {

  out <- lapply(object$params, function(x) x$vr_text) %>%
    stats::setNames(c("")) %>%
    lapply(function(x)
           if(any(is.na(x) | is.null(x) | rlang::is_empty(x))) {
             return(NULL)
           }  else {
             return(x)
           }
    ) %>%
    Filter(Negate(is.null), x = .) %>%
    Filter(Negate(rlang::is_empty), x = .) %>%
    .flatten_to_depth(1L) %>%
    lapply(rlang::parse_expr)

  out <- out[!duplicated(names(out))]

  class(out) <- c("ipmr_vital_rate_exprs", "list")

  attr(out, "proto") <- object

  return(out)
}

#' @rdname accessors
#' @export

vital_rate_exprs.default <- function(object) {

  vital_rate_exprs(object$proto_ipm)

}

#' @rdname accessors
#' @export

vital_rate_funs <- function(ipm) {

  proto    <- .initialize_kernels(ipm$proto_ipm, TRUE, "right")$others

  env_list <- ipm$env_list


  if(length(env_list) < 2) {

    stop("Cannot find sub-kernel evaluation environments in 'ipm'.",
         " Re-run 'make_ipm()' with 'return_all_envs = TRUE'.")

  }

  out <- switch(as.character(grepl("stoch_param|_dd_", class(ipm)[1])),
                "TRUE"  = .vr_funs_stoch_param(proto, env_list, ipm),
                "FALSE" = .vr_funs_det_kerns(proto, env_list, ipm))


  out <- lapply(out,
                function(x) {
                  class(x) <- c("ipmr_vital_rate_funs",
                                class(x))

                  return(x)
                })

  return(out)

}

.vr_funs_stoch_param <- function(proto, env_list, ipm) {

  kern_nms <- proto$kernel_id

  if(any(proto$uses_par_sets)) {

    kern_seq <- ipm$env_seq

    if(inherits(kern_seq, "data.frame")) {
      kern_seq <- as.character(kern_seq$kernel_seq)
    }

  } else {

    kern_seq <- NULL

  }

  n_its <- ncol(ipm$pop_state[[1]]) - 1

  out <- list()

  for(i in seq_len(n_its)) {

    use_it   <- kern_seq[i]
    kern_ind <- paste(kern_seq[i], "it", i, sep = "_")
    kern_ind <- paste("*_", kern_ind, "$", sep = "")

    if(all(is.null(kern_seq))) {
      kern_ind <- gsub("__","_", kern_ind)
    }

    use_nm    <- names(ipm$sub_kernels)[grepl(kern_ind, names(ipm$sub_kernels))]
    use_env   <- env_list[use_nm]
    kern_ind  <- gsub(paste("_it_", i, sep =""), "", use_nm)
    pro_ind   <- which(proto$kernel_id %in% kern_ind)

    kern_vrs  <- lapply(proto$params[pro_ind],
                        function(x) names(x$vr_text))

    out[[i]] <- lapply(seq_along(pro_ind),
                       function(x, use_env, kern_vrs, proto, env_list, pro_ind)
                         .get_vr_fun(use_env[[x]],
                                     kern_vrs[[x]], proto[pro_ind[x], ],
                                     env_list),
                       use_env = use_env,
                       kern_vrs = kern_vrs,
                       proto = proto,
                       env_list = env_list,
                       pro_ind = pro_ind)

    names(out[[i]]) <- use_nm

  }

  out <- .flatten_to_depth(out, 2L)

  return(out)

}


.vr_funs_det_kerns <- function(proto, env_list, ipm) {

  out <- list()

  for(i in seq_len(nrow(proto))) {

    kern_nm <- proto$kernel_id[i]

    kern_vrs <- names(proto$params[[i]]$vr_text)

    if(rlang::is_empty(kern_vrs)) {

      out[[i]] <- "No vital rate functions specified"
      names(out)[i] <- kern_nm
      next

    }

    use_env <- env_list[[kern_nm]]

    out[[i]]      <- .get_vr_fun(use_env, kern_vrs, proto[i, ], env_list)

    names(out)[i] <- kern_nm

  }

  return(out)
}

.get_vr_fun <- function(use_env, kern_vrs, proto, env_list) {

  kern_nm <- proto$kernel_id

  out <- rlang::env_get_list(env = use_env,
                             nms = kern_vrs,
                             inherit = FALSE)

  out <- lapply(out,
                function(x, kern_cls, start, end, main_env, nm) {

                  class(x) <- c(kern_cls, class(x))

                  .fun_to_iteration_mat(x, start, end, main_env, nm)

                },
                kern_cls = proto$params[[1]]$family,
                start    = names(proto$domain[[1]])[1],
                end      = names(proto$domain[[1]])[2],
                main_env = env_list$main_env,
                nm       = kern_nm
  )

  return(out)

}

#' @rdname accessors
#'
#' @param kernel The name of the kernel to insert the new vital rate expression
#' into.
#' @param vital_rate The name of the vital rate to replace. If the vital rate
#' doesn't already exist in the \code{object}, a new one with this name will be
#' created.
#'
#' @export

`vital_rate_exprs<-` <- function(object, kernel, vital_rate, value) {

  UseMethod("vital_rate_exprs<-")

}

#' @rdname accessors
#' @export

`vital_rate_exprs<-.proto_ipm` <- function(object, kernel, vital_rate, value) {


  object$params[[kernel]]$vr_text[[vital_rate]] <- rlang::quo_text(value)

  return(object)

}

#' @rdname accessors
#'
#' @param form An expression representing the new vital rate or kernel formula
#' to insert.
#'
#' @export

new_fun_form <- function(form) {

  rlang::enquo(form)

}

#' @rdname accessors
#' @export

kernel_formulae <- function(object) {

  UseMethod("kernel_formulae")

}

#' @rdname accessors
#' @export

kernel_formulae.proto_ipm <- function(object) {

  out <- lapply(object$params, function(x) x$formula) %>%
    .flatten_to_depth(1L) %>%
    Filter(f = Negate(is.na), x = .) %>%
    lapply(rlang::parse_expr)

  class(out) <- c("ipmr_kernel_exprs", "list")
  attr(out, "proto") <- object

  return(out)
}

#' @rdname accessors
#' @export

kernel_formulae.default <- function(object) {

  kernel_formulae(object$proto_ipm)

}

#' @rdname accessors
#' @export

`kernel_formulae<-` <- function(object, kernel, value) {

  UseMethod("kernel_formulae<-")

}

#' @rdname accessors
#' @export

`kernel_formulae<-.proto_ipm` <- function(object, kernel, value) {

  object$params[[kernel]]$formula <- rlang::quo_text(value)

  return(object)

}

#' @export
#' @rdname accessors

parameters <- function(object) {

  UseMethod("parameters")
}

#' @rdname accessors
#' @export

parameters.proto_ipm <- function(object) {

  out <- lapply(object$params, function(x) x$params) %>%
    stats::setNames(c("")) %>%
    purrr::flatten() %>%
    .[!duplicated(names(.))]

  class(out) <- c("ipmr_parameters", "numeric")
  attr(out, "proto") <- object

  return(out)

}

#' @rdname accessors
#' @export

parameters.default <- function(object) {

  parameters(object$proto_ipm)

}



#' @rdname accessors
#' @export

`parameters<-` <- function(object, value) {

  UseMethod("parameters<-")

}

#' @rdname accessors
#' @param value For \code{parameters<-}, a named list of new parameters. The new list does not need
#' to contain all of the parameters, just the ones to update/append. For
#' \code{vital_rate_exprs<-} and \code{kernel_formulae<-}, a new functional form.
#' The new functional form must be wrapped in a call to \code{new_fun_form}.
#' @export

`parameters<-.proto_ipm` <- function(object, value) {

  value <- lapply(value, .protect_model)

  for(i in seq_len(nrow(object))) {

    old_pars <- object$params[[i]]$params

    old_par_ind <- names(old_pars) %in% names(value)

    if(all(old_par_ind)) {

      object$params[[i]]$params <- value

    } else {

      temp_pars <- c(old_pars[!old_par_ind], value)

      temp_pars <- temp_pars[!duplicated(names(temp_pars))]

      object$params[[i]]$params <- temp_pars

    }
  }

  return(object)

}

#' @rdname accessors
#' @param ipm An object created by \code{make_ipm()}. This argument only applies to
#' \code{int_mesh()} and \code{vital_rate_funs()} (because these are not built
#' until \code{make_ipm()} is called).
#' @export

int_mesh <- function(ipm) {

  if(!"main_env" %in% names(ipm$env_list)) {
    stop("Cannot find the 'main_env' object in the IPM. Do you need to re-run\n",
         "make_ipm() with 'return_main_env = TRUE'?",
         call. = FALSE)
  }

  states <- lapply(ipm$proto_ipm$domain,
                   function(x) {
                     all_nms <- names(x)
                     temp    <- all_nms[!grepl("not_applicable", all_nms)]
                     out     <- lapply(temp,
                                       function(y) {
                                         gsub("_[0-9+]", "", y)
                                       }) %>%
                       unlist() %>%
                       unique()

                     return(out)

                   }
  ) %>%
    Filter(f = Negate(is.null), x = .) %>%
    unlist() %>%
    unique()

  main_env <- ipm$env_list$main_env

  all_nms <- c(paste("d_", states, sep = ""),
               vapply(states,
                      function(x) paste(x, c("_1", "_2"), sep = ""),
                      character(2L)))

  out <- rlang::env_get_list(main_env,
                             all_nms,
                             default = NULL,
                             inherit = FALSE)

  class(out) <- c("ipmr_mp_mesh", "list")

  return(out)

}

#' @rdname accessors
#' @export

pop_state <- function(object) {

  UseMethod("pop_state")

}

#' @rdname accessors
#' @export

pop_state.proto_ipm <- function(object) {

  ps <- object$pop_state %>%
    .flatten_to_depth(1L)

  if(!rlang::is_named(ps)) {
    out <- list("No population state defined.")
    class(out) <- c("ipmr_pop_state", "list")
    return(out)
  }

  doms <- domains(object) %>%
    lapply(function(x) x[3])

  out <- lapply(ps, .pretty_pop_state)

  out <- out[!duplicated(names(out))]

  class(out) <- c("ipmr_pop_state", "list")
  attr(out, "proto") <- object

  return(out)

}

#' @noRd

.pretty_pop_state <- function(pop_state) {

  if(rlang::is_quosure(pop_state)) {
    out <- rlang::quo_squash(pop_state)

  } else if(rlang::is_atomic(pop_state)) {

    out <- "Pre-defined population state."

  } else if(is.na(pop_state)) {

    out <- "Population state not defined."

  }

  return(out)

}

#' @rdname accessors
#' @export

pop_state.default <- function(object) {

  ps <- object$pop_state

  out <- ps[!grepl('lambda', names(ps))]

  return(out)

}

#' @title Extract threshold based population size information
#'
#' @description Given a model object, this function computes population sizes
#' given thresholds for a state variable of interest. For example,
#' the number (or proportion) of individuals shorter than 60 cm tall at the 20th
#' time step of the model.
#'
#' @param ipm An object created by \code{make_ipm}
#' @param time_step the time step to pull out. Can be a single time step or a
#' vector of multiple time steps. In the latter case, one value is computed for
#' each time step.
#' @param ... Named expressions that provide the threshold information for the
#' desired classes. The expression should be logicals with a state variable name
#' on the left side, and a threshold value on the right side.
#'
#' @return A named list of numeric vectors containing the summed population
#' sizes at each requested time step. Names are taken from \code{...}.
#'
#' @examples
#' data(gen_di_det_ex)
#'
#' # Rebuild the model and return_main_env this time
#'
#' gen_di_det_ex <- gen_di_det_ex$proto_ipm %>%
#'   make_ipm(iterate = TRUE, iterations = 50, return_main_env = TRUE)
#'
#' disc_sizes <- collapse_pop_state(gen_di_det_ex,
#'                                  time_step = 20:25,
#'                                  seedlings = ht <= 10,
#'                                  NRA = ht > 10 & ht <= 200,
#'                                  RA = ht > 200)
#'
#' @export

collapse_pop_state <- function(ipm, time_step, ...) {

  thresh <- rlang::enquos(...)

  args <- lapply(thresh, function(x) {
    rlang::quo_text(x) %>%
      .args_from_txt(.)
  }) %>%
    lapply(function(x) unique(Filter(f = .can_be_number, x = x)))

  pos_states <- unlist(ipm$proto_ipm$state_var) %>%
    unique()

  for(i in seq_along(thresh)) {

    if(any(args[[i]] %in% pos_states)) {

      use_state <- args[[i]][which(args[[i]] %in% pos_states)]
      use_state <- paste("n_", use_state, sep = "")

      attr(thresh[[i]], "state_var") <- use_state

    } else {

      stop("could not infer state variable from expression: ",
           rlang::quo_text(thresh[[i]]), ".\n",
           "Is this state variable spelled correctly?")
    }

  }

  mesh   <- int_mesh(ipm) %>%
    lapply(unique) %>%
    setNames(gsub("_[0-9]", "", names(.))) %>%
    .[!duplicated(names(.))]

  ev_env <- rlang::env(!!! mesh)
  thresh <- lapply(thresh,
                   function(x, set_env) rlang::quo_set_env(x, set_env),
                   set_env = ev_env)

  inds <- lapply(thresh,
                 function(x) {
                   temp <- rlang::eval_tidy(x)
                   which(temp)
                 })

  pops <- ipm$pop_state[!grepl("lambda", names(ipm$pop_state))]
  out  <- list()

  for(i in seq_along(thresh)) {

    use_state <- attr(thresh[[i]], "state_var")

    use_ind <- inds[[i]]

    temp    <- matrix(ncol = length(time_step),
                      nrow = length(use_ind))

    for(j in seq_along(time_step)) {
      temp[ , j] <- pops[[use_state]][use_ind , time_step[j]]

    }
    out[[i]] <- apply(temp, 2, sum)
    names(out)[i] <- names(thresh)[i]
  }

  return(out)

}

#' @title Convert ipmr matrix to long data frame
#'
#' @description Converts IPM kernels into long data frames. These are useful for
#' creating plots using \code{ggplot2}.
#'
#' @inheritParams make_iter_kernel
#'
#' @return A data frame with 3 columns named \code{"t"}, \code{"t_1"}, and
#' \code{"value"}.
#'
#' @examples
#'
#' data(gen_di_det_ex)
#'
#' big_mat_df <- ipm_to_df(gen_di_det_ex,
#'                         mega_mat = c(stay_discrete, go_discrete,
#'                                      leave_discrete, P))
#'
#' @export

ipm_to_df <- function(ipm, ...) {

  UseMethod("ipm_to_df")

}

#' @rdname ipm_to_df
#' @export

ipm_to_df.array <- function(ipm, ...) {

  mat_to_df_impl(ipm)

}

#' @rdname ipm_to_df
#' @export

ipm_to_df.default <- function(ipm,
                              ...,
                              mega_mat,
                              name_ps = NULL,
                              f_forms = NULL) {

  mega_mat <- rlang::enquo(mega_mat)

  megas <- make_iter_kernel(ipm,
                            mega_mat = !! mega_mat,
                            name_ps = name_ps,
                            f_forms = f_forms)

  out <- lapply(megas, mat_to_df_impl)

  return(out)

}

#' @title Create iteration kernels from an IPM object
#'
#' @description Creates iteration kernels for IPMs. \code{ipmr} does not create
#' these to iterate models, but they may be useful for further analyses.
#'
#' @param ipm Output from \code{make_ipm}.
#'
#' @param ... Other arguments passed to methods.
#'
#' @param mega_mat A vector with symbols, I's, and/or 0s representing the matrix blocks.
#' They should be specified in ROW MAJOR order! Can also be a character
#' string specifying the call. Parameter set index syntax is supported. When used,
#' \code{format_mega_kernel} will produce as many mega-matrices as there are
#' combinations of \code{par_set_indices} in the \code{proto_ipm}.
#'
#' @param name_ps The prefix(es) for the kernel name that correspond to survival
#' and growth/maturation of existing individuals. For the model
#' \code{K = P_age + F_age}, this would be \code{"P"}. Only applies to
#' age X size models. The \code{"_age"} suffix is appended automatically, so
#' does not need to be supplied.
#'
#' @param f_forms The names of the kernels that correspond to production of new
#' individuals, and possibly, how they are combined. For example, a model that
#' includes sexual (with an "F" kernel) and asexual reproduction (with a "C" kernel),
#' this would be \code{"F + C"}. If data come from multiple sites or years,
#' then this information is supplied using the index syntax (i.e.
#' \code{f_forms = "F_yr + C_yr"}). Only applies to age X size models. The
#' \code{"_age"} index is appended automatically, so does not need to be
#' supplied.
#'
#' @return A list containing a large matrix or many large matrices (when used with
#' suffix syntax). The names in the former case will be \code{"mega_matrix"}
#' and in the latter case, \code{"mega_matrix_<par_sets>"} with the levels of the
#' grouping effects substituted in.
#'
#' @details \code{ipmr} does not generate complete iteration kernels, and uses
#' sub-kernels to iterate models. However, some further analyses are just easier
#' to code with a complete iteration kernel. This handles constructing those for
#' simple and general models of all forms. \code{format_mega_kernel} is used
#' internally by \code{make_iter_kernel} for general IPMs. The difference
#' between these two functions is that \code{make_iter_kernel} always returns
#' a list of objects with \code{c(ipmr_matrix, array, matrix)} classes,
#' whereas \code{format_mega_kernel} always returns a list of objects with
#' \code{c(array, matrix)} classes. The former has \code{plot()} methods while
#' the latter does not.
#'
#' \code{I} and \code{0} represent identity matrices and 0 matrices,
#' respectively. They can be used to fill in blocks that represent either, without
#' having to create those separately and append them to the model object. The function
#' will work out the correct dimensions for both internally, and there is no
#' restriction on the number that may be used in a given call.
#'
#' For \code{age_size_ipm}s, the correct form of \code{mega_mat} is generated
#' internally by creating sub-diagonal matrices for the \code{name_ps} kernels,
#' and a top row using the \code{f_forms}. If parameter set indices are part of the
#' model, the indices should be attached to the \code{name_ps, f_forms} in the
#' function arguments, and the correct block matrices will be generated internally.
#'
#' @examples
#' data(gen_di_det_ex)
#'
#' big_k <- make_iter_kernel(gen_di_det_ex,
#'                             mega_mat = c(0, go_discrete,
#'                                          leave_discrete, P))
#'
#' char_call <- c(0, "go_discrete", "leave_discrete", "P")
#'
#' big_k_c <- make_iter_kernel(gen_di_det_ex, mega_mat = char_call)
#'
#' # Now, with an Identity matrix instead of a 0
#'
#' big_k <- make_iter_kernel(gen_di_det_ex,
#'                             mega_mat = c(I, go_discrete,
#'                                          leave_discrete, P))
#'
#' # For simple IPMs with no grouping effects, this computes the sum of
#' # the sub-kernels (i.e. K = P + F)
#'
#' data(sim_di_det_ex)
#'
#' simple_k <- make_iter_kernel(sim_di_det_ex)
#'
#' @export

make_iter_kernel <- function(ipm,
                             mega_mat = NULL,
                             name_ps = NULL,
                             f_forms = NULL) {

  cls_switch <- as.character(grepl("simple", class(ipm)[1]))

  k_cls <- switch(cls_switch,
                  "TRUE"  = "simple_ipm",
                  "FALSE" = "general_ipm")

  class(ipm) <- c(k_cls, class(ipm))

  mega_mat <- rlang::enquo(mega_mat)

  out <- .make_iter_kernel(ipm,
                           mega_mat = !! mega_mat,
                           name_ps = name_ps,
                           f_forms = f_forms) %>%
    set_ipmr_classes()

  return(out)

}

.make_iter_kernel <- function(ipm, ...) {

  UseMethod(".make_iter_kernel")

}

#' @noRd

.make_iter_kernel.simple_ipm <- function(ipm, ...) {

  proto     <- ipm$proto_ipm

  base_nms  <- proto$kernel_id

  all_kerns <- ipm$sub_kernels

  # Check for grouping effects
  if(any(proto$uses_par_sets)) {

    # First, generate the base_expr for evaluation. In the simple_ipm case, this
    # always just the sum of all sub-kernels. Once we have that, we can sub the
    # par_sets levels, then use args_from_txt to get exact names for evaluating
    # the subbed expression.

    base_expr <- paste(base_nms, collapse = " + ")

    all_effs  <- .flatten_to_depth(proto$par_set_indices, 1L) %>%
      .[!duplicated(names(.))] %>%
      .[!names(.) %in% c("levels", "to_drop")]

    levs      <- .make_par_set_indices(all_effs)

    to_sub    <- paste(names(all_effs), collapse = "_")

    all_args  <- lapply(levs,
                        function(x, base_expr, to_sub) {
                          temp <- gsub(to_sub, x, base_expr)

                          .args_from_txt(temp)

                        },
                        base_expr = base_expr,
                        to_sub = to_sub)

    out <- list()

    for(i in seq_along(all_args)) {

      kern_nms <- all_args[[i]]

      out[[i]] <- do.call(`+`, all_kerns[kern_nms])

      names(out)[i] <- paste("mega_matrix_", levs[i], sep = "")
    }

  } else {

    use_kerns <- ipm$sub_kernels

    out <- list(mega_matrix = do.call(`+`, all_kerns))
  }
  return(out)

}


#' @noRd

.make_iter_kernel.general_ipm <- function(ipm,
                                          ...,
                                          mega_mat,
                                          name_ps = NULL,
                                          f_forms = NULL) {

  mega_mat <- rlang::enquo(mega_mat)

  out <- format_mega_kernel(ipm,
                            mega_mat = !! mega_mat,
                            name_ps  = name_ps,
                            f_forms  = f_forms)

  return(out)

}
