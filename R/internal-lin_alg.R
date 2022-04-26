# Internal linear algebra stuff

# ev helper functions -------------

#' @noRd
# We need the total pop size to standardize our vectors, but we don't really
# know the proper order in which to combine them to create a single output.
# thus, we keep them in a list with entries corresponding to states, and
# use Reduce to compute the total size. The final output will be the list
# of states standardized by this total population size.
# TODO General so that par set models can use this as well.

.extract_conv_ev_general <- function(pop_state, proto) {

  has_ps <- ifelse(any(proto$uses_par_sets), "yes", "no")
  pop_state <- pop_state[!grepl("lambda", names(pop_state))]
  final_it <- dim(pop_state[[1]])[2]

  temp     <- lapply(pop_state,
                     function(x, final_it) {

                       x[ , final_it]

                     },
                     final_it = final_it)

  out <- switch(has_ps,
                "yes" = .standard_pop_ps(temp, proto),
                "no"  = .standard_pop(temp))


  return(out)
}

.standard_pop_ps <- function(pop_state, proto) {

  ps_inds <- .make_par_set_indices(proto$par_set_indices) %>%
    paste0("(_", ., ")$")

  out <- list()

  for(i in seq_along(ps_inds)) {
    use_regex <- ps_inds[i]

    use_pops <- pop_state[grepl(use_regex, names(pop_state))]

    temp <- .standard_pop(use_pops)

    out <- c(out, temp)
  }

  return(out)

}

#' @noRd

.standard_pop <- function(pop_state) {

  pop_std <- Reduce('sum', unlist(pop_state), init = 0)

  out     <- lapply(pop_state,
                    function(x, pop_std) x / pop_std,
                    pop_std = pop_std)

  return(out)
}

#' @noRd
#
.extract_conv_ev_simple <- function(pop_state) {

  pop_state <- pop_state[!grepl("lambda", names(pop_state))]
  final_it <- dim(pop_state[[1]])[2]

  out      <- lapply(pop_state,
                     function(x, final_it) {

                       temp <- x[ , final_it]

                       return(temp / sum(temp))

                     },
                     final_it = final_it)

  return(out)
}


#' @noRd

.get_pop_nm_simple <- function(ipm) {

  pop_nm <- ipm$proto_ipm$state_var %>%
    unlist() %>%
    unique() %>%
    .[1]

  return(pop_nm)
}

#' @noRd

.make_mega_mat <- function(ipm, mega_mat) {

  # Extract names of sub_kernels. We shouldn't need mega-mat after this
  # because call_args() returns the mega_mat call's arguments in order.

  sub_mats    <- rlang::call_args(mega_mat)

  test_ind    <- .args_from_txt(rlang::expr_text(mega_mat))

  sub_mat_nms <- Filter(Negate(function(x) .id_or_0(x)), test_ind)

  sub_kernels <- ipm$sub_kernels

  if(!all(sub_mat_nms %in% names(sub_kernels))) {

    stop("names in 'mega_mat' are not all present in 'ipm$sub_kernels'")

  }

  # Do identity mats first - we need these dimensions to work out
  # the 0s afterwards

  id_test <- vapply(sub_mats,
                    function(x) x == "I",
                    logical(1L))

  if(any(id_test)) {

    lists       <- .fill_Is(sub_mats, sub_kernels)
    sub_mats    <- lists$sub_mats
    sub_kernels <- lists$sub_kernels

  }

  # Next, we need to work out the dimensions of each row + column so that any 0s
  # can be appropriately duplicated. If there aren't any 0s in 'mega_mat', then
  # just skip straight to the creation step.

  zero_test <- vapply(sub_mats,
                      function(x) x == 0,
                      logical(1L))

  if(any(zero_test)) {

    # This actually just creates calls to rep(0, times = dim_1 * dim_2)
    # and makes sure those calls are substituted in for the actual 0s.
    # These get evaluated in .make_mega_mat_impl in the next step.

    sub_mats <- .fill_0s(sub_mats, sub_kernels)

  }

    # *_impl = implementation. This actually generates the full kernel

  out <- .make_mega_mat_impl(sub_mats, sub_kernels)

  return(out)

}

#' @noRd

.make_mega_mat_impl <- function(sub_mats, sub_kernels) {

  # Get number of blocks and then bind

  block_dim <- sqrt(length(sub_mats))

  kern_env  <- rlang::env(!!! sub_kernels)

  sub_mats  <- lapply(sub_mats,
                      function(sub_kern, eval_env) {
                        rlang::eval_tidy(sub_kern, env = eval_env)
                      },
                      eval_env = kern_env)

  # now, we just need to generate an index to subset the sub_mats for
  # generating each row of the block matrix

  blocks <- vector('list', length = block_dim)

  for(i in seq_len(block_dim)) {

    if(i == 1) {
      block_ind <- seq(1, block_dim, by = 1)
    } else {
      block_ind <- block_ind + block_dim
    }

    blocks[[i]] <- do.call(cbind, sub_mats[block_ind])

  }

  out <- do.call(rbind, blocks)

  return(out)
}

#' @noRd

.id_or_0 <- function(x) {

  x == "`0`" || x == "I" || x == "0"

}

#' @noRd

.fill_Is <- function(sub_mats, sub_kernels) {

  mega_block_dim <- sqrt(length(sub_mats))

  dim_mat <- .make_dim_mat(sub_kernels, sub_mats)

  it        <- 1
  id_nm_ind <- 1

  for(ro in seq_len(mega_block_dim)) {

    for(co in seq_len(mega_block_dim)) {

      if(dim_mat[ro, co] == "I") {

        temp <- dim_mat
        temp[ro, co] <- gsub("I", NA_character_, temp[ro, co])

        ro_dim <- vapply(temp[ro, ],
                         function(x) {
                           strsplit(x, ", ")[[1]][1]
                         },
                         character(1L)) %>%
          as.integer()

        if(all(is.na(ro_dim))) {
          ro_dim <- NA_integer_
        } else {
          ro_dim <- max(ro_dim, na.rm = TRUE)
        }

        co_dim <- vapply(temp[ , co],
                         function(x) {
                           strsplit(x, ', ')[[1]][2]
                         },
                         character(1L)) %>%
          as.integer()

        if(all(is.na(co_dim))) {
          co_dim <- NA_integer_
        } else {
          co_dim <- max(co_dim, na.rm = TRUE)
        }

        if(isTRUE(co_dim != ro_dim) || is.na(co_dim != ro_dim)) {

          use_dim <- max(co_dim, ro_dim, na.rm = TRUE)

        } else {

          use_dim <- co_dim

        }

        id_call <- paste("diag(", use_dim, ")", sep = "") %>%
          rlang::parse_expr()

        id_nm <- paste("id_", id_nm_ind, sep = "")

        sub_mats[[it]]      <- rlang::parse_expr(id_nm)
        names(sub_mats)[it] <- id_nm
        id_nm_ind           <- id_nm_ind + 1

        sub_kernels <- c(sub_kernels,
                         rlang::list2(!! id_nm := eval(id_call)))

      } else {

        names(sub_mats)[it] <- as.character(sub_mats[[it]])

      }

      it <- it + 1

    }
  }

  out <- list(
    sub_mats    = sub_mats,
    sub_kernels = sub_kernels
  )

  return(out)
}

#' @noRd

.make_dim_mat <- function(sub_kernels, sub_mats) {

  kern_dims <- lapply(sub_kernels,
                      function(x) {
                        temp <- dim(x)
                        out <- paste(temp, collapse = ', ')
                      })

  dim_env        <- rlang::env()
  mega_block_dim <- sqrt(length(sub_mats))

  rlang::env_bind(!!! kern_dims,
                  .env = dim_env)

  dim_mat        <- lapply(sub_mats,
                           function(x, env_, sub_kernels) {

                             if(as.character(x) %in% names(sub_kernels)) {

                               rlang::eval_tidy(x, env = env_)

                             } else {

                               as.character(x)

                             }
                           },
                           env_ = dim_env,
                           sub_kernels = sub_kernels) %>%
    unlist() %>%
    matrix(ncol = mega_block_dim,
           nrow = mega_block_dim,
           byrow = TRUE)

  return(dim_mat)

}

#' @noRd

.fill_0s <- function(sub_mats, sub_kernels) {

  mega_block_dim <- sqrt(length(sub_mats))

  # Create a character matrix that contains the dimensions of each
  # block in the mega-matrix

  dim_mat <- .make_dim_mat(sub_kernels, sub_mats)

  # Now, loop over and determine the required dimensions for each 0 entry

  it          <- 1
  zero_nm_ind <- 1

  for(ro in seq_len(mega_block_dim)) {

    for(co in seq_len(mega_block_dim)) {

      if(dim_mat[ro, co] == "0") {

        ro_dim <- vapply(dim_mat[ro, ],
                         function(x) {
                           strsplit(x, ', ')[[1]][1]
                         },
                         character(1L)) %>%
          as.integer() %>%
          max(na.rm = TRUE) # Need to remove the NAs generated by entries of 0!

        co_dim <- vapply(dim_mat[ , co],
                         function(x) {
                           strsplit(x, ', ')[[1]][2]
                         },
                         character(1L)) %>%
          as.integer() %>%
          max(na.rm = TRUE)

        # Generate an expression and insert it into sub_mats. These get
        # evaluated and c/rbinded in the next stage
        zero_call <- paste('matrix(rep(0, times = ',
                           ro_dim * co_dim,
                           '), nrow = ',
                           ro_dim,
                           ', ncol = ',
                           co_dim,
                           ')', sep = "") %>%
          rlang::parse_expr()

        sub_mats[[it]]      <- zero_call
        names(sub_mats)[it] <- paste('zero_', zero_nm_ind, sep = "")
        zero_nm_ind         <- zero_nm_ind + 1
      } else {

        names(sub_mats)[it] <- as.character(sub_mats[[it]])

      }
      it <- it + 1

    }   # columns
  }     # end rows


  return(sub_mats)

}

# Lambda helpers----------

#' @noRd

.lambda_pop_size <- function(x, all_lambdas = TRUE) {

  pops    <- x$pop_state
  lam_ind <- grepl("lambda", names(pops))

  if(sum(lam_ind) > 1) {

    lams <- pops[lam_ind]

  } else {

    # Don't drop the list
    lams <- list(lambda = pops$lambda)

  }

  if(!all(is.na(unlist(lams)))) {

    out <- do.call('rbind', lams) %>%
      t()

    if(!all_lambdas) {

      n_its <- dim(out)[1]

      out <- out[n_its, ]
    }

    if(inherits(out, c('matrix', 'array'))) {

      dimnames(out) <- list(NULL, names(lams))

    } else {
      names(out) <- names(lams)
    }
    return(out)

  } else {

    warning("NA's detected in lambda slots - returning NA")

    return(NA_real_)

  }


}

#' @noRd

is_square <- function(x) {

  dim(x)[1] == dim(x)[2]

}

#' @noRd
# Checks for convergence to asymptotic dynamics. Supports either
# lambdas computed by iteration, which

.is_conv_to_asymptotic <- function(x, tol = 1e-6) {

  # Standardize columns in case of dealing w/ population state
  # within right/left_ev. lambdas won't be affected

  x <- apply(x, 2, function(y) y / sum(y))

  end_ind   <- dim(x)[2]
  start_ind <- end_ind - 1

  start_val <- x[ , start_ind]
  end_val   <- x[ , end_ind]



  return(
    isTRUE(
      all.equal(
        start_val, end_val, tolerance = tol
      )
    )
  )

}


#' @rdname check_convergence
#' @title Check for model convergence to asymptotic dynamics
#'
#' @description Checks for convergence to asymptotic dynamics numerically and
#' visually. \code{is_conv_to_asymptotic} checks whether
#' \code{lambda[iterations - 1]} equals \code{lambda[iterations]} within the
#' specified tolerance, \code{tolerance}. \code{conv_plot} plots the time series of
#' \code{lambda} (or \code{log(lambda)}). For stochastic models, a cumulative mean of
#' log(lambda) is used to check for convergence.
#'
#' @param ipm An object returned by \code{make_ipm()}.
#' @param tolerance The tolerance for convergence in lambda or, in the case of stochastic models, the cumulative mean of log(lambda).
#' @param burn_in The proportion of iterations to discard. Default is 0.1
#' (i.e. first 10\% of iterations in the simulation). Ignored for deterministic models.
#'
#' @return \code{is_conv_to_asymptotic}: Either \code{TRUE} or \code{FALSE}.
#' \code{conv_plot}: code{ipm} invisibly.
#' @export
#'

is_conv_to_asymptotic <- function(ipm, tolerance, burn_in) {
  UseMethod("is_conv_to_asymptotic")
}

#' @rdname check_convergence
#' @export

is_conv_to_asymptotic.ipmr_ipm <- function(ipm, tolerance = 1e-6, burn_in = 0.1) {

  pop_state_test <- vapply(ipm$pop_state,
                           function(x) ! any(is.na(x)),
                           logical(1L))

  if(! any(pop_state_test)) {

    stop("pop_state in IPM contains NAs - cannot check for convergence!")
  }

  # If lambda exists, we can just use that and exit early. otherwise, drop
  # lambda entry and proceed with the population state vectors

  lambdas <- ipm$pop_state[grepl("lambda", names(ipm$pop_state))]

  if(!all(is.na(unlist(lambdas)))) {

    #for stochastic models, convert to cummean(log(lambda)) after removing burn_in
    if(grepl("_stoch_", class(ipm)[1])) {
      lambdas <- lapply(lambdas, function(x, burn_in) {

        burn_ind <- seq_len(round(length(x) * burn_in))
        temp <- x[-burn_ind]
        #cummulative mean
        cumsum(log(temp))/1:length(temp)

      },
      burn_in = burn_in)
    }

    convs <- vapply(lambdas, function(x, tol) {

      end <- length(x)
      start <- end - 1

      abs(x[start] - x[end]) <= tol

    },
    logical(1L),
    tol = tolerance)

    return(convs)

  } else {

    warning("Lambda and population state vectors contain NAs. Cannot check for",
            " convergence.")

    return(NA)

  }

}


#' @noRd

# Internal generic to check arguments in lambda() for validity

.check_lambda_args <- function(ipm, type_lambda) {
  UseMethod(".check_lambda_args")
}


#' @noRd

.check_lambda_args.simple_di_det_ipm <- function(ipm, type_lambda) {

  if(!type_lambda %in% c("stochastic", "all", "last")) {
    stop("'type_lambda' must be one of 'all' or 'last'.",
         call. = FALSE)
  }

  if(!attr(ipm, "iterated")) {
    stop("ipmr cannot compute lambda for a model that is not yet",
         " iterated.\n",
         "Re-run with make_ipm(iterate = TRUE).",
         call. = FALSE)
  }

  if(type_lambda == 'stochastic') {
    stop("Cannot compute stochastic lambda for a deterministic IPM.", call. = FALSE)
  }

  invisible(TRUE)

}

#' @noRd

.check_lambda_args.simple_di_stoch_kern_ipm <- function(ipm, type_lambda) {

  if(!attr(ipm, 'iterated')) {
    stop("ipmr cannot compute lambda for a model that is not iterated!",
         call. = FALSE)
  }

  if(!type_lambda %in% c("stochastic", "all", "last")) {
    stop("'type_lambda' must be one of 'stochastic', 'all', or 'last'.",
         call. = FALSE)
  }

  invisible(TRUE)
}

#' @noRd

.check_lambda_args.simple_di_stoch_param_ipm <- function(ipm, type_lambda) {

  if(!attr(ipm, 'iterated')) {
    stop("ipmr cannot compute lambda for a model that is not iterated!",
         call. = FALSE)
  }

  if(!type_lambda %in% c("stochastic", "all", "last")) {
    stop("'type_lambda' must be one of 'stochastic', 'all', or 'last'.",
         call. = FALSE)
  }

  invisible(TRUE)

}

#' @noRd

.check_lambda_args.general_di_det_ipm <- function(ipm, type_lambda) {

  if(!type_lambda %in% c("stochastic", "all", "last")) {

    stop("'type_lambda' must be one of 'all' or 'last'.",
         call. = FALSE)

  } else if(type_lambda == 'stochastic') {

    stop("ipmr cannot compute stochastic lambda for a deterministic IPM.",
         call. = FALSE)
  }

  if(!attr(ipm, "iterated")) {

    stop("Cannot compute lambda for a model that is not yet",
         " iterated.\n",
         "Re-run with make_ipm(iterate = TRUE).",
         call. = FALSE)
  }

  invisible(TRUE)

}

#' @noRd

.check_lambda_args.general_di_stoch_kern_ipm <- function(ipm, type_lambda) {

    if(!type_lambda %in% c("stochastic", "all", "last")) {

    stop("'type_lambda' must be one of 'stochastic', 'all', or 'last'.",
         call. = FALSE)

  }

  if(!attr(ipm, "iterated")) {

    stop("Cannot compute lambda by population size for a model that is not yet",
         " iterated.\n",
         "Re-run with make_ipm(iterate = TRUE).",
         call. = FALSE)
  }

  invisible(TRUE)

}

#' @noRd

.check_lambda_args.general_di_stoch_param_ipm <- function(ipm, type_lambda) {

  if(!type_lambda %in% c("stochastic", "all", "last")) {

    stop("'type_lambda' must be one of 'stochastic', 'all', or 'last'.",
         call. = FALSE)

  }

  if(!attr(ipm, "iterated")) {

    stop("Cannot compute lambda for a model that is not yet",
         " iterated.\n",
         "Re-run with make_ipm(iterate = TRUE).",
         call. = FALSE)
  }

  invisible(TRUE)

}

# check_lambda_args.dd-----------

#' @noRd

.check_lambda_args.simple_dd_det_ipm <- function(ipm, type_lambda) {

  if(!type_lambda %in% c("stochastic", "all", "last")) {
    stop("'type_lambda' must be one of 'all' or 'last'.",
         call. = FALSE)
  }

  if(!attr(ipm, "iterated")) {
    stop("ipmr cannot compute lambda for a model that is not yet",
         " iterated.\n",
         "Re-run with make_ipm(iterate = TRUE).",
         call. = FALSE)
  }

  if(type_lambda == 'stochastic') {
    stop("Cannot compute stochastic lambda for a deterministic IPM.", call. = FALSE)
  }

  invisible(TRUE)

}

#' @noRd

.check_lambda_args.simple_dd_stoch_kern_ipm <- function(ipm, type_lambda) {

  if(!attr(ipm, 'iterated')) {
    stop("ipmr cannot compute lambda for a model that is not iterated!",
         call. = FALSE)
  }

  if(!type_lambda %in% c("stochastic", "all", "last")) {
    stop("'type_lambda' must be one of 'stochastic', 'all', or 'last'.",
         call. = FALSE)
  }

  invisible(TRUE)
}

#' @noRd

.check_lambda_args.simple_dd_stoch_param_ipm <- function(ipm, type_lambda) {

  if(!attr(ipm, 'iterated')) {
    stop("ipmr cannot compute lambda for a model that is not iterated!",
         call. = FALSE)
  }

  if(!type_lambda %in% c("stochastic", "all", "last")) {
    stop("'type_lambda' must be one of 'stochastic', 'all', or 'last'.",
         call. = FALSE)
  }

  invisible(TRUE)

}

#' @noRd

.check_lambda_args.general_dd_det_ipm <- function(ipm, type_lambda) {

  if(!type_lambda %in% c("stochastic", "all", "last")) {

    stop("'type_lambda' must be one of 'all' or 'last'.",
         call. = FALSE)

  } else if(type_lambda == 'stochastic') {

    stop("ipmr cannot compute stochastic lambda for a deterministic IPM.",
         call. = FALSE)
  }

  if(!attr(ipm, "iterated")) {

    stop("Cannot compute lambda for a model that is not yet",
         " iterated.\n",
         "Re-run with make_ipm(iterate = TRUE).",
         call. = FALSE)
  }

  invisible(TRUE)

}

#' @noRd

.check_lambda_args.general_dd_stoch_kern_ipm <- function(ipm, type_lambda) {

  if(!type_lambda %in% c("stochastic", "all", "last")) {

    stop("'type_lambda' must be one of 'stochastic', 'all', or 'last'.",
         call. = FALSE)

  }

  if(!attr(ipm, "iterated")) {

    stop("Cannot compute lambda by population size for a model that is not yet",
         " iterated.\n",
         "Re-run with make_ipm(iterate = TRUE).",
         call. = FALSE)
  }

  invisible(TRUE)

}

#' @noRd

.check_lambda_args.general_dd_stoch_param_ipm <- function(ipm, type_lambda) {

  if(!type_lambda %in% c("stochastic", "all", "last")) {

    stop("'type_lambda' must be one of 'stochastic', 'all', or 'last'.",
         call. = FALSE)

  }

  if(!attr(ipm, "iterated")) {

    stop("Cannot compute lambda for a model that is not yet",
         " iterated.\n",
         "Re-run with make_ipm(iterate = TRUE).",
         call. = FALSE)
  }

  invisible(TRUE)

}



# Helper function to handle when burn_in == 0
#' @noRd

.thin_stoch_lambda <- function(lambdas, burn_ind, log) {

  lambdas <- log(lambdas)

  if(length(burn_ind > 0)) {
    out <- apply(lambdas, 2,
                 function(x, burn_ind) mean(x[-c(burn_ind)]),
                 burn_ind = burn_ind)
  } else {
    out <- apply(lambdas, 2, mean)
  }

  if(!log) out <- exp(out)

  names(out) <- dimnames(lambdas)[[2]]

  return(out)
}

