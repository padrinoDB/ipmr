#' @noRd

.iterate_model <- function(proto_ipm, ...) {

  UseMethod('.iterate_model')
}

#' @noRd
#' @importFrom purrr map2

.iterate_model.simple_di_det <- function(proto_ipm,
                                         sub_kern_list,
                                         iterations,
                                         pop_state,
                                         main_env,
                                         k_row,
                                         normal) {

  # For par_setarchical deterministic models, we need to create a determinstic
  # simulation for each level of the model

  if(any(proto_ipm$uses_par_sets)) {

    par_sets <- proto_ipm$par_set_indices[proto_ipm$uses_par_sets]
    par_sets <- par_sets[!duplicated(par_sets)]

    levs       <- .make_par_set_indices(par_sets)
    lambda_nms <- paste("lambda", levs, sep = "_")

    pop_state$lambda <- NULL

    lambdas        <- vector('list', length = length(levs))
    names(lambdas) <- lambda_nms

    lambdas <- lapply(lambdas,
                      function(x, iterations) {
                        x <- array(NA_real_, dim = c(1, (iterations + 1)))
                      },
                      iterations = iterations)

    pop_state <- c(pop_state, lambdas)


    for(i in seq_along(levs)) {

      use_lev   <- paste("(_", levs[i], ")$", sep = "")

      use_kerns <- sub_kern_list[grepl(use_lev, names(sub_kern_list))]

      if(!all(proto_ipm$uses_par_sets)) {

        append    <- proto_ipm$kernel_id[!proto_ipm$uses_par_sets]
        use_kerns <- c(use_kerns, sub_kern_list[append])

      }

      use_pop <- pop_state[grepl(use_lev, names(pop_state))]
      use_k   <- k_row[grepl(use_lev, k_row$kernel_id), ]

      for(j in seq_len(iterations)) {

        # Bind a helper for the model iteration number. Users can use this
        # to specify lagged effects.

        .bind_iter_var(main_env, j)

        pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                          proto_ipm     = proto_ipm,
                                          sub_kern_list = use_kerns,
                                          main_env    = main_env)

        names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

        pop_list_t_1        <- pop_list_t_1[names(use_pop)]

        # Drop the NA entry created by the lambda slot.

        pop_list_t_1        <- pop_list_t_1[!is.na(names(pop_list_t_1))]

        # Next, compute population sizes and time step growth rates.
        # this is required because if normalize_pop_size = TRUE,
        # there isn't any way to compute growth rates at each
        # time step after the fact

        pop_size_t_1 <- vapply(
          pop_list_t_1,
          function(x) sum(x),
          numeric(1L)
        ) %>%
          Reduce(f = '+', x = ., init = 0)

        pop_size_t <- .pop_size(use_pop, j)

        if(normal) {

          pop_list_t_1 <- lapply(pop_list_t_1,
                                 function(x, pop_size) {
                                   x / pop_size
                                 },
                                 pop_size = pop_size_t_1)

        }

        pop_list_t_1$lambda <- pop_size_t_1 / pop_size_t
        use_l_nm            <- names(use_pop)[grepl("lambda", names(use_pop))]

        names(pop_list_t_1)[names(pop_list_t_1) == 'lambda'] <- use_l_nm

        pop_list_t_1 <- pop_list_t_1[names(use_pop)]


        use_pop <- update_pop_state(use_pop, pop_list_t_1, iter = j)


        # Now, update the names so that we can bind the new n_*_t to
        # main_env so that the next iteration is evaluated correctly

        names(pop_list_t_1) <- names(pop_list_t_1) %>%
          paste(., '_t', sep = "")


        rlang::env_bind(.env = main_env,
                        !!! pop_list_t_1)

      } # end single level simulation

      pop_state[grepl(use_lev, names(pop_state))] <- use_pop

    } # end par_setarchical levels loop


  } else {


    for(i in seq_len(iterations)) {

      # Simple model, no par_set stuff going on. We can just
      # stick the kernels and their expressions right into eval_general_det

      .bind_iter_var(main_env, i)

      pop_list_t_1 <- .eval_general_det(k_row         = k_row,
                                        proto_ipm     = proto_ipm,
                                        sub_kern_list = sub_kern_list,
                                        main_env      = main_env)

      # make names match pop_state and then reorder the list for easy
      # insertion
      names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

      pop_list_t_1        <- pop_list_t_1[names(pop_state)]

      # Drop the NA entry created by the lambda slot.

      pop_list_t_1        <- pop_list_t_1[!is.na(names(pop_list_t_1))]

      # Next, compute population sizes and time step growth rates.
      # this is required because if normalize_pop_size = TRUE,
      # there isn't any way to compute growth rates at each
      # time step after the fact

      pop_size_t_1 <- vapply(
        pop_list_t_1,
        function(x) sum(x),
        numeric(1L)
      ) %>%
        Reduce(f = '+', x = ., init = 0)

      pop_size_t <- .pop_size(pop_state, i)

      if(normal) {

        pop_list_t_1 <- lapply(pop_list_t_1,
                               function(x, pop_size) {
                                 x / pop_size
                               },
                               pop_size = pop_size_t_1)

      }

      pop_list_t_1$lambda <- pop_size_t_1 / pop_size_t

      pop_list_t_1        <- pop_list_t_1[names(pop_state)]

      pop_state <- update_pop_state(pop_state, pop_list_t_1, iter = i)

      # Now, update the names so that we can bind the new n_*_t to
      # main_env so that the next iteration is evaluated correctly

      names(pop_list_t_1) <- names(pop_list_t_1) %>%
        paste(., '_t', sep = "")


      rlang::env_bind(.env = main_env,
                      !!! pop_list_t_1)

    } # End no par_set stuff

  } # End if(par_set){} else {}

  # Remove the NA at 1 - it is a dummy so that purrr::map2 can work
  # properly

  pop_state[grepl('lambda', names(pop_state))] <- lapply(
    pop_state[grepl("lambda", names(pop_state))],
    function(x) {
      out <- x[ , -1, drop = FALSE]
      return(out)
    }
  )

  return(pop_state)

}


#' @noRd

.iterate_model.simple_di_stoch_kern <- function(proto_ipm,
                                                sub_kern_list,
                                                iterations,
                                                kern_seq,
                                                pop_state,
                                                main_env,
                                                k_row,
                                                normal,
                                                report_progress) {

  if(rlang::is_quosure(pop_state)) {
    pop_state <- rlang::eval_tidy(pop_state)
  }


  for(i in seq_len(iterations)) {

    # Select kernels from kern_seq. Similarly, we need to subset the k_rows
    # object so that .eval_general_det doesn't loop over ones which include
    # expressions for other levels of 'kern_seq'

    if(!is.null(kern_seq)) {

      kern_ind <- grepl(kern_seq[i], names(sub_kern_list))

      k_row_ind <- which(grepl(kern_seq[i], k_row$kernel_id))


      use_kerns <- sub_kern_list[kern_ind]
      use_k     <- k_row[k_row_ind, ]

      # Deals with the case where only a subset of kernels have suffixes.
      # In that case, we need to include the ones that don't have them
      # in the use_kerns list every single time!

      if(any(!proto_ipm$uses_par_sets) && any(proto_ipm$uses_par_sets)) {

        nm_ind <- proto_ipm$kernel_id[!proto_ipm$uses_par_sets]
        to_add <- sub_kern_list[nm_ind]

        use_kerns <- c(use_kerns, to_add)

      }


    }

    .bind_iter_var(main_env, i)

    .stoch_progress_message(report_progress, iterations, i)


    # Need to return an iterated pop_state object

    pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                      proto_ipm     = proto_ipm,
                                      sub_kern_list = use_kerns,
                                      main_env    = main_env)

    # make names match pop_state and then reorder the list for easy
    # insertion
    names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

    pop_list_t_1        <- pop_list_t_1[names(pop_state)]

    # Drop the NA entry created by the lambda slot.

    pop_list_t_1        <- pop_list_t_1[!is.na(names(pop_list_t_1))]

    # Next, compute population sizes and time step growth rates.
    # this is required because if normalize_pop_size = TRUE,
    # there isn't any way to compute growth rates at each
    # time step after the fact

    pop_size_t_1 <- vapply(
      pop_list_t_1,
      function(x) sum(x),
      numeric(1L)
    ) %>%
      Reduce(f = '+', x = ., init = 0)

    pop_size_t <- .pop_size(pop_state, i)

    if(normal) {

      pop_list_t_1 <- lapply(pop_list_t_1,
                             function(x, pop_size) {
                               x / pop_size
                             },
                             pop_size = pop_size_t_1)

    }

    pop_list_t_1$lambda <- pop_size_t_1 / pop_size_t

    pop_list_t_1        <- pop_list_t_1[names(pop_state)]

    pop_state <- update_pop_state(pop_state, pop_list_t_1, iter = i)


    # Now, update the names so that we can bind the new n_*_t to
    # main_env so that the next iteration is evaluated correctly

    names(pop_list_t_1) <- names(pop_list_t_1) %>%
      paste(., '_t', sep = "")


    rlang::env_bind(.env = main_env,
                    !!! pop_list_t_1)

  }

  return(pop_state)
}


#' @noRd

.iterate_model.simple_di_stoch_param <- function(proto_ipm,
                                                sub_kern_list,
                                                current_iteration,
                                                kern_seq,
                                                pop_state,
                                                main_env,
                                                k_row,
                                                normal) {

  # Select kernels from kern_seq, or just use all of them if it's
  # a deterministic simulation. Similarly, we need to subset the k_rows
  # object so that .eval_general_det doesn't loop over ones which include
  # expressions for other levels of 'kern_seq'

  if(!is.null(kern_seq)) {

    kern_ind <- grepl(kern_seq[current_iteration], names(sub_kern_list))

    k_row_ind <- which(grepl(kern_seq[current_iteration], k_row$kernel_id))


    use_kerns <- sub_kern_list[kern_ind]
    use_k     <- k_row[k_row_ind, ]

    # Deals with the case where only a subset of kernels have suffixes.
    # In that case, we need to include the ones that don't have them
    # in the use_kerns list every single time!

    if(any(!proto_ipm$uses_par_sets) && any(proto_ipm$uses_par_sets)) {

      nm_ind <- proto_ipm$kernel_id[!proto_ipm$uses_par_sets]
      to_add <- sub_kern_list[nm_ind]

      use_kerns <- c(use_kerns, to_add)

    }


  } else {

    use_kerns <- sub_kern_list
    use_k     <- k_row

  }

  # Need to return an iterated pop_state object

  pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                    proto_ipm     = proto_ipm,
                                    sub_kern_list = use_kerns,
                                    main_env      = main_env)

  # make names match pop_state and then reorder the list for easy
  # insertion
  names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

  pop_list_t_1        <- pop_list_t_1[names(pop_state)]

  # Drop the NA entry created by the lambda slot.

  pop_list_t_1        <- pop_list_t_1[!is.na(names(pop_list_t_1))]

  # Next, compute population sizes and time step growth rates.
  # this is required because if normalize_pop_size = TRUE,
  # there isn't any way to compute growth rates at each
  # time step after the fact

  pop_size_t_1 <- vapply(
    pop_list_t_1,
    function(x) sum(x),
    numeric(1L)
  ) %>%
    Reduce(f = '+', x = ., init = 0)

  pop_size_t <- .pop_size(pop_state, current_iteration)

  if(normal) {

    pop_list_t_1 <- lapply(pop_list_t_1,
                           function(x, pop_size) {
                             x / pop_size
                           },
                           pop_size = pop_size_t_1)

  }

  pop_list_t_1$lambda <- pop_size_t_1 / pop_size_t

  pop_list_t_1        <- pop_list_t_1[names(pop_state)]

  pop_state <- update_pop_state(pop_state,
                                pop_list_t_1,
                                iter = current_iteration)

  # Now, update the names so that we can bind the new n_*_t to
  # main_env so that the next iteration is evaluated correctly

  names(pop_list_t_1) <- names(pop_list_t_1) %>%
    paste(., '_t', sep = "")


  rlang::env_bind(.env = main_env,
                  !!! pop_list_t_1)

  return(pop_state)
}


#' @noRd

.iterate_model.general_di_det <- function(proto_ipm,
                                          k_row,
                                          sub_kern_list,
                                          iterations,
                                          pop_state,
                                          main_env,
                                          normal)  {


  if(any(proto_ipm$uses_par_sets)) {

    par_sets <- proto_ipm$par_set_indices[proto_ipm$uses_par_sets]
    par_sets <- par_sets[!duplicated(par_sets)]

    levs       <- .make_par_set_indices(par_sets)
    lambda_nms <- paste("lambda", levs, sep = "_")

    pop_state$lambda <- NULL

    lambdas    <- vector('list', length = length(levs))
    names(lambdas) <- lambda_nms

    lambdas <- lapply(lambdas,
                      function(x, iterations) {
                        x <- array(NA_real_, dim = c(1, (iterations + 1)))
                      },
                      iterations = iterations)

    pop_state <- c(pop_state, lambdas)


    for(i in seq_along(levs)) {

      use_lev <- paste("(_", levs[i], ")$", sep = "")

      use_kerns <- sub_kern_list[grepl(use_lev, names(sub_kern_list))]

      if(!all(proto_ipm$uses_par_sets)) {

        append <- proto_ipm$kernel_id[!proto_ipm$uses_par_sets]
        use_kerns <- c(use_kerns, sub_kern_list[append])

      }

      use_pop   <- pop_state[grepl(use_lev, names(pop_state))]
      use_k     <- k_row[grepl(use_lev, k_row$kernel_id), ]

      for(j in seq_len(iterations)) {

        .bind_iter_var(main_env, j)

        pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                          proto_ipm     = proto_ipm,
                                          sub_kern_list = use_kerns,
                                          main_env    = main_env)

        names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

        pop_list_t_1        <- pop_list_t_1[names(use_pop)]

        # Drop the NA entry created by the lambda slot.

        pop_list_t_1        <- pop_list_t_1[!is.na(names(pop_list_t_1))]

        # Next, compute population sizes and time step growth rates.
        # this is required because if normalize_pop_size = TRUE,
        # there isn't any way to compute growth rates at each
        # time step after the fact

        pop_size_t_1 <- vapply(
          pop_list_t_1,
          function(x) sum(x),
          numeric(1L)
        ) %>%
          Reduce(f = '+', x = ., init = 0)

        pop_size_t <- .pop_size(use_pop, j)

        if(normal) {

          pop_list_t_1 <- lapply(pop_list_t_1,
                                 function(x, pop_size) {
                                   x / pop_size
                                 },
                                 pop_size = pop_size_t_1)

        }

        pop_list_t_1$lambda <- pop_size_t_1 / pop_size_t
        use_l_nm            <- names(use_pop)[grepl("lambda", names(use_pop))]

        names(pop_list_t_1)[names(pop_list_t_1) == 'lambda'] <- use_l_nm

        pop_list_t_1 <- pop_list_t_1[names(use_pop)]


        use_pop <- update_pop_state(use_pop, pop_list_t_1, j)

        # Now, update the names so that we can bind the new n_*_t to
        # main_env so that the next iteration is evaluated correctly

        names(pop_list_t_1) <- names(pop_list_t_1) %>%
          paste(., '_t', sep = "")


        rlang::env_bind(.env = main_env,
                        !!! pop_list_t_1)

      } # end single level simulation

      pop_state[grepl(use_lev, names(pop_state))] <- use_pop

    } # end par_setarchical levels loop


  } else {

    for(i in seq_len(iterations)) {

      # Select kernels from kern_seq, or just use all of them if it's
      # a deterministic simulation. Similarly, we need to subset the k_rows
      # object so that .eval_general_det doesn't loop over ones which include
      # expressions for other levels of 'kern_seq'

      use_kerns <- sub_kern_list
      use_k     <- k_row

      .bind_iter_var(main_env, i)

      pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                        proto_ipm     = proto_ipm,
                                        sub_kern_list = use_kerns,
                                        main_env    = main_env)


      # make names match pop_state and then reorder the list for easy
      # insertion
      names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

      pop_list_t_1        <- pop_list_t_1[names(pop_state)]

      # Drop the NA entry created by the lambda slot.

      pop_list_t_1        <- pop_list_t_1[!is.na(names(pop_list_t_1))]

      # Next, compute population sizes and time step growth rates.
      # this is required because if normalize_pop_size = TRUE,
      # there isn't any way to compute growth rates at each
      # time step after the fact

      pop_size_t_1 <- vapply(
        pop_list_t_1,
        function(x) sum(x),
        numeric(1L)
      ) %>%
        Reduce(f = '+', x = ., init = 0)

      pop_size_t <- .pop_size(pop_state, i)

      if(normal) {

        pop_list_t_1 <- lapply(pop_list_t_1,
                               function(x, pop_size) {
                                 x / pop_size
                               },
                               pop_size = pop_size_t_1)

      }

      pop_list_t_1$lambda <- pop_size_t_1 / pop_size_t

      pop_list_t_1 <- pop_list_t_1[names(pop_state)]

      pop_state <- update_pop_state(pop_state, pop_list_t_1, i)

      # Now, update the names so that we can bind the new n_*_t to
      # main_env so that the next iteration is evaluated correctly

      names(pop_list_t_1) <- names(pop_list_t_1) %>%
        paste(., '_t', sep = "")


      rlang::env_bind(.env = main_env,
                      !!! pop_list_t_1)

    }
  }

  # Remove the NA at 1 - it is a dummy so that purrr::map2 can work
  # properly

  pop_state[grepl('lambda', names(pop_state))] <- lapply(
    pop_state[grepl("lambda", names(pop_state))],
    function(x) {
      out <- x[ , -1, drop = FALSE]
      return(out)
    }
  )

  return(pop_state)

}

#' @noRd

.iterate_model.general_di_stoch_kern <- function(proto_ipm,
                                                 k_row,
                                                 sub_kern_list,
                                                 iterations,
                                                 kern_seq,
                                                 pop_state,
                                                 main_env,
                                                 normal,
                                                 report_progress) {


  for(i in seq_len(iterations)) {

    # Select kernels from kern_seq, or just use all of them if it's
    # a deterministic simulation. Similarly, we need to subset the k_rows
    # object so that .eval_general_det doesn't loop over ones which include
    # expressions for other levels of 'kern_seq'

    if(!is.null(kern_seq)) {

      if(is.character(kern_seq)) {

        kern_ind <- grepl(kern_seq[i], names(sub_kern_list))

        k_row_ind <- which(grepl(kern_seq[i], k_row$kernel_id))


        use_kerns <- sub_kern_list[kern_ind]
        use_k     <- k_row[k_row_ind, ]

        # Deals with the case where only a subset of kernels have suffixes.
        # In that case, we need to include the ones that don't have them
        # in the use_kerns list every single time!

        if(any(!proto_ipm$uses_par_sets) && any(proto_ipm$uses_par_sets)) {

          nm_ind <- proto_ipm$kernel_id[!proto_ipm$uses_par_sets]
          to_add <- sub_kern_list[nm_ind]

          use_kerns <- c(use_kerns, to_add)

        }

      }

    }

    .bind_iter_var(main_env, i)
    .stoch_progress_message(report_progress, iterations, i)


    pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                      proto_ipm     = proto_ipm,
                                      sub_kern_list = use_kerns,
                                      main_env    = main_env)


    # make names match pop_state and then reorder the list for easy
    # insertion
    names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

    pop_list_t_1        <- pop_list_t_1[names(pop_state)]

    # Drop the NA entry created by the lambda slot.

    pop_list_t_1        <- pop_list_t_1[!is.na(names(pop_list_t_1))]

    # Next, compute population sizes and time step growth rates.
    # this is required because if normalize_pop_size = TRUE,
    # there isn't any way to compute growth rates at each
    # time step after the fact

    pop_size_t_1 <- vapply(
      pop_list_t_1,
      function(x) sum(x),
      numeric(1L)
    ) %>%
      Reduce(f = '+', x = ., init = 0)

    pop_size_t <- .pop_size(pop_state, i)

    if(normal) {

      pop_list_t_1 <- lapply(pop_list_t_1,
                             function(x, pop_size) {
                               x / pop_size
                             },
                             pop_size = pop_size_t_1)

    }

    pop_list_t_1$lambda <- pop_size_t_1 / pop_size_t

    pop_list_t_1 <- pop_list_t_1[names(pop_state)]

    pop_state <- update_pop_state(pop_state, pop_list_t_1, i)

    # Now, update the names so that we can bind the new n_*_t to
    # main_env so that the next iteration is evaluated correctly

    names(pop_list_t_1) <- names(pop_list_t_1) %>%
      paste(., '_t', sep = "")


    rlang::env_bind(.env = main_env,
                    !!! pop_list_t_1)

  }

  return(pop_state)


}


#' @noRd

.iterate_model.general_di_stoch_param <- function(proto_ipm,
                                                  k_row,
                                                  sub_kern_list,
                                                  current_iteration,
                                                  kern_seq,
                                                  pop_state,
                                                  main_env,
                                                  normal) {

  # Select kernels from kern_seq, or just use all of them if it's
  # a deterministic simulation. Similarly, we need to subset the k_rows
  # object so that .eval_general_det doesn't loop over ones which include
  # expressions for other levels of 'kern_seq'

  if(!is.null(kern_seq)) {

    if(is.character(kern_seq)) {

      kern_ind <- grepl(kern_seq[current_iteration], names(sub_kern_list))

      k_row_ind <- which(grepl(kern_seq[current_iteration], k_row$kernel_id))


      use_kerns <- sub_kern_list[kern_ind]
      use_k     <- k_row[k_row_ind, ]

      # Deals with the case where only a subset of kernels have suffixes.
      # In that case, we need to include the ones that don't have them
      # in the use_kerns list every single time!

      if(any(!proto_ipm$uses_par_sets) && any(proto_ipm$uses_par_sets)) {

        nm_ind <- proto_ipm$kernel_id[!proto_ipm$uses_par_sets]
        to_add <- sub_kern_list[nm_ind]

        use_kerns <- c(use_kerns, to_add)

      }

    }

  } else {

    use_kerns <- sub_kern_list
    use_k     <- k_row

  }

  pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                    proto_ipm     = proto_ipm,
                                    sub_kern_list = use_kerns,
                                    main_env    = main_env)


  # make names match pop_state and then reorder the list for easy
  # insertion
  names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

  pop_list_t_1        <- pop_list_t_1[names(pop_state)]

  # Drop the NA entry created by the lambda slot.

  pop_list_t_1        <- pop_list_t_1[!is.na(names(pop_list_t_1))]

  # Next, compute population sizes and time step growth rates.
  # this is required because if normalize_pop_size = TRUE,
  # there isn't any way to compute growth rates at each
  # time step after the fact

  pop_size_t_1 <- vapply(
    pop_list_t_1,
    function(x) sum(x),
    numeric(1L)
  ) %>%
    Reduce(f = '+', x = ., init = 0)

  pop_size_t <- .pop_size(pop_state, current_iteration)

  if(normal) {

    pop_list_t_1 <- lapply(pop_list_t_1,
                           function(x, pop_size) {
                             x / pop_size
                           },
                           pop_size = pop_size_t_1)

  }

  pop_list_t_1$lambda <- pop_size_t_1 / pop_size_t

  pop_list_t_1 <- pop_list_t_1[names(pop_state)]

  pop_state <- update_pop_state(pop_state, pop_list_t_1, current_iteration)

  # Now, update the names so that we can bind the new n_*_t to
  # main_env so that the next iteration is evaluated correctly

  names(pop_list_t_1) <- names(pop_list_t_1) %>%
    paste(., '_t', sep = "")


  rlang::env_bind(.env = main_env,
                  !!! pop_list_t_1)



  return(pop_state)


}


#' @noRd
.iterate_model.simple_dd_det <- function(proto_ipm,
                                         k_row,
                                         sub_kern_list,
                                         current_iteration,
                                         iterations,
                                         pop_state,
                                         main_env,
                                         normal) {


  i <- current_iteration

  # For par_setarchical deterministic models, we need to create a determinstic
  # simulation for each level of the model

  if(any(proto_ipm$uses_par_sets)) {

      par_sets <- proto_ipm$par_set_indices[proto_ipm$uses_par_sets]
      par_sets <- par_sets[!duplicated(par_sets)]

      levs       <- .make_par_set_indices(par_sets)

      # Set up altered lambda list names if first iteration

      if(i == 1) {

        lambda_nms <- paste("lambda", levs, sep = "_")

        pop_state$lambda <- NULL

        lambdas        <- vector('list', length = length(levs))
        names(lambdas) <- lambda_nms

        lambdas <- lapply(lambdas,
                          function(x, iterations) {
                            x <- array(NA_real_, dim = c(1, (iterations + 1)))
                          },
                          iterations = iterations)

        pop_state <- c(pop_state, lambdas)

    }

    # otherwise, just iterate the model

    for(j in seq_along(levs)) {

      use_lev   <- paste("(_", levs[j], ")$", sep = "")

      use_kerns <- sub_kern_list[grepl(use_lev, names(sub_kern_list))]

      if(!all(proto_ipm$uses_par_sets)) {

        append    <- proto_ipm$kernel_id[!proto_ipm$uses_par_sets]
        use_kerns <- c(use_kerns, sub_kern_list[append])

      }

      use_pop <- pop_state[grepl(use_lev, names(pop_state))]
      use_k   <- k_row[grepl(use_lev, k_row$kernel_id), ]


      # Bind a helper for the model iteration number. Users can use this
      # to specify lagged effects.

      .bind_iter_var(main_env, i)

      pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                        proto_ipm     = proto_ipm,
                                        sub_kern_list = use_kerns,
                                        main_env      = main_env)

      names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

      pop_list_t_1        <- pop_list_t_1[names(use_pop)]

      # Drop the NA entry created by the lambda slot.

      pop_list_t_1        <- pop_list_t_1[!is.na(names(pop_list_t_1))]

      # Next, compute population sizes and time step growth rates.
      # this is required because if normalize_pop_size = TRUE,
      # there isn't any way to compute growth rates at each
      # time step after the fact

      pop_size_t_1 <- vapply(
        pop_list_t_1,
        function(x) sum(x),
        numeric(1L)
      ) %>%
        Reduce(f = '+', x = ., init = 0)

      pop_size_t <- .pop_size(use_pop, i)

      if(normal) {

        pop_list_t_1 <- lapply(pop_list_t_1,
                               function(x, pop_size) {
                                 x / pop_size
                               },
                               pop_size = pop_size_t_1)

      }

      pop_list_t_1$lambda <- pop_size_t_1 / pop_size_t
      use_l_nm            <- names(use_pop)[grepl("lambda", names(use_pop))]

      names(pop_list_t_1)[names(pop_list_t_1) == 'lambda'] <- use_l_nm

      pop_list_t_1 <- pop_list_t_1[names(use_pop)]

      use_pop <- update_pop_state(use_pop, pop_list_t_1, i)

      # Now, update the names so that we can bind the new n_*_t to
      # main_env so that the next iteration is evaluated correctly

      names(pop_list_t_1) <- names(pop_list_t_1) %>%
        paste(., '_t', sep = "")


      rlang::env_bind(.env = main_env,
                      !!! pop_list_t_1)


      pop_state[grepl(use_lev, names(pop_state))] <- use_pop

    } # end par_setarchical levels loop


  } else {

    .bind_iter_var(main_env, i)

    pop_list_t_1 <- .eval_general_det(k_row         = k_row,
                                      proto_ipm     = proto_ipm,
                                      sub_kern_list = sub_kern_list,
                                      main_env      = main_env)

    # make names match pop_state and then reorder the list for easy
    # insertion
    names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

    pop_list_t_1        <- pop_list_t_1[names(pop_state)]

    # Drop the NA entry created by the lambda slot.

    pop_list_t_1        <- pop_list_t_1[!is.na(names(pop_list_t_1))]

    # Next, compute population sizes and time step growth rates.
    # this is required because if normalize_pop_size = TRUE,
    # there isn't any way to compute growth rates at each
    # time step after the fact

    pop_size_t_1 <- vapply(
      pop_list_t_1,
      function(x) sum(x),
      numeric(1L)
    ) %>%
      Reduce(f = '+', x = ., init = 0)

    pop_size_t <- .pop_size(pop_state, i)

    if(normal) {

      pop_list_t_1 <- lapply(pop_list_t_1,
                             function(x, pop_size) {
                               x / pop_size
                             },
                             pop_size = pop_size_t_1)

    }

    pop_list_t_1$lambda <- pop_size_t_1 / pop_size_t

    pop_list_t_1        <- pop_list_t_1[names(pop_state)]

    pop_state <- update_pop_state(pop_state, pop_list_t_1, i)

    # Now, update the names so that we can bind the new n_*_t to
    # main_env so that the next iteration is evaluated correctly

    names(pop_list_t_1) <- names(pop_list_t_1) %>%
      paste(., '_t', sep = "")


    rlang::env_bind(.env = main_env,
                    !!! pop_list_t_1)


  } # End if(par_setarchical){} else {}

  return(pop_state)


}



#' @noRd

.iterate_model.simple_dd_stoch_kern <- function(proto_ipm,
                                                k_row,
                                                sub_kern_list,
                                                current_iteration,
                                                iterations,
                                                kern_seq,
                                                pop_state,
                                                main_env,
                                                normal,
                                                report_progress) {


  i <- current_iteration


  # Select kernels from kern_seq, or just use all of them if it's
  # a deterministic simulation. Similarly, we need to subset the k_rows
  # object so that .eval_general_det doesn't loop over ones which include
  # expressions for other levels of 'kern_seq'

  if(!is.null(kern_seq)) {

    if(is.character(kern_seq)) {

      kern_ind <- grepl(kern_seq[i], names(sub_kern_list))

      k_row_ind <- which(grepl(kern_seq[i], k_row$kernel_id))


      use_kerns <- sub_kern_list[kern_ind]
      use_k     <- k_row[k_row_ind, ]

      # Deals with the case where only a subset of kernels have suffixes.
      # In that case, we need to include the ones that don't have them
      # in the use_kerns list every single time!

      if(any(!proto_ipm$uses_par_sets) && any(proto_ipm$uses_par_sets)) {

        nm_ind <- proto_ipm$kernel_id[!proto_ipm$uses_par_sets]
        to_add <- sub_kern_list[nm_ind]

        use_kerns <- c(use_kerns, to_add)

      }

    }

  } else {
    stop("something went wrong!")
  }

  # Need to return an iterated pop_state object

  pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                    proto_ipm     = proto_ipm,
                                    sub_kern_list = use_kerns,
                                    main_env    = main_env)

  # make names match pop_state and then reorder the list for easy
  # insertion
  names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

  pop_list_t_1        <- pop_list_t_1[names(pop_state)]

  # Drop the NA entry created by the lambda slot.

  pop_list_t_1        <- pop_list_t_1[!is.na(names(pop_list_t_1))]

  # Next, compute population sizes and time step growth rates.
  # this is required because if normalize_pop_size = TRUE,
  # there isn't any way to compute growth rates at each
  # time step after the fact

  pop_size_t_1 <- vapply(
    pop_list_t_1,
    function(x) sum(x),
    numeric(1L)
  ) %>%
    Reduce(f = '+', x = ., init = 0)

  pop_size_t <- .pop_size(pop_state, i)

  if(normal) {

    pop_list_t_1 <- lapply(pop_list_t_1,
                           function(x, pop_size) {
                             x / pop_size
                           },
                           pop_size = pop_size_t_1)

  }

  pop_list_t_1$lambda <- pop_size_t_1 / pop_size_t

  pop_list_t_1        <- pop_list_t_1[names(pop_state)]

  pop_state <- update_pop_state(pop_state, pop_list_t_1, i)

  # Now, update the names so that we can bind the new n_*_t to
  # main_env so that the next iteration is evaluated correctly

  names(pop_list_t_1) <- names(pop_list_t_1) %>%
    paste(., '_t', sep = "")


  rlang::env_bind(.env = main_env,
                  !!! pop_list_t_1)



  return(pop_state)
}


.iterate_model.simple_dd_stoch_param <- function(proto_ipm,
                                                 k_row,
                                                 sub_kern_list,
                                                 current_iteration,
                                                 iterations,
                                                 kern_seq,
                                                 pop_state,
                                                 main_env,
                                                 normal,
                                                 report_progress) {


  # Select kernels from kern_seq, or just use all of them if it's
  # a deterministic simulation. Similarly, we need to subset the k_rows
  # object so that .eval_general_det doesn't loop over ones which include
  # expressions for other levels of 'kern_seq'

  if(!is.null(kern_seq)) {

    if(is.character(kern_seq)) {

      kern_ind <- grepl(kern_seq[current_iteration], names(sub_kern_list))

      k_row_ind <- which(grepl(kern_seq[current_iteration], k_row$kernel_id))


      use_kerns <- sub_kern_list[kern_ind]
      use_k     <- k_row[k_row_ind, ]

      # Deals with the case where only a subset of kernels have suffixes.
      # In that case, we need to include the ones that don't have them
      # in the use_kerns list every single time!

      if(any(!proto_ipm$uses_par_sets) && any(proto_ipm$uses_par_sets)) {

        nm_ind <- proto_ipm$kernel_id[!proto_ipm$uses_par_sets]
        to_add <- sub_kern_list[nm_ind]

        use_kerns <- c(use_kerns, to_add)

      }
    }

  } else {

    use_kerns <- sub_kern_list
    use_k     <- k_row

  }

  # Need to return an iterated pop_state object

  pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                    proto_ipm     = proto_ipm,
                                    sub_kern_list = use_kerns,
                                    main_env    = main_env)

  # make names match pop_state and then reorder the list for easy
  # insertion
  names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

  pop_list_t_1        <- pop_list_t_1[names(pop_state)]

  # Drop the NA entry created by the lambda slot.

  pop_list_t_1        <- pop_list_t_1[!is.na(names(pop_list_t_1))]

  # Next, compute population sizes and time step growth rates.
  # this is required because if normalize_pop_size = TRUE,
  # there isn't any way to compute growth rates at each
  # time step after the fact

  pop_size_t_1 <- vapply(
    pop_list_t_1,
    function(x) sum(x),
    numeric(1L)
  ) %>%
    Reduce(f = '+', x = ., init = 0)

  pop_size_t <- .pop_size(pop_state, current_iteration)

  if(normal) {

    pop_list_t_1 <- lapply(pop_list_t_1,
                           function(x, pop_size) {
                             x / pop_size
                           },
                           pop_size = pop_size_t_1)

  }

  pop_list_t_1$lambda <- pop_size_t_1 / pop_size_t

  pop_list_t_1        <- pop_list_t_1[names(pop_state)]

  pop_state <- update_pop_state(pop_state, pop_list_t_1, current_iteration)

  # Now, update the names so that we can bind the new n_*_t to
  # main_env so that the next iteration is evaluated correctly

  names(pop_list_t_1) <- names(pop_list_t_1) %>%
    paste(., '_t', sep = "")


  rlang::env_bind(.env = main_env,
                  !!! pop_list_t_1)



  return(pop_state)

}


#' @noRd
.iterate_model.general_dd_det <- function(proto_ipm,
                                          k_row,
                                          sub_kern_list,
                                          current_iteration,
                                          iterations,
                                          pop_state,
                                          main_env,
                                          normal) {


  i <- current_iteration

  # For par_setarchical deterministic models, we need to create a determinstic
  # simulation for each level of the model

  if(any(proto_ipm$uses_par_sets)) {

    par_sets <- proto_ipm$par_set_indices[proto_ipm$uses_par_sets]
    par_sets <- par_sets[!duplicated(par_sets)]

    levs       <- .make_par_set_indices(par_sets)

    # Set up altered lambda list names if first iteration

    if(i == 1) {

      lambda_nms <- paste("lambda", levs, sep = "_")

      pop_state$lambda <- NULL

      lambdas        <- vector('list', length = length(levs))
      names(lambdas) <- lambda_nms

      lambdas <- lapply(lambdas,
                        function(x, iterations) {
                          x <- array(NA_real_, dim = c(1, (iterations + 1)))
                        },
                        iterations = iterations)

      pop_state <- c(pop_state, lambdas)

    }

    # otherwise, just iterate the model

    for(j in seq_along(levs)) {

      use_lev   <- paste("(_", levs[j], ")$", sep = "")

      use_kerns <- sub_kern_list[grepl(use_lev, names(sub_kern_list))]

      if(!all(proto_ipm$uses_par_sets)) {

        append    <- proto_ipm$kernel_id[!proto_ipm$uses_par_sets]
        use_kerns <- c(use_kerns, sub_kern_list[append])

      }

      use_pop <- pop_state[grepl(use_lev, names(pop_state))]
      use_k   <- k_row[grepl(use_lev, k_row$kernel_id), ]


      # Bind a helper for the model iteration number. Users can use this
      # to specify lagged effects.

      .bind_iter_var(main_env, i)

      pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                        proto_ipm     = proto_ipm,
                                        sub_kern_list = use_kerns,
                                        main_env      = main_env)

      names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

      pop_list_t_1        <- pop_list_t_1[names(use_pop)]

      # Drop the NA entry created by the lambda slot.

      pop_list_t_1        <- pop_list_t_1[!is.na(names(pop_list_t_1))]

      # Next, compute population sizes and time step growth rates.
      # this is required because if normalize_pop_size = TRUE,
      # there isn't any way to compute growth rates at each
      # time step after the fact

      pop_size_t_1 <- vapply(
        pop_list_t_1,
        function(x) sum(x),
        numeric(1L)
      ) %>%
        Reduce(f = '+', x = ., init = 0)

      pop_size_t <- .pop_size(use_pop, i)

      if(normal) {

        pop_list_t_1 <- lapply(pop_list_t_1,
                               function(x, pop_size) {
                                 x / pop_size
                               },
                               pop_size = pop_size_t_1)

      }

      pop_list_t_1$lambda <- pop_size_t_1 / pop_size_t
      use_l_nm            <- names(use_pop)[grepl("lambda", names(use_pop))]

      names(pop_list_t_1)[names(pop_list_t_1) == 'lambda'] <- use_l_nm

      pop_list_t_1 <- pop_list_t_1[names(use_pop)]

      use_pop <- update_pop_state(use_pop, pop_list_t_1, i)

      # Now, update the names so that we can bind the new n_*_t to
      # main_env so that the next iteration is evaluated correctly

      names(pop_list_t_1) <- names(pop_list_t_1) %>%
        paste(., '_t', sep = "")


      rlang::env_bind(.env = main_env,
                      !!! pop_list_t_1)


      pop_state[grepl(use_lev, names(pop_state))] <- use_pop

    } # end par_setarchical levels loop


  } else {

    .bind_iter_var(main_env, i)

    pop_list_t_1 <- .eval_general_det(k_row         = k_row,
                                      proto_ipm     = proto_ipm,
                                      sub_kern_list = sub_kern_list,
                                      main_env      = main_env)

    # make names match pop_state and then reorder the list for easy
    # insertion
    names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

    pop_list_t_1        <- pop_list_t_1[names(pop_state)]

    # Drop the NA entry created by the lambda slot.

    pop_list_t_1        <- pop_list_t_1[!is.na(names(pop_list_t_1))]

    # Next, compute population sizes and time step growth rates.
    # this is required because if normalize_pop_size = TRUE,
    # there isn't any way to compute growth rates at each
    # time step after the fact

    pop_size_t_1 <- vapply(
      pop_list_t_1,
      function(x) sum(x),
      numeric(1L)
    ) %>%
      Reduce(f = '+', x = ., init = 0)

    pop_size_t <- .pop_size(pop_state, i)

    if(normal) {

      pop_list_t_1 <- lapply(pop_list_t_1,
                             function(x, pop_size) {
                               x / pop_size
                             },
                             pop_size = pop_size_t_1)

    }

    pop_list_t_1$lambda <- pop_size_t_1 / pop_size_t

    pop_list_t_1        <- pop_list_t_1[names(pop_state)]

    pop_state           <- update_pop_state(pop_state, pop_list_t_1, i)

    # Now, update the names so that we can bind the new n_*_t to
    # main_env so that the next iteration is evaluated correctly

    names(pop_list_t_1) <- names(pop_list_t_1) %>%
      paste(., '_t', sep = "")


    rlang::env_bind(.env = main_env,
                    !!! pop_list_t_1)


  } # End if(par_setarchical){} else {}

  return(pop_state)


}

#' @noRd

.iterate_model.general_dd_stoch_kern <- function(proto_ipm,
                                                 k_row,
                                                 sub_kern_list,
                                                 current_iteration,
                                                 iterations,
                                                 kern_seq,
                                                 pop_state,
                                                 main_env,
                                                 normal,
                                                 report_progress) {


  i <- current_iteration


  # Select kernels from kern_seq, or just use all of them if it's
  # a deterministic simulation. Similarly, we need to subset the k_rows
  # object so that .eval_general_det doesn't loop over ones which include
  # expressions for other levels of 'kern_seq'

  if(!is.null(kern_seq)) {

    if(is.character(kern_seq)) {

      kern_ind <- grepl(kern_seq[i], names(sub_kern_list))

      k_row_ind <- which(grepl(kern_seq[i], k_row$kernel_id))


      use_kerns <- sub_kern_list[kern_ind]
      use_k     <- k_row[k_row_ind, ]

      # Deals with the case where only a subset of kernels have suffixes.
      # In that case, we need to include the ones that don't have them
      # in the use_kerns list every single time!

      if(any(!proto_ipm$uses_par_sets) && any(proto_ipm$uses_par_sets)) {

        nm_ind <- proto_ipm$kernel_id[!proto_ipm$uses_par_sets]
        to_add <- sub_kern_list[nm_ind]

        use_kerns <- c(use_kerns, to_add)

      }

    }

  } else {
    stop("something went wrong!")
  }

  # Need to return an iterated pop_state object

  pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                    proto_ipm     = proto_ipm,
                                    sub_kern_list = use_kerns,
                                    main_env    = main_env)

  # make names match pop_state and then reorder the list for easy
  # insertion
  names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

  pop_list_t_1        <- pop_list_t_1[names(pop_state)]

  # Drop the NA entry created by the lambda slot.

  pop_list_t_1        <- pop_list_t_1[!is.na(names(pop_list_t_1))]

  # Next, compute population sizes and time step growth rates.
  # this is required because if normalize_pop_size = TRUE,
  # there isn't any way to compute growth rates at each
  # time step after the fact

  pop_size_t_1 <- vapply(
    pop_list_t_1,
    function(x) sum(x),
    numeric(1L)
  ) %>%
    Reduce(f = '+', x = ., init = 0)

  pop_size_t <- .pop_size(pop_state, i)

  if(normal) {

    pop_list_t_1 <- lapply(pop_list_t_1,
                           function(x, pop_size) {
                             x / pop_size
                           },
                           pop_size = pop_size_t_1)

  }

  pop_list_t_1$lambda <- pop_size_t_1 / pop_size_t

  pop_list_t_1        <- pop_list_t_1[names(pop_state)]

  pop_state           <- update_pop_state(pop_state, pop_list_t_1, i)

  # Now, update the names so that we can bind the new n_*_t to
  # main_env so that the next iteration is evaluated correctly

  names(pop_list_t_1) <- names(pop_list_t_1) %>%
    paste(., '_t', sep = "")


  rlang::env_bind(.env = main_env,
                  !!! pop_list_t_1)



  return(pop_state)
}

#' @noRd

.iterate_model.general_dd_stoch_param <- function(proto_ipm,
                                                 k_row,
                                                 sub_kern_list,
                                                 current_iteration,
                                                 iterations,
                                                 kern_seq,
                                                 pop_state,
                                                 main_env,
                                                 normal,
                                                 report_progress) {


  # Select kernels from kern_seq, or just use all of them if it's
  # a deterministic simulation. Similarly, we need to subset the k_rows
  # object so that .eval_general_det doesn't loop over ones which include
  # expressions for other levels of 'kern_seq'

  if(!is.null(kern_seq)) {

    if(is.character(kern_seq)) {

      kern_ind <- grepl(kern_seq[current_iteration], names(sub_kern_list))

      k_row_ind <- which(grepl(kern_seq[current_iteration], k_row$kernel_id))


      use_kerns <- sub_kern_list[kern_ind]
      use_k     <- k_row[k_row_ind, ]

      # Deals with the case where only a subset of kernels have suffixes.
      # In that case, we need to include the ones that don't have them
      # in the use_kerns list every single time!

      if(any(!proto_ipm$uses_par_sets) && any(proto_ipm$uses_par_sets)) {

        nm_ind <- proto_ipm$kernel_id[!proto_ipm$uses_par_sets]
        to_add <- sub_kern_list[nm_ind]

        use_kerns <- c(use_kerns, to_add)

      }
    }

  } else {

    use_kerns <- sub_kern_list
    use_k     <- k_row

  }

  # Need to return an iterated pop_state object

  pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                    proto_ipm     = proto_ipm,
                                    sub_kern_list = use_kerns,
                                    main_env    = main_env)

  # make names match pop_state and then reorder the list for easy
  # insertion
  names(pop_list_t_1) <- gsub('_t_1', '', names(pop_list_t_1))

  pop_list_t_1        <- pop_list_t_1[names(pop_state)]

  # Drop the NA entry created by the lambda slot.

  pop_list_t_1        <- pop_list_t_1[!is.na(names(pop_list_t_1))]

  # Next, compute population sizes and time step growth rates.
  # this is required because if normalize_pop_size = TRUE,
  # there isn't any way to compute growth rates at each
  # time step after the fact

  pop_size_t_1 <- vapply(
    pop_list_t_1,
    function(x) sum(x),
    numeric(1L)
  ) %>%
    Reduce(f = '+', x = ., init = 0)

  pop_size_t <- .pop_size(pop_state, current_iteration)

  if(normal) {

    pop_list_t_1 <- lapply(pop_list_t_1,
                           function(x, pop_size) {
                             x / pop_size
                           },
                           pop_size = pop_size_t_1)

  }

  pop_list_t_1$lambda <- pop_size_t_1 / pop_size_t

  pop_list_t_1        <- pop_list_t_1[names(pop_state)]

  pop_state           <- update_pop_state(pop_state,
                                          pop_list_t_1,
                                          current_iteration)

  # Now, update the names so that we can bind the new n_*_t to
  # main_env so that the next iteration is evaluated correctly

  names(pop_list_t_1) <- names(pop_list_t_1) %>%
    paste(., '_t', sep = "")


  rlang::env_bind(.env = main_env,
                  !!! pop_list_t_1)



  return(pop_state)

}


# Helpers --------------

#' @noRd

.bind_iter_var <- function(env, it) {

  temp <- list(t = it)

  rlang::env_bind(env,
                  !!! temp)

  invisible(TRUE)

}
