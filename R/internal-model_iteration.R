#' @noRd

.iterate_model <- function(proto_ipm, ...) {

  UseMethod('.iterate_model')
}


.iterate_model.simple_di_det <- function(proto_ipm,
                                         iterators,
                                         sub_kern_list,
                                         iterations,
                                         pop_state,
                                         master_env,
                                         k_row,
                                         normal) {

  # For hierarchical deterministic models, we need to create a determinstic
  # simulation for each level of the model

  if(any(proto_ipm$has_hier_effs)) {

    hier_effs <- proto_ipm$levels_hier_effs[proto_ipm$has_hier_effs]
    hier_effs <- hier_effs[!duplicated(hier_effs)]

    levs       <- .make_hier_levels(hier_effs)
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

      use_lev <- levs[i]

      use_kerns <- sub_kern_list[grepl(use_lev, names(sub_kern_list))]

      if(!all(proto_ipm$has_hier_effs)) {

        append <- proto_ipm$kernel_id[!proto_ipm$has_hier_effs]
        use_kerns <- c(use_kerns, sub_kern_list[append])

      }

      use_pop   <- pop_state[grepl(use_lev, names(pop_state))]
      use_k     <- k_row[grepl(use_lev, k_row$kernel_id), ]

      for(j in seq_len(iterations)) {

        pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                          proto_ipm     = proto_ipm,
                                          sub_kern_list = use_kerns,
                                          pop_state     = use_pop,
                                          master_env    = master_env)

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

        use_pop           <- purrr::map2(.x = use_pop,
                                         .y = pop_list_t_1,
                                         .f = function(.x, .y, iteration) {

                                           .x[ , (iteration + 1)] <- .y

                                           return(.x)
                                         },
                                         iteration = j
        )

        # Now, update the names so that we can bind the new n_*_t to
        # master_env so that the next iteration is evaluated correctly

        names(pop_list_t_1) <- names(pop_list_t_1) %>%
          paste(., '_t', sep = "")


        rlang::env_bind(.env = master_env,
                        !!! pop_list_t_1)

      } # end single level simulation

      pop_state[grepl(use_lev, names(pop_state))] <- use_pop

    } # end hierarchical levels loop


  } else {


    for(i in seq_len(iterations)) {

      # Simple model, no hierarchical stuff going on. We can just
      # stick the kernels and their expressions right into eval_general_det

      pop_list_t_1 <- .eval_general_det(k_row         = k_row,
                                        proto_ipm     = proto_ipm,
                                        sub_kern_list = sub_kern_list,
                                        pop_state     = pop_state,
                                        master_env    = master_env)

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

      pop_state           <- purrr::map2(.x = pop_state,
                                         .y = pop_list_t_1,
                                         .f = function(.x, .y, iteration) {

                                           .x[ , (iteration + 1)] <- .y

                                           return(.x)
                                         },
                                         iteration = i
      )

      # Now, update the names so that we can bind the new n_*_t to
      # master_env so that the next iteration is evaluated correctly

      names(pop_list_t_1) <- names(pop_list_t_1) %>%
        paste(., '_t', sep = "")


      rlang::env_bind(.env = master_env,
                      !!! pop_list_t_1)

    } # End no hierarchical stuff

  } # End if(hierarchical){} else {}

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
                                                iterators,
                                                sub_kern_list,
                                                iterations,
                                                kern_seq,
                                                pop_state,
                                                master_env,
                                                k_row,
                                                normal) {

  if(rlang::is_quosure(pop_state)) {
    pop_state <- rlang::eval_tidy(pop_state)
  }


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

        if(any(!proto_ipm$has_hier_effs) && any(proto_ipm$has_hier_effs)) {

          nm_ind <- proto_ipm$kernel_id[!proto_ipm$has_hier_effs]
          to_add <- sub_kern_list[nm_ind]

          use_kerns <- c(use_kerns, to_add)

        }
      } else {

        kern_ind <- vapply(names(sub_kern_list),
                           function(x) as.integer(strsplit(x, '_')[[1]][2]),
                           integer(1L)) %>%
          which(. == kern_seq[i])

        k_row_ind <- which(grepl(kern_seq[i], k_row$kernel_id))


        use_kerns <- sub_kern_list[kern_ind]

        use_k     <- k_row[k_row_ind, ]
      }

    } else {

      use_kerns <- sub_kern_list
      use_k     <- k_row

    }

    # Need to return an iterated pop_state object

    pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                      proto_ipm     = proto_ipm,
                                      sub_kern_list = use_kerns,
                                      pop_state     = pop_state,
                                      master_env    = master_env)

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

    pop_state           <- purrr::map2(.x = pop_state,
                                       .y = pop_list_t_1,
                                       .f = function(.x, .y, iteration) {

                                         .x[ , (iteration + 1)] <- .y

                                         return(.x)
                                       },
                                       iteration = i
    )

    # Now, update the names so that we can bind the new n_*_t to
    # master_env so that the next iteration is evaluated correctly

    names(pop_list_t_1) <- names(pop_list_t_1) %>%
      paste(., '_t', sep = "")


    rlang::env_bind(.env = master_env,
                    !!! pop_list_t_1)

  }

  return(pop_state)
}


#' @noRd

.iterate_model.simple_di_stoch_param <- function(proto_ipm,
                                                iterators,
                                                sub_kern_list,
                                                current_iteration,
                                                kern_seq,
                                                pop_state,
                                                master_env,
                                                k_row,
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

      if(any(!proto_ipm$has_hier_effs) && any(proto_ipm$has_hier_effs)) {

        nm_ind <- proto_ipm$kernel_id[!proto_ipm$has_hier_effs]
        to_add <- sub_kern_list[nm_ind]

        use_kerns <- c(use_kerns, to_add)

      }
    } else {

      kern_ind <- vapply(names(sub_kern_list),
                         function(x) as.integer(strsplit(x, '_')[[1]][2]),
                         integer(1L)) %>%
        which(. == kern_seq[current_iteration])

      k_row_ind <- which(grepl(kern_seq[current_iteration], k_row$kernel_id))


      use_kerns <- sub_kern_list[kern_ind]

      use_k     <- k_row[k_row_ind, ]
    }

  } else {

    use_kerns <- sub_kern_list
    use_k     <- k_row

  }

  # Need to return an iterated pop_state object

  pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                    proto_ipm     = proto_ipm,
                                    sub_kern_list = use_kerns,
                                    pop_state     = pop_state,
                                    master_env    = master_env)

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

  pop_state           <- purrr::map2(.x = pop_state,
                                     .y = pop_list_t_1,
                                     .f = function(.x, .y, iteration) {

                                       .x[ , (iteration + 1)] <- .y

                                       return(.x)
                                     },
                                     iteration = current_iteration
  )

  # Now, update the names so that we can bind the new n_*_t to
  # master_env so that the next iteration is evaluated correctly

  names(pop_list_t_1) <- names(pop_list_t_1) %>%
    paste(., '_t', sep = "")


  rlang::env_bind(.env = master_env,
                  !!! pop_list_t_1)



  return(pop_state)
}


#' @noRd

.iterate_model.general_di_det <- function(proto_ipm,
                                          k_row,
                                          sub_kern_list,
                                          iterations,
                                          pop_state,
                                          master_env,
                                          normal)  {


  if(any(proto_ipm$has_hier_effs)) {

    hier_effs <- proto_ipm$levels_hier_effs[proto_ipm$has_hier_effs]
    hier_effs <- hier_effs[!duplicated(hier_effs)]

    levs       <- .make_hier_levels(hier_effs)
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

      use_lev <- levs[i]

      use_kerns <- sub_kern_list[grepl(use_lev, names(sub_kern_list))]

      if(!all(proto_ipm$has_hier_effs)) {

        append <- proto_ipm$kernel_id[!proto_ipm$has_hier_effs]
        use_kerns <- c(use_kerns, sub_kern_list[append])

      }

      use_pop   <- pop_state[grepl(use_lev, names(pop_state))]
      use_k     <- k_row[grepl(use_lev, k_row$kernel_id), ]

      for(j in seq_len(iterations)) {

        pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                          proto_ipm     = proto_ipm,
                                          sub_kern_list = use_kerns,
                                          pop_state     = use_pop,
                                          master_env    = master_env)

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

        use_pop           <- purrr::map2(.x = use_pop,
                                         .y = pop_list_t_1,
                                         .f = function(.x, .y, iteration) {

                                           .x[ , (iteration + 1)] <- .y

                                           return(.x)
                                         },
                                         iteration = j
        )

        # Now, update the names so that we can bind the new n_*_t to
        # master_env so that the next iteration is evaluated correctly

        names(pop_list_t_1) <- names(pop_list_t_1) %>%
          paste(., '_t', sep = "")


        rlang::env_bind(.env = master_env,
                        !!! pop_list_t_1)

      } # end single level simulation

      pop_state[grepl(use_lev, names(pop_state))] <- use_pop

    } # end hierarchical levels loop


  } else {

    for(i in seq_len(iterations)) {

      # Select kernels from kern_seq, or just use all of them if it's
      # a deterministic simulation. Similarly, we need to subset the k_rows
      # object so that .eval_general_det doesn't loop over ones which include
      # expressions for other levels of 'kern_seq'

      use_kerns <- sub_kern_list
      use_k     <- k_row



      pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                        proto_ipm     = proto_ipm,
                                        sub_kern_list = use_kerns,
                                        pop_state     = pop_state,
                                        master_env    = master_env)


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

      pop_state           <- purrr::map2(.x = pop_state,
                                         .y = pop_list_t_1,
                                         .f = function(.x, .y, iteration) {

                                           .x[ , (iteration + 1)] <- .y

                                           return(.x)
                                         },
                                         iteration = i
      )

      # Now, update the names so that we can bind the new n_*_t to
      # master_env so that the next iteration is evaluated correctly

      names(pop_list_t_1) <- names(pop_list_t_1) %>%
        paste(., '_t', sep = "")


      rlang::env_bind(.env = master_env,
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
                                                 master_env,
                                                 normal) {


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

        if(any(!proto_ipm$has_hier_effs) && any(proto_ipm$has_hier_effs)) {

          nm_ind <- proto_ipm$kernel_id[!proto_ipm$has_hier_effs]
          to_add <- sub_kern_list[nm_ind]

          use_kerns <- c(use_kerns, to_add)

        }

      } else {

        use_kerns <- sub_kern_list[kern_seq[i]]

        use_k     <- k_row[kern_seq[i] , ]
      }

    } else {

      use_kerns <- sub_kern_list
      use_k     <- k_row

    }

    pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                      proto_ipm     = proto_ipm,
                                      sub_kern_list = use_kerns,
                                      pop_state     = pop_state,
                                      master_env    = master_env)


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

    pop_state           <- purrr::map2(.x = pop_state,
                                       .y = pop_list_t_1,
                                       .f = function(.x, .y, iteration) {

                                         .x[ , (iteration + 1)] <- .y

                                         return(.x)
                                       },
                                       iteration = i
    )

    # Now, update the names so that we can bind the new n_*_t to
    # master_env so that the next iteration is evaluated correctly

    names(pop_list_t_1) <- names(pop_list_t_1) %>%
      paste(., '_t', sep = "")


    rlang::env_bind(.env = master_env,
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
                                                  master_env,
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

      if(any(!proto_ipm$has_hier_effs) && any(proto_ipm$has_hier_effs)) {

        nm_ind <- proto_ipm$kernel_id[!proto_ipm$has_hier_effs]
        to_add <- sub_kern_list[nm_ind]

        use_kerns <- c(use_kerns, to_add)

      }

    } else {

      use_kerns <- sub_kern_list[kern_seq[current_iteration]]

      use_k     <- k_row[kern_seq[current_iteration] , ]
    }

  } else {

    use_kerns <- sub_kern_list
    use_k     <- k_row

  }

  pop_list_t_1 <- .eval_general_det(k_row         = use_k,
                                    proto_ipm     = proto_ipm,
                                    sub_kern_list = use_kerns,
                                    pop_state     = pop_state,
                                    master_env    = master_env)


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

  pop_state           <- purrr::map2(.x = pop_state,
                                     .y = pop_list_t_1,
                                     .f = function(.x, .y, iteration) {

                                       .x[ , (iteration + 1)] <- .y

                                       return(.x)
                                     },
                                     iteration = current_iteration
  )

  # Now, update the names so that we can bind the new n_*_t to
  # master_env so that the next iteration is evaluated correctly

  names(pop_list_t_1) <- names(pop_list_t_1) %>%
    paste(., '_t', sep = "")


  rlang::env_bind(.env = master_env,
                  !!! pop_list_t_1)



  return(pop_state)


}


