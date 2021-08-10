

test_mat   <- matrix(c(1, 2, 3,
                       2, 4, 6,
                       3, 6, 9),
                     nrow = 3, byrow = TRUE)

test_mat_2 <- matrix(c(test_mat, test_mat, test_mat),
                   nrow = 3)

test_df <- cbind(expand.grid(t = 1:3, t_1 = 1:3), value = as.vector(test_mat))

test_that('.conv_to_asymptotic and is_square work/fail properly', {

  expect_true(is_square(test_mat))
  expect_true(.is_conv_to_asymptotic(test_mat_2))


  expect_false(is_square(test_mat_2))


})

test_that("ipm_to_df works properly", {

  expect_equal(test_df, mat_to_df_impl(test_mat))

  data(gen_di_det_ex)

  temp_big <- format_mega_kernel(gen_di_det_ex,
                                 c(stay_discrete, go_discrete,
                                   leave_discrete, P))[[1]]

  target_df <- expand.grid(t = 1:201, t_1 = 1:201, value = NA_real_) %>%
    as.list()
  it <- 1

  for(i in seq_len(nrow(temp_big))) {
    for(j in seq_len(ncol(temp_big))) {

      target_df$value[it] <- temp_big[i, j]
      it <- it + 1

    }
  }

  target_df <- as.data.frame(target_df)

  ipmr_df <- ipm_to_df(gen_di_det_ex, mega_mat = c(stay_discrete, go_discrete,
                                                    leave_discrete, P))

  expect_equal(target_df, ipmr_df$mega_matrix)

})

test_that('exported is_conv_to_asymptotic works as well', {
  init_pop_vec   <- runif(500)
  init_seed_bank <- 20

  data_list <- list(
    g_int     = 5.781,
    g_slope   = 0.988,
    g_sd      = 20.55699,
    s_int     = -0.352,
    s_slope   = 0.122,
    s_slope_2 = -0.000213,
    f_r_int   = -11.46,
    f_r_slope = 0.0835,
    f_s_int   = 2.6204,
    f_s_slope = 0.01256,
    f_d_mu    = 5.6655,
    f_d_sd    = 2.0734,
    e_p       = 0.15,
    g_i       = 0.5067
  )
  L <- 1.02
  U <- 624
  n <- 500
  states <- list(c('ht', 'b'))

  inv_logit <- function(int, slope, sv) {
    1/(1 + exp(-(int + slope * sv)))
  }

  inv_logit_2 <- function(int, slope, slope_2, sv) {
    1/(1 + exp(-(int + slope * sv + slope_2 * sv ^ 2)))
  }

  general_ipm <- init_ipm(sim_gen    = "general",
                          di_dd      = "di",
                          det_stoch  = "det") %>%
    define_kernel(
      name          = "P",
      formula       = s * g * d_ht,
      family        = "CC",
      g             = dnorm(ht_2, g_mu, g_sd),
      g_mu          = g_int + g_slope * ht_1,
      s             = inv_logit_2(s_int, s_slope, s_slope_2, ht_1),
      data_list     = data_list,
      states        = states,
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'g')
    ) %>%
    define_kernel(
      name          = "go_discrete",
      formula       = f_r * f_s * g_i,
      family        = 'CD',
      f_r           = inv_logit(f_r_int, f_r_slope, ht_1),
      f_s           = exp(f_s_int + f_s_slope * ht_1),
      data_list     = data_list,
      states        = states,
      uses_par_sets = FALSE
    ) %>%
    define_kernel(
      name    = 'stay_discrete',
      formula = 0,
      family  = "DD",
      states  = states,
      evict_cor = FALSE
    ) %>%
    define_kernel(
      name          = 'leave_discrete',
      formula       = e_p * f_d * d_ht,
      f_d           = dnorm(ht_2, f_d_mu, f_d_sd),
      family        = 'DC',
      data_list     = data_list,
      states        = states,
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'f_d')
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "go_discrete", "stay_discrete", "leave_discrete"),
        int_rule     = c(rep("midpoint", 4)),
        state_start    = c('ht', "ht", "b", "b"),
        state_end      = c('ht', 'b', 'b', 'ht')
      )
    ) %>%
    define_domains(
      ht = c(L, U, n)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_ht = init_pop_vec,
        n_b  = init_seed_bank
      )
    ) %>%
    make_ipm(iterations = 100,
             usr_funs = list(inv_logit   = inv_logit,
                             inv_logit_2 = inv_logit_2),
             normalize_pop_size = FALSE)

  expect_true(is_conv_to_asymptotic(general_ipm))

  few_iterates <- general_ipm$proto_ipm %>%
    make_ipm(iterations = 3,
             normalize_pop_size = FALSE)

  expect_false(is_conv_to_asymptotic(few_iterates))

  my_ipm <- init_ipm(sim_gen    = "simple",
                     di_dd      = "di",
                     det_stoch  = "det") %>%
    define_kernel(
      name      = "P",
      formula   = s * g,
      family    = "CC",

      s         = 1/(1 + exp(-(s_int + s_slope * dbh_1))),
      g         = dnorm(dbh_2, g_mu, g_sd),
      g_mu      = g_int + g_slope * dbh_1,

      data_list = list(
        s_int     = 0.2,
        s_slope   = 0.5,
        g_int     = 0.1,
        g_slope   = 1.033,
        g_sd      = 2.2
      ),
      states        = list(c('dbh')),
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", "g")
    ) %>%
    define_kernel(
      name      = "F",
      formula   = f_r * f_s * f_d,
      family    = "CC",
      f_r       = 1/(1 + exp(-(f_r_int + f_r_slope * dbh_1))),
      f_s       = exp(f_s_int + f_s_slope * dbh_1),
      f_d       = dnorm(dbh_2, f_d_mu, f_d_sd),
      data_list = list(
        f_r_int   = 5,
        f_r_slope = 0.1,
        f_s_int   = 10,
        f_s_slope = 0.03,
        f_d_mu    = 1.2,
        f_d_sd    = 0.7
      ),
      states        = list(c('dbh')),
      uses_par_sets = FALSE,
      evict_cor        = TRUE,
      evict_fun     = truncated_distributions("norm", "f_d")
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "F"),
        int_rule     = rep("midpoint", 2),
        state_start    = rep("dbh", 2),
        state_end      = rep("dbh", 2)
      )
    ) %>%
    define_domains(
      dbh = c(1, 30, 200)
    ) %>%
    define_pop_state(n_dbh = runif(200)) %>%
    make_ipm(normalize_pop_size = FALSE,
             iterations = 100)

  expect_error(is_conv_to_asymptotic(my_ipm),
               regexp = 'pop_state in IPM contains NAs - cannot check for convergence!')

})
