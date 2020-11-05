

test_mat   <- matrix(c(1, 2, 3,
                       2, 4, 6,
                       3, 6, 9),
                     nrow = 3, byrow = TRUE)

test_mat_2 <- matrix(c(test_mat, test_mat, test_mat),
                   nrow = 3)

test_that('.conv_to_asymptotic and is_square work/fail properly', {

  expect_true(is_square(test_mat))
  expect_true(.is_conv_to_asymptotic(test_mat_2))


  expect_false(is_square(test_mat_2))


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
  states <- list(c('ht', 'sb'))

  inv_logit <- function(int, slope, sv) {
    1/(1 + exp(-(int + slope * sv)))
  }

  inv_logit_2 <- function(int, slope, slope_2, sv) {
    1/(1 + exp(-(int + slope * sv + slope_2 * sv ^ 2)))
  }

  general_ipm <- init_ipm("general_di_det") %>%
    define_kernel(
      name          = "P",
      formula       = s * g * d_ht,
      family        = "CC",
      g             = dnorm(ht_2, g_mu, g_sd),
      g_mu          = g_int + g_slope * ht_1,
      s             = inv_logit_2(s_int, s_slope, s_slope_2, ht_1),
      data_list     = data_list,
      states        = states,
      has_hier_effs = FALSE,
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
      has_hier_effs = FALSE
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
      has_hier_effs = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'f_d')
    ) %>%
    define_k(
      name          = "K",
      family        = "IPM",
      n_b_t_1       = stay_discrete %*% n_b_t  + go_discrete %*% n_ht_t,
      n_ht_t_1      = leave_discrete %*% n_b_t + P %*% n_ht_t,
      data_list     = data_list,
      states        = states,
      has_hier_effs = FALSE
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "go_discrete", "stay_discrete", "leave_discrete", "K"),
        int_rule     = c(rep("midpoint", 5)),
        dom_start    = c('ht', "ht", NA_character_, NA_character_, "ht"),
        dom_end      = c('ht', NA_character_, NA_character_, 'ht', 'ht')
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

  my_ipm <- init_ipm('simple_di_det') %>%
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
      has_hier_effs = FALSE,
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
        f_r_int   = 0.5,
        f_r_slope = 0.1,
        f_s_int   = 1.2,
        f_s_slope = 0.03,
        f_d_mu    = 1.2,
        f_d_sd    = 0.7
      ),
      states        = list(c('dbh')),
      has_hier_effs = FALSE,
      evict_cor        = TRUE,
      evict_fun     = truncated_distributions("norm", "f_d")
    ) %>%
    define_k(
      name          = "K",
      family        = "IPM",
      K             = P + F,
      data_list     = list(),
      states        = list(c('dbh')),
      has_hier_effs = FALSE,
      evict_cor        = FALSE
    )  %>%
    # Alternative 2, put the call to make_impl_args_list() inside of define_impl().
    define_impl(
      make_impl_args_list(
        kernel_names = c("K", "P", "F"),
        int_rule     = rep("midpoint", 3),
        dom_start    = rep("dbh", 3),
        dom_end      = rep("dbh",3)
      )
    ) %>%
    define_domains(
      dbh = c(1, 30, 200)
    ) %>%
    make_ipm(normalize_pop_size = FALSE)

  expect_error(is_conv_to_asymptotic(my_ipm),
               regexp = 'pop_state in IPM contains NAs - cannot check for convergence!')

})
