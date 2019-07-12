context('test define_* functions')


states_single <- c('dbh')

impl_args_single <- make_impl_args_list(c('P', 'F', 'K'),
                                        int_rule = rep('midpoint', 3),
                                        dom_start = rep('dbh', 3),
                                        dom_end = rep('dbh', 3))

single_state <- init_ipm('simple_di_det') %>%
  define_kernel("P",
                formula = s_g_mult(s, g),
                family = "CC",
                s = inv_logit(s_int, s_slope, dbh_1),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * dbh_1,
                data_list = list(s_int = 2.2,
                                 s_slope = 0.25,
                                 g_int = 0.2,
                                 g_slope = 1.02,
                                 sd_g = 0.7),
                states = states_single,
                evict = TRUE,
                evict_fun = truncated_distributions(g,
                                                    n_mesh_p = 100)) %>%
  define_kernel('F',
                formula = f_r * f_s * f_d,
                family = 'CC',
                f_r = inv_logit(f_r_int, f_r_slope, dbh_1),
                f_s = exp(f_s_int + f_s_slope * dbh_1),
                f_d = dnorm(dbh_2, mu_fd, sd_fd),
                data_list = list(f_r_int = 0.03,
                                 f_r_slope = 0.015,
                                 f_s_int = 1.3,
                                 f_s_slope = 0.075,
                                 mu_fd = 0.5,
                                 sd_fd = 0.2),
                states = states_single,
                evict = FALSE) %>%#,
  # evict_fun = truncated_distributions(f_d,
  #                                     n_mesh_p = 100)) %>%
  define_k(
    name = 'K',
    K = P + F,
    family = 'IPM',
    data_list = list(),
    states = states_single,
    evict = FALSE) %>%
  define_impl(impl_args_single) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  define_pop_state(dbh = rep(1:50, 2))


states_2 <- c('dbh', 'ht')
impl_args_2 <- make_impl_args_list(c('P_1',"P_2", 'F', 'K'),
                                   int_rule = rep('midpoint', 4),
                                   dom_start = c("ht", "dbh", 'dbh', "dbh"),
                                   dom_end = c("dbh", "dbh", 'ht', "dbh"))

two_state <- init_ipm('simple_di_det') %>%
  define_kernel("P_1",
                formula = s_g_mult(s, g), # make helper for double transposes
                family = "CC",
                s = inv_logit(s_int, s_slope, ht_1),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * ht_1,
                data_list = list(s_int = 0.3,
                                 s_slope = 0.2,
                                 g_int = 1.3,
                                 g_slope = 0.99,
                                 sd_g = 0.6),
                states = states_2,
                evict = TRUE,
                evict_fun = truncated_distributions(g,
                                                    n_mesh_p = 100)) %>%
  define_kernel("P_2",
                formula = s_g_mult(s, g), # make helper for double transposes
                family = "CC",
                s = inv_logit(s_int, s_slope, dbh_1),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * dbh_1,
                data_list = list(s_int = 2.2,
                                 s_slope = 0.25,
                                 g_int = 0.2,
                                 g_slope = 1.02,
                                 sd_g = 0.7),
                states = states_2,
                evict = TRUE,
                evict_fun = truncated_distributions(g,
                                                    n_mesh_p = 100)) %>%
  define_kernel('F',
                formula = f_r * f_s * f_d,
                family = 'CC',
                f_r = inv_logit(f_r_int, f_r_slope, dbh_1),
                f_s = exp(f_s_int + f_s_slope * dbh_1),
                f_d = dnorm(ht_2, mu_fd, sd_fd),
                data_list = list(f_r_int = 0.03,
                                 f_r_slope = 0.015,
                                 f_s_int = 1.3,
                                 f_s_slope = 0.075,
                                 mu_fd = 0.5,
                                 sd_fd = 0.2),
                states = states_2,
                evict = FALSE) %>%#,
  # evict_fun = truncated_distributions(f_d,
  #                                     n_mesh_p = 100)) %>%
  define_k(
    name = 'K',
    K = P + F,
    family = 'IPM',
    data_list = list(),
    states = states_2,
    evict = FALSE) %>%
  define_impl(impl_args_2) %>%
  define_domains(dbh = c(0, 50, 100),
                 ht = c(0, 10, 100)) %>%
  define_pop_state(dbh = rep(1:50, 2),
                   ht = runif(50, 0, 10))


test_that('define_pop_state produces expected outputs', {

  pop_state <- two_state$pop_state

  expect_true(is.list(pop_state))

  classes <- vapply(pop_state, function(x) inherits(x[[1]], 'quosures'), logical(1L))

  expect_true((all(classes)))

  nms <- vapply(pop_state, function(x) names(x[[1]]), character(2L)) %>%
    as.vector() %>%
    unique()

  expect_true(all(nms %in% c('ht', 'dbh')))

  text_exprs <- lapply(unlist(pop_state), rlang::quo_text) %>%
    unlist() %>%
    as.vector() %>%
    unique()

  expect_true(grepl('runif', text_exprs[2]))
  expect_true(grepl('rep', text_exprs[1]))

  # define_pop_state should not produce errors related to evaluation as no
  # evaluation occurs - .check_pop_state needs separate unit tests
})

