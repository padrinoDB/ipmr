
states_single <- c('dbh')

impl_args_single <- make_impl_args_list(c('P', 'F'),
                                        int_rule = rep('midpoint', 2),
                                        state_start = rep('dbh', 2),
                                        state_end = rep('dbh', 2))

L <- 100
U <- 500
n_mesh_p <- 1000

single_state <- init_ipm(sim_gen    = "simple",
                         di_dd      = "di",
                         det_stoch  = "det") %>%
  define_kernel("P",
                formula = s * g,
                family = "CC",
                s = plogis(s_int, s_slope, dbh_1),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * dbh_1,
                data_list = list(s_int = 2.2,
                                 s_slope = 0.25,
                                 g_int = 0.2,
                                 g_slope = 1.02,
                                 sd_g = 0.7),
                states = states_single,
                evict_cor = TRUE,
                evict_fun = truncated_distributions(g,
                                                    n_mesh_p = 100)) %>%
  define_kernel('F',
                formula = f_r * f_s * f_d,
                family = 'CC',
                f_r = plogis(f_r_int, f_r_slope, dbh_1),
                f_s = exp(f_s_int + f_s_slope * dbh_1),
                f_d = dnorm(dbh_2, mu_fd, sd_fd),
                data_list = list(f_r_int = 0.03,
                                 f_r_slope = 0.015,
                                 f_s_int = 1.3,
                                 f_s_slope = 0.075,
                                 mu_fd = 0.5,
                                 sd_fd = 0.2),
                states = states_single,
                evict_cor = FALSE) %>%#,
  # evict_fun = truncated_distributions(f_d,
  #                                     n_mesh_p = 100)) %>%
  define_impl(impl_args_single) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  define_pop_state(n_dbh = rep(1:50, 2))


states_2 <- c('dbh', 'ht')
impl_args_2 <- make_impl_args_list(c('P_1',"P_2", 'F'),
                                   int_rule = rep('midpoint', 3),
                                   state_start = c("ht", "dbh", "dbh"),
                                   state_end = c("dbh", "dbh", 'ht'))

two_state <- init_ipm(sim_gen    = "simple",
                      di_dd      = "di",
                      det_stoch  = "det") %>%
  define_kernel("P_1",
                formula = s * g, # make helper for double transposes
                family = "CC",
                s = plogis(s_int, s_slope, ht_1),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * ht_1,
                data_list = list(s_int = 0.3,
                                 s_slope = 0.2,
                                 g_int = 1.3,
                                 g_slope = 0.99,
                                 sd_g = 0.6),
                states = states_2,
                evict_cor = TRUE,
                evict_fun = truncated_distributions("norm", "g")) %>%
  define_kernel("P_2",
                formula = s * g,
                family = "CC",
                s = plogis(s_int, s_slope, dbh_1),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * dbh_1,
                data_list = list(s_int = 2.2,
                                 s_slope = 0.25,
                                 g_int = 0.2,
                                 g_slope = 1.02,
                                 sd_g = 0.7),
                states = states_2,
                evict_cor = TRUE,
                evict_fun = truncated_distributions("norm", "g")) %>%
  define_kernel('F',
                formula = f_r * f_s * f_d,
                family = 'CC',
                f_r = plogis(f_r_int, f_r_slope, dbh_1),
                f_s = exp(f_s_int + f_s_slope * dbh_1),
                f_d = dnorm(ht_2, mu_fd, sd_fd),
                data_list = list(f_r_int = 0.03,
                                 f_r_slope = 0.015,
                                 f_s_int = 1.3,
                                 f_s_slope = 0.075,
                                 mu_fd = 0.5,
                                 sd_fd = 0.2),
                states = states_2,
                evict_cor = FALSE)  %>%
  define_impl(impl_args_2) %>%
  define_domains(dbh = c(0, 50, 100),
                 ht = c(0, 10, 150)) %>%
  define_pop_state(n_dbh = rep(1:50, 2),
                   n_ht = runif(150, 0, 10))

test_that('define_pop_state produces expected outputs', {

  pop_state <- two_state$pop_state

  expect_true(is.list(pop_state))

  classes <- vapply(pop_state, function(x) inherits(x[[1]], 'quosure'), logical(1L))

  expect_true((all(classes)))

  nms <- vapply(pop_state, function(x) names(x), character(2L)) %>%
    as.vector() %>%
    unique()

  expect_true(all(nms %in% c('pop_state_ht', 'pop_state_dbh')))

  text_exprs <- lapply(unlist(pop_state), rlang::quo_text) %>%
    unlist() %>%
    as.vector() %>%
    unique()

  expect_true(grepl('runif', text_exprs[2]))
  expect_true(grepl('rep', text_exprs[1]))

  expect_error(define_pop_state(two_state,
                                dbh = runif(100)),
               regexp = "All population state names must start with 'n_'")
})

test_that('define_domains can use global variables', {


  two_state <- define_domains(two_state,
                               dbh = c(L, U, n_mesh_p))

  new_dom <- two_state$domain[[1]]$dbh

  expect_equal(new_dom, c(L, U, n_mesh_p))
})


test_that("NA warnings are generated correctly", {

  expect_warning(parameters(single_state) <- list (s_int = NA))
  expect_error(make_ipm(single_state))
  expect_warning(
    init_ipm("simple", "di", "det") %>%
      define_kernel("P_1",
                    formula = s * g, # make helper for double transposes
                    family = "CC",
                    s = inv_logit(s_int, s_slope, ht_1),
                    g = dnorm(dbh_2, mu_g, sd_g),
                    mu_g = g_int + g_slope * ht_1,
                    data_list = list(s_int = 0.3,
                                     s_slope = 0.2,
                                     g_int = 1.3,
                                     g_slope = 0.99,
                                     sd_g = NA),
                    states = states_2))

  expect_warning(
    init_ipm("simple", "di", "stoch", "param") %>%
      define_kernel("P_1",
                    formula = s * g,
                    family = "CC",
                    s = inv_logit(s_int, s_slope, ht_1),
                    g = dnorm(dbh_2, mu_g, sd_g),
                    mu_g = g_int + g_slope * ht_1,
                    data_list = list(s_int = 0.3,
                                     s_slope = 0.2,
                                     g_int = 1.3,
                                     g_slope = 0.99,
                                     sd_g = 0.2),
                    states = states_2) %>%
      define_env_state(env_params = inv_logit(2 * 2),
                       data_list = list(inv_logit = NA))
  )

})
