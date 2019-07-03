context('test define_* functions')


states_single <- c('dbh')

impl_args_single <- make_impl_args_list(c('P', 'F', 'K'),
                                        int_rule = rep('midpoint', 3),
                                        dom_start = rep('dbh', 3),
                                        dom_end = rep('dbh', 3))

single_state <- init_ipm('simple_di_det') %>%
  define_kernel("P",
                formula = t(s * t(g)), # make helper for double transposes
                family = "CC",
                s = inv_logit(s_int, s_slope, dbh_1),
                g = cell_size_dbh * dnorm(dbh_2, mu_g, sd_g),
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
                f_d = cell_size_dbh * dnorm(dbh_2, mu_fd, sd_fd),
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
    formula = P + F,
    family = 'IPM',
    data_list = list(),
    states = states_single,
    evict = FALSE) %>%
  define_impl(impl_args_single) %>%
  #define_domains(dbh = c(0, 50, 100)) %>%
  define_pop_state(dbh = rep(1:50, 2))


states_2 <- c('dbh', 'ht')
impl_args_2 <- make_impl_args_list(c('P', 'F', 'K'),
                                   int_rule = rep('midpoint', 3),
                                   dom_start = c("dbh", "ht", "dbh"),
                                   dom_end = c("dbh", "dbh", "dbh"))

two_state <- init_ipm('simple_di_det') %>%
  define_kernel("P",
                formula = t(s * t(g)), # make helper for double transposes
                family = "CC",
                s = inv_logit(s_int, s_slope, dbh_1),
                g = cell_size_dbh * dnorm(dbh_2, mu_g, sd_g),
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
                f_r = inv_logit(f_r_int, f_r_slope, ht_1),
                f_s = exp(f_s_int + f_s_slope * ht_1),
                f_d = cell_size_ht * dnorm(ht_2, mu_fd, sd_fd),
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
    formula = P + F,
    family = 'IPM',
    data_list = list(),
    states = states_2,
    evict = FALSE) %>%
  define_impl(impl_args_2) %>%
  define_domains(dbh = c(0, 50, 100),
                 ht = c(0, 10, 50)) #%>%
  # define_pop_state(dbh = rep(1:50, 2),
  #                  ht = runif(50, 0, 10))
