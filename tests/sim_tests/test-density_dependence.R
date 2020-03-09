# This implementation won't work because we need to build sub-kernels lazily
# at each iteration. This means that it may work for stoch_param methods right now,
# but dd_det/dd_stoch_kern methods need to use .make_sub_kern_lazy instead of
# their current set up.

context('density dependent IPMs')


inv_logit <- function(int, slope, sv) {
  return(1/(1 + exp(-(int + slope * sv))))
}

impl_args <- make_impl_args_list(c('P', 'F', 'K'),
                                 int_rule = rep('midpoint', 3),
                                 dom_start = rep('dbh', 3),
                                 dom_end = rep('dbh', 3))

data_list = list(s_int = -19,
                 s_slope = 0.25,
                 g_int = 0.2,
                 g_slope = 1.02,
                 sd_g = 0.7,
                 f_r_int = 0.003,
                 f_r_slope = 0.015,
                 f_s_int = 1.3,
                 f_s_slope = 0.075,
                 mu_fd = 2,
                 sd_fd = 0.3)

states <- c('dbh', 'dbh')

x <- init_ipm('simple_dd_det') %>%
  define_kernel("P",
                formula = s_g_mult(s, g),
                family = "CC",
                s = inv_logit(s_int, s_slope, dbh_1),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * dbh_1,
                data_list = data_list,
                states = states,
                evict_cor =  TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'g')
  ) %>%
  define_kernel('F',
                formula = f_r * f_s * f_d,
                family = 'CC',
                f_r = inv_logit(f_r_int, f_r_slope, dbh_1),
                f_s = exp(f_s_int + (f_s_slope * exp(-32 * sum(n_dbh_t))) * dbh_1),
                f_d = dnorm(dbh_2, mu_fd, sd_fd),
                data_list = data_list,
                states = states,
                evict_cor =  TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'f_d')
  ) %>%
  define_k('K',
           K = P + F,
           n_dbh_t_1 = K %*% n_dbh_t,
           family = 'IPM',
           data_list = list(),
           states = states,
           evict_cor =  FALSE) %>%
  define_impl(impl_args) %>%
  define_pop_state(
    n_dbh = rpois(100, 50:100)
  ) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit),
           iterate = TRUE,
           iterations = 50)
