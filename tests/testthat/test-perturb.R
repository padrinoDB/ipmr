# Test prospective perturbation functions

context('sensitivity.simple_di_det_ipm works as expected')

data_list = list(s_int = 2.2,
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

impl_args <- make_impl_args_list(c('P', 'F', 'K'),
                                 int_rule = rep('midpoint', 3),
                                 dom_start = rep('dbh', 3),
                                 dom_end = rep('dbh', 3))

inv_logit <- function(int, slope, sv) {
  return(1/(1 + exp(-(int + slope * sv))))
}

impl_args <- make_impl_args_list(c('P', 'F', 'K'),
                                 int_rule = rep('midpoint', 3),
                                 dom_start = rep('dbh', 3),
                                 dom_end = rep('dbh', 3))

states <- c('dbh', 'dbh')

sim_di_det_1 <- init_ipm('simple_di_det') %>%
  define_kernel("P",
                formula = s_g_mult(s, g),
                family = "CC",
                s = inv_logit(s_int, s_slope, dbh_1),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * dbh_1,
                data_list = data_list,
                states = states,
                evict = TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'g')
  ) %>%
  define_kernel('F',
                formula = f_r * f_s * f_d,
                family = 'CC',
                f_r = inv_logit(f_r_int, f_r_slope, dbh_1),
                f_s = exp(f_s_int + f_s_slope * dbh_1),
                f_d = dnorm(dbh_2, mu_fd, sd_fd),
                data_list = data_list,
                states = states,
                evict = TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'f_d')
  ) %>%
  define_k('K',
           K         = P + F,
           family    = 'IPM',
           data_list = list(),
           states    = states,
           evict     = FALSE) %>%
  define_impl(impl_args) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit))

K <- sim_di_det_1$iterators$K
w <- Re(eigen(K)$vectors[ , 1])
w <- w / sum(w)
v <- Re(eigen(t(K))$vectors[ , 1])
v <- v / sum(v)
l <- Re(eigen(K)$values[1])
h <- 0.5

sens_hand <- outer(v, w) / sum(v * w * h)
elas_hand <- sens_hand * (K / h) / l

sens_ipmr <- sensitivity(sim_di_det_1, what = 'lambda', level = 'kernel')
elas_ipmr <- elasticity(sim_di_det_1, what = 'lambda', level = 'kernel')

test_that('sensitivity/elasticity.simple_di_det_ipm are working', {

  expect_equal(sens_ipmr, sens_hand, tol = 1e-10, check.attributes = FALSE)
  expect_equal(elas_ipmr, elas_hand, tol = 1e-10, check.attributes = FALSE)

})

test_that('elasticity is doing math correctly', {

  expect_equal(sum(elas_ipmr * .get_d_z_vars(sim_di_det_1)^2),
               1L,
               tol = 1e-13)

})

sim_di_det_2 <- init_ipm('simple_di_det') %>%
  define_kernel("P",
                formula = s_g_mult(s, g),
                family = "CC",
                s = inv_logit(s_int, s_slope, dbh_1),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * dbh_1,
                data_list = data_list,
                states = states,
                evict = TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'g')
  ) %>%
  define_kernel('F',
                formula = f_r * f_s * f_d,
                family = 'CC',
                f_r = inv_logit(f_r_int, f_r_slope, dbh_1),
                f_s = exp(f_s_int + f_s_slope * dbh_1),
                f_d = dnorm(dbh_2, mu_fd, sd_fd),
                data_list = data_list,
                states = states,
                evict = TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'f_d')
  ) %>%
  define_k('K',
           K         = P + F,
           family    = 'IPM',
           data_list = list(),
           states    = states,
           evict     = FALSE) %>%
  define_impl(impl_args) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit),
           return_all = TRUE)

test_that('.get_d_z_vars can find master_env in env_list correctly', {

  master_d_z <- .get_d_z_vars(sim_di_det_2)

  expect_equal(master_d_z, h, check.attributes = FALSE)

})


context('sensitivity/elasticity general_di_det_ipm is working')


