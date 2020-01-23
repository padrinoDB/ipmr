# Test prospective perturbation functions

context('sensitivity/elasticity.simple_di_det_ipm works as expected')

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

  expect_equal(sens_ipmr, sens_hand, tol = 1e-9, check.attributes = FALSE)
  expect_equal(elas_ipmr, elas_hand, tol = 1e-9, check.attributes = FALSE)

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

# The lower and upper bounds for the continuous state variable and the number
# of meshpoints for the midpoint rule integration.

L <- 1.02
U <- 624
n <- 500

# Initialize the state list and add some helper functions. The survival function
# in this model is a quadratic function.

states <- list(c('ht', 'b'))

inv_logit <- function(int, slope, sv) {
  1/(1 + exp(-(int + slope * sv)))
}

inv_logit_2 <- function(int, slope, slope_2, sv) {
  1/(1 + exp(-(int + slope * sv + slope_2 * sv ^ 2)))
}

gen_di_det_1 <- init_ipm("general_di_det") %>%
  define_kernel(
    name          = "P",
    formula       = s_g_mult(s, g) * d_ht,
    family        = "CC",
    g             = dnorm(ht_2, g_mu, g_sd),
    g_mu          = g_int + g_slope * ht_1,
    s             = inv_logit_2(s_int, s_slope, s_slope_2, ht_1),
    data_list     = data_list,
    states        = states,
    has_hier_effs = FALSE,
    evict         = TRUE,
    evict_fun     = truncated_distributions('norm',
                                            'g')
  ) %>%
  define_kernel(
    name          = "go_discrete",
    formula       = f_r * f_s * g_i * d_ht,
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
    evict = FALSE
  ) %>%
  define_kernel(
    name          = 'leave_discrete',
    formula       = e_p * f_d * d_ht,
    f_d           = dnorm(ht_2, f_d_mu, f_d_sd),
    family        = 'DC',
    data_list     = data_list,
    states        = states,
    has_hier_effs = FALSE,
    evict         = TRUE,
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
                           inv_logit_2 = inv_logit_2))

mega_k <- rbind(
  cbind(gen_di_det_1$sub_kernels$stay_discrete,
        gen_di_det_1$sub_kernels$go_discrete),
  cbind(gen_di_det_1$sub_kernels$leave_discrete,
        gen_di_det_1$sub_kernels$P)
)

w <- Re(eigen(mega_k)$vectors[ , 1])
w <- w / sum(w)
v <- Re(eigen(t(mega_k))$vectors[ , 1])
v <- v / sum(v)
l <- Re(eigen(mega_k)$values[1])

bounds <- seq(L, U, length.out = (n + 1))
h      <- bounds[2] - bounds[1]


mega_k_fun <- rbind(
  cbind(gen_di_det_1$sub_kernels$stay_discrete,
        gen_di_det_1$sub_kernels$go_discrete),
  cbind(gen_di_det_1$sub_kernels$leave_discrete,
        gen_di_det_1$sub_kernels$P)
) / h

sens_mega <- outer(v, w) / sum(v * w * h)
elas_mega <- sens_mega * (mega_k_fun) / l

sens_ipmr <- sensitivity(gen_di_det_1,
                         mega_mat = c(stay_discrete, go_discrete,
                                      leave_discrete, P),
                         mega_vec = c(b, ht))


elas_ipmr <- elasticity(gen_di_det_1,
                        mega_mat = c(stay_discrete, go_discrete,
                                     leave_discrete, P),
                        mega_vec = c(b, ht))


test_that('sensitivity/elasticity.general_di_det_ipm are working', {

  expect_equal(sens_ipmr, sens_mega, tol = 1e-9, check.attributes = FALSE)
  expect_equal(elas_ipmr, elas_mega, tol = 1e-9, check.attributes = FALSE)

  expect_equal(sum(elas_ipmr * .get_d_z_vars(gen_di_det_1) ^ 2),
               1L)

})



