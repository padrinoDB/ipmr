context('simple_di_det plot and print methods')
library(rlang)
library(purrr)
library(mvtnorm)

# lambda methods are tested in each specific classes unit tests

# simple_di_det methods ---------
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

sim_di_det <- init_ipm('simple_di_det') %>%
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
           K = P + F,
           family = 'IPM',
           data_list = list(),
           states = states,
           evict = FALSE) %>%
  define_impl(impl_args) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit))

test_that('print.simple_di_det returns correctly', {

  print_str <- print(sim_di_det)

  expect_s3_class(print_str,
                  'simple_di_det_ipm')


})

test_that('plot.simple_di_det returns correctly', {

  plot_str <- plot(sim_di_det,
                   sub_kernels = TRUE,
                   exponent = 0.05)

  expect_s3_class(plot_str,
                  'simple_di_det_ipm')

})


# simple_di_stoch_kern methods --------------
context('simple_di_stoch_kern print and plot methods')


hier_levels <- list(yr = 1:5)

# additional usr_funs to be passed into make_ipm()

inv_logit <- function(sv, int, slope) {
  return(
    1/(1 + exp(-(int + slope * sv)))
  )
}

inv_logit_r <- function(sv, int, slope, r_eff) {
  return(
    1/(1 + exp(-(int + slope * sv + r_eff)))
  )
}

pois_r <- function(sv, int, slope, r_eff) {
  return(
    exp(
      int + slope * sv + r_eff
    )
  )
}

# Define some fixed parameters
data_list = list(
  s_int     = 1.03,
  s_slope   = 2.2,
  g_int     = 8,
  g_slope   = 0.92,
  sd_g      = 0.9,
  f_r_int   = 0.09,
  f_r_slope = 0.05,
  f_s_int   = 0.1,
  f_s_slope = 0.005,
  mu_fd     = 9,
  sd_fd     = 2
)

# Now, simulate some random intercepts for growth, survival, and offspring production

g_r_int   <- rnorm(5, 0, 0.3)
s_r_int   <- rnorm(5, 0, 0.7)
f_s_r_int <- rnorm(5, 0, 0.2)

nms <- paste("r_", 1:5, sep = "")

names(g_r_int)   <- paste('g_', nms, sep = "")
names(s_r_int)   <- paste('s_', nms, sep = "")
names(f_s_r_int) <- paste('f_s_', nms, sep = "")

g_params   <- list2(!!! g_r_int)
s_params   <- list2(!!! s_r_int)
f_s_params <- list2(!!! f_s_r_int)

params     <- splice(data_list, g_params, s_params, f_s_params)

# The "model_class" argument is now changed to reflect a different model type

sim_di_stoch_kern <- init_ipm('simple_di_stoch_kern') %>%
  define_kernel(
    name             = 'P_yr',
    formula          = s_g_mult(s_yr, g_yr) ,
    family           = "CC",
    s_yr             = inv_logit_r(ht_1, s_int, s_slope, s_r_yr) *
      (1 - inv_logit(ht_1, f_r_int, f_r_slope)),
    g_yr             = dnorm(ht_2, mu_g_yr, sd_g),
    mu_g_yr          = g_int + g_slope * ht_1 + g_r_yr,
    data_list        = params,
    states           = list(c('ht')),
    has_hier_effs    = TRUE,
    levels_hier_effs = hier_levels,
    evict = TRUE,
    evict_fun        = truncated_distributions('norm', 'g_yr')
  ) %>%
  define_kernel(
    name             = "F_yr",
    formula          = f_r * f_s_yr * f_d,
    family           = "CC",
    f_r              = inv_logit(ht_1, f_r_int, f_r_slope),
    f_s_yr           = pois_r(ht_1, f_s_int, f_s_slope, f_s_r_yr),
    f_d              = dnorm(ht_2, mu_fd, sd_fd),
    data_list        = params,
    states           = list(c('ht')),
    has_hier_effs    = TRUE,
    levels_hier_effs = hier_levels,
    evict            = TRUE,
    evict_fun        = truncated_distributions('norm', 'f_d')
  ) %>%
  define_k(
    name             = 'K_yr',
    K_yr             = P_yr + F_yr,
    family           = "IPM",
    data_list        = params,
    states           = list(c("ht")),
    has_hier_effs    = TRUE,
    levels_hier_effs = hier_levels
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("K_yr", "P_yr", "F_yr"),
      int_rule     = rep("midpoint", 3),
      dom_start    = rep("ht", 3),
      dom_end      = rep("ht", 3)
    )
  ) %>%
  define_domains(ht = c(0.2, 40, 100)) %>%
  make_ipm(usr_funs = list(inv_logit   = inv_logit,
                           inv_logit_r = inv_logit_r,
                           pois_r      = pois_r))


test_that('print.simple_di_stoch_kern returns correctly', {

  print_str <- print(sim_di_stoch_kern,
                     lambda_type = 'deterministic',
                     compute_type = 'eigen')

  expect_s3_class(print_str,
                  'simple_di_stoch_kern_ipm')


})

test_that('plot.simple_di_stoch_kern returns correctly', {

  plot_str <- plot(sim_di_stoch_kern,
                   sub_kernels = TRUE,
                   exponent = 0.05)

  expect_s3_class(plot_str,
                  'simple_di_stoch_kern_ipm')

})

# simple_di_stoch_param methods -----

context('simple_di_stoch_param print and plot methods')

data_list <- list(s_slope = 0.2,
                  g_slope = 0.99,
                  g_sd = 0.2,
                  f_r_slope = 0.003,
                  f_s_slope = 0.01,
                  f_d_mu = 2,
                  f_d_sd = 0.75)

r_means <- c(s_int_yr = 0.8,
             g_int_yr = 0.1,
             f_r_int_yr = 0.3,
             f_s_int_yr = 0.01)

# These are likely nonsense values - using for testing purposes only!
r_sigma <- runif(16)
r_sigma <- matrix(r_sigma, nrow = 4)

r_sigma <- r_sigma %*% t(r_sigma) # make symmetrical


inv_logit <- function(int, slope, sv1) {
  1/(1 + exp(-(int + slope * sv1)))
}

mvt_wrapper <- function(r_means, r_sigma, nms) {
  out <- rmvnorm(1, r_means, r_sigma) %>%
    as.list()

  names(out) <- nms
  return(out)
}

init_pop_vec <- runif(100, 0, 10)

sim_di_stoch_param <- init_ipm('simple_di_stoch_param') %>%
  define_kernel(
    'P',
    formula = s_g_mult(s, g),
    family = 'CC',
    g_mu = env_params$g_int_yr + g_slope * surf_area_1,
    s = inv_logit(env_params$s_int_yr, s_slope, surf_area_1),
    g = dnorm(surf_area_2, g_mu, g_sd),
    data_list = data_list,
    states = list(c('surf_area')),
    has_hier_effs = FALSE,
    evict = TRUE,
    evict_fun = truncated_distributions('norm', 'g')
  ) %>%
  define_kernel(
    'F',
    formula = f_r * f_s * f_d,
    family = 'CC',
    f_r = inv_logit(env_params$f_r_int_yr, f_r_slope, surf_area_1),
    f_s = exp(env_params$f_s_int_yr + f_s_slope * surf_area_1),
    f_d = dnorm(surf_area_2, f_d_mu, f_d_sd),
    data_list = data_list,
    states = list(c('surf_area')),
    has_hier_effs = FALSE,
    evict = TRUE,
    evict_fun = truncated_distributions('norm', 'f_d')
  ) %>%
  define_k(
    'K',
    K = P + F,
    n_surf_area_t_1 = right_mult(K, n_surf_area_t),
    family = 'IPM',
    data_list = data_list,
    states = list(c('surf_area')),
    has_hier_effs = FALSE,
    evict = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c('P', "F", "K"),
      int_rule = rep('midpoint', 3),
      dom_start = rep('surf_area',3),
      dom_end = rep('surf_area', 3)
    )
  ) %>%
  define_domains(surf_area = c(0, 10, 100)) %>%
  define_env_state(
    env_params = mvt_wrapper(r_means, r_sigma, nms = c('s_int_yr',
                                                       'g_int_yr',
                                                       'f_r_int_yr',
                                                       'f_s_int_yr')),
    data_list = list(
      r_means = r_means,
      r_sigma = r_sigma
    )
  ) %>%
  define_pop_state(
    pop_vectors = list(n_surf_area_t = init_pop_vec),
  ) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit,
                           mvt_wrapper = mvt_wrapper),
           iterate = TRUE,
           iterations = 3)


test_that('print.simple_di_stoch_param returns correctly', {

  print_str <- print(sim_di_stoch_param)

  expect_s3_class(print_str,
                  'simple_di_stoch_param_ipm')


})

test_that('plot.simple_di_stoch_param returns correctly', {

  plot_str <- plot(sim_di_stoch_param,
                   sub_kernels = TRUE,
                   exponent = 0.05)

  expect_s3_class(plot_str,
                  'simple_di_stoch_param_ipm')

})
