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

impl_args <- make_impl_args_list(c('P', 'F'),
                                 int_rule = rep('midpoint', 2),
                                 state_start = rep('dbh', 2),
                                 state_end = rep('dbh', 2))

inv_logit <- function(int, slope, sv) {
  return(1/(1 + exp(-(int + slope * sv))))
}

states <- c('dbh', 'dbh')

sim_di_det_1 <- init_ipm(sim_gen    = "simple",
                         di_dd      = "di",
                         det_stoch  = "det") %>%
  define_kernel("P",
                formula = s * g,
                family = "CC",
                s = inv_logit(s_int, s_slope, dbh_1),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * dbh_1,
                data_list = data_list,
                states = states,
                evict_cor = TRUE,
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
                evict_cor = TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'f_d')
  ) %>%
  define_impl(impl_args) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  define_pop_state(n_dbh = runif(100)) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit),
           normalize_pop_size = FALSE)


sim_di_det_2 <- init_ipm(sim_gen    = "simple",
                         di_dd      = "di",
                         det_stoch  = "det") %>%
  define_kernel("P",
                formula = s * g,
                family = "CC",
                s = inv_logit(s_int, s_slope, dbh_1),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * dbh_1,
                data_list = data_list,
                states = states,
                evict_cor = TRUE,
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
                evict_cor = TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'f_d')
  ) %>%
  define_impl(impl_args) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  define_pop_state(n_dbh = runif(100)) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit),
           iterate = TRUE,
           iterations = 100,
           normalize_pop_size = FALSE)


sim_di_det_3 <- init_ipm(sim_gen    = "simple",
                         di_dd      = "di",
                         det_stoch  = "det") %>%
  define_kernel("P",
                formula = s * g,
                family = "CC",
                s = inv_logit(s_int, s_slope, dbh_1),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * dbh_1,
                data_list = data_list,
                states = states,
                evict_cor = TRUE,
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
                evict_cor = TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'f_d')
  ) %>%
  define_impl(impl_args) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  define_pop_state(n_dbh = runif(100)) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit),
           iterate = TRUE,
           iterations = 5,
           normalize_pop_size = FALSE)

test_that('print.simple_di_det returns correctly', {

  p <- capture_output(print_str <- print(sim_di_det_1))

  print_msg <- capture_output_lines(print(sim_di_det_2,
                                          type_lambda = 'all'))

  expect_s3_class(print_str,
                  'simple_di_det_ipm')

  expect_true(
    any(
      grepl('A simple, density independent, deterministic IPM',
            print_msg)
    )
  )

  wrn_msg    <- catch_cnd(print(sim_di_det_3))$message

  expect_true(
    grepl(
      'sim_di_det_3 has has not converged to asymptotic dynamics!',
      wrn_msg
    )
  )
})

test_that('plot.simple_di_det returns correctly', {

  plot_str <- plot(sim_di_det_1,
                   # sub_kernels = TRUE,
                   exponent = 0.025,
                   do_contour = TRUE)

  expect_s3_class(plot_str,
                  'simple_di_det_ipm')

})


# simple_di_stoch_kern methods --------------


par_set_indices <- list(yr = 1:5)

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

params     <- c(data_list, g_params, s_params, f_s_params)


sim_di_stoch_kern <- init_ipm(sim_gen    = "simple",
                              di_dd      = "di",
                              det_stoch  = "stoch",
                              "kern") %>%
  define_kernel(
    name             = 'P_yr',
    formula          = s_yr * g_yr ,
    family           = "CC",
    s_yr             = inv_logit_r(ht_1, s_int, s_slope, s_r_yr) *
      (1 - inv_logit(ht_1, f_r_int, f_r_slope)),
    g_yr             = dnorm(ht_2, mu_g_yr, sd_g),
    mu_g_yr          = g_int + g_slope * ht_1 + g_r_yr,
    data_list        = params,
    states           = list(c('ht')),
    uses_par_sets    = TRUE,
    par_set_indices = par_set_indices,
    evict_cor        = TRUE,
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
    uses_par_sets    = TRUE,
    par_set_indices = par_set_indices,
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions('norm', 'f_d')
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P_yr", "F_yr"),
      int_rule     = rep("midpoint", 2),
      state_start    = rep("ht", 2),
      state_end      = rep("ht", 2)
    )
  ) %>%
  define_domains(ht = c(0.2, 40, 100)) %>%
  define_pop_state(n_ht = runif(100)) %>%
  make_ipm(usr_funs = list(inv_logit   = inv_logit,
                           inv_logit_r = inv_logit_r,
                           pois_r      = pois_r),
           normalize_pop_size = FALSE,
           iterate = TRUE,
           kernel_seq = sample(1:5, size = 50, replace = TRUE))


test_that('print.simple_di_stoch_kern returns correctly', {

  p <- capture_output(print_str <- print(sim_di_stoch_kern))

  expect_s3_class(print_str,
                  'simple_di_stoch_kern_ipm')

  expect_true(
    any(
      grepl('A simple, density independent, stochastic, kernel-resampled IPM',
            p)
    )
  )

})

test_that('plot.simple_di_stoch_kern returns correctly', {

  plot_str <- plot(sim_di_stoch_kern,
                   sub_kernels = TRUE,
                   exponent = 0.05)

  expect_s3_class(plot_str,
                  'simple_di_stoch_kern_ipm')

})

# simple_di_stoch_param methods -----


data_list <- list(s_slope = 0.2,
                  g_slope = 0.99,
                  g_sd = 0.2,
                  f_r_slope = 0.003,
                  f_s_slope = 0.01,
                  f_d_mu = 2,
                  f_d_sd = 0.75)

r_means <- c(s_int_yr = 0.8,
             g_int_yr = 0.8,
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

sim_di_stoch_param <- init_ipm(sim_gen    = "simple",
                               di_dd      = "di",
                               det_stoch  = "stoch",
                               "param") %>%
  define_kernel(
    'P',
    formula = s * g,
    family = 'CC',
    g_mu = g_int_yr + g_slope * surf_area_1,
    s = inv_logit(s_int_yr, s_slope, surf_area_1),
    g = dnorm(surf_area_2, g_mu, g_sd),
    data_list = data_list,
    states = list(c('surf_area')),
    uses_par_sets = FALSE,
    evict_cor = TRUE,
    evict_fun = truncated_distributions('norm', 'g')
  ) %>%
  define_kernel(
    'F',
    formula = f_r * f_s * f_d,
    family = 'CC',
    f_r = inv_logit(f_r_int_yr, f_r_slope, surf_area_1),
    f_s = exp(f_s_int_yr + f_s_slope * surf_area_1),
    f_d = dnorm(surf_area_2, f_d_mu, f_d_sd),
    data_list = data_list,
    states = list(c('surf_area')),
    uses_par_sets = FALSE,
    evict_cor = TRUE,
    evict_fun = truncated_distributions('norm', 'f_d')
  )  %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c('P', "F"),
      int_rule = rep('midpoint', 2),
      state_start = rep('surf_area', 2),
      state_end = rep('surf_area', 2)
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
    pop_vectors = list(n_surf_area = init_pop_vec),
  ) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit,
                           mvt_wrapper = mvt_wrapper),
           iterate = TRUE,
           iterations = 3,
           normalize_pop_size = FALSE)


test_that('print.simple_di_stoch_param returns correctly', {

  p <- capture_output(print_str <- print(sim_di_stoch_param,
                                         type_lambda = "all"))

  expect_s3_class(print_str,
                  'simple_di_stoch_param_ipm')

  print_msg <- capture_output_lines(
    print(sim_di_stoch_param)
  )

  expect_true(
    any(
      grepl('A simple, density independent, stochastic, parameter-resampled IPM',
            print_msg)
    )
  )


})

test_that('plot.simple_di_stoch_param returns correctly', {

  # There is a bug in do_legend = TRUE. Return later...

  plot_str <- plot(sim_di_stoch_param,
                   exponent = 0.05,
                   do_legend = FALSE)

  expect_s3_class(plot_str,
                  'simple_di_stoch_param_ipm')

})


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

states <- list(c('ht', "b"))

inv_logit <- function(int, slope, sv) {
  1/(1 + exp(-(int + slope * sv)))
}

inv_logit_2 <- function(int, slope, slope_2, sv) {
  1/(1 + exp(-(int + slope * sv + slope_2 * sv ^ 2)))
}

gen_di_det_1 <- init_ipm(sim_gen    = "general",
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
      state_end      = c('ht', "b", "b", 'ht')
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

# Test unconverged warnings

gen_di_det_2 <- init_ipm(sim_gen    = "general",
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
      state_end      = c('ht', "b", "b", 'ht')
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
  make_ipm(iterations = 10,
           usr_funs = list(inv_logit   = inv_logit,
                           inv_logit_2 = inv_logit_2),
           normalize_pop_size = FALSE)

test_that('print.general_di_det returns correctly', {

  p <- capture_output(print_str <- print(gen_di_det_1))

  print_msg <- capture_output_lines(print(gen_di_det_1,
                                          type_lambda = 'all'))

  expect_s3_class(print_str,
                  'general_di_det_ipm')

  expect_true(
    any(
      grepl('A general, density independent, deterministic IPM',
            print_msg)
    )
  )

  wrn_msg    <- catch_cnd(print(gen_di_det_2))$message

  expect_true(
    grepl(
      'gen_di_det_2 has has not converged to asymptotic dynamics!',
      wrn_msg
    )
  )
})



flatten_to_depth <- ipmr:::.flatten_to_depth

par_sets <- list(
  site = c(
    'whitetop',
    'mt_rogers',
    'roan',
    'big_bald',
    'bob_bald',
    'oak_knob'
  )
)

inv_logit <- function(int, slope, sv) {

  1/(1 + exp(-(int + slope * sv)))

}

pois      <- function(int, slope, sv) {

  exp(int + slope * sv)

}

lin_prob <- function(int, slope, sigma, sv1, sv2, L, U) {

  mus <- int + slope * sv1

  ev <- pnorm(U, mus, sigma) - pnorm(L, mus, sigma)

  dnorm(sv2, mus, sigma) / ev

}

# pass to usr_funs in make_ipm
vr_funs <- list(
  pois      = pois,
  inv_logit = inv_logit
)

# non-reproductive transitions + vrs
nr_data_list <- list(
  nr_s_z_int = c(-1.1032, -1.0789, -1.3436,
                 -1.3188, -0.9635, -1.3243), # survival
  nr_s_z_b   = c(1.6641, 0.7961, 1.5443,
                 1.4561, 1.2375, 1.2168),    # survival
  nr_d_z_int = rep(-1.009, 6),               # dormancy prob
  nr_d_z_b   = rep(-1.3180, 6),              # dormancy prob
  nr_f_z_int = c(-7.3702, -4.3854, -9.0733,
                 -6.7287, -9.0416, -7.4223), # flowering prob
  nr_f_z_b   = c(4.0272, 2.2895, 4.7366,
                 3.1670, 4.7410, 4.0834),    # flowering prob
  nr_nr_int  = c(0.4160, 0.7812, 0.5697,
                 0.5426, 0.3744, 0.6548),    # growth w/ stasis in NR stage
  nr_nr_b    = c(0.6254, 0.3742, 0.5325,
                 0.6442, 0.6596, 0.4809),    # growth w/ stasis in NR stage
  nr_nr_sd   = rep(0.3707, 6),               # growth w/ stasis in NR stage
  nr_ra_int  = c(0.8695, 0.7554, 0.9789,
                 1.1824, 0.8918, 1.0027),    # growth w/ move to RA stage
  nr_ra_b    = rep(0.5929, 6),               # growth w/ move to RA stage
  nr_ra_sd   = rep(0.4656, 6)                # growth w/ move to RA stage

)

# reproductive output parameters. not sure if this applies only to reproductive ones,
# or if they include it for non-reproductive ones that transition to flowering
# as well!

fec_data_list <- list(
  f_s_int   = c(4.1798, 3.6093, 4.0945,
                3.8689, 3.4776, 2.4253), # flower/seed number intercept
  f_s_slope = c(0.9103, 0.5606, 1.0898,
                0.9101, 1.2352, 1.6022), # flower/seed number slope
  tau_int   = c(3.2428, 0.4263, 2.6831,
                1.5475, 1.3831, 2.5715), # Survival through summer for flowering plants
  tau_b     = rep(0, 6)
)

# reproductive vital rates

ra_data_list <- list(
  ra_s_z_int = c(-0.6624, -0.9061, -0.0038,
                 -1.3844, -0.1084, -0.6477), # survival
  ra_s_z_b   = rep(0, 6),                    # survival
  ra_d_z_int = rep(-1.6094, 6),              # dormancy prob
  ra_d_z_b   = rep(0, 6),                    # dormancy prob
  ra_n_z_int = rep(0.9463, 6),               # transition to non-reproductive status
  ra_n_z_b   = rep(0, 6),                    # transition to non-reproductive status
  ra_n_z_sd  = rep(0.4017, 6)                # transition to non-reproductive status
)

# dormany - nra + recruitment parameters (e.g. discrete -> continuous parameters)

dc_data_list <- list(
  dc_nr_int = rep(0.9687, 6),            # dormant -> NR transition prob
  dc_nr_sd  = rep(0.4846, 6),            # dormant -> NR transition sd
  sdl_z_int = c(0.4992, 0.5678, 0.6379,
                0.6066, 0.5147, 0.5872), # recruit size mu
  sdl_z_b   = rep(0, 6),
  sdl_z_sd  = rep(0.1631, 6),            # recruit size sd
  sdl_es_r  = c(0.1233, 0.1217, 0.0817,
                0.1800, 0.1517, 0.0533)  # establishment probabilities
)

# Extracted domains from figure 4F in main text - these are approximations
# but should be pretty damn close to the actual limits of integration.

domains <- list(sqrt_area = c(0.63 * 0.9, 3.87 * 1.1, 50),
                ln_leaf_l = c(0.26 * 0.9, 2.70 * 1.1, 50))


# helper to get data_list's into a format that ipmr can actually use

rename_data_list <- function(data_list, nms) {

  for(j in seq_along(data_list)) {

    # Isolate entries in big list
    temp <- data_list[[j]]

    for(i in seq_along(temp)) {

      # loop over parameter entries - there should be 6 for every one
      # e.g. 1 for every population

      nm_i <- names(temp)[i]

      # make new names

      x <- temp[[i]]
      names(x) <- paste(nm_i, '_', nms, sep = "")

      temp[[i]] <- as.list(x)
    }

    # replace unnamed with freshly named!

    data_list[[j]] <- temp

  }

  # finito

  return(data_list)

}

all_params <- list(nr_data_list,
                   fec_data_list,
                   ra_data_list,
                   dc_data_list)

full_data_list <- rename_data_list(data_list = all_params,
                                   nms       = unlist(par_sets)) %>%
  flatten_to_depth(1)

init_pop_vec <- list(
  ln_leaf_l = runif(50),
  sqrt_area = runif(50)
)

gen_di_stoch_kern_1 <- init_ipm(sim_gen    = "general",
                                di_dd      = "di",
                                det_stoch  = "stoch",
                                "kern") %>%
  define_kernel(
    name             = 'k_xx_site',
    family           = "CC",
    formula          = sig_n_site          *
      (1 - mu_n_site)   *
      (1 - beta_n_site) *
      gamma_nn_site     *
      d_ln_leaf_l,

    sig_n_site        = inv_logit(nr_s_z_int_site, nr_s_z_b_site, ln_leaf_l_1),
    gamma_nn_site     = dnorm(ln_leaf_l_2, nr_nr_mu_site, nr_nr_sd_site),
    nr_nr_mu_site     = nr_nr_int_site + nr_nr_b_site * ln_leaf_l_1,
    mu_n_site         = inv_logit(nr_d_z_int_site, nr_d_z_b_site, ln_leaf_l_1),
    beta_n_site       = inv_logit(nr_f_z_int_site, nr_f_z_b_site, ln_leaf_l_1),

    data_list        = full_data_list,
    states           = list(c('ln_leaf_l')),
    uses_par_sets    = TRUE,
    par_set_indices = par_sets,
    evict_cor        =  TRUE,
    evict_fun        = truncated_distributions('norm',
                                               'gamma_nn_site')
  ) %>%
  define_kernel(
    name             = 'k_zx_site',
    family           = 'CC',

    formula          = (
      phi_site        *
        nu_site         *
        gamma_sr_site   +
        sig_r_site      *
        (1 - mu_r_site) *
        gamma_nr_site
    )                *
      d_sqrt_area,

    phi_site      = pois(f_s_int_site, f_s_slope_site, sqrt_area_1),
    nu_site       = sdl_es_r_site,
    gamma_sr_site = dnorm(ln_leaf_l_2, sdl_z_int_site, sdl_z_sd_site),
    sig_r_site    = inv_logit(ra_s_z_int_site, ra_s_z_b_site, sqrt_area_1),
    mu_r_site     = inv_logit(ra_d_z_int_site, ra_d_z_b_site, sqrt_area_1),
    gamma_nr_site = dnorm(ln_leaf_l_2, mu_ra_nr_site, ra_n_z_sd_site),
    mu_ra_nr_site = ra_n_z_int_site + ra_n_z_b_site * sqrt_area_1,

    data_list    = full_data_list,
    states       = list(c('sqrt_area', 'ln_leaf_l')),
    uses_par_sets = TRUE,
    par_set_indices = par_sets,
    evict_cor = TRUE,
    evict_fun = truncated_distributions(c('norm', 'norm'),
                                        c('gamma_nr_site',
                                          'gamma_sr_site'))
  ) %>%
  define_kernel(
    name = 'k_dx_site',
    family = 'DC',

    formula = gamma_nd_site * d_ln_leaf_l,

    gamma_nd_site = dnorm(ln_leaf_l_2, dc_nr_int_site, dc_nr_sd_site),

    data_list = full_data_list,
    states = list(c('ln_leaf_l', "d")),
    uses_par_sets = TRUE,
    par_set_indices = par_sets,
    evict_cor = TRUE,
    evict_fun = truncated_distributions('norm',
                                        'gamma_nd_site')
  ) %>%
  define_kernel(
    name             = 'k_xz_site',
    family           = 'CC',
    formula          = sig_n_site         *
      (1 - mu_n_site)    *
      beta_n_site       *
      gamma_rn_site     *
      tau               *
      d_ln_leaf_l,

    sig_n_site        = inv_logit(nr_s_z_int_site, nr_s_z_b_site, ln_leaf_l_1),
    mu_n_site         = inv_logit(nr_d_z_int_site, nr_d_z_b_site, ln_leaf_l_1),
    beta_n_site       = inv_logit(nr_f_z_int_site, nr_f_z_b_site, ln_leaf_l_1),
    gamma_rn_site     = dnorm(sqrt_area_2, mu_nr_ra_site, nr_ra_sd_site),
    mu_nr_ra_site     = nr_ra_int_site + nr_ra_b_site * ln_leaf_l_1,
    tau               = inv_logit(tau_int_site, tau_b_site, ln_leaf_l_1),

    data_list        = full_data_list,
    states           = list(c('sqrt_area', 'ln_leaf_l')),
    uses_par_sets    = TRUE,
    par_set_indices = par_sets,
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions('norm',
                                               'gamma_rn_site')

  ) %>%
  define_kernel(
    name             = 'k_xd_site',
    family           = 'CD',
    formula          = sig_n_site * mu_n_site * d_ln_leaf_l,
    sig_n_site       = inv_logit(nr_s_z_int_site, nr_s_z_b_site, ln_leaf_l_1),
    mu_n_site        = inv_logit(nr_d_z_int_site, nr_d_z_b_site, ln_leaf_l_1),
    data_list        = full_data_list,
    states           = list(c('ln_leaf_l', "d")),
    uses_par_sets    = TRUE,
    par_set_indices = par_sets,
    evict_cor        = FALSE
  ) %>%
  define_kernel(
    name             = 'k_zd_site',
    family           = 'CD',
    formula          = sig_r_site * mu_r_site * d_sqrt_area,
    sig_r_site        = inv_logit(ra_s_z_int_site, ra_s_z_b_site, sqrt_area_1),
    mu_r_site         = inv_logit(ra_d_z_int_site, ra_d_z_b_site, sqrt_area_1),
    data_list        = full_data_list,
    states           = list(c('sqrt_area', "d")),
    uses_par_sets    = TRUE,
    par_set_indices = par_sets,
    evict_cor        = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c(paste('k_',
                             c('xx',
                               'zx', 'dx', 'xz', 'xd', 'zd'),
                             '_site',
                             sep = "")),
      int_rule     = rep('midpoint', 6),
      state_start    = c('ln_leaf_l',
                       'sqrt_area',
                       "d",
                       'ln_leaf_l',
                       'ln_leaf_l',
                       'sqrt_area'),
      state_end      = c('ln_leaf_l',
                       'ln_leaf_l',
                       'ln_leaf_l',
                       'sqrt_area',
                       "d",
                       "d")
    )
  ) %>%
  define_domains(
    sqrt_area = c(0.63 * 0.9, 3.87 * 1.1, 50),
    ln_leaf_l = c(0.26 * 0.9, 2.70 * 1.1, 50)
  ) %>%
  define_pop_state(
    pop_vectors = list(
      n_ln_leaf_l = init_pop_vec$ln_leaf_l,
      n_sqrt_area = init_pop_vec$sqrt_area,
      n_d         = 10
    )
  ) %>%
  make_ipm(
    return_all = TRUE,
    usr_funs   = list(
      inv_logit = inv_logit,
      pois      = pois
    ),
    iterations = 100,
    normalize_pop_size = FALSE
  )

test_that('print.general_di_stoch_kern works', {

  p <- capture_output(print_str <- print(gen_di_stoch_kern_1))

  print_msg <- capture_output_lines(print(gen_di_stoch_kern_1,
                                          type_lambda = 'stochastic'))

  expect_s3_class(print_str,
                  'general_di_stoch_kern_ipm')

  expect_true(
    any(
      grepl('A general, density independent, stochastic, kernel-resampled IPM',
            print_msg)
    )
  )

  expect_length(print_msg, 4L)

})


to_mu_sd <- function(x) {
  lapply(x,
         function(y){
           list(
             mu    = mean(y),
             sigma = sd(y)
           )
         })
}

flatten_to_depth <- ipmr:::.flatten_to_depth

# Hacky estimates of parameter means and variances. Those without any variance
# will be converted to fixed point estimates, those with variance will get
# converted into a multivariate joint distribution with simulated values for the
# variance-covariance matrix

nr_data_list <- list(
  nr_s_z_int = c(-1.1032, -1.0789, -1.3436,
                 -1.3188, -0.9635, -1.3243), # survival
  nr_s_z_b   = c(1.6641, 0.7961, 1.5443,
                 1.4561, 1.2375, 1.2168),    # survival
  nr_d_z_int = rep(-1.009, 6),               # dormancy prob
  nr_d_z_b   = rep(-1.3180, 6),              # dormancy prob
  nr_f_z_int = c(-7.3702, -4.3854, -9.0733,
                 -6.7287, -9.0416, -7.4223), # flowering prob
  nr_f_z_b   = c(4.0272, 2.2895, 4.7366,
                 3.1670, 4.7410, 4.0834),    # flowering prob
  nr_nr_int  = c(0.4160, 0.7812, 0.5697,
                 0.5426, 0.3744, 0.6548),    # growth w/ stasis in NR stage
  nr_nr_b    = c(0.6254, 0.3742, 0.5325,
                 0.6442, 0.6596, 0.4809),    # growth w/ stasis in NR stage
  nr_nr_sd   = rep(0.3707, 6),               # growth w/ stasis in NR stage
  nr_ra_int  = c(0.8695, 0.7554, 0.9789,
                 1.1824, 0.8918, 1.0027),    # growth w/ move to RA stage
  nr_ra_b    = rep(0.5929, 6),               # growth w/ move to RA stage
  nr_ra_sd   = rep(0.4656, 6)                # growth w/ move to RA stage

) %>%
  to_mu_sd()

fec_data_list <- list(
  f_s_int   = c(4.1798, 3.6093, 4.0945,
                3.8689, 3.4776, 2.4253), # flower/seed number intercept
  f_s_slope = c(0.9103, 0.5606, 1.0898,
                0.9101, 1.2352, 1.6022), # flower/seed number slope
  tau_int   = c(3.2428, 0.4263, 2.6831,
                1.5475, 1.3831, 2.5715), # Survival through summer for flowering plants
  tau_b     = rep(0, 6)
) %>%
  to_mu_sd()

# reproductive vital rates

ra_data_list <- list(
  ra_s_z_int = c(-0.6624, -0.9061, -0.0038,
                 -1.3844, -0.1084, -0.6477), # survival
  ra_s_z_b   = rep(0, 6),                    # survival
  ra_d_z_int = rep(-1.6094, 6),              # dormancy prob
  ra_d_z_b   = rep(0, 6),                    # dormancy prob
  ra_n_z_int = rep(0.9463, 6),               # transition to non-reproductive status
  ra_n_z_b   = rep(0, 6),                    # transition to non-reproductive status
  ra_n_z_sd  = rep(0.4017, 6)                # transition to non-reproductive status
) %>%
  to_mu_sd()

# dormany - nra + recruitment parameters (e.g. discrete -> continuous parameters)

dc_data_list <- list(
  dc_nr_int = rep(0.9687, 6),            # dormant -> NR transition prob
  dc_nr_sd  = rep(0.4846, 6),            # dormant -> NR transition sd
  sdl_z_int = c(0.4992, 0.5678, 0.6379,
                0.6066, 0.5147, 0.5872), # recruit size mu
  sdl_z_b   = rep(0, 6),
  sdl_z_sd  = rep(0.1631, 6),            # recruit size sd
  sdl_es_r  = c(0.1233, 0.1217, 0.0817,
                0.1800, 0.1517, 0.0533)  # establishment probabilities
) %>%
  to_mu_sd()


# Compile full list, then split out fixed parameters into list form (passed to
# define_kernel) and random parameters into vector form (passed to mvtnorm)

all_params <- c(nr_data_list,
                fec_data_list,
                ra_data_list,
                dc_data_list)

ind_rand <- map_lgl(all_params,
                    .f = function(.x) {
                      .x$sigma != 0
                    })

ind_fixed <- ! ind_rand

# Next, split fixed and random variables. Flatten fixed ones out since
# we are only using their point estimates

fixed_params <- all_params[ind_fixed]

fixed_params <- lapply(fixed_params,
                       function(x) x$mu)

rando_means  <- vapply(all_params[ind_rand],
                       function(x) x$mu + 1,
                       numeric(1L))

rando_sigs  <- vapply(all_params[ind_rand],
                      function(x) x$sigma,
                      numeric(1L))

# Check names really quickly
stopifnot(names(rando_means) == names(rando_sigs))

# Var-covar matrix - this makes it symmetrical. perhaps not the "correct" way to
# generate it, but this is a unit test for something else entirely.

rando_sigmas <- rando_sigs %*% t(rando_sigs)

rando_names  <- names(rando_means)

# Define a wrapper for mvtnorm that generates a list of parameter draws and
# names them correctly. This gets passed to define_env_state

mvt_wrapper <- function(r_means, r_sigmas, nms) {

  out <- rmvnorm(1, r_means, r_sigmas) %>%
    as.list()

  names(out) <- nms

  return(out)

}

inv_logit <- function(int, slope, sv) {

  1 / (1 + exp(-(int + slope * sv)))
}

pois <- function(int, slope, sv) {

  exp(int + slope * sv)

}

init_pop_vec <- list(
  ln_leaf_l = runif(50),
  sqrt_area = runif(50)
)

# Implement the ipmr version. We'll use the env_seq slot from the output
# to iterate a hand coded version with identical parameter estimates to ensure
# we're creating the same model. No seeds set, so this test changes slightly
# every time

gen_di_stoch_param_1 <- init_ipm(sim_gen    = "general",
                                 di_dd      = "di",
                                 det_stoch  = "stoch",
                                 "param") %>%
  define_kernel(
    name             = 'k_xx',
    family           = "CC",
    formula          = sig_n        *
      (1 - mu_n)   *
      (1 - beta_n) *
      gamma_nn     *
      d_ln_leaf_l,

    sig_n        = inv_logit(nr_s_z_int, nr_s_z_b, ln_leaf_l_1),
    gamma_nn     = dnorm(ln_leaf_l_2, nr_nr_mu, nr_nr_sd),
    nr_nr_mu     = nr_nr_int + nr_nr_b * ln_leaf_l_1,
    mu_n         = inv_logit(nr_d_z_int, nr_d_z_b, ln_leaf_l_1),
    beta_n       = inv_logit(nr_f_z_int, nr_f_z_b, ln_leaf_l_1),

    data_list        = fixed_params,
    states           = list(c('ln_leaf_l')),
    uses_par_sets    = FALSE,
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions('norm',
                                               'gamma_nn')
  ) %>%
  define_kernel(
    name             = 'k_zx',
    family           = 'CC',

    formula          = (
      phi        *
        nu         *
        gamma_sr   +
        sig_r      *
        (1 - mu_r) *
        gamma_nr
    )            *
      d_sqrt_area,

    phi           = pois(f_s_int, f_s_slope, sqrt_area_1),
    nu            = sdl_es_r,
    gamma_sr      = dnorm(ln_leaf_l_2, sdl_z_int, sdl_z_sd),
    sig_r         = inv_logit(ra_s_z_int, ra_s_z_b, sqrt_area_1),
    mu_r          = inv_logit(ra_d_z_int, ra_d_z_b, sqrt_area_1),
    gamma_nr      = dnorm(ln_leaf_l_2, mu_ra_nr, ra_n_z_sd),
    mu_ra_nr      = ra_n_z_int + ra_n_z_b * sqrt_area_1,

    data_list     = fixed_params,
    states        = list(c('sqrt_area', 'ln_leaf_l')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(c('norm', 'norm'),
                                            c('gamma_nr',
                                              'gamma_sr'))
  ) %>%
  define_kernel(
    name      = 'k_dx',
    family    = 'DC',

    formula   = gamma_nd * d_ln_leaf_l,

    gamma_nd  = dnorm(ln_leaf_l_2, dc_nr_int, dc_nr_sd),

    data_list = fixed_params,
    states    = list(c('ln_leaf_l', "d")),

    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm',
                                            'gamma_nd')
  ) %>%
  define_kernel(
    name             = 'k_xz',
    family           = 'CC',
    formula          = sig_n         *
      (1 - mu_n)    *
      beta_n        *
      gamma_rn      *
      tau           *
      d_ln_leaf_l,

    sig_n            = inv_logit(nr_s_z_int, nr_s_z_b, ln_leaf_l_1),
    mu_n             = inv_logit(nr_d_z_int, nr_d_z_b, ln_leaf_l_1),
    beta_n           = inv_logit(nr_f_z_int, nr_f_z_b, ln_leaf_l_1),
    gamma_rn         = dnorm(sqrt_area_2, mu_nr_ra, nr_ra_sd),
    mu_nr_ra         = nr_ra_int + nr_ra_b * ln_leaf_l_1,
    tau              = inv_logit(tau_int, tau_b, ln_leaf_l_1),

    data_list        = fixed_params,
    states           = list(c('sqrt_area', 'ln_leaf_l')),
    uses_par_sets    = FALSE,
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions('norm',
                                               'gamma_rn')

  ) %>%
  define_kernel(
    name             = 'k_xd',
    family           = 'CD',
    formula          = sig_n * mu_n * d_ln_leaf_l,
    sig_n            = inv_logit(nr_s_z_int, nr_s_z_b, ln_leaf_l_1),
    mu_n             = inv_logit(nr_d_z_int, nr_d_z_b, ln_leaf_l_1),
    data_list        = fixed_params,
    states           = list(c('ln_leaf_l', "d")),
    uses_par_sets    = FALSE,
    evict_cor        = FALSE
  ) %>%
  define_kernel(
    name             = 'k_zd',
    family           = 'CD',
    formula          = sig_r * mu_r * d_sqrt_area,
    sig_r            = inv_logit(ra_s_z_int, ra_s_z_b, sqrt_area_1),
    mu_r             = inv_logit(ra_d_z_int, ra_d_z_b, sqrt_area_1),
    data_list        = fixed_params,
    states           = list(c('sqrt_area', "d")),
    uses_par_sets    = FALSE,
    evict_cor        = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c(paste('k_',
                             c('xx',
                               'zx', 'dx', 'xz', 'xd', 'zd'),
                             sep = "")),
      int_rule     = rep('midpoint', 6),
      state_start    = c('ln_leaf_l',
                       'sqrt_area',
                       "d",
                       'ln_leaf_l',
                       'ln_leaf_l',
                       'sqrt_area'),
      state_end      = c('ln_leaf_l',
                       'ln_leaf_l',
                       'ln_leaf_l',
                       'sqrt_area',
                       "d",
                       "d")
    )
  ) %>%
  define_domains(
    sqrt_area = c(0.63 * 0.9, 3.87 * 1.1, 50),
    ln_leaf_l = c(0.26 * 0.9, 2.70 * 1.1, 50)
  ) %>%
  define_pop_state(
    pop_vectors = list(
      n_ln_leaf_l = init_pop_vec$ln_leaf_l,
      n_sqrt_area = init_pop_vec$sqrt_area,
      n_d         = 10
    )
  ) %>%
  define_env_state(
    env_params = mvt_wrapper(rando_means,
                             rand_sigs,
                             nms = rando_names),
    data_list  = list(
      rando_means = rando_means,
      rand_sigs   = rando_sigmas,
      rando_names = rando_names
    )
  ) %>%
  make_ipm(
    return_all = TRUE,
    usr_funs   = list(
      inv_logit   = inv_logit,
      pois        = pois,
      mvt_wrapper = mvt_wrapper
    ),
    iterations = 100,
    normalize_pop_size = FALSE
  )

test_that('print.general_di_stoch_param works', {

  p <- capture_output(print_str <- print(gen_di_stoch_param_1))

  print_msg <- capture_output_lines(print(gen_di_stoch_param_1,
                                          type_lambda = 'stochastic'))

  expect_s3_class(print_str,
                  'general_di_stoch_param_ipm')

  expect_true(
    any(
      grepl('A general, density independent, stochastic, parameter-resampled IPM',
            print_msg)
    )
  )

  expect_length(print_msg, 5L)


})


test_that('print.proto_ipm is working correctly', {

  test_proto <- sim_di_det_1$proto_ipm

  print_msg <- capture_output_lines(
    print(test_proto)
  )

  expect_true(
    any(
      grepl(
        'A simple, density independent, deterministic proto_ipm with 2 kernels defined',
        print_msg
      )
    )
  )

})


test_that("`%^%` is working correctly", {

  # Test from expm::`%^%` helper file
  test_mat <- cbind(1, 2 * diag(3)[,-1])

  target   <- matrix(c(1, 0, 0,
                       3, 4, 0,
                       3, 0, 4),
                     nrow = 3, byrow = TRUE)

  test_exp <- test_mat %^% 2L

  expect_equal(target, test_exp)

  test_mat_2 <- runif(25) %>%
    matrix(nrow = 5, ncol = 5)

  target_2   <- test_mat_2 %*% test_mat_2 %*% test_mat_2 %*% test_mat_2

  test_exp_2 <- test_mat_2 %^% 4L

  expect_equal(test_exp_2, target_2)

  err_mat <- rnorm(6) %>%
    matrix(ncol = 3, nrow = 2)

  errs <- catch_cnd(err_mat %^% 2L)
  expect_true(grepl("not implemented for non-square matrices", errs))

  expect_equal(mat_power(target, 3L), target %^% 3L)

})


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

impl_args <- make_impl_args_list(c('P', 'F'),
                                 int_rule = rep('midpoint', 2),
                                 state_start = rep('dbh', 2),
                                 state_end = rep('dbh', 2))

inv_logit <- function(int, slope, sv) {
  return(1/(1 + exp(-(int + slope * sv))))
}


states <- c('dbh', 'dbh')

sim_di_det_1 <- init_ipm(sim_gen    = "simple",
                         di_dd      = "di",
                         det_stoch  = "det") %>%
  define_kernel("P",
                formula = s * g,
                family = "CC",
                s = inv_logit(s_int, s_slope, dbh_1),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * dbh_1,
                data_list = data_list,
                states = states,
                evict_cor = TRUE,
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
                evict_cor = TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'f_d')
  )  %>%
  define_impl(impl_args) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit),
           normalize_pop_size = FALSE,
           iterate = FALSE)

hand_k <- sim_di_det_1$sub_kernels$P + sim_di_det_1$sub_kernels$F

target <- Re(eigen(hand_k)$vectors[ , 1])
target <- target / sum(target)

test_that('right_ev.simple_di_det iterates un-iterated model correctly', {

  ipmr_w <- right_ev(sim_di_det_1)[[1]]

  expect_equal(target, ipmr_w, tolerance = 1e-10)

})

init_pop <- runif(100)

sim_di_det_2 <- init_ipm(sim_gen    = "simple",
                         di_dd      = "di",
                         det_stoch  = "det") %>%
  define_kernel("P",
                formula = s * g,
                family = "CC",
                s = inv_logit(s_int, s_slope, dbh_1),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * dbh_1,
                data_list = data_list,
                states = states,
                evict_cor = TRUE,
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
                evict_cor = TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'f_d')
  )  %>%
  define_impl(impl_args) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  define_pop_state(n_dbh = init_pop) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit),
           iterate = TRUE,
           iterations = 200,
           normalize_pop_size = FALSE)

sim_di_det_3 <- init_ipm(sim_gen    = "simple",
                         di_dd      = "di",
                         det_stoch  = "det") %>%
  define_kernel("P",
                formula = s * g,
                family = "CC",
                s = inv_logit(s_int, s_slope, dbh_1),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * dbh_1,
                data_list = data_list,
                states = states,
                evict_cor = TRUE,
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
                evict_cor = TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'f_d')
  ) %>%
  define_impl(impl_args) %>%
  define_pop_state(n_dbh = init_pop) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit),
           iterate = TRUE,
           iterations = 5,
           normalize_pop_size = FALSE)

test_that('right_ev can handle iterated models', {

  ipmr_w <- right_ev(sim_di_det_2)[[1]]

  expect_equal(target, ipmr_w, tolerance = 1e-10)

  ipmr_w <- right_ev(sim_di_det_3, iterations = 200)[[1]]

  expect_equal(target, ipmr_w, tolerance = 1e-10)

})

test_that('right_ev messages are being produced correctly', {

  test_message <- catch_cnd(
    right_ev(sim_di_det_1)
  )

  expect_true(
    grepl(
      "Generating a population vector using runif\\(\\)",
      test_message$message
    )
  )

  # If the model is already iterated, there will not be any messages!
  # skip onward to sim_di_det_3

  test_message <- catch_cnd(
    right_ev(sim_di_det_3)
  )

  expect_true(
    grepl(
      "did not converge to asymptotic dynamics after 6 iterations",
      test_message$message
    )
  )

  # Now, test out failure to converge warnings from right_ev

  test_message <- capture_warning(
    right_ev(sim_di_det_3,
             iterations = 8)
  )

  expect_true(
    grepl(
      'did not converge after 14 iterations\\. Returning NA',
      test_message$message
    )
  )

  test_message <- capture_warning(
    right_ev(sim_di_det_1,
             iterations = 3)
  )

  expect_true(
    grepl(
      'did not converge after 3 iterations\\. Returning NA',
      test_message$message
    )
  )

})


target_v <- Re(eigen(t(hand_k))$vectors[ , 1])
target_v <- target_v / sum(target_v)

test_that('left_ev.simple_di_det_ipm is working correctly', {

  ipmr_v <- left_ev(sim_di_det_1)[[1]]
  expect_equal(target_v, ipmr_v, tolerance = 1e-10)

})

test_that('left_ev.simple_di_det_ipm can re-iterate models', {

  ipmr_v <- left_ev(sim_di_det_3)[[1]]
  expect_equal(target_v, ipmr_v, tolerance = 1e-10)

})

test_that('warnings are produced correctly', {

  msg <- catch_cnd(left_ev(sim_di_det_3,
                           iterations = 3))

  expect_true(
    grepl(
      'did not converge after 3 iterations\\. Returning NA',
      msg$message
    )
  )


  msg <- capture_warning(
    left_ev(sim_di_det_1,
            iterations = 2
    )
  )

  expect_true(
    grepl(
      'did not converge after 2 iterations\\. Returning NA',
      msg$message
    )
  )


})


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

gen_di_det_1 <- init_ipm(sim_gen    = "general",
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
      state_end      = c('ht', "b", "b", 'ht')
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

mega_k <- rbind(
  cbind(gen_di_det_1$sub_kernels$stay_discrete,
        gen_di_det_1$sub_kernels$go_discrete),
  cbind(gen_di_det_1$sub_kernels$leave_discrete,
        gen_di_det_1$sub_kernels$P)
)

target <- Re(eigen(mega_k)$vectors[ , 1])
target <- target / sum(target)

test_that('right_ev.general_di_det returns the same as mega-matrix methods', {

  ipmr_temp <- right_ev(gen_di_det_1)
  ipmr_w    <- c(ipmr_temp$b_w,
                 ipmr_temp$ht_w)

  expect_equal(target, ipmr_w, tolerance = 1e-9)

})

gen_di_det_2 <- init_ipm(sim_gen    = "general",
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
      state_end      = c('ht', "b", "b", 'ht')
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
  make_ipm(iterations = 5,
           usr_funs = list(inv_logit   = inv_logit,
                           inv_logit_2 = inv_logit_2),
           normalize_pop_size = FALSE)



test_that('right_ev.general_di_det can re-iterate models', {

  test_reiterate <- right_ev(gen_di_det_2,
                             iterations = 100)

  ipmr_w <- c(test_reiterate$b_w,
              test_reiterate$ht_w)

  expect_equal(target, ipmr_w, tolerance = 1e-9)

})


test_that('right_ev.general_di_det returns NAs and warnings properly', {

  test_message <- capture_warning(
    right_ev(gen_di_det_2,
             iterations = 3)
  )

  expect_true(
    grepl(
      'did not converge after 3 iterations\\. Returning NA',
      test_message$message
    )
  )

})

target_v <- Re(eigen(t(mega_k))$vectors[ , 1])
target_v <- target_v / sum(target_v)

test_that('left_ev.gen_di_det returns the same as mega-matrix models', {

  ipmr_v <- left_ev(gen_di_det_1,
                    iterations = 200)

  ipmr_v <- c(ipmr_v$b_v, ipmr_v$ht_v)

  expect_equal(target_v, ipmr_v, tolerance = 1e-9)

})

test_that('left_ev.general_di_det can re-iterate models', {

  # THis really shouldn't make a difference for the left eigenvector,
  # We don't actually use make_ipm() to re-iterate things

  ipmr_v <- left_ev(gen_di_det_2, iterations = 200)

  ipmr_v <- c(ipmr_v$b_v, ipmr_v$ht_v)

  expect_equal(target_v, ipmr_v, tolerance = 1e-9)


})


test_that('left_ev.general_di_det returns warnings properly', {

  test_message <- capture_warning(
    left_ev(gen_di_det_2,
            iterations = 3)
  )

  expect_true(
    grepl(
      'did not converge after 3 iterations\\. Returning NA',
      test_message$message
    )
  )

})

test_that("left/right_ev work with par set deterministic models", {

  fixed_list <- list(s_int = 2.2,
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

  s_r_ints <- rnorm(5, -2.2, sd = 1) %>%
    as.list()

  g_r_ints <- rnorm(5, - 0.2, 1) %>%
    as.list

  names(s_r_ints) <- paste('s_r_', 1:5, sep = "")
  names(g_r_ints) <- paste('g_r_', 1:5, sep = "")

  data_list <- c(fixed_list, s_r_ints, g_r_ints)

  states <- list(c("dbh"))

  impl_args <-  make_impl_args_list(c('P_yr', 'F'),
                                    int_rule = rep('midpoint', 2),
                                    state_start = rep('dbh', 2),
                                    state_end = rep('dbh', 2))

  inv_logit <- function(int, slope, sv) {

    1 / (1 + exp(-(int + slope * sv)))
  }

  par_set_mod <- init_ipm(sim_gen    = "simple",
                          di_dd      = "di",
                          det_stoch  = "det") %>%
    define_kernel("P_yr",
                  formula = s_yr * g_yr,
                  family = "CC",
                  s_yr = inv_logit(s_int + s_r_yr, s_slope, dbh_1),
                  g_yr = dnorm(dbh_2, mu_g_yr, sd_g),
                  mu_g_yr = g_int + g_r_yr + g_slope * dbh_1,
                  data_list = data_list,
                  states = states,
                  uses_par_sets = TRUE,
                  par_set_indices = list(yr = 1:5),
                  evict_cor = TRUE,
                  evict_fun = truncated_distributions('norm',
                                                      'g_yr')
    ) %>%
    define_kernel('F',
                  formula = f_r * f_s * f_d,
                  family = 'CC',
                  f_r = inv_logit(f_r_int, f_r_slope, dbh_1),
                  f_s = exp(f_s_int + f_s_slope * dbh_1),
                  f_d = dnorm(dbh_2, mu_fd, sd_fd),
                  data_list = data_list,
                  states = states,
                  evict_cor = TRUE,
                  evict_fun = truncated_distributions('norm',
                                                      'f_d')
    ) %>%
    define_impl(impl_args) %>%
    define_pop_state(
      n_dbh_yr = runif(100)
    ) %>%
    define_domains(dbh = c(0, 50, 100)) %>%
    make_ipm(usr_funs = list(inv_logit = inv_logit),
             iterate = TRUE,
             iterations = 100,
             normalize_pop_size = FALSE)

  w <- right_ev(par_set_mod)
  v <- left_ev(par_set_mod)

  nms_w <- paste0("dbh_", 1:5, "_w")
  nms_v <- paste0("dbh_", 1:5, "_v")

  expect_equal(names(w), nms_w)
  expect_equal(names(v), nms_v)
  expect_s3_class(w, "ipmr_w")
  expect_s3_class(v, "ipmr_v")
  expect_type(w[[1]], "double")
  expect_type(v[[1]], "double")


  # General models
  set.seed(2312)
  init_pop_vec <- runif(500)
  init_b <- 20
  g_ints <- rnorm(3) %>% as.list()
  s_ints <- rnorm(3) %>% as.list()

  names(g_ints) <- paste("g_int_", 1:3, sep = "")
  names(s_ints) <- paste("s_int_", 1:3, sep = "")
  data_list_control <- list(
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
  data_list_control <- c(data_list_control, g_ints, s_ints)

  L <- 1.02
  U <- 624
  n <- 500

  inv_logit <- function(int, slope, z) {
    lin_p <- int + slope * z
    return(1/(1+exp(-lin_p)))
  }

  inv_logit_2 <- function(int, slope_1, slope_2, z) {
    lin_p <- int + slope_1 * z + slope_2 * z^2
    return(1/(1+exp(-lin_p)))
  }

  states <- list(c("ht", "b"))

  ipmr_control <- init_ipm(sim_gen    = "general",
                           di_dd      = "di",
                           det_stoch  = "det") %>%
    define_kernel(
      name             = "P_site",
      formula          = s_site * g_site * d_ht,
      family           = "CC",
      g_site           = dnorm(ht_2, g_mu_site, g_sd),
      g_mu_site        = g_int + g_int_site + g_slope * ht_1,
      s_site           = inv_logit_2(s_int + s_int_site, s_slope, s_slope_2, ht_1),
      data_list        = data_list_control,
      states           = states,
      uses_par_sets    = TRUE,
      par_set_indices = list(site = 1:3),
      evict_cor        = TRUE,
      evict_fun        = truncated_distributions('norm',
                                                 'g_site')
    ) %>%
    define_kernel(
      name          = "go_discrete",
      formula       = f_r * f_s * g_i,
      family        = 'CD',
      f_r           = inv_logit(f_r_int, f_r_slope, ht_1),
      f_s           = exp(f_s_int + f_s_slope * ht_1),
      data_list     = data_list_control,
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
      data_list     = data_list_control,
      states        = states,
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'f_d')
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P_site", "go_discrete", "stay_discrete", "leave_discrete"),
        int_rule     = c(rep("midpoint", 4)),
        state_start    = c('ht', "ht", "b", "b"),
        state_end      = c('ht', "b", "b", 'ht')
      )
    ) %>%
    define_domains(
      ht = c(1.02, 624, 500)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_ht_site = init_pop_vec,
        n_b_site  = init_b
      )
    ) %>%
    make_ipm(iterations = 200,
             usr_funs = list(inv_logit   = inv_logit,
                             inv_logit_2 = inv_logit_2),
             return_all_envs = FALSE,
             normalize_pop_size = TRUE)

  w <- right_ev(ipmr_control)
  v <- left_ev(ipmr_control, iterations = 200)

  nms_w <- paste0(c("ht_", "b_"), c(1:3, 1:3), "_w")
  nms_v <- paste0(c("ht_", "b_"), c(1:3, 1:3), "_v")

  expect_true(all(names(w) %in% nms_w))
  expect_true(all(names(v) %in% nms_v))
  expect_s3_class(w, "ipmr_w")
  expect_s3_class(v, "ipmr_v")
  expect_type(w[[1]], "double")
  expect_type(v[[1]], "double")

})

test_that("age_size models left/right_ev", {

  m.par.true <- c(## survival
    surv.int  = -1.70e+1,
    surv.z    =  6.68e+0,
    surv.a    = -3.34e-1,
    ## growth
    grow.int  =  1.27e+0,
    grow.z    =  6.12e-1,
    grow.a    = -7.24e-3,
    grow.sd   =  7.87e-2,
    ## reproduce or not
    repr.int  = -7.88e+0,
    repr.z    =  3.11e+0,
    repr.a    = -7.80e-2,
    ## recruit or not
    recr.int  =  1.11e+0,
    recr.a    =  1.84e-1,
    ## recruit size
    rcsz.int  =  3.62e-1,
    rcsz.z    =  7.09e-1,
    rcsz.sd   =  1.59e-1)

  param_list <- as.list(m.par.true) %>%
    setNames(gsub(pattern = "\\.", replacement = "_", x = names(.)))

  inv_logit <- function(x) {

    return( 1 / (1 + exp(-x)) )
  }

  f_fun <- function(age, s_age, pb_age, pr_age, recr) {

    if(age == 0) return(0)

    s_age * pb_age * pr_age * recr * 0.5

  }

  a_s_ipm <- init_ipm(sim_gen    = "general",
                      di_dd      = "di",
                      det_stoch  = "det",
                      uses_age = TRUE) %>%
    define_kernel(
      name          = "P_age",
      family        = "CC",
      formula       = s_age * g_age * d_wt,
      s_age         = inv_logit(surv_int + surv_z * wt_1 + surv_a * age),
      g_age         = dnorm(wt_2, mu_g_age, grow_sd),
      mu_g_age      = grow_int + grow_z * wt_1 + grow_a * age,
      data_list     = param_list,
      states        = list(c("wt_age")),
      uses_par_sets = FALSE,
      age_indices   = list(age = c(0:20), max_age = 21),
      evict_cor     = FALSE
    ) %>%
    define_kernel(
      name          = "F_age",
      family        = "CC",
      formula       = f_fun(age, s_age, pb_age, pr_age, recr) * d_wt,
      s_age         = inv_logit(surv_int + surv_z * wt_1 + surv_a * age),
      pb_age        = inv_logit(repr_int + repr_z * wt_1 + repr_a * age),
      pr_age        = inv_logit(recr_int + recr_a * age),
      recr          = dnorm(wt_2, rcsz_mu, rcsz_sd),
      rcsz_mu       = rcsz_int + rcsz_z * wt_1,
      data_list     = param_list,
      states        = list(c("wt")),
      uses_par_sets = FALSE,
      age_indices   = list(age = c(0:20), max_age = 21),
      evict_cor     = FALSE
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P_age", "F_age"),
        int_rule     = rep("midpoint", 2),
        state_start    = c("wt_age", "wt_age"),
        state_end      = c("wt_age", "wt_0")
      )
    ) %>%
    define_domains(
      wt = c(1.6, 3.7, 100)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_wt_age = runif(100))
    ) %>%
    make_ipm(
      usr_funs = list(inv_logit = plogis,
                      f_fun     = f_fun),
      iterate  = TRUE,
      iterations = 100,
      return_all_envs = TRUE
    )

  mega <- format_mega_kernel(a_s_ipm,
                             name_ps = "P",
                             f_forms = "F")$mega_matrix

  target_w <- target_v <- runif(2200)

  for(i in 1:100) {

    target_w <- mega %*% target_w
    target_v <- t(mega) %*% target_v

  }

  ipmr_w   <- right_ev(a_s_ipm) %>%
    unlist() %>%
    unname()
  ipmr_v   <- left_ev(a_s_ipm) %>%
    unlist() %>%
    unname()

  target_w <- (target_w / sum(target_w)) %>% as.vector
  target_v <- (target_v / sum(target_v)) %>% as.vector

  expect_equal(target_w, ipmr_w)
  expect_equal(target_v, ipmr_v)

  vr_texts <- vital_rate_exprs(a_s_ipm$proto_ipm)

  cat_out  <- capture_output_lines(print(vr_texts))[1:5]

  test_ind <- grepl("<age>", cat_out)

  expect_true(all(test_ind))

})


test_that('fill_0s is working as it should', {

  mega_mat <- rlang::quo(
    c(
      x, 0, y,
      z, a, 0,
      b, c, d
    )
  )

  sub_kernels <- list(
    x = matrix(rnorm(9, 10), 3, 3),
    y = matrix(rnorm(12), 3, 4),
    z = matrix(rnorm(15, 25), 5, 3),
    a = matrix(rnorm(25, 50), 5, 5),
    b = matrix(runif(18), 6, 3),
    c = matrix(rnorm(30, -20), 6, 5),
    d = matrix(rnorm(24, -70), 6, 4)
  )

  ex_ipm <- list(
    iterators = NA_real_,
    sub_kernels = sub_kernels,
    proto_ipm = data.frame(uses_par_sets = FALSE)
  )

  class(ex_ipm) <- c("simple_di_det_ipm", "list")

  mega_k   <- .make_mega_mat(ex_ipm, mega_mat)
  mega_dim <- dim(mega_k)

  expect_equal(14, mega_dim[1])
  expect_equal(12, mega_dim[2])

  expect_equal(mega_k[1:3, 1:3], sub_kernels$x)
  expect_equal(mega_k[4:8, 4:8], sub_kernels$a)
  expect_equal(matrix(rep(0, 15), 3, 5), mega_k[1:3, 4:8])
  expect_equal(matrix(rep(0, 20), 5, 4), mega_k[4:8, 9:12])
  expect_equal(mega_k[9:14, 9:12], sub_kernels$d)

})

test_that('format_mega_mat works as advertised', {

  test_mat <- format_mega_kernel(gen_di_det_2,
                                 mega_mat = c(
                                   stay_discrete, go_discrete,
                                   leave_discrete, P
                                 ))

  expect_equal(unclass(gen_di_det_2$sub_kernels$P),
               test_mat[[1]][2:501, 2:501])

  disc <- unclass(gen_di_det_2$sub_kernels$go_discrete) %>%
    as.vector()

  expect_equal(disc,
               test_mat[[1]][1, 2:501])


})


test_that("format_mega_kernel can handle par_sets", {

  Ps <- lapply(1:5,
               function(x) {
                 matrix(runif(2500),
                        nrow = 50,
                        ncol = 50)
               })

  names(Ps) <- paste("P", 1:5, sep = "_")

  go_discs  <- lapply(1:5,
                      function(x) {
                        matrix(rpois(50, 3), nrow = 1)
                      })

  names(go_discs) <- paste('go_disc', 1:5, sep = "_")

  leave_discs <- lapply(1:5,
                        function(x) {
                          matrix(runif(50), ncol = 1)
                        })

  names(leave_discs) <- paste('leave_disc', 1:5, sep = "_")

  mat_expr <- rlang::expr(c(0, go_disc_site, leave_disc_site, P_site))

  actuals <- lapply(1:5,
                    function(x, leaves, gos, Ps) {

                      tr <- cbind(0, go_discs[[x]])
                      br <- cbind(leaves[[x]], Ps[[x]])

                      return(rbind(tr, br))

                    },
                    leaves = leave_discs,
                    gos    = go_discs,
                    Ps     = Ps) %>%
    setNames(paste("mega_matrix_", names(.), sep = ""))

  ex_ipm <- list(
    iterators = NA_integer_,
    sub_kernels = c(
      leave_discs,
      go_discs,
      Ps
    ),
    env_list = list(),
    env_seq  = sample(1:5, 100, replace = TRUE),
    pop_state = list(),
    proto_ipm = data.frame(uses_par_sets = TRUE,
                           par_set_indices = I(list(site = 1:5)))
  )

  class(ex_ipm)           <- c("general_di_det_ipm", 'list')
  class(ex_ipm$proto_ipm) <- c("general_di_det", "proto_ipm", "data.frame")

  ipmr_megas <- format_mega_kernel(ex_ipm,
                                   c(0, go_disc_site,
                                     leave_disc_site, P_site))

  test_vals  <- vapply(1:5,
                       function(x, actual, pkg_val) {

                         isTRUE(
                           all.equal(actual[[x]],
                                     pkg_val[[x]],
                                     tolerance = 1e-10)
                         )

                       },
                       logical(1L),
                       actual = actuals,
                       pkg_val = ipmr_megas)

  expect_true(all(test_vals))

})

test_that("format_mega_kernel works w/ multiple par_sets", {


  nms <- expand.grid(list(site = c("A", "B"),
                          yr   = 1:3),
                     stringsAsFactors = FALSE)

  out_nms <- character(6L)

  for(i in seq_len(6)) {
    out_nms[i] <- paste(nms[i, ], collapse = "_")
  }

  Ps <- lapply(1:6,
               function(x) {
                 matrix(runif(2500),
                        nrow = 50,
                        ncol = 50)
               })

  names(Ps) <- paste("P_", out_nms, sep = "")
  go_discs  <- lapply(1:6,
                      function(x) {
                        matrix(rpois(50, 3), nrow = 1)
                      })

  names(go_discs) <- paste("go_disc_", out_nms, sep = "")

  leave_discs <- lapply(1:6,
                        function(x) {
                          matrix(runif(50), ncol = 1)
                        })
  names(leave_discs) <- paste("leave_disc_", out_nms, sep = "")


  mat_expr <- rlang::expr(c(0, go_disc_site_yr, leave_disc_site_yr, P_site_yr))

  actuals <- lapply(1:6,
                    function(x, leaves, gos, Ps) {

                      tr <- cbind(0, go_discs[[x]])
                      br <- cbind(leaves[[x]], Ps[[x]])

                      return(rbind(tr, br))

                    },
                    leaves = leave_discs,
                    gos    = go_discs,
                    Ps     = Ps)

  names(actuals) <- paste("mega_matrix", out_nms, sep = "_")

  ex_ipm <- list(
    iterators = NA_integer_,
    sub_kernels = c(
      leave_discs,
      go_discs,
      Ps
    ),
    env_list = list(),
    env_seq  = sample(out_nms, 100, replace = TRUE),
    pop_state = list(),
    proto_ipm = data.frame(uses_par_sets = TRUE,
                           par_set_indices = I(list(site = c("A","B"),
                                                     yr   = 1:3)))
  )

  class(ex_ipm)           <- c("general_di_det_ipm", 'list')
  class(ex_ipm$proto_ipm) <- c("general_di_det", "proto_ipm", "data.frame")

  ipmr_megas <- format_mega_kernel(ex_ipm,
                                   c(0, go_disc_site_yr,
                                     leave_disc_site_yr, P_site_yr))

  test_vals  <- vapply(1:6,
                       function(x, actual, pkg_val) {

                         isTRUE(
                           all.equal(actual[[x]],
                                     pkg_val[[x]],
                                     tolerance = 1e-10)
                         )

                       },
                       logical(1L),
                       actual = actuals,
                       pkg_val = ipmr_megas)

  expect_true(all(test_vals))
  expect_equal(names(ipmr_megas), names(actuals))

})


test_that('format_mega_kernel can handle character vectors', {

  x <- matrix(rnorm(25), ncol = 5)
  y <- matrix(rnorm(25), ncol = 5)
  z <- matrix(runif(25), ncol = 5)

  target <- rbind(
    cbind(x, y),
    cbind(z, matrix(0, nrow = 5, ncol = 5))
  )

  ipm <- list(sub_kernels = list(x = x, y = y, z = z))


  sym_mat <- format_mega_kernel(ipm, mega_mat = c(x , y , z , 0))

  # First, make sure we've fooled the ipm argument
  expect_equal(sym_mat[[1]], target)

  # now, try with a character vector
  vec_text <- 'c(x , y, z, 0)'

  test_mat <- format_mega_kernel(ipm, vec_text)

  expect_equal(sym_mat, test_mat)

  test_vec <- c("x", "y", "z", 0)

  test_mat <- format_mega_kernel(ipm, test_vec)

  expect_equal(sym_mat, test_mat)

})

test_that("format_mega_kernel can handles identity matrices as advertized", {

  x <- matrix(rnorm(25), ncol = 5)
  y <- matrix(rnorm(25), ncol = 5)
  z <- matrix(runif(25), ncol = 5)
  I <- diag(5)

  target <- rbind(
    cbind(x, y),
    cbind(z, I)
  )

  ipm <- list(sub_kernels = list(x = x, y = y, z = z))

  test_mat <- format_mega_kernel(ipm, c(x, y, z, I))

  expect_equal(test_mat[[1]], target)

})

test_that("We can fill 0s and Identity matrices", {

  x <- matrix(rnorm(25), ncol = 5)
  y <- matrix(rnorm(25), ncol = 5)
  z <- matrix(runif(25), ncol = 5)
  I <- diag(5)

  target <- rbind(
    cbind(x, matrix(0, ncol = 5, nrow = 5)),
    cbind(z, I)
  )

  ipm <- list(sub_kernels = list(x = x, z = z))

  test_mat <- format_mega_kernel(ipm, c(x, 0, z, I))

  expect_equal(test_mat[[1]], target)



})


test_that("format_mega_kernel works w/ drop_levels", {


  nms <- expand.grid(list(site = c("A", "B"),
                          yr   = 1:3),
                     stringsAsFactors = FALSE)

  out_nms <- character(6L)

  for(i in seq_len(6)) {
    out_nms[i] <- paste(nms[i, ], collapse = "_")
  }

  to_drop <- c("A_2", "B_1")

  out_nms <- out_nms[!out_nms %in% to_drop]

  Ps <- lapply(1:4,
               function(x) {
                 matrix(runif(2500),
                        nrow = 50,
                        ncol = 50)
               })

  names(Ps) <- paste("P_", out_nms, sep = "")
  go_discs  <- lapply(1:4,
                      function(x) {
                        matrix(rpois(50, 3), nrow = 1)
                      })

  names(go_discs) <- paste("go_disc_", out_nms, sep = "")

  leave_discs <- lapply(1:4,
                        function(x) {
                          matrix(runif(50), ncol = 1)
                        })
  names(leave_discs) <- paste("leave_disc_", out_nms, sep = "")


  mat_expr <- rlang::expr(c(0, go_disc_site_yr, leave_disc_site_yr, P_site_yr))

  actuals <- lapply(1:4,
                    function(x, leaves, gos, Ps) {

                      tr <- cbind(0, go_discs[[x]])
                      br <- cbind(leaves[[x]], Ps[[x]])

                      return(rbind(tr, br))

                    },
                    leaves = leave_discs,
                    gos    = go_discs,
                    Ps     = Ps)

  names(actuals) <- paste("mega_matrix", out_nms, sep = "_")

  ex_ipm <- list(
    iterators = NA_integer_,
    sub_kernels = c(
      leave_discs,
      go_discs,
      Ps
    ),
    env_list = list(),
    env_seq  = sample(out_nms, 100, replace = TRUE),
    pop_state = list(),
    proto_ipm = data.frame(uses_par_sets = TRUE,
                           par_set_indices = I(list(site = c("A","B"),
                                                     yr   = 1:3,
                                                     drop_levels = to_drop)))
  )

  class(ex_ipm)           <- c("general_di_det_ipm", 'list')
  class(ex_ipm$proto_ipm) <- c("general_di_det", "proto_ipm", "data.frame")

  ipmr_megas <- format_mega_kernel(ex_ipm,
                                   c(0, go_disc_site_yr,
                                     leave_disc_site_yr, P_site_yr))

  test_vals  <- vapply(1:4,
                       function(x, actual, pkg_val) {

                         isTRUE(
                           all.equal(actual[[x]],
                                     pkg_val[[x]],
                                     tolerance = 1e-10)
                         )

                       },
                       logical(1L),
                       actual = actuals,
                       pkg_val = ipmr_megas)

  expect_true(all(test_vals))
  expect_equal(names(ipmr_megas), names(actuals))

})
