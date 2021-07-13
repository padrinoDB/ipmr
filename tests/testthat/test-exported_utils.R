
test_that("exported utils return expected values w general IPMs", {

  data(gen_di_det_ex)

  prot <- gen_di_det_ex$proto_ipm

  dom <- domains(prot)

  vr_exprs <- vital_rate_exprs(prot)

  forms <- kernel_formulae(prot)

  vr_types <- vapply(vr_exprs,
                     rlang::is_call,
                     logical(1L))

  kern_types <- vapply(forms,
                       function(x)
                         rlang::is_call(x) || rlang::is_bare_atomic(x),
                       logical(1L))

  expect_true(all(vr_types))
  expect_true(all(kern_types))

  dom_types <- vapply(dom, is.numeric, logical(1L))

  expect_true(all(dom_types))

  env_ipm <- make_ipm(prot,
                      return_all_envs = TRUE)

  vr_funs <- vital_rate_funs(env_ipm)

  expect_s3_class(vr_funs$P$g, "CC")
  expect_s3_class(vr_funs$P, "ipmr_vital_rate_funs")

  g <- matrix(env_ipm$env_list$P$g,
              nrow = 200,
              ncol = 200,
              byrow = TRUE)


  g_vr <- unclass(vr_funs$P$g)

  expect_equal(g, g_vr)

})

test_that("exported utils return expected values w simple IPMs", {

  data(sim_di_det_ex)
  data(gen_di_det_ex)

  prot <- sim_di_det_ex$proto_ipm

  dom <- domains(prot)

  vr_exprs <- vital_rate_exprs(prot)

  forms <- kernel_formulae(prot)

  vr_types <- vapply(vr_exprs,
                     rlang::is_call,
                     logical(1L))

  kern_types <- vapply(forms,
                       rlang::is_call,
                       logical(1L))

  expect_true(all(vr_types))
  expect_true(all(kern_types))

  dom_types <- vapply(dom, is.numeric, logical(1L))

  expect_true(all(dom_types))

  vr_ipm <- vital_rate_exprs(sim_di_det_ex)

  expect_identical(vr_exprs, vr_ipm)

  kern_ipm <- kernel_formulae(sim_di_det_ex)

  expect_identical(kern_ipm, forms)

  dom_proto <- domains(prot)
  dom_ipm   <- domains(sim_di_det_ex)

  expect_identical(dom_proto, dom_ipm)


  # pop_state shouldn't actually return the same thing here,
  # so these expectations need to be adjusted
  prot    <- gen_di_det_ex$proto_ipm

  ps_prot <- pop_state(prot)

  prot_tst <- vapply(ps_prot,
                     function(x) x == "Pre-defined population state.",
                     logical(1L))

  expect_true(all(prot_tst))
  expect_s3_class(ps_prot, "ipmr_pop_state")

  ps_ipm  <- pop_state(gen_di_det_ex)

  expect_true(all(is.array(ps_ipm$n_ht), is.array(ps_ipm$n_b)))
  expect_type(ps_ipm, "list")

  env_ipm <- make_ipm(prot,
                      return_all_envs = TRUE)

  vr_funs <- vital_rate_funs(env_ipm)

  expect_s3_class(vr_funs$P$g, "CC")
  expect_s3_class(vr_funs$P, "ipmr_vital_rate_funs")


})


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

states <- list(c("dbh"))

sim_di_det_2 <- init_ipm(sim_gen    = "simple",
                         di_dd      = "di",
                         det_stoch  = "det") %>%
  define_kernel("P",
                formula = s * g,
                family = "CC",
                s = plogis(s_int, s_slope, dbh_1),
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
                f_r = plogis(f_r_int, f_r_slope, dbh_1),
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
  define_pop_state(n_dbh = runif(100)) %>%
  make_ipm(iterate = TRUE,
           iterations = 100,
           normalize_pop_size = FALSE,
           return_all_envs = TRUE,
           return_main_env = TRUE)

test_that("int_mesh works as expected", {

  mesh_ps <- int_mesh(sim_di_det_2)

  expect_equal(names(mesh_ps), c("d_dbh", "dbh_1", "dbh_2"))

  expect_equal(length(mesh_ps$dbh_1), 10000L)
  expect_equal(length(mesh_ps$d_dbh), 1L)


  mesh_ps <- int_mesh(sim_di_det_2, full_mesh = FALSE)
  dummy_seq <- seq(0, 50, length.out = 101)
  dummy_z   <- 0.5 * (dummy_seq[2:101] + dummy_seq[1:100])

  expect_equal(mesh_ps$dbh_1, dummy_z)
  expect_equal(mesh_ps$dbh_2, dummy_z)
  expect_equal(mesh_ps$d_dbh, dummy_z[2] - dummy_z[1])

})

test_that("parameters gets and sets correctly", {

  pars <- parameters(sim_di_det_2$proto_ipm)
  expect_s3_class(pars, "ipmr_parameters")

  pars <- unlist(pars)

  expect_equal(unlist(data_list), pars)

  new_pars <- lapply(data_list,
                     function(x) x + rnorm(1, 0, 0.2))

  new_proto <- sim_di_det_2$proto_ipm
  parameters(new_proto) <- new_pars

  expect_equal(unlist(parameters(new_proto)), unlist(new_pars))

  # Now, test subsetted assignment

  newer_pars <- new_pars[1:5]
  newer_pars <- lapply(newer_pars, function(x) x  + rnorm(1, 0, 0.1))

  parameters(new_proto) <- newer_pars

  test_pars <- parameters(new_proto)
  test_pars <- test_pars[names(new_pars)]

  expect_equal(unlist(test_pars[names(newer_pars)]), unlist(newer_pars))
  expect_equal(unlist(test_pars[6:11]), unlist(new_pars)[6:11])


  data(gen_di_det_ex)

  ipm_params <- parameters(gen_di_det_ex)

  proto_params <- parameters(gen_di_det_ex$proto_ipm)

  expect_s3_class(ipm_params, "ipmr_parameters")
  expect_s3_class(proto_params, "ipmr_parameters")

  expect_identical(ipm_params, proto_params)


})

test_that("collapse_pop_state works", {

  data("gen_di_det_ex")

  gen_di_det_ex <- gen_di_det_ex$proto_ipm %>%
    make_ipm(iterate = TRUE,
             iterations = 100,
             return_main_env = TRUE)

  temp <- collapse_pop_state(gen_di_det_ex,
                             time_step = 100,
                             seedlings = ht <= 10,
                             NRA = ht > 10 & ht <= 200,
                             RA = ht > 200) %>%
    unlist() %>%
    round(digits = 6)

  target <- c(seedlings = 0.104932,
              NRA       = 0.113471,
              RA        = 0.000762)
  expect_equal(temp, target)

})


# Test vr_funs for stoch_param/dd models


library(mvtnorm)
library(rlang)
library(purrr)

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

gen_di_stoch_param <- init_ipm(sim_gen    = "general",
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
    return_all_envs = TRUE,
    usr_funs   = list(
      inv_logit   = inv_logit,
      pois        = pois,
      mvt_wrapper = mvt_wrapper
    ),
    iterations = 100,
    normalize_pop_size = FALSE
  )

test_that("vital_rate_funs returns correctly for stoch_param/dd models", {

  vrs <- vital_rate_funs(gen_di_stoch_param)

  cls <- vapply(vrs, function(x) inherits(x, "ipmr_vital_rate_funs"), logical(1L))

  expect_true(all(cls))

  kern_nms <- gen_di_stoch_param$proto_ipm$kernel_id

  nms <- expand.grid(kern_nms, paste("it", 1:100, sep = "_"))
  nms <- paste(nms[,1], nms[,2], sep = "_")

  expect_true(all(nms %in% names(vrs)))

})

test_that("conv_plot works correctly", {

  x <- conv_plot(gen_di_stoch_param)
  expect_s3_class(x, "general_di_stoch_param_ipm")

  x <- conv_plot(gen_di_stoch_param, log = TRUE)

  expect_s3_class(x, "general_di_stoch_param_ipm")


})


test_that("discretize_pop_vec works correctly", {

  data(iceplant_ex)
  zs <- iceplant_ex$log_size
  pv_ipmr <- discretize_pop_vector(zs,
                                   100,
                                   1.2,
                                   1.2,
                                   normalize = TRUE)

  expect_equal(names(pv_ipmr), c("n_zs", "midpoints_zs"))

  expect_warning(discretize_pop_vector(zs,
                                       100,
                                       1.2,
                                       1.2,
                                       na.rm = FALSE))

  temp <- suppressWarnings(discretize_pop_vector(zs,
                                                 100,
                                                 1.2,
                                                 1.2,
                                                 na.rm = FALSE))

  expect_equal(temp$n_zs, NA_real_)
  expect_equal(temp$midpoints_zs, NA_real_)

})
