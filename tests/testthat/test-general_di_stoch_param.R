
# test general_di_stoch_param models. This will use a simulation with some
# parameters from Aikens & Roach, with the random effects part being stuff I
# tweak.

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

all_params <- purrr::splice(nr_data_list,
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
    normalize_pop_size = TRUE
  )


# Store for comparison!
ipmr_lambdas <- lambda(gen_di_stoch_param,
                       type_lambda = 'all') %>%
  as.vector()

# Now, get the sequence of environmental parameters. These will get inserted into
# the test case to make sure we're recapturing the population dynamics
use_param_seq <- as.data.frame(gen_di_stoch_param$env_seq)
names(use_param_seq) <- gsub('env_params\\.', '', names(use_param_seq))


k_xx <- function(fixed_params, env_params, dom_1, dom_2) {

  d_z <- dom_1[2] - dom_1[1]
  l   <- min(dom_1) - d_z / 2
  u   <- max(dom_1) + d_z / 2

  nr_nr_mu <- env_params$nr_nr_int + env_params$nr_nr_b * dom_1

  sig_n    <- inv_logit(env_params$nr_s_z_int, env_params$nr_s_z_b, dom_1)

  ev <- pnorm(u, nr_nr_mu, fixed_params$nr_nr_sd) -
        pnorm(l, nr_nr_mu, fixed_params$nr_nr_sd)

  gamma_nn <- dnorm(dom_2, nr_nr_mu, fixed_params$nr_nr_sd) /
              ev

  mu_n     <- inv_logit(fixed_params$nr_d_z_int, fixed_params$nr_d_z_b, dom_1)

  beta_n   <- inv_logit(env_params$nr_f_z_int, env_params$nr_f_z_b, dom_1)

  disc_fun <- (
    sig_n        *
    (1 - mu_n)   *
    (1 - beta_n) *
    gamma_nn     *
    d_z
  ) %>%
    matrix(nrow = 50, ncol = 50, byrow = TRUE)

  return(disc_fun)
}

k_zx <- function(fixed_params, env_params, dom_1, dom_2, dom_3, dom_4) {

  d_z  <- dom_3[2] - dom_3[1]
  d_z1 <- dom_1[2] - dom_1[1]
  l   <- min(dom_2) - d_z1 / 2
  u   <- max(dom_2) + d_z1 / 2

  phi           <- pois(env_params$f_s_int, env_params$f_s_slope, dom_3)
  nu            <- env_params$sdl_es_r

  ev_sr         <- pnorm(u, env_params$sdl_z_int, fixed_params$sdl_z_sd) -
                   pnorm(l, env_params$sdl_z_int, fixed_params$sdl_z_sd)

  gamma_sr      <- dnorm(dom_2, env_params$sdl_z_int, fixed_params$sdl_z_sd) /
                   ev_sr
  sig_r         <- inv_logit(env_params$ra_s_z_int, fixed_params$ra_s_z_b, dom_3)
  mu_r          <- inv_logit(fixed_params$ra_d_z_int, fixed_params$ra_d_z_b, dom_3)
  mu_ra_nr      <- fixed_params$ra_n_z_int + fixed_params$ra_n_z_b * dom_3

  ev_nr         <- pnorm(u, mu_ra_nr, fixed_params$ra_n_z_sd) -
                   pnorm(l, mu_ra_nr, fixed_params$ra_n_z_sd)

  gamma_nr      <- dnorm(dom_2, mu_ra_nr, fixed_params$ra_n_z_sd) /
                    ev_nr
  disc_fun <- (
    (
      phi        *
      nu         *
      gamma_sr   +

      sig_r      *
      (1 - mu_r) *
      gamma_nr
    )            *
      d_z

  ) %>%
    matrix(nrow = 50, ncol = 50, byrow = TRUE)

  return(disc_fun)
}

k_dx <- function(fixed_params, dom_1, dom_2) {

  d_z <- dom_1[2] - dom_1[1]
  l   <- min(dom_2) - d_z / 2
  u   <- max(dom_2) + d_z / 2

  ev        <- pnorm(u, fixed_params$dc_nr_int, fixed_params$dc_nr_sd) -
               pnorm(l, fixed_params$dc_nr_int, fixed_params$dc_nr_sd)

  gamma_nd  <-  dnorm(dom_2, fixed_params$dc_nr_int, fixed_params$dc_nr_sd) /
    ev

  ind <- seq(1, length(gamma_nd), by = 50)

  disc_fun <- (gamma_nd * d_z) %>%
    .[ind] %>%
    matrix(nrow = 50, ncol = 1, byrow = TRUE)

  return(disc_fun)

}

k_xz <- function(fixed_params, env_params, dom_1, dom_2, dom_3, dom_4) {

  d_z  <- dom_1[2] - dom_1[1]
  d_z1 <- dom_3[2] - dom_3[1]
  l   <- min(dom_3) - d_z1 / 2
  u   <- max(dom_3) + d_z1 / 2

  sig_n    <- inv_logit(env_params$nr_s_z_int, env_params$nr_s_z_b, dom_1)
  mu_n     <- inv_logit(fixed_params$nr_d_z_int, fixed_params$nr_d_z_b, dom_1)
  beta_n   <- inv_logit(env_params$nr_f_z_int, env_params$nr_f_z_b, dom_1)
  mu_nr_ra <- env_params$nr_ra_int + fixed_params$nr_ra_b * dom_1
  tau      <- inv_logit(env_params$tau_int, fixed_params$tau_b, dom_1)

  ev       <- pnorm(u, mu_nr_ra, fixed_params$nr_ra_sd) -
              pnorm(l, mu_nr_ra, fixed_params$nr_ra_sd)

  gamma_rn <- dnorm(dom_4, mu_nr_ra, fixed_params$nr_ra_sd) /
    ev


  disc_fun <- (
    sig_n         *
    (1 - mu_n)    *
    beta_n        *
    gamma_rn      *
    tau           *
    d_z

  ) %>%
    matrix(nrow = 50, ncol = 50, byrow = TRUE)

  return(disc_fun)
}

k_xd <- function(fixed_params, env_params, dom_1, dom_2) {

  d_z <- dom_1[2] - dom_1[1]
  sig_n <- inv_logit(env_params$nr_s_z_int, env_params$nr_s_z_b, dom_1)
  mu_n  <- inv_logit(fixed_params$nr_d_z_int, fixed_params$nr_d_z_b, dom_1)

  disc_fun <- (sig_n * mu_n * d_z) %>%
    unique() %>%
    matrix(nrow = 1, ncol = 50, byrow = TRUE)

  return(disc_fun)
}

k_zd <- function(fixed_params, env_params, dom_1, dom_2) {

  d_z <- dom_1[2] - dom_1[1]

  sig_r <- inv_logit(env_params$ra_s_z_int, fixed_params$ra_s_z_b, dom_1)
  mu_r  <- inv_logit(fixed_params$ra_d_z_int, fixed_params$ra_d_z_b, dom_1)

  disc_fun <- (sig_r * mu_r * d_z) %>%
    unique() %>%
    matrix(nrow = 1, ncol = 50, byrow = TRUE)

  return(disc_fun)
}

iterate_model <- function(env_params,
                          fixed_params,
                          dom_1,
                          dom_2,
                          dom_3,
                          dom_4,
                          pop_list,
                          v_holder,
                          iteration) {

  k_xx_temp <- k_xx(fixed_params,
                    env_params,
                    dom_1, dom_2)
  k_zx_temp <- k_zx(fixed_params,
                    env_params,
                    dom_1, dom_2,
                    dom_3, dom_4)
  k_dx_temp <- k_dx(fixed_params,
                    dom_1, dom_2)
  k_xz_temp <- k_xz(fixed_params,
                    env_params,
                    dom_1,
                    dom_2,
                    dom_3,
                    dom_4)
  k_xd_temp <- k_xd(fixed_params,
                    env_params,
                    dom_1,
                    dom_2)
  k_zd_temp <- k_zd(fixed_params,
                    env_params,
                    dom_3,
                    dom_4)


  n_ln_leaf_l_t_1 <- k_xx_temp %*% pop_list$n_ln_leaf_l[ , iteration] +
                     k_zx_temp %*% pop_list$n_sqrt_area[ , iteration] +
                     k_dx_temp %*% pop_list$n_d[ , iteration]

  n_sqrt_area_t_1 <- k_xz_temp %*% pop_list$n_ln_leaf_l[ , iteration]

  n_d_t_1         <- k_xd_temp %*% pop_list$n_ln_leaf_l[ , iteration] +
                     k_zd_temp %*% pop_list$n_sqrt_area[ , iteration]

  v_x_t1 <- left_mult(k_xx_temp, v_holder$ln_leaf_l[iteration , ]) +
            left_mult(k_xd_temp, v_holder$d[iteration , ]) +
            left_mult(k_xz_temp, v_holder$sqrt_area[iteration , ])

  v_z_t1 <- left_mult(k_zx_temp, v_holder$ln_leaf_l[iteration , ]) +
            left_mult(k_zd_temp, v_holder$d[iteration , ])

  v_d_t1 <- left_mult(k_dx_temp, v_holder$ln_leaf_l[iteration , ])

  pop_size <- sum(n_ln_leaf_l_t_1, n_sqrt_area_t_1, n_d_t_1)

  pop_list$n_ln_leaf_l[ , (iteration + 1)] <- n_ln_leaf_l_t_1 / pop_size
  pop_list$n_sqrt_area[ , (iteration + 1)] <- n_sqrt_area_t_1 / pop_size
  pop_list$n_d[ , (iteration + 1)]         <- n_d_t_1 / pop_size

  pop_list$lambda[ , iteration] <- sum(n_ln_leaf_l_t_1,
                                       n_sqrt_area_t_1,
                                       n_d_t_1)

  tot_v <- sum(v_x_t1, v_z_t1, v_d_t1)

  v_holder$ln_leaf_l[(iteration + 1), ] <- v_x_t1 / tot_v
  v_holder$sqrt_area[(iteration + 1), ] <- v_z_t1 / tot_v
  v_holder$d[(iteration + 1) , ]        <- v_d_t1 / tot_v

  kerns_temp <- list(k_xx = k_xx_temp,
                     k_zx = k_zx_temp,
                     k_dx = k_dx_temp,
                     k_xz = k_xz_temp,
                     k_xd = k_xd_temp,
                     k_zd = k_zd_temp)

  out_list = list(kernels = kerns_temp,
                  pop_list = pop_list,
                  v_holder = v_holder)

  return(out_list)

}


# Now, we can iterate the model 100 times and check to see if we get identical
# values. We absolutely should be if we're using the same parameter values as
# in the ipmr version.

pop_list <- list(
  n_ln_leaf_l = matrix(NA_real_, nrow = 50, ncol = 101),
  n_sqrt_area = matrix(NA_real_, nrow = 50, ncol = 101),
  n_d         = matrix(NA_real_, nrow = 1,  ncol = 101),
  lambda      = matrix(NA_real_, nrow = 1,  ncol = 100)
)

init_size <- sum(unlist(init_pop_vec), 10)

pop_list$n_ln_leaf_l[ , 1] <- init_pop_vec$ln_leaf_l / init_size
pop_list$n_sqrt_area[ , 1] <- init_pop_vec$sqrt_area / init_size
pop_list$n_d[ , 1]         <- 10 / init_size

v_holder <- lapply(pop_list, t) %>%
  setNames(c("ln_leaf_l", "sqrt_area", "d", "lambda"))

v_holder$lambda <- NULL

sqrt_area_bounds <- seq(0.63 * 0.9,
                        3.87 * 1.1,
                        length.out = 51)

ln_leaf_l_bounds <- seq(0.26 * 0.9,
                        2.70 * 1.1,
                        length.out = 51)

sqrt_area_mids   <- (sqrt_area_bounds[2:51] + sqrt_area_bounds[1:50]) * 0.5
ln_leaf_l_mids   <- (ln_leaf_l_bounds[2:51] + ln_leaf_l_bounds[1:50]) * 0.5

domain_list <- list(expand.grid(sqrt_area_1 = sqrt_area_mids,
                                sqrt_area_2 = sqrt_area_mids),
                    expand.grid(ln_leaf_l_1 = ln_leaf_l_mids,
                                ln_leaf_l_2 = ln_leaf_l_mids)) %>%
  flatten_to_depth(1)



kernel_holder <- list()

for(i in seq(1, 100, 1)) {

  env_params_it <- as.list(use_param_seq[i , ])

  temp <- iterate_model(env_params_it,
                        fixed_params,
                        domain_list$ln_leaf_l_1,
                        domain_list$ln_leaf_l_2,
                        domain_list$sqrt_area_1,
                        domain_list$sqrt_area_2,
                        pop_list,
                        v_holder,
                        iteration = i)

  pop_list          <- temp$pop_list
  kerns_temp        <- temp$kernels
  v_holder          <- temp$v_holder
  names(kerns_temp) <- paste(names(kerns_temp), i, sep = '_')
  kernel_holder     <- c(kernel_holder, kerns_temp)

}


usr_lambdas <- as.vector(pop_list$lambda)

ipmr_sub_kernels <- gen_di_stoch_param$sub_kernels

test_that('ipmr lambdas match user generated ones', {

  expect_equal(ipmr_lambdas, usr_lambdas, tolerance = 1e-10)

  kern_test <- map2_lgl(
    ipmr_sub_kernels,
    kernel_holder,
    .f = ~isTRUE(
      all.equal(
        unclass(.x),
                .y
      )
    )
  )

  expect_true(all(kern_test))


})


test_that("other outputs are of the expected form", {

  expect_equal(names(gen_di_stoch_param$env_seq),
               rando_names)

  ipmr_w <- right_ev(gen_di_stoch_param)
  hand_w <- lapply(pop_list[1:3], function(x) x[ , 26:101, drop = FALSE])

  expect_equal(ipmr_w, hand_w, ignore_attr = TRUE)
  expect_s3_class(ipmr_w, "ipmr_w")

  ipmr_v <- left_ev(gen_di_stoch_param, iterations = 100)
  expect_s3_class(ipmr_v, "ipmr_v")

  hand_v <- lapply(v_holder, function(x) t(x[26:101, ]))
  expect_equal(ipmr_v, hand_v, ignore_attr = TRUE)


})

mean_kernel_r <- function(kernel_holder, n_unique) {

  kern_nms <- gsub("_[0-9]$", "", names(kernel_holder)[1:n_unique]) %>%
    unique()

  out <- list()

  for(i in seq_along(kern_nms)) {

    use_kerns <- kernel_holder[grepl(kern_nms[i], names(kernel_holder))]

    dims <- c(dim(use_kerns[[1]])[1], dim(use_kerns[[1]])[2], length(use_kerns))

    holder <- array(0, dim = dims)

    for(j in seq_along(use_kerns)) {

      holder[, , j] <- use_kerns[[j]]

    }

    out[[i]] <- apply(holder, 1:2, mean)

    names(out)[i] <- paste("mean_", kern_nms[i], sep = "")

  }

  return(out)

}

test_that("mean_kernel works correctly for non-par_setarchical models", {


  mean_hand_kerns <- mean_kernel_r(kernel_holder, 6) %>%
    set_ipmr_classes()

  mean_ipmr_kerns <- mean_kernel(gen_di_stoch_param)

  expect_equal(mean_hand_kerns, mean_ipmr_kerns)


})


iterate_model_norm <- function(env_params,
                               fixed_params,
                               dom_1,
                               dom_2,
                               dom_3,
                               dom_4,
                               pop_list,
                               lambdas,
                               iteration) {

  k_xx_temp <- k_xx(fixed_params,
                    env_params,
                    dom_1, dom_2)
  k_zx_temp <- k_zx(fixed_params,
                    env_params,
                    dom_1, dom_2,
                    dom_3, dom_4)
  k_dx_temp <- k_dx(fixed_params,
                    dom_1, dom_2)
  k_xz_temp <- k_xz(fixed_params,
                    env_params,
                    dom_1,
                    dom_2,
                    dom_3,
                    dom_4)
  k_xd_temp <- k_xd(fixed_params,
                    env_params,
                    dom_1,
                    dom_2)
  k_zd_temp <- k_zd(fixed_params,
                    env_params,
                    dom_3,
                    dom_4)

  n_ln_leaf_l_t_1 <- k_xx_temp %*% pop_list$n_ln_leaf_l[ , iteration] +
    k_zx_temp %*% pop_list$n_sqrt_area[ , iteration] +
    k_dx_temp %*% pop_list$n_d[ , iteration]

  n_sqrt_area_t_1 <- k_xz_temp %*% pop_list$n_ln_leaf_l[ , iteration]

  n_d_t_1         <- k_xd_temp %*% pop_list$n_ln_leaf_l[ , iteration] +
    k_zd_temp %*% pop_list$n_sqrt_area[ , iteration]

  tot_size <- Reduce('sum',
                     c(n_ln_leaf_l_t_1,
                       n_sqrt_area_t_1,
                       n_d_t_1),
                     init = 0)

  pop_list$n_ln_leaf_l[ , (iteration + 1)] <- n_ln_leaf_l_t_1 / tot_size
  pop_list$n_sqrt_area[ , (iteration + 1)] <- n_sqrt_area_t_1 / tot_size
  pop_list$n_d[ , (iteration + 1)]         <- n_d_t_1 / tot_size

  lambdas[iteration]                       <- tot_size

  kerns_temp <- list(k_xx = k_xx_temp,
                     k_zx = k_zx_temp,
                     k_dx = k_dx_temp,
                     k_xz = k_xz_temp,
                     k_xd = k_xd_temp,
                     k_zd = k_zd_temp)

  out_list = list(kernels = kerns_temp,
                  pop_list = pop_list,
                  lambdas = lambdas)

  return(out_list)

}


test_that('normalize_pop_vec works', {


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
      normalize_pop_size = TRUE
    )

  lambdas_ipmr_norm <- lambda(gen_di_stoch_param,
                              type_lambda = 'all') %>%
    as.vector()


  param_seq <- gen_di_stoch_param$env_seq

  pop_list <- list(
    n_ln_leaf_l = matrix(NA_real_, nrow = 50, ncol = 101),
    n_sqrt_area = matrix(NA_real_, nrow = 50, ncol = 101),
    n_d         = matrix(NA_real_, nrow = 1,  ncol = 101)
  )

  tot_size <- Reduce('sum', unlist(init_pop_vec), init = 10)

  pop_list$n_ln_leaf_l[ , 1] <- init_pop_vec$ln_leaf_l / tot_size
  pop_list$n_sqrt_area[ , 1] <- init_pop_vec$sqrt_area / tot_size
  pop_list$n_d[ , 1]         <- 10 / tot_size
  lambdas <- numeric(100L)

  for(i in seq(1, 100, 1)) {

    env_params_it <- as.list(param_seq[i , ])
    names(env_params_it) <- gsub('env_params\\.', '', names(env_params_it))

    temp <- iterate_model_norm(env_params_it,
                               fixed_params,
                               domain_list$ln_leaf_l_1,
                               domain_list$ln_leaf_l_2,
                               domain_list$sqrt_area_1,
                               domain_list$sqrt_area_2,
                               pop_list,
                               lambdas,
                               iteration = i)

    pop_list          <- temp$pop_list
    kerns_temp        <- temp$kernels
    lambdas           <- temp$lambdas
    names(kerns_temp) <- paste(names(kerns_temp), i, sep = '_')
    kernel_holder     <- c(kernel_holder, kerns_temp)

  }


  expect_equal(lambdas, lambdas_ipmr_norm, tolerance = 1e-10)

  pop_sizes_ipmr <- lapply(gen_di_stoch_param$pop_state[1:3],
                           function(x) {
                             colSums(x)
                           }) %>%
    lapply(X = 1:101,
           function(x, sizes) {
             sum(sizes[[1]][x], sizes[[2]][x], sizes[[3]][x])
           },
           sizes = .) %>%
    unlist(use.names = FALSE)

  expect_equal(pop_sizes_ipmr, rep(1, 101), tolerance = 1e-15)

})


test_that("t variable works as advertised", {

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
      normalize_pop_size = TRUE
    )


  env_state <- gen_di_stoch_param$env_seq

  env_sampler <- function(environ_seq, iteration, rando_names) {

      out <- as.list(environ_seq[iteration, ]) %>%
        setNames(rando_names)

      return(out)

  }

  test_det_seq <- init_ipm(sim_gen    = "general",
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
      states    = list(c('ln_leaf_l', 'd')),

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
      states           = list(c('ln_leaf_l', 'd')),
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
      states           = list(c('sqrt_area', 'd')),
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
                         'd',
                         'ln_leaf_l',
                         'ln_leaf_l',
                         'sqrt_area'),
        state_end      = c('ln_leaf_l',
                         'ln_leaf_l',
                         'ln_leaf_l',
                         'sqrt_area',
                         'd',
                         'd')
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
      env_params = env_sampler(env_state, t, rando_names),
      data_list  = list(
        env_state   = env_state,
        rando_names = rando_names,
        env_sampler = env_sampler
      )
    ) %>%
    make_ipm(
      return_all_envs = TRUE,
      usr_funs   = list(
        inv_logit   = inv_logit,
        pois        = pois
      ),
      iterations = 100,
      normalize_pop_size = TRUE
    )

  det_l <- lambda(test_det_seq)
  ran_l <- lambda(gen_di_stoch_param)

  expect_equal(det_l, ran_l, tolerance = 1e-10)

})


test_that("Parameter sets work in parameter re-sampled model", {

  par_set_vals           <- rnorm(5) %>% as.list()
  names(par_set_vals)    <- paste("nr_s_z_int_", 1:5, sep = "")
  rando_names         <- gsub("nr_s_z_int",
                              "nr_s_z_int_fixed",
                              rando_names)
  fixed_params        <- c(fixed_params, par_set_vals)


  gen_di_stoch_param <- init_ipm(sim_gen    = "general",
                                 di_dd      = "di",
                                 det_stoch  = "stoch",
                                 "param") %>%
    define_kernel(
      name             = 'k_xx_site',
      family           = "CC",
      formula          = sig_n        *
        (1 - mu_n)   *
        (1 - beta_n) *
        gamma_nn     *
        d_ln_leaf_l,

      sig_n        = inv_logit(nr_s_z_int, nr_s_z_b, ln_leaf_l_1),
      nr_s_z_int   = nr_s_z_int_site + nr_s_z_int_fixed,
      gamma_nn     = dnorm(ln_leaf_l_2, nr_nr_mu, nr_nr_sd),
      nr_nr_mu     = nr_nr_int + nr_nr_b * ln_leaf_l_1,
      mu_n         = inv_logit(nr_d_z_int, nr_d_z_b, ln_leaf_l_1),
      beta_n       = inv_logit(nr_f_z_int, nr_f_z_b, ln_leaf_l_1),

      data_list        = fixed_params,
      states           = list(c('ln_leaf_l')),
      uses_par_sets    = TRUE,
      par_set_indices = list(site = 1:5),
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
      uses_par_sets    = FALSE,
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
      states    = list(c('ln_leaf_l', 'd')),

      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'gamma_nd')
    ) %>%
    define_kernel(
      name             = 'k_xz_site',
      family           = 'CC',
      formula          = sig_n         *
        (1 - mu_n)    *
        beta_n        *
        gamma_rn      *
        tau           *
        d_ln_leaf_l,

      sig_n            = inv_logit(nr_s_z_int, nr_s_z_b, ln_leaf_l_1),
      nr_s_z_int       = nr_s_z_int_site + nr_s_z_int_fixed,
      mu_n             = inv_logit(nr_d_z_int, nr_d_z_b, ln_leaf_l_1),
      beta_n           = inv_logit(nr_f_z_int, nr_f_z_b, ln_leaf_l_1),
      gamma_rn         = dnorm(sqrt_area_2, mu_nr_ra, nr_ra_sd),
      mu_nr_ra         = nr_ra_int + nr_ra_b * ln_leaf_l_1,
      tau              = inv_logit(tau_int, tau_b, ln_leaf_l_1),

      data_list        = fixed_params,
      states           = list(c('sqrt_area', 'ln_leaf_l')),
      uses_par_sets    = TRUE,
      par_set_indices = list(site = 1:5),
      evict_cor        = TRUE,
      evict_fun        = truncated_distributions('norm',
                                                 'gamma_rn')

    ) %>%
    define_kernel(
      name             = 'k_xd_site',
      family           = 'CD',
      formula          = sig_n * mu_n * d_ln_leaf_l,
      sig_n            = inv_logit(nr_s_z_int, nr_s_z_b, ln_leaf_l_1),
      nr_s_z_int       = nr_s_z_int_site + nr_s_z_int_fixed,
      mu_n             = inv_logit(nr_d_z_int, nr_d_z_b, ln_leaf_l_1),
      data_list        = fixed_params,
      states           = list(c('ln_leaf_l', 'd')),
      uses_par_sets    = TRUE,
      par_set_indices = list(site = 1:5),
      evict_cor        = FALSE
    ) %>%
    define_kernel(
      name             = 'k_zd',
      family           = 'CD',
      formula          = sig_r * mu_r * d_sqrt_area,
      sig_r            = inv_logit(ra_s_z_int, ra_s_z_b, sqrt_area_1),
      mu_r             = inv_logit(ra_d_z_int, ra_d_z_b, sqrt_area_1),
      data_list        = fixed_params,
      states           = list(c('sqrt_area', 'd')),
      uses_par_sets    = FALSE,
      evict_cor        = FALSE
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c(paste('k_',
                               c('xx_site',
                                 'zx', 'dx', 'xz_site', 'xd_site', 'zd'),
                               sep = "")),
        int_rule     = rep('midpoint', 6),
        state_start    = c('ln_leaf_l',
                         'sqrt_area',
                         'd',
                         'ln_leaf_l',
                         'ln_leaf_l',
                         'sqrt_area'),
        state_end      = c('ln_leaf_l',
                         'ln_leaf_l',
                         'ln_leaf_l',
                         'sqrt_area',
                         'd',
                         'd')
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
      normalize_pop_size = TRUE,
      kernel_seq = sample(1:5, 100, TRUE)
    )

  use_param_seq <- gen_di_stoch_param$env_seq

  k_xx <- function(fixed_params, env_params, use_kern, dom_1, dom_2) {

    d_z <- dom_1[2] - dom_1[1]
    l   <- min(dom_1) - d_z / 2
    u   <- max(dom_1) + d_z / 2

    int_nm <- paste("nr_s_z_int_", use_kern, sep = "")

    nr_nr_mu <- env_params$nr_nr_int + env_params$nr_nr_b * dom_1

    sig_n_int <- env_params$nr_s_z_int_fixed + fixed_params[[int_nm]]
    sig_n    <- inv_logit(sig_n_int, env_params$nr_s_z_b, dom_1)

    ev <- pnorm(u, nr_nr_mu, fixed_params$nr_nr_sd) -
      pnorm(l, nr_nr_mu, fixed_params$nr_nr_sd)

    gamma_nn <- dnorm(dom_2, nr_nr_mu, fixed_params$nr_nr_sd) /
      ev

    mu_n     <- inv_logit(fixed_params$nr_d_z_int, fixed_params$nr_d_z_b, dom_1)

    beta_n   <- inv_logit(env_params$nr_f_z_int, env_params$nr_f_z_b, dom_1)

    disc_fun <- (
      sig_n        *
        (1 - mu_n)   *
        (1 - beta_n) *
        gamma_nn     *
        d_z
    ) %>%
      matrix(nrow = 50, ncol = 50, byrow = TRUE)

    return(disc_fun)
  }

  k_xz <- function(fixed_params, env_params, use_kern, dom_1, dom_2, dom_3, dom_4) {

    d_z  <- dom_1[2] - dom_1[1]
    d_z1 <- dom_3[2] - dom_3[1]
    l   <- min(dom_3) - d_z1 / 2
    u   <- max(dom_3) + d_z1 / 2

    int_nm <- paste("nr_s_z_int_", use_kern, sep = "")

    sig_n_int <- env_params$nr_s_z_int_fixed + fixed_params[[int_nm]]
    sig_n    <- inv_logit(sig_n_int, env_params$nr_s_z_b, dom_1)
    mu_n     <- inv_logit(fixed_params$nr_d_z_int, fixed_params$nr_d_z_b, dom_1)
    beta_n   <- inv_logit(env_params$nr_f_z_int, env_params$nr_f_z_b, dom_1)
    mu_nr_ra <- env_params$nr_ra_int + fixed_params$nr_ra_b * dom_1
    tau      <- inv_logit(env_params$tau_int, fixed_params$tau_b, dom_1)

    ev       <- pnorm(u, mu_nr_ra, fixed_params$nr_ra_sd) -
      pnorm(l, mu_nr_ra, fixed_params$nr_ra_sd)

    gamma_rn <- dnorm(dom_4, mu_nr_ra, fixed_params$nr_ra_sd) /
      ev


    disc_fun <- (
      sig_n         *
        (1 - mu_n)    *
        beta_n        *
        gamma_rn      *
        tau           *
        d_z

    ) %>%
      matrix(nrow = 50, ncol = 50, byrow = TRUE)

    return(disc_fun)
  }

  k_xd <- function(fixed_params, env_params, use_kern, dom_1, dom_2) {

    d_z <- dom_1[2] - dom_1[1]
    int_nm <- paste("nr_s_z_int_", use_kern, sep = "")

    sig_n_int <- env_params$nr_s_z_int_fixed + fixed_params[[int_nm]]
    sig_n    <- inv_logit(sig_n_int, env_params$nr_s_z_b, dom_1)
    mu_n  <- inv_logit(fixed_params$nr_d_z_int, fixed_params$nr_d_z_b, dom_1)

    disc_fun <- (sig_n * mu_n * d_z) %>%
      unique() %>%
      matrix(nrow = 1, ncol = 50, byrow = TRUE)

    return(disc_fun)
  }

  iterate_model_norm <- function(env_params,
                                 use_kern,
                                 fixed_params,
                                 dom_1,
                                 dom_2,
                                 dom_3,
                                 dom_4,
                                 pop_list,
                                 lambdas,
                                 iteration) {

    k_xx_temp <- k_xx(fixed_params,
                      env_params,
                      use_kern,
                      dom_1, dom_2)
    k_zx_temp <- k_zx(fixed_params,
                      env_params,
                      dom_1, dom_2,
                      dom_3, dom_4)
    k_dx_temp <- k_dx(fixed_params,
                      dom_1, dom_2)
    k_xz_temp <- k_xz(fixed_params,
                      env_params,
                      use_kern,
                      dom_1,
                      dom_2,
                      dom_3,
                      dom_4)
    k_xd_temp <- k_xd(fixed_params,
                      env_params,
                      use_kern,
                      dom_1,
                      dom_2)
    k_zd_temp <- k_zd(fixed_params,
                      env_params,
                      dom_3,
                      dom_4)

    n_ln_leaf_l_t_1 <- k_xx_temp %*% pop_list$n_ln_leaf_l[ , iteration] +
      k_zx_temp %*% pop_list$n_sqrt_area[ , iteration] +
      k_dx_temp %*% pop_list$n_d[ , iteration]

    n_sqrt_area_t_1 <- k_xz_temp %*% pop_list$n_ln_leaf_l[ , iteration]

    n_d_t_1         <- k_xd_temp %*% pop_list$n_ln_leaf_l[ , iteration] +
      k_zd_temp %*% pop_list$n_sqrt_area[ , iteration]

    tot_size <- Reduce('sum',
                       c(n_ln_leaf_l_t_1,
                         n_sqrt_area_t_1,
                         n_d_t_1),
                       init = 0)

    pop_list$n_ln_leaf_l[ , (iteration + 1)] <- n_ln_leaf_l_t_1 / tot_size
    pop_list$n_sqrt_area[ , (iteration + 1)] <- n_sqrt_area_t_1 / tot_size
    pop_list$n_d[ , (iteration + 1)]         <- n_d_t_1 / tot_size

    lambdas[iteration]                       <- tot_size

    kerns_temp <- list(k_xx = k_xx_temp,
                       k_zx = k_zx_temp,
                       k_dx = k_dx_temp,
                       k_xz = k_xz_temp,
                       k_xd = k_xd_temp,
                       k_zd = k_zd_temp)

    out_list = list(kernels = kerns_temp,
                    pop_list = pop_list,
                    lambdas = lambdas)

    return(out_list)

  }

  kernel_holder <- list()
  lambdas <- numeric
  pop_list <- list(
    n_ln_leaf_l = matrix(NA_real_, nrow = 50, ncol = 101),
    n_sqrt_area = matrix(NA_real_, nrow = 50, ncol = 101),
    n_d         = matrix(NA_real_, nrow = 1,  ncol = 101)
  )

  tot_size <- Reduce('sum', unlist(init_pop_vec), init = 10)

  pop_list$n_ln_leaf_l[ , 1] <- init_pop_vec$ln_leaf_l / tot_size
  pop_list$n_sqrt_area[ , 1] <- init_pop_vec$sqrt_area / tot_size
  pop_list$n_d[ , 1]         <- 10 / tot_size
  lambdas <- numeric(100L)


  for(i in seq(1, 100, 1)) {

    env_params_it <- as.list(use_param_seq[i , ])

    use_kerns     <- env_params_it[[14]]
    env_params_it <- env_params_it[-c(14)]



    temp <- iterate_model_norm(env_params_it,
                               use_kerns,
                               fixed_params,
                               domain_list$ln_leaf_l_1,
                               domain_list$ln_leaf_l_2,
                               domain_list$sqrt_area_1,
                               domain_list$sqrt_area_2,
                               pop_list,
                               lambdas,
                               iteration = i)

    pop_list          <- temp$pop_list
    kerns_temp        <- temp$kernels
    lambdas           <- temp$lambdas
    names(kerns_temp) <- paste(names(kerns_temp), i, sep = '_')
    kernel_holder     <- splice(kernel_holder, kerns_temp)

  }

  hand_lam <- mean(log(lambdas[11:100]))

  ipmr_lam <- lambda(gen_di_stoch_param)

  expect_equal(hand_lam, ipmr_lam, tolerance = 1e-9)

  mean_hand_kerns <- mean_kernel_r(kernel_holder, 6L) %>%
    set_ipmr_classes()
  mean_ipmr_kerns <- mean_kernel(gen_di_stoch_param)
  names(mean_ipmr_kerns) <- gsub("_site", "", names(mean_ipmr_kerns))

  expect_equal(mean_hand_kerns, mean_ipmr_kerns)

})
