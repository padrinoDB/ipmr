context('general density independent stochastic parameter resampled models')

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

fixed_params <- all_params[ind_fixed]
rando_means  <- vapply(all_params[ind_rand],
                       function(x) x$mu,
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







