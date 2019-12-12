# print(search())
context('general density independent stochastic kernel models')

# Test with parameters from Aikens & Roach 2014 Population dynamics in central
# and edge populations of a narrowly endemic plant. They only report
# lambdas on a figure, and so I won't try to hit those perfectly.
# Rather, just going to simulate a model using their parameters and then
# try to recover it using ipmr.

library(rlang)
library(purrr)

flatten_to_depth <- ipmr:::.flatten_to_depth

hier_effs <- list(
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
                                   nms       = unlist(hier_effs)) %>%
  flatten_to_depth(1)


# Target lambdas extracted from Figure 5 using Data Thief. Not precise, but
# probably good to within 1 or 2 decimal points

target_lambdas <- c(
  whitetop  = 1.193,
  mt_rogers = 0.752,
  roan      = 1.147,
  big_bald  = 1.262,
  bob_bald  = 1.058,
  oak_knob  = 0.921
)

# Now, compose functions for vital rates and kernels!
# Functions are similar between vital rates and stages, but I've split them
# out here to clarify when each part is being used

# x_1/2 = ln_leaf_l t, t+1
# z_1/2 = sqrt_area t, t+1
# *_n   = non-reproductive plants
# *_r   = reproductive plants

# survival functions
sig_n <- function(x_1, params) {

  inv_logit(params[[1]], params[[2]], x_1)

}

sig_r <- function(z_1, params) {

  inv_logit(params[[1]], params[[2]], z_1)

}

# probability of dormancy
mu_n <- function(x_1, params) {

  inv_logit(params[[1]], params[[2]], x_1)

}

mu_r <- function(z_1, params) {

  inv_logit(params[[1]], params[[2]], z_1)

}

# probability of flowering
beta_n <- function(x_1, params) {

  inv_logit(params[[1]], params[[2]], x_1)

}

# Growth functions

gamma_nn <- function(x_2, x_1, params, L, U) {


  lin_prob(params[[1]], params[[2]], params[[3]], x_1, x_2, L, U)

}

gamma_nr <- function(x_2, z_1, params, L, U) {

  lin_prob(params[[1]], params[[2]], params[[3]], z_1, x_2, L, U)


}

gamma_nd <- function(x_2, x_1, params, L, U) {

  params <- list(params[[1]], 0, params[[2]])

  lin_prob(params[[1]], params[[2]], params[[3]], x_1, x_2, L, U)

  # ev <- pnorm(U,
  #             mean = x_2,
  #             sd   = params[[2]]) - pnorm(L,
  #                                         mean = x_2,
  #                                         sd   = params[[2]])
  #
  # out <- dnorm(x_2, params[[1]], params[[2]]) / ev

  # return(out)

}

gamma_rn <- function(z_2, x_1, params, L, U) {

  lin_prob(params[[1]], params[[2]], params[[3]], x_1, z_2, L, U)

}

gamma_sr <- function(x_2, x_1, params, L, U) {

  lin_prob(params[[1]], params[[2]], params[[3]], x_1, x_2, L, U)

}

# Flower production

phi <- function(z_1, params) {

  pois(params[[1]], params[[2]], z_1)

}

tau <- function(z_1, params) {

  inv_logit(params[[1]], params[[2]], z_1)

}

# Now, kernel functions. named as K_s(t)->s(t+1) where s = states %in% c(d, x, z)

# eq 1a

k_xx <- function(x_2, x_1, params) {

  d_x <- x_1[2] - x_1[1]
  L_x <- 0.234 # size bounds for eviction correction
  U_x <- 2.97

  sig_n(x_1, params[1:2])                      *
    (1 - mu_n(x_1, params[3:4]))               *
    (1 - beta_n(x_1, params[5:6]))             *
    gamma_nn(x_2, x_1, params[7:9], L_x, U_x)  *
    d_x

}

# eq 1b + eq 1d

k_zx <- function(x_2, z_2, x_1, z_1, params){

  d_z <- z_1[2] - z_1[1]
  L_x <- 0.234 # size bounds for eviction correction
  U_x <- 2.97



  nu <- params[[length(params)]]

  # Two types of movement - survival and avoidance of dormancy AND
  # creation of new plants through seeds

  sig_r(z_1, params[17:18])                     *
    (1 - mu_r(z_1, params[19:20]))              *
    gamma_nr(x_2, z_1, params[21:23], L_x, U_x) *
    d_z                                         + # And now... new plants!

    phi(z_1, params[13:14])                     *
    nu                                          *
    gamma_sr(x_2, x_1, params[26:28], L_x, U_x) *
    d_z
}

# eq 1c
k_dx <- function(x_2, x_1, params) {

  d_x <- x_1[2] - x_1[1]
  L_x <- 0.234 # size bounds for eviction correction
  U_x <- 2.97

  gamma_nd(x_2, x_1, params[24:25], L_x, U_x) * d_x

}

# kernel for reproductive individuals

# eq 2

k_xz <- function(x_2, z_2, x_1, z_1, params) {

  d_x <- x_1[2] - x_1[1]
  L_z <- 0.567
  U_z <- 4.257

  sig_n(x_1, params[1:2])                       *
    (1 - mu_n(x_1, params[3:4]))                *
    beta_n(x_1, params[5:6])                    *
    gamma_rn(z_2, x_1, params[10:12], L_z, U_z) *
    tau(z_1, params[15:16])                     *
    d_x

}

# eq 3a

k_xd <- function(x_1, params) {

  d_x <- x_1[2] - x_1[1]

  sig_n(x_1, params[1:2])  *
    mu_n(x_1, params[3:4]) *
    d_x

}


# eq 3b

k_zd <- function(z_1, params) {

  d_z <- z_1[2] - z_1[1]

  sig_r(z_1, params[17:18])  *
    mu_r(z_1, params[19:20]) *
    d_z

}

# Now initialize kernel lists and population vectors
# Use 100 iterations of each model to compute lambdas

n_iterations   <- 100

models         <- vector('list', length(hier_effs$site))
names(models)  <- hier_effs$site

pop_holders        <- vector('list', length(hier_effs$site))
names(pop_holders) <- hier_effs$site

# All populations will get the same initial population vector, and it's
# just a random set of numbers in their range of possible values

set.seed(214512)

init_pop_vec      <- list(ln_leaf_l = runif(domains$sqrt_area[3]),
                          sqrt_area = runif(domains$sqrt_area[3]))

pop_holders <- lapply(seq_along(pop_holders),
                      function(x, pop_holders, init_pop_vecs, iterations) {

                        pop_t_1 <- init_pop_vec
                        sqrt_area <- array(NA_real_,
                                           dim = c(length(pop_t_1$sqrt_area),
                                                   (iterations + 1)))
                        ln_leaf_l <- array(NA_real_,
                                           dim = c(length(pop_t_1$ln_leaf_l),
                                                   (iterations + 1)))

                        d <- array(NA_real_, dim = c(1, (iterations + 1)))

                        sqrt_area[ , 1] <- pop_t_1$sqrt_area
                        ln_leaf_l[ , 1] <- pop_t_1$ln_leaf_l
                        d[ , 1]         <- 10

                        nm <- names(pop_holders)[x]

                        pop_holders[[x]] <- rlang::list2(!! nm := list(
                          n_sqrt_area = sqrt_area,
                          n_ln_leaf_l = ln_leaf_l,
                          n_d         = d
                        ))

                        return(pop_holders[[x]])
                      },
                      pop_holders   = pop_holders,
                      init_pop_vecs = init_pop_vec,
                      iterations    = n_iterations) %>%
  flatten_to_depth(2L)


# Now, we're all set to begin simulating! The first step is to see if
# the simulation I code below reproduces what happens in the paper.
# If that works, then I'll create an ipmr version of the model and
# see if I can match the simulation. The reason I don't try to match the ipmr
# version to the paper version is that coefficients are only reported to 4 decimals,
# and we only have deterministic lambdas for each kernel (which themselves were
# obtained using dataThief - not the most precise method!!), and so our expected
# numerical accuracy would be much lower for that comparison.
# This could potentially obscure small, but important flaws in the logic of the
# make_ipm.general_di_stoch_kern method!

# Make mesh bounds and mid points. length.out = length + 1 because we are making
# bounds. when the meshpoints are generated, the resulting length is length(bounds - 1)

bounds <- list(sqrt_area = seq(domains$sqrt_area[1],
                               domains$sqrt_area[2],
                               length.out = (domains$sqrt_area[3] + 1)),
               ln_leaf_l = seq(domains$ln_leaf_l[1],
                               domains$ln_leaf_l[2],
                               length.out = (domains$ln_leaf_l[3] + 1)))


mesh_p <- lapply(bounds,
                 function(x) {
                   l <- length(x) - 1

                   out <- 0.5 * (x[1:l] + x[2:(l + 1)])
                   return(out)

                 })

all_mesh_p <- list(sqrt_area = expand.grid(z_1 = mesh_p$sqrt_area,
                                           z_2 = mesh_p$sqrt_area),
                   ln_leaf_l = expand.grid(x_1 = mesh_p$ln_leaf_l,
                                           x_2 = mesh_p$ln_leaf_l)) %>%
  flatten_to_depth(1L)


# Everything is initialized. We'll loop over the populations and iterate
# the model for each one individually, compute lambdas, and store those
# for comparison using test_that()

for(i in seq_along(models)) {

  pop <- names(pop_holders)[i]

  # each parameter has 6 entries in full_data_list - 1 for each population.
  # this generates an index to grab parameters unique to each population.
  # populations are always in the same order for each parameter, so this grabs
  # one for each.

  param_ind   <- seq(i, length(full_data_list), by = 6)

  site_params <- full_data_list[param_ind]

  # just to make sure we haven't lost our minds

  stopifnot(all(grepl(pop, names(site_params))))

  # Next, use each kernel function to generate a (discretized!) vector.
  # transform these into iteration matrices, and then iterate using a loop over
  # n_iterations. We'll also store the complete set of kernels in the `models`
  # object in case we need to debug later.

  kernels <- list(
    kern_xx = k_xx(all_mesh_p$x_2,
                   all_mesh_p$x_1,
                   site_params) %>%
      matrix(nrow = 50, ncol = 50, byrow = TRUE),
    kern_zx = k_zx(all_mesh_p$x_2,
                   all_mesh_p$z_2,
                   all_mesh_p$x_1,
                   all_mesh_p$z_1,
                   site_params) %>%
      matrix(nrow = 50, ncol = 50, byrow = TRUE),

    # This requires the subsetting because we only want a single
    # column vector (e.g. the first entry of each row).
    kern_dx = k_dx(all_mesh_p$x_2,
                   all_mesh_p$x_1,
                   site_params)[seq(1, 2500, by = 50)] %>%
      matrix(nrow = 50, ncol = 1, byrow = TRUE),
    kern_xz = k_xz(all_mesh_p$x_2,
                   all_mesh_p$z_2,
                   all_mesh_p$x_1,
                   all_mesh_p$z_1,
                   site_params) %>%
      matrix(nrow = 50, ncol = 50, byrow = TRUE),

    # These are simpler to subset because we just want the first row of the
    # matrix.
    kern_xd = k_xd(all_mesh_p$x_1,
                   site_params)[1:50] %>%
      matrix(nrow = 1, ncol = 50, byrow = TRUE),
    kern_zd = k_zd(all_mesh_p$z_1,
                   site_params)[1:50]  %>%
      matrix(nrow = 1, ncol = 50, byrow = TRUE)

  )

  temp_pop <- pop_holders[pop]

  for(j in seq_len(n_iterations)) {

    # These are the full kernel equations from pages 1853-1854.

    x_t_1 <- kernels$kern_xx %*% temp_pop[[1]]$n_ln_leaf_l[ , j] +
      kernels$kern_zx %*% temp_pop[[1]]$n_sqrt_area[ , j] +
      kernels$kern_dx %*% temp_pop[[1]]$n_d[j]

    z_t_1 <- kernels$kern_xz %*% temp_pop[[1]]$n_ln_leaf_l[ , j]

    d_t_1 <- kernels$kern_xd %*% temp_pop[[1]]$n_ln_leaf_l[ , j] +
      kernels$kern_zd %*% temp_pop[[1]]$n_sqrt_area[ , j]

    temp_pop[[1]]$n_ln_leaf_l[ , (j + 1)] <- x_t_1
    temp_pop[[1]]$n_sqrt_area[ , (j + 1)] <- z_t_1
    temp_pop[[1]]$n_d[ , (j + 1)]         <- d_t_1

  }

  # Some house keeping so we can use ipmr convenience functions to compute
  # lambda for populations which have converged to stable growth/decline

  temp_pop <- list(pop_state = flatten_to_depth(temp_pop, 1))

  pop_holders[[i]] <- temp_pop
  models[pop]      <- list(flatten_to_depth(kernels, 1))

}

actual_lambdas <- vapply(pop_holders,
                         function(x) ipmr:::.stoch_lambda_pop_size(x),
                         numeric(1L))

# And now, for the ipmr version.


gen_di_stoch_kern <- init_ipm('general_di_stoch_kern') %>%
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
    has_hier_effs    = TRUE,
    levels_hier_effs = hier_effs,
    evict            =  TRUE,
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
    has_hier_effs = TRUE,
    levels_hier_effs = hier_effs,
    evict = TRUE,
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
    states = list(c('ln_leaf_l')),
    has_hier_effs = TRUE,
    levels_hier_effs = hier_effs,
    evict = TRUE,
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
    has_hier_effs    = TRUE,
    levels_hier_effs = hier_effs,
    evict            = TRUE,
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
    states           = list(c('ln_leaf_l')),
    has_hier_effs    = TRUE,
    levels_hier_effs = hier_effs,
    evict            = FALSE
  ) %>%
  define_kernel(
    name             = 'k_zd_site',
    family           = 'CD',
    formula          = sig_r_site * mu_r_site * d_sqrt_area,
    sig_r_site        = inv_logit(ra_s_z_int_site, ra_s_z_b_site, sqrt_area_1),
    mu_r_site         = inv_logit(ra_d_z_int_site, ra_d_z_b_site, sqrt_area_1),
    data_list        = full_data_list,
    states           = list(c('sqrt_area')),
    has_hier_effs    = TRUE,
    levels_hier_effs = hier_effs,
    evict            = FALSE
  ) %>%
  define_k(
    name = 'K_site',
    n_ln_leaf_l_site_t_1 = k_xx_site %*% n_ln_leaf_l_site_t +
      k_zx_site %*% n_sqrt_area_site_t +
      k_dx_site %*% n_d_site_t,
    n_sqrt_area_site_t_1 = k_xz_site %*% n_ln_leaf_l_site_t,
    n_d_site_t_1 =         k_xd_site %*% n_ln_leaf_l_site_t +
      k_zd_site %*% n_sqrt_area_site_t,
    family = 'IPM',
    data_list = full_data_list,
    states    = list(c('sqrt_area', 'ln_leaf_l')),
    has_hier_effs    = TRUE,
    levels_hier_effs = hier_effs
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c(paste('k_',
                             c('xx',
                               'zx', 'dx', 'xz', 'xd', 'zd'),
                             '_site',
                             sep = ""),
                       'K_site'),
      int_rule     = rep('midpoint', 7),
      dom_start    = c('ln_leaf_l',
                       'sqrt_area',
                       NA_character_,
                       'ln_leaf_l',
                       'ln_leaf_l',
                       'sqrt_area',
                       NA_character_),
      dom_end      = c('ln_leaf_l',
                       'ln_leaf_l',
                       'ln_leaf_l',
                       'sqrt_area',
                       NA_character_,
                       NA_character_,
                       NA_character_)
    )
  ) %>%
  define_domains(
    sqrt_area = c(0.63 * 0.9, 3.87 * 1.1, 50),
    ln_leaf_l = c(0.26 * 0.9, 2.70 * 1.1, 50)
  ) %>%
  define_pop_state(
    pop_vectors = list(
      n_ln_leaf_l_site = init_pop_vec$ln_leaf_l,
      n_sqrt_area_site = init_pop_vec$sqrt_area,
      n_d_site         = 10
    )
  ) %>%
  make_ipm(
    return_all = TRUE,
    usr_funs   = list(
      inv_logit = inv_logit,
      pois      = pois
    ),
    iterations = 100
  )

# need to rethink how output is generated - no reason I shouldn't be able to
# plug model object straight into lambda generic!

ipmr_pop_sizes <- list(
  whitetop  = list(pop_state = gen_di_stoch_kern$pop_state[1:3]),
  mt_rogers = list(pop_state = gen_di_stoch_kern$pop_state[4:6]),
  roan      = list(pop_state = gen_di_stoch_kern$pop_state[7:9]),
  big_bald  = list(pop_state = gen_di_stoch_kern$pop_state[10:12]),
  bob_bald  = list(pop_state = gen_di_stoch_kern$pop_state[13:15]),
  oak_knob  = list(pop_state = gen_di_stoch_kern$pop_state[16:18])
)

ipmr_lambdas <- vapply(ipmr_pop_sizes,
                       function(x) ipmr:::.stoch_lambda_pop_size(x),
                       numeric(1L))

test_that('ipmr version matches simulation', {

  expect_equal(ipmr_lambdas, actual_lambdas, tol = 1e-10)

  # Standardized right eigenvectors

  ipmr_w <- ipmr_pop_sizes$whitetop$pop_state$pop_state_ln_leaf_l_whitetop[ , 101] /
    sum(ipmr_pop_sizes$whitetop$pop_state$pop_state_ln_leaf_l_whitetop[ , 101])

  hand_w <- pop_holders$whitetop$pop_state$n_ln_leaf_l[ , 101] /
    sum(pop_holders$whitetop$pop_state$n_ln_leaf_l [ , 101])

  expect_equal(ipmr_w, hand_w, tol = 1e-10)

  kern_tests <- logical(length(gen_di_stoch_kern$sub_kernels))

  kern_suffs <- c('xx', 'zx', 'xd', 'xz', 'dx', 'dz')

  for(i in seq_along(hier_effs$site)) {

    ind <- seq(i,
               length(gen_di_stoch_kern$sub_kernels),
               by = length(hier_effs$site))

    kern_tests[ind] <- compare_kernels('gen_di_stoch_kern',
                                       'models',
                                       nms_ipmr = paste('k_',
                                                        paste(kern_suffs,
                                                              hier_effs$site[i],
                                                              sep = "_"),
                                                        sep = ""),
                                       nms_hand = paste(hier_effs$site[i],
                                                        '$kern_',
                                                        kern_suffs,
                                                        sep = ""))

  }

  if(any(! kern_tests)) warning('kernel indices: ', which(! kern_tests),
                                ' are not identical')

  expect_true(all(kern_tests))



})

# Next, we need to test the internal iteration machinery - particularly the
# kern_seq properties. This tests whether or not we can recover stochastic simulations
# given the same sequence of kernels (i.e. did I write the user-defined sequence
# machinery correctly). Since we've already verified that the kernels are all identical,
# there should never be any issues here.

usr_seq <- sample(hier_effs$site,
                  n_iterations,
                  replace = TRUE)


pop_holder_stoch <- list(

  n_ln_leaf_l = array(NA_real_,
                      dim = c(domains$ln_leaf_l[3],
                              (n_iterations + 1))),
  n_sqrt_area = array(NA_real_,
                      dim = c(domains$sqrt_area[3],
                              (n_iterations + 1))),
  n_d         = array(NA_real_,
                      dim = c(1,
                              (n_iterations + 1)))
)

pop_holder_stoch$n_ln_leaf_l[ , 1] <- init_pop_vec$ln_leaf_l
pop_holder_stoch$n_sqrt_area[ , 1] <- init_pop_vec$sqrt_area
pop_holder_stoch$n_d[ , 1]         <- 10

for(i in seq_len(n_iterations)) {

  temp <- models[[usr_seq[i]]]

  x_t_1 <- temp$kern_xx %*% pop_holder_stoch$n_ln_leaf_l[ , i] +
           temp$kern_zx %*% pop_holder_stoch$n_sqrt_area[ , i] +
           temp$kern_dx %*% pop_holder_stoch$n_d[ , i]

  z_t_1 <- temp$kern_xz %*% pop_holder_stoch$n_ln_leaf_l[ , i]

  d_t_1 <- temp$kern_xd %*% pop_holder_stoch$n_ln_leaf_l[ , i] +
           temp$kern_zd %*% pop_holder_stoch$n_sqrt_area[ , i]


  pop_holder_stoch$n_ln_leaf_l[ , (i + 1)] <- x_t_1
  pop_holder_stoch$n_sqrt_area[ , (i + 1)] <- z_t_1
  pop_holder_stoch$n_d[ , (i + 1)]         <- d_t_1

}

stoch_mod_hand <- list(pop_state = pop_holder_stoch)


# ipmr_version

temp_proto <- gen_di_stoch_kern$proto_ipm[!grepl('K_site',
                                                 gen_di_stoch_kern$proto_ipm$kernel_id), ]

stoch_mod_ipmr <- temp_proto %>%
  define_k(
    name = 'K_site',
    n_ln_leaf_l_t_1 = k_xx_site %*% n_ln_leaf_l_t +
                      k_zx_site %*% n_sqrt_area_t +
                      k_dx_site %*% n_d_t,
    n_sqrt_area_t_1 = k_xz_site %*% n_ln_leaf_l_t,
    n_d_t_1         = k_xd_site %*% n_ln_leaf_l_t +
                      k_zd_site %*% n_sqrt_area_t,
    family = 'IPM',
    data_list = full_data_list,
    states    = list(c('sqrt_area', 'ln_leaf_l')),
    has_hier_effs    = TRUE,
    levels_hier_effs = hier_effs
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c(paste('k_',
                             c('xx',
                               'zx', 'dx', 'xz', 'xd', 'zd'),
                             '_site',
                             sep = ""),
                       'K_site'),
      int_rule     = rep('midpoint', 7),
      dom_start    = c('ln_leaf_l',
                       'sqrt_area', NA_character_,
                       'ln_leaf_l', 'ln_leaf_l', 'sqrt_area', NA_character_),
      dom_end      = c('ln_leaf_l',
                       'ln_leaf_l', 'ln_leaf_l',
                       'sqrt_area', NA_character_, NA_character_, NA_character_)
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
    kernel_seq = usr_seq
  )


pop_ipmr    <- stoch_mod_ipmr$pop_state
all_lambdas <- list(
  hand = numeric(100L),
  ipmr = numeric(100L)
)

for(i in seq_len(n_iterations)) {

  hand_t <- lapply(pop_holder_stoch,
                   function(x, it) {
                     sum(x[ , it])

                   },
                   it = i) %>%
    unlist() %>%
    sum()

  hand_t_1 <- lapply(pop_holder_stoch,
                     function(x, it) {
                       sum(x[ , (it + 1)])
                     },
                     it = i) %>%
    unlist() %>%
    sum()


  ipmr_t <- lapply(pop_ipmr,
                   function(x, it) {
                     sum(x[ , it])

                   },
                   it = i) %>%
    unlist() %>%
    sum()

  ipmr_t_1 <- lapply(pop_ipmr,
                     function(x, it) {
                       sum(x[ , (it + 1)])
                     },
                     it = i) %>%
    unlist() %>%
    sum()

  all_lambdas$hand[i] <- hand_t_1 / hand_t
  all_lambdas$ipmr[i] <- ipmr_t_1 / ipmr_t
}

test_that('general stochastic simulations match hand generated ones', {

  expect_equal(all_lambdas$hand, all_lambdas$ipmr, tol = 1e-4)

})


test_that('evict_fun warnings are correctly generated', {

  test_evict_fun <-
    gen_di_stoch_kern$proto_ipm[!grepl('k_zx_site', gen_di_stoch_kern$proto_ipm$kernel_id), ]

  wrngs <- capture_warnings(
    test_evict_fun_warning <- test_evict_fun %>%
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
        has_hier_effs = TRUE,
        levels_hier_effs = hier_effs,
        evict = TRUE,
        evict_fun = truncated_distributions('norm',
                                            c('gamma_nr_site',
                                              'gamma_sr_site'))
      ) %>%
      define_impl(
        make_impl_args_list(
          kernel_names = c(
            'k_xx_site',
            'k_dx_site',
            'k_xz_site',
            'k_xd_site',
            'k_zd_site',
            'K_site',
            'k_zx_site'
          ),
          int_rule     = rep('midpoint', 7),
          dom_start    = c('ln_leaf_l',
                           NA_character_,
                           'ln_leaf_l',
                           'ln_leaf_l',
                           'sqrt_area',
                           NA_character_,
                           'sqrt_area'),
          dom_end      = c('ln_leaf_l',
                           'ln_leaf_l',
                           'sqrt_area',
                           NA_character_,
                           NA_character_,
                           NA_character_,
                           'ln_leaf_l')
        )
      ) %>%
      define_domains(
        sqrt_area = c(0.63 * 0.9, 3.87 * 1.1, 50),
        ln_leaf_l = c(0.26 * 0.9, 2.70 * 1.1, 50)
      ) %>%
      define_pop_state(
        pop_vectors = list(
          n_ln_leaf_l_site = init_pop_vec$ln_leaf_l,
          n_sqrt_area_site = init_pop_vec$sqrt_area,
          n_d_site         = 10
        )
      ) %>%
      make_ipm(
        return_all = TRUE,
        usr_funs   = list(
          inv_logit = inv_logit,
          pois      = pois
        ),
        iterations = 100
      )
    )

  test_text <-
    "length of 'fun' in 'truncated_distributions()' is not equal to length of 'param'. Recycling 'fun'."

  expect_equal(wrngs[1], test_text)

})
