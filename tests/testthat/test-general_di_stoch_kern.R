
# Test with parameters from Aikens & Roach 2014 Population dynamics in central
# and edge populations of a narrowly endemic plant. They only report
# lambdas on a figure, and so I won't try to hit those perfectly.
# Rather, just going to simulate the kernels using their parameters and then
# try to recover it using ipmr.

# These will get fed into a stochastic simulation afterwards.

library(rlang)
library(purrr)

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

PF_xx <- function(x_2, x_1, params) {

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

PF_zx <- function(x_2, z_2, x_1, z_1, params){

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
PF_dx <- function(x_2, x_1, params) {

  d_x <- x_1[2] - x_1[1]
  L_x <- 0.234 # size bounds for eviction correction
  U_x <- 2.97

  gamma_nd(x_2, x_1, params[24:25], L_x, U_x) * d_x

}

# kernel for reproductive individuals

# eq 2

PF_xz <- function(x_2, z_2, x_1, z_1, params) {

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

PF_xd <- function(x_1, params) {

  d_x <- x_1[2] - x_1[1]

  sig_n(x_1, params[1:2])  *
    mu_n(x_1, params[3:4]) *
    d_x

}


# eq 3b

PF_zd <- function(z_1, params) {

  d_z <- z_1[2] - z_1[1]

  sig_r(z_1, params[17:18])  *
    mu_r(z_1, params[19:20]) *
    d_z

}

# All populations will get the same initial population vector, and it's
# just a random set of numbers in their range of possible values

set.seed(214512)

init_pop_vec      <- list(ln_leaf_l = runif(domains$sqrt_area[3]),
                          sqrt_area = runif(domains$sqrt_area[3]))

# Now initialize kernel lists and population vectors
# Use 100 iterations of each model to compute lambdas

n_iterations   <- 100

models         <- vector('list', length(par_sets$site))
names(models)  <- par_sets$site

pop_holder     <- list(
  sqrt_area = array(NA_real_, dim = c(50, 101)),
  ln_leaf_l = array(NA_real_, dim = c(50, 101)),
  d         = array(NA_real_, dim = c(1, 101))
)

pop_holder$sqrt_area[ , 1] <- init_pop_vec$sqrt_area
pop_holder$ln_leaf_l[ , 1] <- init_pop_vec$ln_leaf_l
pop_holder$d[ , 1]         <- 20
pop_holder$lambda          <- array(NA_real_, dim = c(1, 101))

tot_size <- sum(unlist(pop_holder), na.rm = TRUE)

v_holder <- pop_holder <- lapply(
  pop_holder,
  function(x, size) {

    x[ , 1] <- x[ , 1] / size

    return(x)
  },
  size = tot_size
)

v_holder <- lapply(v_holder, t)
v_holder$lambda <- NULL

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

k_seq <- sample(par_sets$site, 100, replace = TRUE)

# Everything is initialized. We'll loop over the populations and iterate
# the model for each one individually, compute lambdas, and store those
# for comparison using test_that()

for(i in seq_len(n_iterations)) {

  pop <- k_seq[i]

  param_name_ind <- grepl(pop, names(full_data_list))

  site_params <- full_data_list[param_name_ind]

  # Next, use each kernel function to generate a (discretized!) vector.
  # transform these into iteration matrices, and then iterate using a loop over
  # n_iterations. We'll also store the complete set of kernels in the `models`
  # object in case we need to debug later.

  kernels <- list(
    kern_xx = PF_xx(all_mesh_p$x_2,
                   all_mesh_p$x_1,
                   site_params) %>%
      matrix(nrow = 50, ncol = 50, byrow = TRUE),
    kern_zx = PF_zx(all_mesh_p$x_2,
                   all_mesh_p$z_2,
                   all_mesh_p$x_1,
                   all_mesh_p$z_1,
                   site_params) %>%
      matrix(nrow = 50, ncol = 50, byrow = TRUE),

    # This requires the subsetting because we only want a single
    # column vector (e.g. the first entry of each row).
    kern_dx = PF_dx(all_mesh_p$x_2,
                   all_mesh_p$x_1,
                   site_params)[seq(1, 2500, by = 50)] %>%
      matrix(nrow = 50, ncol = 1, byrow = TRUE),
    kern_xz = PF_xz(all_mesh_p$x_2,
                   all_mesh_p$z_2,
                   all_mesh_p$x_1,
                   all_mesh_p$z_1,
                   site_params) %>%
      matrix(nrow = 50, ncol = 50, byrow = TRUE),

    # These are simpler to subset because we just want the first row of the
    # matrix.
    kern_xd = PF_xd(all_mesh_p$x_1,
                   site_params)[1:50] %>%
      matrix(nrow = 1, ncol = 50, byrow = TRUE),
    kern_zd = PF_zd(all_mesh_p$z_1,
                   site_params)[1:50]  %>%
      matrix(nrow = 1, ncol = 50, byrow = TRUE)

  )



  # These are the full kernel equations from pages 1853-1854.

  x_t_1 <- kernels$kern_xx %*% pop_holder$ln_leaf_l[ , i] +
    kernels$kern_zx %*% pop_holder$sqrt_area[ , i] +
    kernels$kern_dx %*% pop_holder$d[ , i]

  z_t_1 <- kernels$kern_xz %*% pop_holder$ln_leaf_l[ , i]

  d_t_1 <- kernels$kern_xd %*% pop_holder$ln_leaf_l[ , i] +
    kernels$kern_zd %*% pop_holder$sqrt_area[ , i]


  v_x_t1 <- left_mult(kernels$kern_xx, v_holder$ln_leaf_l[i , ]) +
    left_mult(kernels$kern_xd, v_holder$d[i , ]) +
    left_mult(kernels$kern_xz, v_holder$sqrt_area[i , ])

  v_z_t1 <- left_mult(kernels$kern_zx, v_holder$ln_leaf_l[i , ]) +
    left_mult(kernels$kern_zd, v_holder$d[i , ])

  v_d_t1 <- left_mult(kernels$kern_dx, v_holder$ln_leaf_l[i , ])

  tot_size <- lambda <- sum(x_t_1, z_t_1, d_t_1)

  tot_v    <- sum(v_x_t1, v_z_t1, v_d_t1)

  pop_holder$ln_leaf_l[ , (i + 1)] <- x_t_1 / tot_size
  pop_holder$sqrt_area[ , (i + 1)] <- z_t_1 / tot_size
  pop_holder$d[ , (i + 1)]         <- d_t_1 / tot_size

  v_holder$ln_leaf_l[(i + 1) , ] <- v_x_t1 / tot_v
  v_holder$sqrt_area[(i + 1) , ] <- v_z_t1 / tot_v
  v_holder$d[(i + 1) , ]         <- v_d_t1 / tot_v

  pop_holder$lambda[ (i + 1) ]   <- lambda

  # Store the kernels for subsequent kernel comparisons

  if(rlang::is_empty(models[[pop]])) {
    models[[pop]] <- kernels
  }
}

actual_lambdas <- pop_holder$lambda[ , -1]

# And now, for the ipmr version.


gen_di_stoch_kern <- init_ipm(sim_gen    = "general",
                              di_dd      = "di",
                              det_stoch  = "stoch",
                              'kern') %>%
  define_kernel(
    name             = 'PF_xx_site',
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
    name             = 'PF_zx_site',
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
    name = 'PF_dx_site',
    family = 'DC',

    formula = gamma_nd_site * d_ln_leaf_l,

    gamma_nd_site = dnorm(ln_leaf_l_2, dc_nr_int_site, dc_nr_sd_site),

    data_list = full_data_list,
    states = list(c('ln_leaf_l', 'd')),
    uses_par_sets = TRUE,
    par_set_indices = par_sets,
    evict_cor = TRUE,
    evict_fun = truncated_distributions('norm',
                                        'gamma_nd_site')
  ) %>%
  define_kernel(
    name             = 'PF_xz_site',
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
    name             = 'PF_xd_site',
    family           = 'CD',
    formula          = sig_n_site * mu_n_site * d_ln_leaf_l,
    sig_n_site       = inv_logit(nr_s_z_int_site, nr_s_z_b_site, ln_leaf_l_1),
    mu_n_site        = inv_logit(nr_d_z_int_site, nr_d_z_b_site, ln_leaf_l_1),
    data_list        = full_data_list,
    states           = list(c('ln_leaf_l', 'd')),
    uses_par_sets    = TRUE,
    par_set_indices = par_sets,
    evict_cor        = FALSE
  ) %>%
  define_kernel(
    name             = 'PF_zd_site',
    family           = 'CD',
    formula          = sig_r_site * mu_r_site * d_sqrt_area,
    sig_r_site        = inv_logit(ra_s_z_int_site, ra_s_z_b_site, sqrt_area_1),
    mu_r_site         = inv_logit(ra_d_z_int_site, ra_d_z_b_site, sqrt_area_1),
    data_list        = full_data_list,
    states           = list(c('sqrt_area', 'd')),
    uses_par_sets    = TRUE,
    par_set_indices = par_sets,
    evict_cor        = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c(paste('PF_',
                             c('xx',
                               'zx', 'dx', 'xz', 'xd', 'zd'),
                             '_site',
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
      n_d         = 20
    )
  ) %>%
  make_ipm(
    return_all_envs = TRUE,
    usr_funs   = list(
      inv_logit = inv_logit,
      pois      = pois
    ),
    iterations = 100,
    normalize_pop_size = TRUE,
    kernel_seq = k_seq
  )

# need to rethink how output is generated - no reason I shouldn't be able to
# plug model object straight into lambda generic!

ipmr_lambdas <- lambda(gen_di_stoch_kern,
                       type_lambda = "all") %>%
  as.vector()

test_that('ipmr version matches simulation', {

  expect_equal(ipmr_lambdas, actual_lambdas, tolerance = 1e-10)

  expect_s3_class(gen_di_stoch_kern, "general_di_stoch_kern_ipm")
  expect_s3_class(gen_di_stoch_kern, "ipmr_ipm")

  # Standardized right eigenvectors

  ipmr_w <- right_ev(gen_di_stoch_kern)
  expect_s3_class(ipmr_w, "ipmr_w")

  hand_w <- lapply(pop_holder, function(x) x[ , 26:101, drop = FALSE]) %>%
    .[sort(names(.))]
  hand_w$lambda <- NULL

  ipmr_w <- ipmr_w[sort(names(ipmr_w))]

  expect_equal(ipmr_w, hand_w, tolerance = 1e-10, ignore_attr = TRUE)

  kern_tests <- logical(length(gen_di_stoch_kern$sub_kernels))

  kern_suffs <- c('xx', 'zx', 'xd', 'xz', 'dx', 'dz')

  for(i in seq_along(par_sets$site)) {

    ind <- seq(i,
               length(gen_di_stoch_kern$sub_kernels),
               by = length(par_sets$site))

    kern_tests[ind] <- compare_kernels('gen_di_stoch_kern',
                                       'models',
                                       nms_ipmr = paste('PF_',
                                                        paste(kern_suffs,
                                                              par_sets$site[i],
                                                              sep = "_"),
                                                        sep = ""),
                                       nms_hand = paste(par_sets$site[i],
                                                        '$kern_',
                                                        kern_suffs,
                                                        sep = ""))

  }

  if(any(! kern_tests)) warning('kernel indices: ', which(! kern_tests),
                                ' are not identical')

  expect_true(all(kern_tests))


})

test_that("left_ev works on general_di_stoch_kern IPMs", {

  v_ipmr <- left_ev(gen_di_stoch_kern,
                    iterations = 100,
                    kernel_seq = k_seq)
  expect_s3_class(v_ipmr, "ipmr_v")

  v_ipmr <- v_ipmr[sort(names(v_ipmr))]

  v_hand <- lapply(v_holder, function(x) t(x[26:101, , drop = FALSE])) %>%
    .[sort(names(.))]

  expect_equal(v_ipmr[[1]], v_hand[[1]], ignore_attr = TRUE)
  expect_equal(v_ipmr[[2]], v_hand[[2]], ignore_attr = TRUE)
  expect_equal(v_ipmr[[3]], v_hand[[3]], ignore_attr = TRUE)


})

# Next, we need to test the internal iteration machinery - particularly the
# kern_seq properties. This tests whether or not we can recover stochastic simulations
# given the same sequence of kernels (i.e. did I write the user-defined sequence
# machinery correctly). Since we've already verified that the kernels are all identical,
# there should never be any issues here.

usr_seq <- sample(par_sets$site,
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

temp_proto <- gen_di_stoch_kern$proto_ipm

stoch_mod_ipmr <- temp_proto %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c(paste('PF_',
                             c('xx',
                               'zx', 'dx', 'xz', 'xd', 'zd'),
                             '_site',
                             sep = "")),
      int_rule     = rep('midpoint', 6),
      state_start    = c('ln_leaf_l',
                       'sqrt_area',
                       'd',
                       'ln_leaf_l', 'ln_leaf_l', 'sqrt_area'),
      state_end      = c('ln_leaf_l',
                       'ln_leaf_l', 'ln_leaf_l',
                       'sqrt_area', 'd', 'd')
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
    return_all_envs = TRUE,
    usr_funs   = list(
      inv_logit = inv_logit,
      pois      = pois
    ),
    iterations = 100,
    kernel_seq = usr_seq,
    normalize_pop_size = FALSE
  )


pop_ipmr    <- stoch_mod_ipmr$pop_state
pop_ipmr    <- pop_ipmr[names(pop_ipmr) != 'lambda']

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

  expect_equal(all_lambdas$hand, all_lambdas$ipmr, tolerance = 1e-7)

})


test_that('evict_fun warnings are correctly generated', {

  test_evict_fun <-
    gen_di_stoch_kern$proto_ipm[!grepl('PF_zx_site',
                                       gen_di_stoch_kern$proto_ipm$kernel_id), ]

  wrngs <- capture_warnings(
    test_evict_fun_warning <- test_evict_fun %>%
      define_kernel(
        name             = 'PF_zx_site',
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
        evict_fun = truncated_distributions('norm',
                                            c('gamma_nr_site',
                                              'gamma_sr_site'))
      ) %>%
      define_impl(
        make_impl_args_list(
          kernel_names = c(
            'PF_xx_site',
            'PF_dx_site',
            'PF_xz_site',
            'PF_xd_site',
            'PF_zd_site',
            'PF_zx_site'
          ),
          int_rule     = rep('midpoint', 6),
          state_start    = c('ln_leaf_l',
                           'd',
                           'ln_leaf_l',
                           'ln_leaf_l',
                           'sqrt_area',
                           'sqrt_area'),
          state_end      = c('ln_leaf_l',
                           'ln_leaf_l',
                           'sqrt_area',
                           'd',
                           'd',
                           'ln_leaf_l')
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
        return_all_envs = TRUE,
        usr_funs   = list(
          inv_logit = inv_logit,
          pois      = pois
        ),
        iterations = 100,
        normalize_pop_size = FALSE
      )
    )

  test_text <-
    "length of 'fun' in 'truncated_distributions()' is not equal to length of 'target'. Recycling 'fun'."

  expect_equal(wrngs[1], test_text)

})


test_that('normalize pop vec works', {

  usr_seq <- sample(par_sets$site, size = 100, replace = TRUE)

  gen_di_stoch_kern <- init_ipm(sim_gen    = "general",
                                di_dd      = "di",
                                det_stoch  = "stoch",
                                'kern') %>%
    define_kernel(
      name             = 'PF_xx_site',
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
      name             = 'PF_zx_site',
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
      name = 'PF_dx_site',
      family = 'DC',

      formula = gamma_nd_site * d_ln_leaf_l,

      gamma_nd_site = dnorm(ln_leaf_l_2, dc_nr_int_site, dc_nr_sd_site),

      data_list = full_data_list,
      states = list(c('ln_leaf_l', 'd')),
      uses_par_sets = TRUE,
      par_set_indices = par_sets,
      evict_cor = TRUE,
      evict_fun = truncated_distributions('norm',
                                          'gamma_nd_site')
    ) %>%
    define_kernel(
      name             = 'PF_xz_site',
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
      name             = 'PF_xd_site',
      family           = 'CD',
      formula          = sig_n_site * mu_n_site * d_ln_leaf_l,
      sig_n_site       = inv_logit(nr_s_z_int_site, nr_s_z_b_site, ln_leaf_l_1),
      mu_n_site        = inv_logit(nr_d_z_int_site, nr_d_z_b_site, ln_leaf_l_1),
      data_list        = full_data_list,
      states           = list(c('ln_leaf_l', 'd')),
      uses_par_sets    = TRUE,
      par_set_indices = par_sets,
      evict_cor        = FALSE
    ) %>%
    define_kernel(
      name             = 'PF_zd_site',
      family           = 'CD',
      formula          = sig_r_site * mu_r_site * d_sqrt_area,
      sig_r_site        = inv_logit(ra_s_z_int_site, ra_s_z_b_site, sqrt_area_1),
      mu_r_site         = inv_logit(ra_d_z_int_site, ra_d_z_b_site, sqrt_area_1),
      data_list        = full_data_list,
      states           = list(c('sqrt_area', 'd')),
      uses_par_sets    = TRUE,
      par_set_indices = par_sets,
      evict_cor        = FALSE
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c(paste('PF_',
                               c('xx',
                                 'zx', 'dx', 'xz', 'xd', 'zd'),
                               '_site',
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
    make_ipm(
      return_all_envs = TRUE,
      usr_funs   = list(
        inv_logit = inv_logit,
        pois      = pois
      ),
      iterations = 100,
      normalize_pop_size = TRUE,
      kernel_seq = usr_seq
    )


  lambdas_ipmr <- lambda(gen_di_stoch_kern,
                         type_lambda = 'all') %>%
    as.vector()

  pop_holder <- list(n_leaf_l = array(NA_real_, dim = c(50, 101)),
                     n_sqt_ar = array(NA_real_, dim = c(50, 101)),
                     n_d      = array(NA_real_, dim = c(1, 101)))

  tot_size <- Reduce('sum', init_pop_vec, init = 10)

  pop_holder[[1]][ , 1] <- init_pop_vec$ln_leaf_l / tot_size
  pop_holder[[2]][ , 1] <- init_pop_vec$sqrt_area / tot_size
  pop_holder[[3]][ , 1] <- 10 / tot_size

  lambdas_hand <- numeric(100L)

  for(i in seq_len(100)) {

    selector <- gen_di_stoch_kern$env_seq[i]
    kerns    <- models[[selector]]

    temp_n_ln_t_1 <- kerns$kern_xx %*% pop_holder[[1]][ , i] +
                     kerns$kern_zx %*% pop_holder[[2]][ , i] +
                     kerns$kern_dx %*% pop_holder[[3]][ , i]

    temp_n_sq_t_1 <- kerns$kern_xz %*% pop_holder[[1]][ , i]

    temp_n_do_t_1 <- kerns$kern_xd %*% pop_holder[[1]][ , i] +
                     kerns$kern_zd %*% pop_holder[[2]][ , i]

    tot_size <- Reduce('sum',
                       c(temp_n_sq_t_1, temp_n_do_t_1, temp_n_ln_t_1),
                       init = 0)

    lambdas_hand[i] <- tot_size

    pop_holder[[1]][ , (i + 1)] <- temp_n_ln_t_1 / tot_size
    pop_holder[[2]][ , (i + 1)] <- temp_n_sq_t_1 / tot_size
    pop_holder[[3]][ , (i + 1)] <- temp_n_do_t_1 / tot_size

  }

  expect_equal(lambdas_ipmr, lambdas_hand, tolerance = 1e-10)

  pop_sizes_ipmr <- lapply(gen_di_stoch_kern$pop_state[1:3],
                           function(x) colSums(x)) %>%
    lapply(X = 1:101,
           function(x, sizes) {
             sum(sizes[[1]][x], sizes[[2]][x], sizes[[3]][x])
           },
           sizes = .) %>%
    unlist()

  expect_equal(pop_sizes_ipmr, rep(1, 101), tolerance = 1e-15)

})

mean_kernel_r <- function(kernel_holder, n_unique) {

  kern_nms <- gsub("whitetop_", "", names(kernel_holder)[1:n_unique]) %>%
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

test_that("mean_kernel works for fully par_set models", {

  mean_ipmr_kernels <- mean_kernel(gen_di_stoch_kern)
  names(mean_ipmr_kernels) <- gsub("_site", "", names(mean_ipmr_kernels))
  names(mean_ipmr_kernels) <- gsub("PF", "kern", names(mean_ipmr_kernels))

  kernel_holder <- unlist(models, recursive = FALSE) %>%
    setNames(gsub("\\.", "_", names(.)))

  mean_hand_kernels <- mean_kernel_r(kernel_holder, 6L) %>%
    set_ipmr_classes()

  expect_equal(mean_hand_kernels, mean_ipmr_kernels)

})


test_that('partially stochastic models also work', {

  all_g_int   <- as.list(rnorm(5, mean = 5.781, sd = 0.9))
  all_f_s_int <- as.list(rnorm(5, mean = 2.6204, sd = 0.3))

  names(all_g_int)   <- paste("g_int_", 1:5, sep = "")
  names(all_f_s_int) <- paste("f_s_int_", 1:5, sep = "")

  constant_list <- list(
    g_slope   = 0.988,
    g_sd      = 20.55699,
    s_int     = -0.352,
    s_slope   = 0.122,
    s_slope_2 = -0.000213,
    f_r_int   = -11.46,
    f_r_slope = 0.0835,
    f_s_slope = 0.01256,
    f_d_mu    = 5.6655,
    f_d_sd    = 2.0734,
    e_p       = 0.15,
    g_i       = 0.5067
  )

  all_params <- c(constant_list, all_g_int, all_f_s_int)

  data_list_cr <- list(
    g_int     = 7.229,
    g_slope   = 0.988,
    g_sd      = 21.72262,
    s_int     = 0.0209,
    s_slope   = 0.0831,
    s_slope_2 = -0.00012999,
    f_r_int   = -11.46,
    f_r_slope = 0.0835,
    f_s_int   = 2.6204,
    f_s_slope = 0.01256,
    f_d_mu    = 5.6655,
    f_d_sd    = 2.0734,
    e_p       = 0.15,
    g_i       = 0.5067
  )

  L <- 1.02
  U <- 624
  n <- 500

  set.seed(2312)

  init_pop_vec   <- runif(500)
  init_seed_bank <- 20

  inv_logit <- function(int, slope, sv) {
    1/(1 + exp(-(int + slope * sv)))
  }

  inv_logit_2 <- function(int, slope, slope_2, sv) {
    1/(1 + exp(-(int + slope * sv + slope_2 * sv ^ 2)))
  }


  general_stoch_kern_ipm <- init_ipm(sim_gen    = "general",
                                     di_dd      = "di",
                                     det_stoch  = "stoch",
                                     'kern') %>%
    define_kernel(
      name             = "P_year",
      formula          = s * g_year * d_ht,
      family           = "CC",
      g_year           = dnorm(ht_2, g_mu_year, g_sd),
      g_mu_year        = g_int_year + g_slope * ht_1,
      s                = inv_logit_2(s_int, s_slope, s_slope_2, ht_1),
      data_list        = all_params,
      states           = list(c('ht')),
      uses_par_sets    = TRUE,
      par_set_indices = list(year = 1:5),
      evict_cor        = TRUE,
      evict_fun        = truncated_distributions('norm',
                                                 'g_year')
    ) %>%
    define_kernel(
      name          = "go_discrete_year",
      formula       = f_r * f_s_year * g_i,
      family        = 'CD',
      f_r           = inv_logit(f_r_int, f_r_slope, ht_1),
      f_s_year      = exp(f_s_int_year + f_s_slope * ht_1),
      data_list     = all_params,
      states        = list(c('ht', 'b')),
      uses_par_sets    = TRUE,
      par_set_indices = list(year = 1:5)
    ) %>%
    define_kernel(
      name    = 'stay_discrete',
      formula = 0,
      family  = "DD",
      states  = list(c('b')),
      uses_par_sets = FALSE,
      evict_cor = FALSE
    ) %>%
    define_kernel(
      name          = 'leave_discrete',
      formula       = e_p * f_d * d_ht,
      f_d           = dnorm(ht_2, f_d_mu, f_d_sd),
      family        = 'DC',
      data_list     = all_params,
      states        = list(c('ht', 'b')),
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'f_d')
    )  %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P_year", "go_discrete_year", "stay_discrete", "leave_discrete"),
        int_rule     = c(rep("midpoint", 4)),
        state_start    = c('ht', "ht",'b', 'b'),
        state_end      = c('ht', 'b', 'b', 'ht')
      )
    ) %>%
    define_domains(
      ht = c(L, U, n)
    ) %>%
    define_pop_state(
      n_ht = init_pop_vec,
      n_b  = init_seed_bank
    ) %>%
    make_ipm(iterations = 3,
             kernel_seq = sample(1:5, size = 3, replace = TRUE),
             usr_funs = list(inv_logit   = inv_logit,
                             inv_logit_2 = inv_logit_2),
             return_all_envs = TRUE)

  states <- list(c('ht', 'b'))

  det_version <- init_ipm(sim_gen    = "general",
                          di_dd      = "di",
                          det_stoch  = "det") %>%
    define_kernel(
      name          = "P",
      formula       = s * g * d_ht,
      family        = "CC",
      g             = dnorm(ht_2, g_mu, g_sd),
      g_mu          = g_int + g_slope * ht_1,
      s             = inv_logit_2(s_int, s_slope, s_slope_2, ht_1),
      data_list     = data_list_cr,
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
      data_list     = data_list_cr,
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
      data_list     = data_list_cr,
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
        state_start    = c('ht', "ht", 'b', 'b'),
        state_end      = c('ht', 'b', 'b', 'ht')
      )
    ) %>%
    define_domains(
      ht = c(1.02, 624, 500)
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
             return_all_envs = TRUE)

  leave_discrete_stoch <- general_stoch_kern_ipm$sub_kernels$leave_discrete
  leave_discrete_det   <- det_version$sub_kernels$leave_discrete

  # if partially par_setarchical models work the way they should,
  # then there should never be an difference between the deterministic model
  # version of this kernel and the stochastic model version of it!

  expect_equal(leave_discrete_det, leave_discrete_stoch)

  mean_kernel_r <- function(kernel_holder) {

    kern_nms <- gsub("_[0-9]", "", names(kernel_holder)) %>%
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

  sub_kerns <- general_stoch_kern_ipm$sub_kernels

  mean_hand_kerns <- mean_kernel_r(sub_kerns)

  mean_ipmr_kerns <- mean_kernel(general_stoch_kern_ipm) %>%
    lapply(unclass)

  names(mean_ipmr_kerns) <- gsub("_year", "", names(mean_ipmr_kerns))

  mean_hand_kerns <- mean_hand_kerns[sort(names(mean_hand_kerns))]
  mean_ipmr_kerns <- mean_ipmr_kerns[sort(names(mean_ipmr_kerns))]

  expect_equal(mean_hand_kerns, mean_ipmr_kerns)


})
