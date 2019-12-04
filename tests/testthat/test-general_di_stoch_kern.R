if(rlang::is_interactive()) {
  library(testthat)
}
context('general density independent stochastic kernel models')

# Test with Aikens & Roach 2014 Population dynamics in central and edge
# populations of a narrowly endemic plant

library(rlang)
library(purrr)

flatten_to_depth <- ipmr:::.flatten_to_depth

hier_effs <- list(
  pop = c(
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

gamma_nd <- function(x_2, params, L, U) {

  ev <- pnorm(U, x_2, params[[2]]) - pnorm(L, x_2, params[[2]])

  dnorm(x_2, params[[1]], params[[2]]) / ev

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

  gamma_nd(x_2, params[24:25], L_x, U_x) * d_x

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

models         <- vector('list', length(hier_effs$pop))
names(models)  <- hier_effs$pop

pop_holders        <- vector('list', length(hier_effs$pop))
names(pop_holders) <- hier_effs$pop

# All populations will get the same initial population vector, and it's
# just a random set of numbers in their range of possible values

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
    kern_dx = k_dx(all_mesh_p$x_2,
                   all_mesh_p$x_1,
                   site_params) %>%
      matrix(nrow = 50, ncol = 1, byrow = TRUE),
    kern_xz = k_xz(all_mesh_p$x_2,
                   all_mesh_p$z_2,
                   all_mesh_p$x_1,
                   all_mesh_p$z_1,
                   site_params) %>%
      matrix(nrow = 50, ncol = 50, byrow = TRUE),
    kern_xd = k_xd(all_mesh_p$x_1,
                   site_params) %>%
      matrix(nrow = 1, ncol = 50, byrow = TRUE),
    kern_zd = k_zd(all_mesh_p$z_1,
                   site_params)  %>%
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


test_that('simulation can somewhat recover the kernels', {

  expect_equal(target_lambdas, actual_lambdas, tol = 3e-2)

})


