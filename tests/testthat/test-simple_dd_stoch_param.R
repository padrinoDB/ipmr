
library(mvtnorm)


set.seed(25129)

s <- function(sv1, params, pop_size) {
  1/(1 + exp(-(params[1] + params[2] * sv1 + params[3] * pop_size)))
}


g <- function(sv1, sv2, params, L, U) {
  mu <- params[1] + params[2] * sv1
  ev <- pnorm(U, mu, params[3]) - pnorm(L, mu, params[3])
  dnorm(sv2, mean = mu, sd = params[3]) / ev
}

f_r <- function(sv1, params) {
  1/(1 + exp(-(params[1] + params[2] * sv1)))
}

f_s <- function(sv1, params, pop_size) {
  exp(params[1] + params[2] * sv1 + params[3] * pop_size)
}

f_d <- function(sv2, params) {
  dnorm(sv2, mean = params[1], sd = params[2])
}

fec <- function(sv1, sv2, params, L, U, pop_size) {
  ev <- pnorm(U, params[6], params[7]) - pnorm(L, params[6], params[7])
  f_r(sv1, params[1:2]) * f_s(sv1, params[3:5], pop_size) * (f_d(sv2, params[6:7]) / ev)
}


data_list <- list(s_slope = 0.2,
                  s_dd    = -0.005,
                  g_slope = 0.99,
                  g_sd = 0.2,
                  f_r_slope = 0.003,
                  f_s_slope = 0.01,
                  f_s_dd    = -0.0002,
                  f_d_mu = 2,
                  f_d_sd = 0.75)

b <- seq(0, 10, length.out = 101)
d1 <- d2 <- (b[2:101] + b[1:100]) * 0.5
h <- d1[3] - d1[2]

r_means <- c(s_int_yr = 0.8,
             g_int_yr = 0.1,
             f_r_int_yr = 0.3,
             f_s_int_yr = 0.01)

# These are likely nonsense values - using for testing purposes only!
r_sigma <- runif(16)
r_sigma <- matrix(r_sigma, nrow = 4)

r_sigma <- r_sigma %*% t(r_sigma) # make symmetrical

iterations <- 10

pop_vec <- matrix(0, ncol = iterations + 1, nrow = length(d1))

pop_vec[ ,1] <- init_pop_vec <- runif(length(d1), 0, 10)

## IPMR Version

inv_logit <- function(int, slope, sv1) {
  1/(1 + exp(-(int + slope * sv1)))
}

mvt_wrapper <- function(r_means, r_sigma, nms) {
  out <- rmvnorm(1, r_means, r_sigma) %>%
    as.list()

  names(out) <- nms
  return(out)
}

test_stoch_param <- init_ipm(sim_gen    = "simple",
                             di_dd      = "dd",
                             det_stoch  = "stoch",
                             "param") %>%
  define_kernel(
    'P',
    formula = s * g,
    family = 'CC',
    g_mu = g_int_yr + g_slope * surf_area_1,
    s = plogis(s_int_yr +  s_slope * surf_area_1 + s_dd * sum(n_surf_area_t)),
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
    f_r = plogis(f_r_int_yr + f_r_slope * surf_area_1),
    f_s = exp(f_s_int_yr + f_s_slope * surf_area_1 + f_s_dd * sum(n_surf_area_t)),
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
  make_ipm(usr_funs = list(mvt_wrapper = mvt_wrapper),
           iterate = TRUE,
           iterations = 10,
           normalize_pop_size = FALSE)

pop_state_ipmr <- test_stoch_param$pop_state$n_surf_area
ipmr_lam       <- lambda(test_stoch_param, type_lambda = "all") %>%
  as.vector()

# Hand version
hand_lam <- vector('numeric', iterations)

ks <- list()
ps <- list()
fs <- list()

domains <- expand.grid(list(d2 = d1, d1 = d1))

for(i in seq_len(iterations)) {

  pop_size <- sum(pop_vec[ , i])

  r_ests <- test_stoch_param$env_seq[i, ] %>%
    as.list() %>%
    setNames(names(r_means))

  temp <- c(data_list,  r_ests)

  g_mat <- g(domains$d2, domains$d1,
             params = c(temp$g_int_yr,
                        temp$g_slope,
                        temp$g_sd),
             L = 0,
             U = 10)

  s_vec <- s(d1, params = c(temp$s_int_yr,
                            temp$s_slope,
                            temp$s_dd),
             pop_size)

  P <- (t(s_vec * t(g_mat)) * h) %>%
    matrix(nrow = 100, ncol = 100, byrow = TRUE)

  F_mat <- h * fec(domains$d2, domains$d1,
                   params = c(temp$f_r_int_yr,
                              temp$f_r_slope,
                              temp$f_s_int_yr,
                              temp$f_s_slope,
                              temp$f_s_dd,
                              temp$f_d_mu,
                              temp$f_d_sd),
                   L = 0,
                   U = 10,
                   pop_size) %>%
    matrix(nrow = 100, ncol = 100, byrow = TRUE)

  K_temp <- P + F_mat

  n_t_1 <- K_temp %*% pop_vec[ , i]

  hand_lam[i] <- sum(n_t_1)/sum(pop_vec[ , i])

  pop_vec[ , (i + 1)] <- n_t_1

}


test_that("math works", {

  expect_equal(ipmr_lam, hand_lam)
  expect_equal(pop_state_ipmr, pop_vec)

})

# Hierarchical models -------------

g_int_effs <- rnorm(3) %>%
  as.list %>%
  setNames(paste("g_int_yr_", 1:3, sep = ''))

f_r_int_effs <- rnorm(3) %>%
  as.list %>%
  setNames(paste("f_r_int_yr_", 1:3, sep = ''))

data_list <- c(data_list, g_int_effs, f_r_int_effs)


test_stoch_param <- init_ipm(sim_gen    = "simple",
                             di_dd      = "dd",
                             det_stoch  = "stoch",
                             "param") %>%
  define_kernel(
    'P_site',
    formula = s * g_site,
    family = 'CC',
    g_mu_site = g_int_yr + g_int_yr_site + g_slope * surf_area_1,
    s = plogis(s_int_yr +  s_slope * surf_area_1 + s_dd * sum(n_surf_area_t)),
    g_site = dnorm(surf_area_2, g_mu_site, g_sd),
    data_list = data_list,
    states = list(c('surf_area')),
    uses_par_sets = TRUE,
    par_set_indices = list(site = 1:3),
    evict_cor = TRUE,
    evict_fun = truncated_distributions('norm', 'g_site')
  ) %>%
  define_kernel(
    'F_site',
    formula = f_r_site * f_s * f_d,
    family = 'CC',
    f_r_site = plogis(f_r_int_yr + f_r_int_yr_site + f_r_slope * surf_area_1),
    f_s = exp(f_s_int_yr + f_s_slope * surf_area_1 + f_s_dd * sum(n_surf_area_t)),
    f_d = dnorm(surf_area_2, f_d_mu, f_d_sd),
    data_list = data_list,
    states = list(c('surf_area')),
    uses_par_sets = TRUE,
    par_set_indices = list(site = 1:3),
    evict_cor = TRUE,
    evict_fun = truncated_distributions('norm', 'f_d')
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c('P_site', "F_site"),
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
  make_ipm(usr_funs = list(mvt_wrapper = mvt_wrapper),
           iterate = TRUE,
           iterations = 10,
           kernel_seq = sample(1:3, 10, TRUE),
           normalize_pop_size = FALSE)


# Hand version

iterations <- 10

pop_vec <- matrix(0, ncol = iterations + 1, nrow = length(d1))

pop_vec[ ,1] <- test_stoch_param$pop_state$n_surf_area[ , 1]

g <- function(sv1, sv2, params, L, U) {
  mu <- params[1] + params[2] * sv1 + params[4]
  ev <- pnorm(U, mu, params[3]) - pnorm(L, mu, params[3])
  dnorm(sv2, mean = mu, sd = params[3]) / ev
}

fec <- function(sv1, sv2, params, L, U, pop_size) {
  ev <- pnorm(U, params[6], params[7]) - pnorm(L, params[6], params[7])
  f_r(sv1, params[c(1:2,8)]) *
    f_s(sv1, params[3:5], pop_size) *
    (f_d(sv2, params[6:7]) / ev)
}

f_r <- function(sv1, params) {
  1/(1 + exp(-(params[1] + params[2] * sv1 + params[3])))
}

hand_lam <- vector('numeric', iterations)

ks <- list()
ps <- list()
fs <- list()

domains <- expand.grid(list(d2 = d1, d1 = d1))

for(i in seq_len(iterations)) {

  pop_size <- sum(pop_vec[ , i])

  r_ests <- test_stoch_param$env_seq[i, 1:4] %>%
    as.list() %>%
    setNames(names(r_means))

  temp <- c(data_list,  r_ests)

  yr <- test_stoch_param$env_seq[i , 5]

  g_yr_eff <- g_int_effs[grepl(yr, names(g_int_effs))] %>%
    unlist()
  f_r_yr_eff <- f_r_int_effs[grepl(yr, names(f_r_int_effs))] %>%
    unlist()

  g_mat <- g(domains$d2, domains$d1,
             params = c(temp$g_int_yr,
                        temp$g_slope,
                        temp$g_sd,
                        g_yr_eff),
             L = 0,
             U = 10)

  s_vec <- s(d1, params = c(temp$s_int_yr,
                            temp$s_slope,
                            temp$s_dd),
             pop_size)

  P <- (t(s_vec * t(g_mat)) * h) %>%
    matrix(nrow = 100, ncol = 100, byrow = TRUE)

  F_mat <- h * fec(domains$d2, domains$d1,
                   params = c(temp$f_r_int_yr,
                              temp$f_r_slope,
                              temp$f_s_int_yr,
                              temp$f_s_slope,
                              temp$f_s_dd,
                              temp$f_d_mu,
                              temp$f_d_sd,
                              f_r_yr_eff),
                   L = 0,
                   U = 10,
                   pop_size) %>%
    matrix(nrow = 100, ncol = 100, byrow = TRUE)

  K_temp <- P + F_mat

  n_t_1 <- K_temp %*% pop_vec[ , i]

  hand_lam[i] <- sum(n_t_1)/sum(pop_vec[ , i])

  pop_vec[ , (i + 1)] <- n_t_1

}

hand_pop_sizes <- colSums(pop_vec)
ipmr_pop_sizes <- colSums(test_stoch_param$pop_state$n_surf_area)

hand_lams <- hand_pop_sizes[2:11] / hand_pop_sizes[1:10]
ipmr_lams <- lambda(test_stoch_param, type_lambda = 'all') %>%
  as.vector()

test_that("parameter resampled par_set models work", {


  expect_equal(hand_pop_sizes, ipmr_pop_sizes)
  expect_equal(hand_lams, ipmr_lams)

  expect_s3_class(test_stoch_param, "simple_dd_stoch_param_ipm")
  expect_s3_class(test_stoch_param, "ipmr_ipm")

})




