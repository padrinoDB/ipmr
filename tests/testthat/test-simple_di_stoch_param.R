# test simple_di_stoch_param
library(mvtnorm)
library(rlang)
library(ipmr)

context('Simple density independent stochastic parameter resampled models')

s <- function(sv1, params) {
  1/(1 + exp(-(params[1] + params[2] * sv1)))
}


g <- function(sv1, sv2, params) {
  mu <- params[1] + params[2] * sv1
  dnorm(sv2, mean = mu, sd = params[3])
}

f_r <- function(sv1, params) {
  1/(1 + exp(-(params[1] + params[2] * sv1)))
}

f_s <- function(sv1, params) {
  exp(params[1] + params[2] * sv1)
}

f_d <- function(sv2, params) {
  dnorm(sv2, mean = params[1], sd = params[2])
}

fec <- function(sv1, sv2, params) {
  f_r(sv1, params[1:2]) * f_s(sv1, params[3:4]) * f_d(sv2, params[5:6])
}
update_r_effs <- function(to_add, r_eff_list) {
  purrr::map2(r_eff_list, to_add, .f = function(.x, .y) c(.y, .x))
}

data_list <- list(s_slope = 0.2,
                  g_slope = 0.99,
                  g_sd = 0.2,
                  f_r_slope = 0.003,
                  f_s_slope = 0.0001,
                  f_d_mu = 0.09,
                  f_d_sd = 0.01)

b <- seq(0, 10, length.out = 101)
d1 <- d2 <- (b[2:101] + b[1:100]) * 0.5
h <- d1[3] - d1[2]

r_means <- c(s_int_yr = 0.8,
             g_int_yr = 0.1,
             f_r_int_yr = 0.3,
             f_s_int_yr = 0.004)

set.seed(5000)

# These are likely nonsense values - using for testing purposes only!
r_sigma <- runif(16)
r_sigma <- matrix(r_sigma, nrow = 4)

r_sigma <- r_sigma %*% t(r_sigma) # make symmetrical

iterations <- 50


seeds <- sample(1:iterations, iterations, replace = TRUE)

pop_vec <- matrix(0, ncol = iterations + 1, nrow = length(d1))

pop_vec[ ,1] <- init_pop_vec <- runif(length(d1), 0, 10)
lambda <- vector('numeric', iterations)

for(i in seq_len(iterations)) {

  set.seed(seeds[i])

  r_ests <- rmvnorm(1, mean = r_means, sigma = r_sigma) %>%
    as.list() %>%
    setNames(names(r_means))

  temp <- purrr::splice(data_list,  r_ests)

  g_mat <- h * outer(d1, d2, FUN = g,
                    params = c(temp$g_int_yr,
                               temp$g_slope,
                               temp$g_sd))

  s_vec <- s(d1, params = c(temp$s_int_yr,
                            temp$s_slope))

  g_mat <- truncated_distributions(g_mat, 100)
  P <- t(s_vec * t(g_mat))

  F_mat <- h * outer(d1, d2, FUN = fec,
                     params = c(temp$f_r_int_yr,
                                temp$f_r_slope,
                                temp$f_s_int_yr,
                                temp$f_s_slope,
                                temp$f_d_mu,
                                temp$f_d_sd))

  K_temp <- P + F_mat

  n_t_1 <- K_temp %*% pop_vec[ , i]

  lambda[i] <- sum(n_t_1)/sum(pop_vec[ , i])

  pop_vec[ , (i + 1)] <- n_t_1

}


## IPMR Version

data_list <- list(s_slope = 0.2,
                  g_slope = 0.99,
                  g_sd = 0.2,
                  f_r_slope = 0.003,
                  f_s_slope = 0.0001,
                  f_d_mu = 0.09,
                  f_d_sd = 0.01)

inv_logit <- function(int, slope, sv1) {
  1/(1 + exp(-(int + slope * sv1)))
}

mvt_wrapper <- function(r_means, r_sigma, nms) {
  out <- rmvnorm(1, r_means, r_sigma) %>%
    as.list()

  names(out) <- nms
  return(out)
}

test_stoch_param <- init_ipm('simple_di_stoch_param') %>%
  define_kernel(
    'P',
    formula = s_g_mult(s, g),
    family = 'CC',
    g_mu = env_params$g_int_yr + g_slope * surf_area_1,
    s = inv_logit(env_params$s_int_yr, s_slope, surf_area_1),
    g = cell_size_surf_area * dnorm(surf_area_2, g_mu, g_sd),
    data_list = data_list,
    states = list(c('surf_area')),
    has_hier_effs = FALSE,
    evict = TRUE,
    evict_fun = truncated_distributions(g,
                                        100)
  ) %>%
  define_kernel(
    'F',
    formula = f_r * f_s * f_d,
    family = 'CC',
    f_r = inv_logit(env_params$f_r_int_yr, f_r_slope, surf_area_1),
    f_s = exp(env_params$f_s_int_yr + f_s_slope * surf_area_1),
    f_d = cell_size_surf_area * dnorm(surf_area_2, f_d_mu, f_d_sd),
    data_list = data_list,
    states = list(c('surf_area')),
    has_hier_effs = FALSE,
    evict = FALSE
  ) %>%
  define_k(
    'K',
    n_surf_area_t_1 = right_mult((P+F), n_surf_area_t),
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
           iterations = 10)


