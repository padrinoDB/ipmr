# test simple_di_stoch_param
library(mvtnorm)
library(rlang)

context('Simple density independent stochastic parameter resampled models')

set.seed(25129)

s <- function(sv1, params) {
  1/(1 + exp(-(params[1] + params[2] * sv1)))
}


g <- function(sv1, sv2, params, L, U) {
  mu <- params[1] + params[2] * sv1
  ev <- pnorm(U, mu, params[3]) - pnorm(L, mu, params[3])
  dnorm(sv2, mean = mu, sd = params[3]) / ev
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

fec <- function(sv1, sv2, params, L, U) {
  ev <- pnorm(U, params[5], params[6]) - pnorm(L, params[5], params[6])
  f_r(sv1, params[1:2]) * f_s(sv1, params[3:4]) * (f_d(sv2, params[5:6]) / ev)
}


data_list <- list(s_slope = 0.2,
                  g_slope = 0.99,
                  g_sd = 0.2,
                  f_r_slope = 0.003,
                  f_s_slope = 0.01,
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

data_list <- list(s_slope = 0.2,
                  g_slope = 0.99,
                  g_sd = 0.2,
                  f_r_slope = 0.003,
                  f_s_slope = 0.01,
                  f_d_mu = 2,
                  f_d_sd = 0.75)

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
           iterations = 10)

pop_state_ipmr <- test_stoch_param$pop_state$pop_state_surf_area
lambda_ipmr <- lambda(test_stoch_param)



# Now, use the env_seq to plug into this loop for each iteration and see if
# lambdas are identical. If not, then find a bridge and jump

lambda <- vector('numeric', iterations)

ks <- list()
ps <- list()
fs <- list()

domains <- expand.grid(list(d2 = d1, d1 = d1))

for(i in seq_len(iterations)) {

  r_ests <- test_stoch_param$env_seq[i, ] %>%
    as.list() %>%
    setNames(names(r_means))

  temp <- purrr::splice(data_list,  r_ests)

  g_mat <- g(domains$d2, domains$d1,
             params = c(temp$g_int_yr,
                        temp$g_slope,
                        temp$g_sd),
             L = 0,
             U = 10)

  s_vec <- s(d1, params = c(temp$s_int_yr,
                            temp$s_slope))

  P <- (t(s_vec * t(g_mat)) * h) %>%
    matrix(nrow = 100, ncol = 100, byrow = TRUE)

  F_mat <- h * fec(domains$d2, domains$d1,
                   params = c(temp$f_r_int_yr,
                              temp$f_r_slope,
                              temp$f_s_int_yr,
                              temp$f_s_slope,
                              temp$f_d_mu,
                              temp$f_d_sd),
                   L = 0,
                   U = 10) %>%
    matrix(nrow = 100, ncol = 100, byrow = TRUE)

  K_temp <- P + F_mat

  nm_k <- paste0('K_', i, sep = "")
  nm_p <- paste0('P_', i, sep = "")
  nm_f <- paste0('F_', i, sep = "")
  ks[[i]] <- K_temp
  ps[[i]] <- P
  fs[[i]] <- F_mat
  names(ks)[i] <- nm_k
  names(ps)[i] <- nm_p
  names(fs)[i] <- nm_f

  n_t_1 <- K_temp %*% pop_vec[ , i]

  lambda[i] <- sum(n_t_1)/sum(pop_vec[ , i])

  pop_vec[ , (i + 1)] <- n_t_1

}

ws <- vapply(ks, function(x) Re(eigen(x)$vectors[ , 1]), numeric(100L))
ws_ipmr <- vapply(test_stoch_param$iterators,
                  function(x) Re(eigen(x)$vectors[ , 1]),
                  numeric(100L))

vs <- vapply(ks, function(x) Re(eigen(t(x))$vectors[ , 1]), numeric(100L))
vs_ipmr <- vapply(test_stoch_param$iterators,
                  function(x) Re(eigen(t(x))$vectors[ , 1]),
                  numeric(100L))

test_that('eigenvectors/values are all good', {

  # Deterministic lambdas
  expect_equal(lambda_ipmr, lambda, tolerance = 1e-10)

  # Stable stage distributions
  expect_equal(ws_ipmr[ ,1], ws[ ,1], tolerance = 1e-13)
  expect_equal(ws_ipmr[ ,2], ws[ ,2], tolerance = 1e-13)
  expect_equal(ws_ipmr[ ,3], ws[ ,3], tolerance = 1e-13)
  expect_equal(ws_ipmr[ ,4], ws[ ,4], tolerance = 1e-13)
  expect_equal(ws_ipmr[ ,5], ws[ ,5], tolerance = 1e-13)
})


test_that('classes are correctly set', {

  k_cls <- vapply(test_stoch_param$sub_kernels,
                  function(x) class(x)[1],
                  character(1L))

  expect_true(all(k_cls == 'ipmr_matrix'))

  sub_cls <- vapply(test_stoch_param$sub_kernels,
                    function(x) class(x)[1],
                    character(1L))

  expect_true(all(sub_cls == 'ipmr_matrix'))
  expect_s3_class(test_stoch_param, 'simple_di_stoch_param_ipm')


})
