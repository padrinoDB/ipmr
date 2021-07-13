# test simple_di_stoch_param
library(mvtnorm)
library(rlang)


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

iterations <- 100

pop_vec <- matrix(0, ncol = iterations + 1, nrow = length(d1))

init_pop_vec <- runif(length(d1), 0, 10)
pop_vec[ ,1] <- init_pop_vec <- init_pop_vec / sum(init_pop_vec)

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

test_stoch_param <- init_ipm(sim_gen    = "simple",
                             di_dd      = "di",
                             det_stoch  = "stoch",
                             kern_param = "param") %>%
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
           iterations = 100,
           normalize_pop_size = TRUE)

pop_state_ipmr <- test_stoch_param$pop_state$n_surf_area
lambda_ipmr <- lambda(test_stoch_param,
                      comp_method = 'pop_size',
                      type_lambda = 'all') %>%
  as.vector

# Now, use the env_seq to plug into this loop for each iteration and see if
# lambdas are identical. If not, then find a bridge and jump

lambda <- vector('numeric', iterations)

ks <- list()
ps <- list()
fs <- list()

domains <- expand.grid(list(d2 = d1, d1 = d1))

v_holder <- t(pop_vec)

for(i in seq_len(iterations)) {

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
  v_t_1 <- left_mult(K_temp, v_holder[i, ])
  v_t_1 <- v_t_1 / sum(v_t_1)

  lambda[i] <- sum(n_t_1)

  pop_vec[ , (i + 1)] <- n_t_1 / sum(n_t_1)

  v_holder[(i + 1), ] <- v_t_1

}

ws <- list(pop_vec[ , 26:101])

ws_ipmr <- right_ev(test_stoch_param)

vs <- list(t(v_holder[26:101 , ]))
vs_ipmr <- left_ev(test_stoch_param, iterations = 100, burn_in = 0.25)

test_that('eigenvectors/values are all good', {

  # Deterministic lambdas
  expect_equal(lambda_ipmr, lambda, tolerance = 1e-10)

  # Stable stage distributions
  expect_equal(ws_ipmr, ws, ignore_attr = TRUE)

  # Reproductive values
  expect_equal(vs_ipmr, vs, ignore_attr = TRUE)

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


test_that('normalize pop_vectors works as it should', {

  test_stoch_param <- init_ipm(sim_gen    = "simple",
                               di_dd      = "di",
                               det_stoch  = "stoch",
                               kern_param = "param") %>%
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
             iterations = 10,
             normalize_pop_size = TRUE)


  lambdas_norm <- lambda(test_stoch_param,
                         comp_method = 'pop_size',
                         type_lambda = 'all') %>%
    as.vector()

  pop_holder       <- array(NA_real_, dim = c(100, 11))
  pop_holder[ , 1] <- init_pop_vec / sum(init_pop_vec)
  lambdas_hand     <- numeric(10L)

  for(i in seq_len(10)) {

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

    n_t_1 <- K_temp %*% pop_holder[ , i]

    lambdas_hand[i] <- sum(n_t_1)
    pop_holder[ , (i + 1)] <- n_t_1 / sum(n_t_1)

  }


  expect_equal(lambdas_norm, lambdas_hand, tolerance = 1e-10)
  expect_equal(test_stoch_param$pop_state$n_surf_area,
               pop_holder)

  pop_sizes <- colSums(test_stoch_param$pop_state$n_surf_area)

  expect_equal(pop_sizes, rep(1, 11), tolerance = 1e-15)

})


test_that("t helper variable works as advertised", {

  env_state <- test_stoch_param$env_seq

  env_sampler <- function(environ_seq, iteration) {

    temp <- as.list(environ_seq[iteration, ]) %>%
      setNames(
        c('s_int_yr',
          'g_int_yr',
          'f_r_int_yr',
          'f_s_int_yr')
      )

    return(temp)
  }

  test_det_seq <- init_ipm(sim_gen    = "simple",
                           di_dd      = "di",
                           det_stoch  = "stoch",
                           kern_param = "param") %>%
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
    )%>%
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
      env_params = env_sampler(environ_seq = env_state,
                               iteration   = t),
      data_list = list(
        env_state   = env_state,
        env_sampler = env_sampler
      )
    ) %>%
    define_pop_state(
      pop_vectors = list(n_surf_area = init_pop_vec),
    ) %>%
    make_ipm(usr_funs = list(inv_logit = inv_logit),
             iterate = TRUE,
             iterations = 100,
             normalize_pop_size = TRUE)

  det_seq_l <- lambda(test_det_seq)
  ran_seq_l <- lambda(test_stoch_param)

  expect_equal(det_seq_l, ran_seq_l, tolerance = 1e-10)

  pop_state_det <- test_det_seq$pop_state$n_surf_area[ , 11]
  pop_state_ran <- test_stoch_param$pop_state$n_surf_area[ , 11]

  expect_equal(pop_state_det, pop_state_ran, tolerance = 1e-10)

})



test_that("stoch_param can handle pararameter sets", {

  par_set_vals <- rnorm(5, 0, 2) %>% as.list
  names(par_set_vals) <- paste("s_int", 2000:2004, sep = "_")


  data_list <- list(s_slope = 0.2,
                    s_int   = -0.8,
                    g_slope = 0.99,
                    g_sd = 0.2,
                    f_r_slope = 0.003,
                    f_s_slope = 0.01,
                    f_d_mu = 2,
                    f_d_sd = 0.75)

  data_list <- c(data_list, par_set_vals)

  b <- seq(0, 10, length.out = 101)
  d1 <- d2 <- (b[2:101] + b[1:100]) * 0.5
  h <- d1[3] - d1[2]

  r_means <- c(g_int_yr = 0.1,
               f_r_int_yr = 0.3,
               f_s_int_yr = 0.01)

  # These are likely nonsense values - using for testing purposes only!
  r_sigma <- runif(9)
  r_sigma <- matrix(r_sigma, nrow = 3)

  r_sigma <- r_sigma %*% t(r_sigma) # make symmetrical

  iterations <- 10

  pop_vec <- matrix(0, ncol = iterations + 1, nrow = length(d1))

  pop_vec[ ,1] <- init_pop_vec <- runif(length(d1), 0, 10)

  ## IPMR Version

  mvt_wrapper <- function(r_means, r_sigma, nms) {
    out <- rmvnorm(1, r_means, r_sigma) %>%
      as.list()

    names(out) <- nms
    return(out)
  }

  test_stoch_param <- init_ipm(sim_gen    = "simple",
                               di_dd      = "di",
                               det_stoch  = "stoch",
                               kern_param = "param") %>%
    define_kernel(
      'P_yr',
      formula = s_yr * g,
      family = 'CC',
      g_mu = g_int_r + g_slope * surf_area_1,
      s_yr = plogis(s_int + s_int_yr + s_slope * surf_area_1),
      g = dnorm(surf_area_2, g_mu, g_sd),
      data_list = data_list,
      states = list(c('surf_area')),
      uses_par_sets = TRUE,
      par_set_indices = list(yr = 2000:2004),
      evict_cor = TRUE,
      evict_fun = truncated_distributions('norm', 'g')
    ) %>%
    define_kernel(
      'F',
      formula = f_r * f_s * f_d,
      family = 'CC',
      f_r = plogis(f_r_int_r + f_r_slope * surf_area_1),
      f_s = exp(f_s_int_r + f_s_slope * surf_area_1),
      f_d = dnorm(surf_area_2, f_d_mu, f_d_sd),
      data_list = data_list,
      states = list(c('surf_area')),
      uses_par_sets = FALSE,
      evict_cor = TRUE,
      evict_fun = truncated_distributions('norm', 'f_d')
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c('P_yr', "F"),
        int_rule = rep('midpoint', 2),
        state_start = rep('surf_area', 2),
        state_end = rep('surf_area', 2)
      )
    ) %>%
    define_domains(surf_area = c(0, 10, 100)) %>%
    define_env_state(
      env_params = mvt_wrapper(r_means, r_sigma, nms = c('g_int_r',
                                                         'f_r_int_r',
                                                         'f_s_int_r')),
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
             normalize_pop_size = FALSE,
             kernel_seq = sample(2000:2004, iterations, TRUE))

  hand_lambda <- vector('numeric', iterations)

  ks <- list()
  ps <- list()
  fs <- list()

  s_r <- function(dom, params) {

    plogis(params[1] + params[2] + params[3] * dom)

  }

  domains <- expand.grid(list(d2 = d1, d1 = d1))

  temp_data <- data_list[!grepl('[0-9]', names(data_list))]

  for(i in seq_len(iterations)) {

    r_ests <- test_stoch_param$env_seq[i , 1:3] %>%
      as.list() %>%
      setNames(names(r_means))

    yr <- test_stoch_param$env_seq[i , 4]

    temp <- c(data_list,  r_ests)
    s_r_yr <- unlist(par_set_vals[grepl(yr, names(par_set_vals))])

    g_mat <- g(domains$d2, domains$d1,
               params = c(temp$g_int_yr,
                          temp$g_slope,
                          temp$g_sd),
               L = 0,
               U = 10)

    s_vec <- s_r(d1, params = c(temp$s_int,
                                s_r_yr,
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

    hand_lambda[i] <- sum(n_t_1)/sum(pop_vec[ , i])

    pop_vec[ , (i + 1)] <- n_t_1

  }

  ipmr_lambda <- lambda(test_stoch_param, type_lambda = "all") %>%
    as.vector()

  expect_equal(hand_lambda, ipmr_lambda)

  ipmr_ps <- test_stoch_param$pop_state$n_surf_area

  expect_equal(pop_vec, ipmr_ps)

})
