
# simple-ndd-det IPM test

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

s <- function(sv1, params) {
  1/(1 + exp(-(params[1] + params[2] * sv1)))
}


g <- function(sv1, sv2, params, L, U) {
  mu <- params[1] + params[2] * sv1
  ev <- (pnorm(U, mu, params[3]) - pnorm(L, mu, params[3]))
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
  ev <- (pnorm(U, params[5], params[6]) - pnorm(L, params[5], params[6]))
  f_r(sv1, params[1:2]) * f_s(sv1, params[3:4]) * (f_d(sv2, params[5:6] / ev))
}

b <- seq(0, 50, length.out = 101)
d1 <- (b[2:101] + b[1:100]) * 0.5
h <- d1[3] - d1[2]

domains <- expand.grid(list(d2 = d1, d1 = d1))

G <- g(domains$d2,
       domains$d1,
       params = c(data_list$g_int,
                  data_list$g_slope,
                  data_list$sd_g),
       L = 0,
       U = 50)



S <- s(d1, c(data_list$s_int, data_list$s_slope))

P <- h * t(S * t(G))

Fm <- h * fec(domains$d2,
              domains$d1,
              params = unlist(data_list[6:11]),
              L = 0, U = 50)

K <- P + Fm

K <- matrix(K, nrow = 100, ncol = 100, byrow = TRUE)

lambda <- Re(eigen(K)$values[1])
w <- Re(eigen(K)$vectors[ , 1])
w <- w / sum(w)


inv_logit <- function(int, slope, sv) {
  return(1/(1 + exp(-(int + slope * sv))))
}

impl_args <- make_impl_args_list(c('P', 'F'),
                                 int_rule = rep('midpoint', 2),
                                 state_start = rep('dbh', 2),
                                 state_end = rep('dbh', 2))

states <- c('dbh', 'dbh')

x <- init_ipm(sim_gen    = "simple",
              di_dd      = "di",
              det_stoch  = "det") %>%
  define_kernel("P",
                formula = s * g,
                family = "CC",
                s = inv_logit(s_int, s_slope, dbh_1),
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
                f_r = inv_logit(f_r_int, f_r_slope, dbh_1),
                f_s = exp(f_s_int + f_s_slope * dbh_1),
                f_d = dnorm(dbh_2, mu_fd, sd_fd),
                data_list = data_list,
                states = states,
                evict_cor = TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'f_d')
  ) %>%
  define_impl(impl_args) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  define_pop_state(n_dbh = runif(100)) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit),
           normalize_pop_size = FALSE,
           iterate = TRUE,
           iterations = 200)


lambda_ipmr <- unname(lambda(x))
w_ipmr <- x$pop_state$n_dbh[ , 201] / sum(x$pop_state$n_dbh[ , 201])

init_pop_vec <- x$pop_state$n_dbh[ , 1]


test_that('lambdas are equal', {
  expect_equal(lambda - lambda_ipmr, 0, tolerance = 1e-6)
  expect_equal(w_ipmr, w, tolerance = 1e-6)
})

test_that('classes are correctly set', {

  expect_s3_class(x$sub_kernels$P, 'ipmr_matrix')
  expect_s3_class(x$sub_kernels$F, 'ipmr_matrix')
  expect_s3_class(x, 'simple_di_det_ipm')
  expect_s3_class(x, "ipmr_ipm")

})


test_that('define_impl() can handle mismatched argument lengths', {

  expect_warning(
    test_impl <- x$proto_ipm %>%
      define_impl(
        make_impl_args_list(
          kernel_names = c('P', "F"),
          int_rule     = 'midpoint',
          state_start    = rep('dbh', 2),
          state_end      = rep('dbh', 2)
        )
      ) %>%
      define_domains(dbh = c(0, 50, 100)) %>%
      make_ipm(usr_funs = list(inv_logit = inv_logit),
               normalize_pop_size = FALSE),

    regexp = "Assuming that all kernels are implemented with the same 'int_rule'."

  )

  expect_warning(
    test_impl <- x$proto_ipm %>%
      define_impl(
        make_impl_args_list(
          kernel_names = c('P', "F"),
          int_rule     = rep('midpoint', 2),
          state_start    = rep('dbh', 2),
          state_end      = 'dbh'
        )
      ) %>%
      define_domains(dbh = c(0, 50, 100)) %>%
      make_ipm(usr_funs = list(inv_logit = inv_logit),
               normalize_pop_size = FALSE),

    regexp = "Assuming that all kernels are implemented with the same 'state_end'."

  )

  expect_warning(
    test_impl <- x$proto_ipm %>%
      define_impl(
        make_impl_args_list(
          kernel_names = c('P', "F"),
          int_rule     = rep('midpoint', 2),
          state_start    = 'dbh',
          state_end      = rep('dbh', 2)
        )
      ) %>%
      define_domains(dbh = c(0, 50, 100)) %>%
      make_ipm(usr_funs = list(inv_logit = inv_logit),
               normalize_pop_size = FALSE,
               iterate = TRUE,
               iterations = 100),

    regexp = "Assuming that all kernels are implemented with the same 'state_start'."

  )

  test_impl_lambda <- unname(lambda(test_impl))

  expect_equal(test_impl_lambda, lambda_ipmr, tolerance = 1e-10)

})


test_that("order of kernel_definition doesn't matter", {

  y <- init_ipm(sim_gen    = "simple",
                di_dd      = "di",
                det_stoch  = "det") %>%
    define_kernel('F',
                  formula = f_r * f_s * f_d,
                  family = 'CC',
                  f_r = inv_logit(f_r_int, f_r_slope, dbh_1),
                  f_s = exp(f_s_int + f_s_slope * dbh_1),
                  f_d = dnorm(dbh_2, mu_fd, sd_fd),
                  data_list = data_list,
                  states = states,
                  evict_cor = TRUE,
                  evict_fun = truncated_distributions('norm',
                                                      'f_d')
    ) %>%
    define_kernel("P",
                  formula = s * g,
                  family = "CC",
                  s = inv_logit(s_int, s_slope, dbh_1),
                  g = dnorm(dbh_2, mu_g, sd_g),
                  mu_g = g_int + g_slope * dbh_1,
                  data_list = data_list,
                  states = states,
                  evict_cor = TRUE,
                  evict_fun = truncated_distributions('norm',
                                                      'g')
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c('F', 'P'),
        int_rule = rep('midpoint', 2),
        state_start = rep('dbh', 2),
        state_end   = rep('dbh', 2)
      )
    ) %>%
    define_domains(dbh = c(0, 50, 100)) %>%
    define_pop_state(n_dbh = init_pop_vec) %>%
    make_ipm(usr_funs = list(inv_logit = inv_logit),
             normalize_pop_size = FALSE,
             iterate = TRUE)

  lambda_out_of_order <- unname(lambda(y))

  expect_equal(lambda_ipmr, lambda_out_of_order)

})


test_that('iteration methods work the same as eigenvalue methods', {

  it <- init_ipm(sim_gen    = "simple",
                 di_dd      = "di",
                 det_stoch  = "det") %>%
    define_kernel("P",
                  formula = s * g,
                  family = "CC",
                  s = inv_logit(s_int, s_slope, dbh_1),
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
                  f_r = inv_logit(f_r_int, f_r_slope, dbh_1),
                  f_s = exp(f_s_int + f_s_slope * dbh_1),
                  f_d = dnorm(dbh_2, mu_fd, sd_fd),
                  data_list = data_list,
                  states = states,
                  evict_cor = TRUE,
                  evict_fun = truncated_distributions('norm',
                                                      'f_d')
    ) %>%
    define_impl(impl_args) %>%
    define_pop_state(
      n_dbh = runif(100)
    ) %>%
    define_domains(dbh = c(0, 50, 100)) %>%
    make_ipm(usr_funs = list(inv_logit = inv_logit),
             iterate = TRUE,
             iterations = 100,
             normalize_pop_size = FALSE)

  K <- it$sub_kernels$F + it$sub_kernels$P

  lambda_it <- unname(lambda(it))
  lambda_eig <- Re(eigen(K)$values[1])

  expect_equal(lambda_it, lambda_ipmr)
  expect_equal(lambda_it, lambda_eig)


})


test_that('normalizing pop vector produces same lambdas as eigen methods', {

  it <- init_ipm(sim_gen    = "simple",
                 di_dd      = "di",
                 det_stoch  = "det") %>%
    define_kernel("P",
                  formula = s * g,
                  family = "CC",
                  s = inv_logit(s_int, s_slope, dbh_1),
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
                  f_r = inv_logit(f_r_int, f_r_slope, dbh_1),
                  f_s = exp(f_s_int + f_s_slope * dbh_1),
                  f_d = dnorm(dbh_2, mu_fd, sd_fd),
                  data_list = data_list,
                  states = states,
                  evict_cor = TRUE,
                  evict_fun = truncated_distributions('norm',
                                                      'f_d')
    ) %>%
    define_impl(impl_args) %>%
    define_pop_state(
      n_dbh = runif(100)
    ) %>%
    define_domains(dbh = c(0, 50, 100)) %>%
    make_ipm(usr_funs = list(inv_logit = inv_logit),
             iterate = TRUE,
             iterations = 100,
             normalize_pop_size = TRUE)

  K <- it$sub_kernels$F + it$sub_kernels$P

  lambda_it <- unname(lambda(it, type_lambda = 'last'))
  lambda_eig <- Re(eigen(K)$values[1])
  names(lambda_eig) <- NULL
  names(lambda_it)  <- NULL

  conv_lamb <- lambda_it

  expect_equal(conv_lamb, lambda_ipmr)
  expect_equal(conv_lamb, lambda_eig)

})



fixed_list <- list(s_int = 2.2,
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

s_r_ints <- rnorm(5, -2.2, sd = 1) %>%
  as.list()

g_r_ints <- rnorm(5, - 0.2, 1) %>%
  as.list

names(s_r_ints) <- paste('s_r_', 1:5, sep = "")
names(g_r_ints) <- paste('g_r_', 1:5, sep = "")

data_list <- c(fixed_list, s_r_ints, g_r_ints)

s <- function(sv1, params, r_params) {
  1/(1 + exp(-(params[1] + params[2] * sv1 + r_params)))
}


g <- function(sv1, sv2, params, r_params, L, U) {
  mu <- params[1] + params[2] * sv1 + r_params
  ev <- (pnorm(U, mu, params[3]) - pnorm(L, mu, params[3]))
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
  ev <- (pnorm(U, params[5], params[6]) - pnorm(L, params[5], params[6]))
  f_r(sv1, params[1:2]) * f_s(sv1, params[3:4]) * (f_d(sv2, params[5:6] / ev))
}

b <- seq(0, 50, length.out = 101)
d1 <- (b[2:101] + b[1:100]) * 0.5
h <- d1[3] - d1[2]

domains <- expand.grid(list(d2 = d1, d1 = d1))

G_1 <- g(domains$d2,
         domains$d1,
         params = c(data_list$g_int,
                    data_list$g_slope,
                    data_list$sd_g),
         r_params = data_list$g_r_1,
         L = 0,
         U = 50)
G_5 <- g(domains$d2,
         domains$d1,
         params = c(data_list$g_int,
                    data_list$g_slope,
                    data_list$sd_g),
         r_params = data_list$g_r_5,
         L = 0,
         U = 50)
G_2 <- g(domains$d2,
         domains$d1,
         params = c(data_list$g_int,
                    data_list$g_slope,
                    data_list$sd_g),
         r_params = data_list$g_r_2,
         L = 0,
         U = 50)
G_3 <- g(domains$d2,
         domains$d1,
         params = c(data_list$g_int,
                    data_list$g_slope,
                    data_list$sd_g),
         r_params = data_list$g_r_3,
         L = 0,
         U = 50)
G_4 <- g(domains$d2,
         domains$d1,
         params = c(data_list$g_int,
                    data_list$g_slope,
                    data_list$sd_g),
         r_params = data_list$g_r_4,
         L = 0,
         U = 50)


S_1 <- s(d1, c(data_list$s_int, data_list$s_slope), data_list$s_r_1)
S_2 <- s(d1, c(data_list$s_int, data_list$s_slope), data_list$s_r_2)
S_3 <- s(d1, c(data_list$s_int, data_list$s_slope), data_list$s_r_3)
S_4 <- s(d1, c(data_list$s_int, data_list$s_slope), data_list$s_r_4)
S_5 <- s(d1, c(data_list$s_int, data_list$s_slope), data_list$s_r_5)

P_1 <- h * S_1 * G_1
P_2 <- h * S_2 * G_2
P_3 <- h * S_3 * G_3
P_4 <- h * S_4 * G_4
P_5 <- h * S_5 * G_5

Fm <- h * fec(domains$d2,
              domains$d1,
              params = unlist(data_list[6:11]),
              L = 0, U = 50)

K_1 <- P_1 + Fm
K_2 <- P_2 + Fm
K_3 <- P_3 + Fm
K_4 <- P_4 + Fm
K_5 <- P_5 + Fm

K_1 <- matrix(K_1, nrow = 100, ncol = 100, byrow = TRUE)
K_2 <- matrix(K_2, nrow = 100, ncol = 100, byrow = TRUE)
K_3 <- matrix(K_3, nrow = 100, ncol = 100, byrow = TRUE)
K_4 <- matrix(K_4, nrow = 100, ncol = 100, byrow = TRUE)
K_5 <- matrix(K_5, nrow = 100, ncol = 100, byrow = TRUE)

lambda_hand_1 <- Re(eigen(K_1)$values[1])
lambda_hand_2 <- Re(eigen(K_2)$values[1])
lambda_hand_3 <- Re(eigen(K_3)$values[1])
lambda_hand_4 <- Re(eigen(K_4)$values[1])
lambda_hand_5 <- Re(eigen(K_5)$values[1])

lambdas_hand <- c(lambda_hand_1, lambda_hand_2, lambda_hand_3,
                  lambda_hand_4, lambda_hand_5)

names(lambdas_hand) <- paste("K_", 1:5, sep = "")

states <- list(c("dbh"))

impl_args <-  make_impl_args_list(c('P_yr', 'F'),
                                  int_rule = rep('midpoint', 2),
                                  state_start = rep('dbh', 2),
                                  state_end = rep('dbh', 2))

inv_logit <- function(int, slope, sv) {

  1 / (1 + exp(-(int + slope * sv)))
}

par_set_mod <- init_ipm(sim_gen    = "simple",
                     di_dd      = "di",
                     det_stoch  = "det") %>%
  define_kernel("P_yr",
                formula = s_yr * g_yr,
                family = "CC",
                s_yr = inv_logit(s_int + s_r_yr, s_slope, dbh_1),
                g_yr = dnorm(dbh_2, mu_g_yr, sd_g),
                mu_g_yr = g_int + g_r_yr + g_slope * dbh_1,
                data_list = data_list,
                states = states,
                uses_par_sets = TRUE,
                par_set_indices = list(yr = 1:5),
                evict_cor = TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'g_yr')
  ) %>%
  define_kernel('F',
                formula = f_r * f_s * f_d,
                family = 'CC',
                f_r = inv_logit(f_r_int, f_r_slope, dbh_1),
                f_s = exp(f_s_int + f_s_slope * dbh_1),
                f_d = dnorm(dbh_2, mu_fd, sd_fd),
                data_list = data_list,
                states = states,
                evict_cor = TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'f_d')
  ) %>%
  define_impl(impl_args) %>%
  define_pop_state(
    n_dbh_yr = runif(100)
  ) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit),
           iterate = TRUE,
           iterations = 100,
           normalize_pop_size = FALSE)

lambdas_ipmr_pop   <- lambda(par_set_mod)

names(lambdas_ipmr_pop) <- paste("K_", 1:5, sep = "")

test_that('par_setarchical deterministic simulations work', {

  expect_equal(lambdas_ipmr_pop, lambdas_hand, tolerance = 1e-10)

})

test_that("make_iter_kernel works", {

  k_list <- make_iter_kernel(par_set_mod) %>%
    lapply(unclass)

  expect_equal(k_list[[1]], K_1)
  expect_equal(k_list[[2]], K_2)
  expect_equal(k_list[[3]], K_3)
  expect_equal(k_list[[4]], K_4)
  expect_equal(k_list[[5]], K_5)


})


test_that("discrete_extrema_works for simple models", {

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

  s <- function(sv1, params) {
    1/(1 + exp(-(params[1] + params[2] * sv1)))
  }


  g <- function(sv1, sv2, params, L, U) {
    mu <- params[1] + params[2] * sv1
    dnorm(sv2, mean = mu, sd = params[3])
  }

  f_r <- function(sv1, params) {
    1/(1 + exp(-(params[1] + params[2] * sv1)))
  }

  f_s <- function(sv1, params) {
    exp(params[1] + params[2] * sv1)
  }

  f_d <- function(sv2, params, h) {
    temp <- matrix(dnorm(sv2, mean = params[1], sd = params[2]),
                   nrow = 100,
                   ncol = 100,
                   byrow = TRUE) * h

    for(i in seq_len(50)) {

      temp[1, i]          <- temp[1, i] + (1 - sum(temp[ , i]))

      temp[100, (50 + i)] <- temp[100, (50 + i)] + (1 - sum(temp[ , (50 + i)]))
    }

    return(temp)
  }

  fec <- function(sv1, sv2, params, L, U, h) {
    m_fr <- matrix(f_r(sv1, params[1:2]), nrow = 100, ncol = 100, byrow = TRUE)
    m_fs <- matrix(f_s(sv1, params[3:4]), nrow = 100, ncol = 100, byrow = TRUE)
    m_fd <- f_d(sv2, params[5:6], h)

    m_fr * m_fs * m_fd
  }

  b <- seq(0, 50, length.out = 101)
  d1 <- (b[2:101] + b[1:100]) * 0.5
  h <- d1[3] - d1[2]

  domains <- expand.grid(list(d1 = d1, d2 = d1))

  G <- g(domains$d1,
         domains$d2,
         params = c(data_list$g_int,
                    data_list$g_slope,
                    data_list$sd_g),
         L = 0,
         U = 50) * h

  G <- matrix(G, nrow = 100, ncol = 100, byrow = TRUE)

  for(i in seq_len(50)) {
    G[1, i] <- G[1, i] + (1 - sum(G[ , i]))
    G[100, (50 + i)] <- G[100, (50 + i)] + (1 - sum(G[ , (50 + i)]))
  }

  S <- s(d1, c(data_list$s_int, data_list$s_slope))

  P <- t(S * t(G))

  Fm <- fec(domains$d1,
            domains$d2,
            params = unlist(data_list[6:11]),
            L = 0, U = 50, h = h)

  K <- P + Fm

  K <- matrix(K, nrow = 100, ncol = 100, byrow = TRUE)

  lambda <- Re(eigen(K)$values[1])

  # now for ipmr

  inv_logit <- function(int, slope, sv) {
    return(1/(1 + exp(-(int + slope * sv))))
  }

  impl_args <- make_impl_args_list(c('P', 'F'),
                                   int_rule = rep('midpoint', 2),
                                   state_start = rep('dbh', 2),
                                   state_end = rep('dbh', 2))

  states <- c('dbh')

  x <- init_ipm(sim_gen    = "simple",
                di_dd      = "di",
                det_stoch  = "det") %>%
    define_kernel("P",
                  formula = s * g,
                  family = "CC",
                  s = inv_logit(s_int, s_slope, dbh_1),
                  g = dnorm(dbh_2, mu_g, sd_g),
                  mu_g = g_int + g_slope * dbh_1,
                  data_list = data_list,
                  states = states,
                  evict_cor = TRUE,
                  evict_fun = discrete_extrema('g', "dbh")
    ) %>%
    define_kernel('F',
                  formula = f_r * f_s * f_d,
                  family = 'CC',
                  f_r = inv_logit(f_r_int, f_r_slope, dbh_1),
                  f_s = exp(f_s_int + f_s_slope * dbh_1),
                  f_d = dnorm(dbh_2, mu_fd, sd_fd),
                  data_list = data_list,
                  states = states,
                  evict_cor = TRUE,
                  evict_fun = discrete_extrema('f_d', "dbh")
    )  %>%
    define_impl(impl_args) %>%
    define_domains(dbh = c(0, 50, 100)) %>%
    make_ipm(usr_funs = list(inv_logit = inv_logit),
             normalize_pop_size = FALSE,
             iterate = FALSE)

  K_ipmr <- do.call(`+`, x$sub_kernels)

  lambda_ipmr <- Re(eigen(K_ipmr)$values[1])

  expect_equal(lambda, lambda_ipmr)

})
