library(rlang)
library(purrr)

# define functions for target ipm

# Survival - logistic regression
s <- function(sv1, params, r_effect) {
  1/(1 + exp(-(params[1] + params[2] * sv1 + r_effect))) *
    (1 - f_r(sv1, params[3:4]))
}

# Growth
g <- function(sv1, sv2, params, r_effect, L, U) {
  mu <- params[1] + params[2] * sv1 + r_effect
  ev <- pnorm(U, mu, params[3]) - pnorm(L, mu, params[3])
  dnorm(sv2, mean = mu, sd = params[3]) / ev
}

# probability of reproducing
f_r <- function(sv1, params) {
  1/(1 + exp(-(params[1] + params[2] * sv1)))
}

# offspring production
f_s <- function(sv1, params, r_effect) {
  exp(params[1] + params[2] * sv1 + r_effect)
}

# offspring size distribution
f_d <- function(sv2, params, L, U) {
  ev <- pnorm(U, params[1], params[2]) - pnorm(L, params[1], params[2])
  dnorm(sv2, mean = params[1], sd = params[2]) / ev
}

# constructor function for the F kernel
fec <- function(sv1, sv2, params, r_effect, L, U) {
  f_r(sv1, params[1:2]) * f_s(sv1, params[3:4], r_effect) * f_d(sv2, params[5:6], L, U)
}


set.seed(50127)


# Define some fixed parameters
data_list = list(
  s_int     = 1.03,
  s_slope   = 2.2,
  g_int     = 8,
  g_slope   = 0.92,
  sd_g      = 0.9,
  f_r_int   = 0.09,
  f_r_slope = 0.05,
  f_s_int   = 0.1,
  f_s_slope = 0.005,
  mu_fd     = 9,
  sd_fd     = 2
)

# Now, simulate some random intercepts for growth, survival, and offspring production

g_r_int   <- rnorm(5, 0, 0.3)
s_r_int   <- rnorm(5, 0, 0.7)
f_s_r_int <- rnorm(5, 0, 0.2)

nms <- paste("r_", 1:5, sep = "")

names(g_r_int)   <- paste('g_', nms, sep = "")
names(s_r_int)   <- paste('s_', nms, sep = "")
names(f_s_r_int) <- paste('f_s_', nms, sep = "")

# The !!! operator used inside of list2 from rlang takes the named vector
# and converts it to a named list. This can be spliced into the data list
# to rapidly make a parameter set suitable for usage in the data_list argument
# of define_kernel

g_params   <- list2(!!! g_r_int)
s_params   <- list2(!!! s_r_int)
f_s_params <- list2(!!! f_s_r_int)

params     <- c(data_list, g_params, s_params, f_s_params)


b   <- seq(0.2, 40, length.out = 101)

sv1 <- (b[2:101] + b[1:100]) * 0.5

domains <- expand.grid(list(d2 = sv1, d1 = sv1))

h   <- sv1[2] - sv1[1]

# repetitive to demonstrate the typical kernel construction process.

g_1 <- g(domains$d2, domains$d1,
         params = c(params$g_int,
                    params$g_slope,
                    params$sd_g),
         r_effect = params$g_r_1,
         L = 0.2,
         U = 40)

g_2 <- g(domains$d2, domains$d1,
         params = c(params$g_int,
                    params$g_slope,
                    params$sd_g),
         r_effect = params$g_r_2,
         L = 0.2,
         U = 40)

g_3 <- g(domains$d2, domains$d1,
         params = c(params$g_int,
                    params$g_slope,
                    params$sd_g),
         r_effect = params$g_r_3,
         L = 0.2,
         U = 40)


g_4 <- g(domains$d2, domains$d1,
         params = c(params$g_int,
                    params$g_slope,
                    params$sd_g),
         r_effect = params$g_r_4,
         L = 0.2,
         U = 40)

g_5 <-g(domains$d2, domains$d1,
        params = c(params$g_int,
                   params$g_slope,
                   params$sd_g),
        r_effect = params$g_r_5,
        L = 0.2,
        U = 40)

s_1 <- s(sv1, c(params$s_int, params$s_slope,
                params$f_r_int, params$f_r_slope), params$s_r_1)
s_2 <- s(sv1, c(params$s_int, params$s_slope,
                params$f_r_int, params$f_r_slope), params$s_r_2)
s_3 <- s(sv1, c(params$s_int, params$s_slope,
                params$f_r_int, params$f_r_slope), params$s_r_3)
s_4 <- s(sv1, c(params$s_int, params$s_slope,
                params$f_r_int, params$f_r_slope), params$s_r_4)
s_5 <- s(sv1, c(params$s_int, params$s_slope,
                params$f_r_int, params$f_r_slope), params$s_r_5)

P_1 <- t(s_1 * t(g_1)) * h
P_2 <- t(s_2 * t(g_2)) * h
P_3 <- t(s_3 * t(g_3)) * h
P_4 <- t(s_4 * t(g_4)) * h
P_5 <- t(s_5 * t(g_5)) * h

# These are not corrected for eviction, but they probably should be

F_1 <- h * fec(domains$d2, domains$d1,
               params = unlist(params[6:11]),
               r_effect = params$f_s_r_1,
               L = 0.2,
               U = 40)

F_2 <- h * fec(domains$d2, domains$d1,
               params = unlist(params[6:11]),
               r_effect = params$f_s_r_2,
               L = 0.2,
               U = 40)

F_3 <- h * fec(domains$d2, domains$d1,
               params = unlist(params[6:11]),
               r_effect = params$f_s_r_3,
               L = 0.2,
               U = 40)
F_4 <- h * fec(domains$d2, domains$d1,
               params = unlist(params[6:11]),
               r_effect = params$f_s_r_4,
               L = 0.2,
               U = 40)
F_5 <- h * fec(domains$d2, domains$d1,
               params = unlist(params[6:11]),
               r_effect = params$f_s_r_5,
               L = 0.2,
               U = 40)

K_1 <- (P_1 + F_1) %>%
  matrix(nrow = 100, ncol = 100, byrow = TRUE)
K_2 <- (P_2 + F_2) %>%
  matrix(nrow = 100, ncol = 100, byrow = TRUE)
K_3 <- (P_3 + F_3) %>%
  matrix(nrow = 100, ncol = 100, byrow = TRUE)
K_4 <- (P_4 + F_4) %>%
  matrix(nrow = 100, ncol = 100, byrow = TRUE)
K_5 <- (P_5 + F_5) %>%
  matrix(nrow = 100, ncol = 100, byrow = TRUE)

sys <- list(K_1 = K_1,
            K_2 = K_2,
            K_3 = K_3,
            K_4 = K_4,
            K_5 = K_5)

eigen_sys <- lapply(sys, eigen)

lambdas <- vapply(eigen_sys, function(x) Re(x$values[1]), numeric(1)) %>%
  unname()
ws      <- vapply(eigen_sys, function(x) Re(x$vectors[ , 1]), numeric(100)) %>%
  unname()

## ipmr version

# define the levels of the par_setarchical variable and save them in a named
# list that corresponds to the suffix in the kernel notation

par_set_indices <- list(yr = 1:5)

# additional usr_funs to be passed into make_ipm()

inv_logit <- function(sv, int, slope) {
  return(
    1/(1 + exp(-(int + slope * sv)))
  )
}

inv_logit_r <- function(sv, int, slope, r_eff) {
  return(
    1/(1 + exp(-(int + slope * sv + r_eff)))
  )
}

pois_r <- function(sv, int, slope, r_eff) {
  return(
    exp(
      int + slope * sv + r_eff
    )
  )
}


monocarp_sys <- init_ipm(sim_gen    = "simple",
                         di_dd      = "di",
                         det_stoch  = "stoch",
                         kern_param = "kern") %>%
  define_kernel(

    name             = 'P_yr',
    formula          = s_yr * g_yr,
    family           = "CC",
    s_yr             = inv_logit_r(ht_1, s_int, s_slope, s_r_yr) *
      (1 - inv_logit(ht_1, f_r_int, f_r_slope)),
    g_yr             = dnorm(ht_2, mu_g_yr, sd_g),
    mu_g_yr          = g_int + g_slope * ht_1 + g_r_yr,
    data_list        = params,
    states           = list(c('ht')),
    uses_par_sets    = TRUE,
    par_set_indices = par_set_indices,
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions('norm', 'g_yr')
  ) %>%
  define_kernel(
    name             = "F_yr",
    formula          = f_r * f_s_yr * f_d,
    family           = "CC",
    f_r              = inv_logit(ht_1, f_r_int, f_r_slope),
    f_s_yr           = pois_r(ht_1, f_s_int, f_s_slope, f_s_r_yr),
    f_d              = dnorm(ht_2, mu_fd, sd_fd),
    data_list        = params,
    states           = list(c('ht')),
    uses_par_sets    = TRUE,
    par_set_indices = par_set_indices,
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions('norm', 'f_d')
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P_yr", "F_yr"),
      int_rule     = rep("midpoint", 2),
      state_start    = rep("ht", 2),
      state_end      = rep("ht", 2)
    )
  ) %>%
  define_domains(ht = c(0.2, 40, 100)) %>%
  define_pop_state(n_ht = runif(100)) %>%
  make_ipm(usr_funs = list(inv_logit   = inv_logit,
                           inv_logit_r = inv_logit_r,
                           pois_r      = pois_r),
           normalize_pop_size = TRUE,
           iterate = TRUE,
           iterations = 100)

Ks <- make_iter_kernel(monocarp_sys) %>%
  lapply(unclass)

lambdas_ipmr <- vapply(Ks,
                       function(x) Re(eigen(x)$values[1]),
                       numeric(1)) %>%
  unname()

ws_ipmr <- vapply(Ks,
                  function(x) Re(eigen(x)$vectors[ , 1]),
                  numeric(100L))

test_that('eigenvectors and values are correct', {

  expect_equal(lambdas_ipmr, lambdas, tolerance = 1e-10)
  expect_equal(ws_ipmr[ ,1], ws[ ,1], tolerance = 1e-13)
  expect_equal(ws_ipmr[ ,2], ws[ ,2], tolerance = 1e-13)
  expect_equal(ws_ipmr[ ,3], ws[ ,3], tolerance = 1e-13)
  expect_equal(ws_ipmr[ ,4], ws[ ,4], tolerance = 1e-13)
  expect_equal(ws_ipmr[ ,5], ws[ ,5], tolerance = 1e-13)

})

test_that("make_iter_kernel works as expected", {

  expect_equal(Ks[[1]], K_1)
  expect_equal(Ks[[2]], K_2)
  expect_equal(Ks[[3]], K_3)
  expect_equal(Ks[[4]], K_4)
  expect_equal(Ks[[5]], K_5)

})


# Test whether .iterate kerns does what it should
kern_seq <- sample(1:5, 50, replace = TRUE)

proto <- monocarp_sys$proto_ipm

init_pop_vec <- runif(100)

init_pop_vec <- init_pop_vec / sum(init_pop_vec)


iterated_sys <- proto %>%
  define_pop_state(
    n_ht = init_pop_vec
  ) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit,
                           inv_logit_r = inv_logit_r,
                           pois_r = pois_r),
           kernel_seq = kern_seq,
           iterate = TRUE,
           iterations = 50,
           normalize_pop_size = TRUE)

ipmr_pop_state <- iterated_sys$pop_state$n_ht

pop_holder <- array(NA_real_, dim = c(100, 51))

pop_holder[ , 1] <- init_pop_vec / sum(init_pop_vec)


v_holder <- array(NA_real_, dim = c(51, 100))
v_holder[1 , ] <- init_pop_vec / sum(init_pop_vec)

pop_size_lambdas <- numeric(50L)

for(i in seq_len(50)) {

  k_selector <- kern_seq[i]
  k_temp <- sys[[k_selector]]

  n_t_1 <- k_temp %*% pop_holder[ , i]

  pop_holder[ , (i + 1)] <- n_t_1 / sum(n_t_1)

  pop_size_lambdas[i] <- sum(n_t_1)

  v_t_1 <- left_mult(k_temp, v_holder[i , ])

  v_holder[(i + 1), ] <- v_t_1 / sum(v_t_1)

}

lambdas_ipmr <- as.vector(lambda(iterated_sys,
                                 type_lambda = 'all'))

test_that('.iterate_kerns is acting correctly', {

  expect_equal(pop_size_lambdas, lambdas_ipmr, tolerance = 1e-10)

})

test_that("burn_in works", {

  burn_lams <- c(lambda = mean(log(iterated_sys$pop_state$lambda[6:50])))
  stoch_lam <- lambda(iterated_sys)

  expect_equal(burn_lams, stoch_lam)

  unburn_lams      <- c(lambda = mean(log(iterated_sys$pop_state$lambda)))
  unburn_stoch_lam <- lambda(iterated_sys, burn_in = 0)

  expect_equal(unburn_lams, unburn_stoch_lam)

})

test_that("log = TRUE works for all type_lambda= 'last' and 'all'", {
  expect_equal(log(lambda(monocarp_sys, type_lambda = "last")),
               lambda(monocarp_sys, type_lambda = "last", log = TRUE))
  expect_equal(log(lambda(monocarp_sys, type_lambda = "all")),
               lambda(monocarp_sys, type_lambda = "all", log = TRUE))
})

test_that('classes are correctly set', {

  sub_cls <- vapply(monocarp_sys$sub_kernels,
                  function(x) class(x)[1],
                  character(1L))

  expect_true(all(sub_cls == 'ipmr_matrix'))
  expect_s3_class(monocarp_sys, 'simple_di_stoch_kern_ipm')
  expect_s3_class(monocarp_sys, 'ipmr_ipm')

})

hand_v <- list(t(v_holder[13:51, ]))
ipmr_v <- left_ev(iterated_sys, iterations = 50)

hand_w <- list(pop_holder[ , 13:51])
ipmr_w <- right_ev(iterated_sys)

test_that("left and right_ev work correctly", {

  expect_s3_class(ipmr_v, "ipmr_v")
  expect_equal(hand_v, ipmr_v, ignore_attr = TRUE)

  expect_s3_class(ipmr_w, "ipmr_w")
  expect_equal(hand_w, ipmr_w, ignore_attr = TRUE)

  expect_equal(colSums(ipmr_w[[1]]), rep(1L, dim(ipmr_w[[1]])[2]))
  expect_equal(colSums(ipmr_v[[1]]), rep(1L, dim(ipmr_v[[1]])[2]))

})


test_that("order of kernel definition doesn't matter", {

  test_order_1 <- init_ipm(sim_gen    = "simple",
                           di_dd      = "di",
                           det_stoch  = "stoch",
                           kern_param = "kern") %>%
    define_kernel(
      name             = "F_yr",
      formula          = f_r * f_s_yr * f_d,
      family           = "CC",
      f_r              = inv_logit(ht_1, f_r_int, f_r_slope),
      f_s_yr           = pois_r(ht_1, f_s_int, f_s_slope, f_s_r_yr),
      f_d              = dnorm(ht_2, mu_fd, sd_fd),
      data_list        = params,
      states           = list(c('ht')),
      uses_par_sets    = TRUE,
      par_set_indices = par_set_indices,
      evict_cor        = TRUE,
      evict_fun        = truncated_distributions('norm', 'f_d')
    ) %>%
    define_kernel(
      name             = 'P_yr',
      formula          = s_yr * g_yr,
      family           = "CC",
      s_yr             = inv_logit_r(ht_1, s_int, s_slope, s_r_yr) *
        (1 - inv_logit(ht_1, f_r_int, f_r_slope)),
      g_yr             = dnorm(ht_2, mu_g_yr, sd_g),
      mu_g_yr          = g_int + g_slope * ht_1 + g_r_yr,
      data_list        = params,
      states           = list(c('ht')),
      uses_par_sets    = TRUE,
      par_set_indices = par_set_indices,
      evict_cor        = TRUE,
      evict_fun        = truncated_distributions('norm', 'g_yr')
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("F_yr", "P_yr"),
        int_rule     = rep("midpoint", 2),
        state_start    = rep("ht", 2),
        state_end      = rep("ht", 2)
      )
    ) %>%
    define_domains(ht = c(0.2, 40, 100)) %>%
    define_pop_state(n_ht = init_pop_vec) %>%
    make_ipm(usr_funs = list(inv_logit   = inv_logit,
                             inv_logit_r = inv_logit_r,
                             pois_r      = pois_r),
             kernel_seq = kern_seq,
             normalize_pop_size = TRUE)

  Ks <- make_iter_kernel(test_order_1) %>%
    lapply(unclass)

  lambdas_test <- as.vector(lambda(test_order_1, type_lambda = "all"))

  w_test <- right_ev(test_order_1)
  v_test <- left_ev(test_order_1, iterations = 50)


  expect_equal(lambdas_ipmr, lambdas_test, tolerance = 1e-10)
  expect_equal(ipmr_w, w_test)
  expect_equal(ipmr_v, v_test)
})

test_that("return_all gets all of the environments back", {

  test_order_1 <- init_ipm(sim_gen    = "simple",
                           di_dd      = "di",
                           det_stoch  = "stoch",
                           kern_param = "kern") %>%
    define_kernel(
      name             = "F_yr",
      formula          = f_r * f_s_yr * f_d,
      family           = "CC",
      f_r              = inv_logit(ht_1, f_r_int, f_r_slope),
      f_s_yr           = pois_r(ht_1, f_s_int, f_s_slope, f_s_r_yr),
      f_d              = dnorm(ht_2, mu_fd, sd_fd),
      data_list        = params,
      states           = list(c('ht')),
      uses_par_sets    = TRUE,
      par_set_indices = par_set_indices,
      evict_cor        = TRUE,
      evict_fun        = truncated_distributions('norm', 'f_d')
    ) %>%
    define_kernel(
      name             = 'P_yr',
      formula          = s_yr * g_yr,
      family           = "CC",
      s_yr             = inv_logit_r(ht_1, s_int, s_slope, s_r_yr) *
        (1 - inv_logit(ht_1, f_r_int, f_r_slope)),
      g_yr             = dnorm(ht_2, mu_g_yr, sd_g),
      mu_g_yr          = g_int + g_slope * ht_1 + g_r_yr,
      data_list        = params,
      states           = list(c('ht')),
      uses_par_sets    = TRUE,
      par_set_indices = par_set_indices,
      evict_cor        = TRUE,
      evict_fun        = truncated_distributions('norm', 'g_yr')
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("F_yr", "P_yr"),
        int_rule     = rep("midpoint", 2),
        state_start    = rep("ht", 2),
        state_end      = rep("ht", 2)
      )
    ) %>%
    define_domains(ht = c(0.2, 40, 100)) %>%
    define_pop_state(
      n_ht = runif(100)
    ) %>%
    make_ipm(usr_funs = list(inv_logit   = inv_logit,
                             inv_logit_r = inv_logit_r,
                             pois_r      = pois_r),
             return_all_envs = TRUE,
             normalize_pop_size = FALSE)

  env_list_nms <- c('main_env', c(paste('F', 1:5, sep = '_'),
                                    paste("P", 1:5, sep = "_")))

  expect_equal(names(test_order_1$env_list), env_list_nms)
})


test_that('normalizing pop vector gets same lambdas as before', {

  init_pop <- runif(100)
  usr_seq  <- sample(1:5, size = 100, replace = TRUE)

  test_norm_1 <- init_ipm(sim_gen    = "simple",
                          di_dd      = "di",
                          det_stoch  = "stoch",
                          kern_param = "kern") %>%
    define_kernel(
      name             = "F_yr",
      formula          = f_r * f_s_yr * f_d,
      family           = "CC",
      f_r              = inv_logit(ht_1, f_r_int, f_r_slope),
      f_s_yr           = pois_r(ht_1, f_s_int, f_s_slope, f_s_r_yr),
      f_d              = dnorm(ht_2, mu_fd, sd_fd),
      data_list        = params,
      states           = list(c('ht')),
      uses_par_sets    = TRUE,
      par_set_indices = par_set_indices,
      evict_cor        = TRUE,
      evict_fun        = truncated_distributions('norm', 'f_d')
    ) %>%
    define_kernel(
      name             = 'P_yr',
      formula          = s_yr * g_yr,
      family           = "CC",
      s_yr             = inv_logit_r(ht_1, s_int, s_slope, s_r_yr) *
        (1 - inv_logit(ht_1, f_r_int, f_r_slope)),
      g_yr             = dnorm(ht_2, mu_g_yr, sd_g),
      mu_g_yr          = g_int + g_slope * ht_1 + g_r_yr,
      data_list        = params,
      states           = list(c('ht')),
      uses_par_sets    = TRUE,
      par_set_indices = par_set_indices,
      evict_cor        = TRUE,
      evict_fun        = truncated_distributions('norm', 'g_yr')
    )  %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("F_yr", "P_yr"),
        int_rule     = rep("midpoint", 2),
        state_start    = rep("ht", 2),
        state_end      = rep("ht", 2)
      )
    ) %>%
    define_domains(
      ht = c(0.2, 40, 100)
    ) %>%
    define_pop_state(
      n_ht = init_pop
    ) %>%
    make_ipm(usr_funs = list(inv_logit   = inv_logit,
                             inv_logit_r = inv_logit_r,
                             pois_r      = pois_r),
             normalize_pop_size = TRUE,
             iterate = TRUE,
             iterations = 100,
             kernel_seq = usr_seq)

  lambdas_test <- lambda(test_norm_1,
                         type_lambda = 'all') %>%
    as.vector()

  pop_holder       <- array(NA_real_, dim = c(100, 101))
  pop_holder[ , 1] <- init_pop / sum(init_pop)
  lambdas_hand     <- numeric(100L)

  for(i in seq_len(100)) {

    k_selector <- as.integer(usr_seq[i])
    use_k      <- sys[[k_selector]]

    n_t_1      <- use_k %*% pop_holder[ , i]

    # Store lambda, normalize pop vec and stick into holder

    lambdas_hand[i] <- sum(n_t_1)

    pop_holder[ , (i + 1)] <- n_t_1 / sum(n_t_1)

  }

  expect_equal(lambdas_hand, lambdas_test, tolerance = 1e-10)
  expect_equal(pop_holder, test_norm_1$pop_state$n_ht)

  pop_sizes <- colSums(test_norm_1$pop_state$n_ht)

  expect_equal(pop_sizes, rep(1, 101), tolerance = 1e-15)

})


test_that("drop_levels works in par_sets", {

  # Define some fixed parameters
  data_list = list(
    s_int     = 1.03,
    s_slope   = 2.2,
    g_int     = 8,
    g_slope   = 0.92,
    sd_g      = 0.9,
    f_r_int   = 0.09,
    f_r_slope = 0.05,
    f_s_int   = 0.1,
    f_s_slope = 0.005,
    mu_fd     = 9,
    sd_fd     = 2
  )

  # Now, simulate some random intercepts for growth, survival, and offspring production

  g_r_int   <- rnorm(8, 0, 0.3)
  s_r_int   <- rnorm(8, 0, 0.7)
  f_s_r_int <- rnorm(8, 0, 0.2)

  par_set_indices <- list(yr = 1:5, site = c("a", "b"))

  levels <- expand.grid(par_set_indices)

  levels <- apply(levels, 1, function(x) paste(x[1], x[2], sep = "_"))

  to_drop <- c("1_a", "4_b")

  levels <- levels[!levels %in% to_drop]

  nms <- paste("r_", levels, sep = "")

  names(g_r_int)   <- paste('g_', nms, sep = "")
  names(s_r_int)   <- paste('s_', nms, sep = "")
  names(f_s_r_int) <- paste('f_s_', nms, sep = "")

  par_set_indices$drop_levels <- to_drop

  g_params   <- list2(!!! g_r_int)
  s_params   <- list2(!!! s_r_int)
  f_s_params <- list2(!!! f_s_r_int)

  params     <- c(data_list, g_params, s_params, f_s_params)


  init_pop <- runif(100)
  usr_seq  <- sample(levels, size = 100, replace = TRUE)


  test_drop <- init_ipm(sim_gen    = "simple",
                        di_dd      = "di",
                        det_stoch  = "stoch",
                        kern_param = "kern") %>%
    define_kernel(
      name             = "F_yr_site",
      formula          = f_r * f_s_yr_site * f_d,
      family           = "CC",
      f_r              = inv_logit(ht_1, f_r_int, f_r_slope),
      f_s_yr_site      = pois_r(ht_1, f_s_int, f_s_slope, f_s_r_yr_site),
      f_d              = dnorm(ht_2, mu_fd, sd_fd),
      data_list        = params,
      states           = list(c('ht')),
      uses_par_sets    = TRUE,
      par_set_indices = par_set_indices,
      evict_cor        = TRUE,
      evict_fun        = truncated_distributions('norm', 'f_d')
    ) %>%
    define_kernel(
      name             = 'P_yr_site',
      formula          = s_yr_site * g_yr_site,
      family           = "CC",
      s_yr_site        = inv_logit_r(ht_1, s_int, s_slope, s_r_yr_site) *
        (1 - inv_logit(ht_1, f_r_int, f_r_slope)),
      g_yr_site        = dnorm(ht_2, mu_g_yr_site, sd_g),
      mu_g_yr_site     = g_int + g_slope * ht_1 + g_r_yr_site,
      data_list        = params,
      states           = list(c('ht')),
      uses_par_sets    = TRUE,
      par_set_indices = par_set_indices,
      evict_cor        = TRUE,
      evict_fun        = truncated_distributions('norm', 'g_yr_site')
    )  %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("F_yr_site", "P_yr_site"),
        int_rule     = rep("midpoint", 2),
        state_start    = rep("ht", 2),
        state_end      = rep("ht", 2)
      )
    ) %>%
    define_domains(
      ht = c(0.2, 40, 100)
    ) %>%
    define_pop_state(
      n_ht = init_pop
    ) %>%
    make_ipm(usr_funs = list(inv_logit   = inv_logit,
                             inv_logit_r = inv_logit_r,
                             pois_r      = pois_r),
             normalize_pop_size = TRUE,
             iterate = TRUE,
             iterations = 100,
             kernel_seq = usr_seq)

  kerns <-test_drop$sub_kernels

  for(i in seq_along(to_drop)) {

    expect_true(all(!grepl(to_drop[i], names(kerns))))

  }

})
