context('simple density indepenent models')

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


inv_logit <- function(int, slope, sv) {
  return(1/(1 + exp(-(int + slope * sv))))
}

impl_args <- make_impl_args_list(c('P', 'F', 'K'),
                                 int_rule = rep('midpoint', 3),
                                 dom_start = rep('dbh', 3),
                                 dom_end = rep('dbh', 3))

states <- c('dbh', 'dbh')

x <- init_ipm('simple_di_det') %>%
  define_kernel("P",
                formula = s_g_mult(s, g),
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
  define_k('K',
           K = P + F,
           family = 'IPM',
           data_list = list(),
           states = states,
           evict_cor = FALSE) %>%
  define_impl(impl_args) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit))

K_ipmr <- x$iterators$K

lambda_ipmr <- Re(eigen(K_ipmr)$values[1])
w_ipmr <- Re(eigen(K_ipmr)$vectors[ , 1])


test_that('lambdas are equal', {
  expect_equal(lambda - lambda_ipmr, 0, tolerance = 1e-6)
  expect_equal(w_ipmr, w, tolerance = 1e-6)
})

test_that('classes are correctly set', {

  expect_s3_class(x$iterators$K, 'ipmr_matrix')
  expect_s3_class(x$sub_kernels$P, 'ipmr_matrix')
  expect_s3_class(x$sub_kernels$F, 'ipmr_matrix')
  expect_s3_class(x, 'simple_di_det_ipm')

})


test_that('define_impl() can handle mismatched argument lengths', {

  expect_warning(
    test_impl <- x$proto_ipm %>%
      define_impl(
        make_impl_args_list(
          kernel_names = c('K', 'P', "F"),
          int_rule     = 'midpoint',
          dom_start    = rep('dbh', 3),
          dom_end      = rep('dbh', 3)
        )
      ) %>%
      define_domains(dbh = c(0, 50, 100)) %>%
      make_ipm(usr_funs = list(inv_logit = inv_logit)),

    regexp = "Assuming that all kernels are implemented with the same 'int_rule'."

  )

  expect_warning(
    test_impl <- x$proto_ipm %>%
      define_impl(
        make_impl_args_list(
          kernel_names = c('K', 'P', "F"),
          int_rule     = 'midpoint',
          dom_start    = rep('dbh', 3),
          dom_end      = 'dbh'
        )
      ) %>%
      define_domains(dbh = c(0, 50, 100)) %>%
      make_ipm(usr_funs = list(inv_logit = inv_logit)),

    regexp = "Assuming that all kernels are implemented with the same 'dom_end'."

  )

  expect_warning(
    test_impl <- x$proto_ipm %>%
      define_impl(
        make_impl_args_list(
          kernel_names = c('K', 'P', "F"),
          int_rule     = rep('midpoint', 3),
          dom_start    = 'dbh',
          dom_end      = rep('dbh', 3)
        )
      ) %>%
      define_domains(dbh = c(0, 50, 100)) %>%
      make_ipm(usr_funs = list(inv_logit = inv_logit)),

    regexp = "Assuming that all kernels are implemented with the same 'dom_start'."

  )

  test_impl_lambda <- Re(eigen(test_impl$iterators$K)$values[1])

  expect_equal(test_impl_lambda, lambda_ipmr, tol = 1e-10)

})


test_that("order of kernel_definition doesn't matter", {

  y <- init_ipm('simple_di_det') %>%
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
                  formula = s_g_mult(s, g),
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
    define_k('K',
             K = P + F,
             family = 'IPM',
             data_list = list(),
             states = states,
             evict_cor = FALSE) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c('F', 'P', "K"),
        int_rule = rep('midpoint', 3),
        dom_start = rep('dbh', 3),
        dom_end   = rep('dbh', 3)
      )) %>%
    define_domains(dbh = c(0, 50, 100)) %>%
    make_ipm(usr_funs = list(inv_logit = inv_logit))

  lambda_out_of_order <- Re(eigen(y$iterators$K)$values[1])

  expect_equal(lambda_ipmr, lambda_out_of_order)

  y <- init_ipm('simple_di_det') %>%
    define_k('K',
             K = P + F,
             family = 'IPM',
             data_list = list(),
             states = states,
             evict_cor = FALSE) %>%
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
                  formula = s_g_mult(s, g),
                  family = "CC",
                  s = inv_logit(s_int, s_slope, dbh_1),
                  g = dnorm(dbh_2, mu_g, sd_g),
                  mu_g = g_int + g_slope * dbh_1,
                  data_list = data_list,
                  states = states,
                  evict_cor = TRUE,
                  evict_fun = truncated_distributions('norm',
                                                      'g')
    )  %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c('K', 'F', 'P'),
        int_rule = rep('midpoint', 3),
        dom_start = rep('dbh', 3),
        dom_end   = rep('dbh', 3)
      )
    ) %>%
    define_domains(dbh = c(0, 50, 100)) %>%
    make_ipm(usr_funs = list(inv_logit = inv_logit))

  lambda_out_of_order_2 <- Re(eigen(y$iterators$K)$values[1])

  expect_equal(lambda_ipmr, lambda_out_of_order_2)

})


test_that('iteration methods work the same as eigenvalue methods', {

  it <- init_ipm('simple_di_det') %>%
    define_kernel("P",
                  formula = s_g_mult(s, g),
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
    define_k('K',
             K = P + F,
             n_dbh_t_1 = K %*% n_dbh_t,
             family = 'IPM',
             data_list = list(),
             states = states,
             evict_cor = FALSE) %>%
    define_impl(impl_args) %>%
    define_pop_state(
      n_dbh = runif(100)
    ) %>%
    define_domains(dbh = c(0, 50, 100)) %>%
    make_ipm(usr_funs = list(inv_logit = inv_logit),
             iterate = TRUE,
             iterations = 100)

  lambda_it <- lambda(it, comp_method = 'pop_size')
  lambda_eig <- lambda(it, comp_method = 'eigen')
  names(lambda_eig) <- NULL

  conv_lamb <- lambda_it[100]

  expect_equal(conv_lamb, lambda_ipmr)
  expect_equal(conv_lamb, lambda_eig)


})
