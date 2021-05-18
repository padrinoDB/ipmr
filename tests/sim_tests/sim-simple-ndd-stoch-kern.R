# sim-simple_di_stoch_kern

library(rlang)
library(ipmr)
library(purrr)

s <- function(sv1, params, r_effect) {
  1/(1 + exp(-(params[1] + params[2] * sv1 + r_effect))) *
    (1 - f_r(sv1, params[3:4]))
}


g <- function(sv1, sv2, params, r_effect) {
  mu <- params[1] + params[2] * sv1 + r_effect
  dnorm(sv2, mean = mu, sd = params[3])
}

f_r <- function(sv1, params) {
  1/(1 + exp(-(params[1] + params[2] * sv1)))
}

f_s <- function(sv1, params, r_effect) {
  exp(params[1] + params[2] * sv1 + r_effect)
}

f_d <- function(sv2, params) {
  dnorm(sv2, mean = params[1], sd = params[2])
}

fec <- function(sv1, sv2, params, r_effect) {
  f_r(sv1, params[1:2]) * f_s(sv1, params[3:4], r_effect) * f_d(sv2, params[5:6])
}



data_list = list(s_int = 1.03,
                 s_slope = 2.2,
                 g_int = 8,
                 g_slope = 0.92,
                 sd_g = 0.9,
                 f_r_int = 0.09,
                 f_r_slope = 0.05,
                 f_s_int = 0.1,
                 f_s_slope = 0.005,
                 mu_fd = 9,
                 sd_fd = 2)

# Now, simulate some random intercepts
g_r_int <- rnorm(5, 0, 0.3)
s_r_int <- rnorm(5, 0, 0.7)
f_s_r_int <- rnorm(5, 0, 0.2)

nms <- paste("r_", 1:5, sep = "")

names(g_r_int) <- paste('g_', nms, sep = "")
names(s_r_int) <- paste('s_', nms, sep = "")
names(f_s_r_int) <- paste('f_s_', nms, sep = "")

g_params <- list2(!!! g_r_int)
s_params <- list2(!!! s_r_int)
f_s_params <- list2(!!! f_s_r_int)

params <- splice(data_list, g_params, s_params, f_s_params)


b <- seq(0.2, 400, length.out = 501)
sv1 <- sv2 <- (b[2:501] + b[1:500]) * 0.5
h <- sv1[2] - sv1[1]

g_1 <- h * outer(sv1, sv2, FUN = g, params = c(params$g_int,
                                               params$g_slope,
                                               params$sd_g),
                 r_effect = params$g_r_1)

g_2 <- h * outer(sv1, sv2, FUN = g, params = c(params$g_int,
                                               params$g_slope,
                                               params$sd_g),
                 r_effect = params$g_r_2)
g_3 <- h * outer(sv1, sv2, FUN = g, params = c(params$g_int,
                                        params$g_slope,
                                        params$sd_g),
          r_effect = params$g_r_3)

g_4 <- h * outer(sv1, sv2, FUN = g, params = c(params$g_int,
                                               params$g_slope,
                                               params$sd_g),
                 r_effect = params$g_r_5)

g_5 <- h * outer(sv1, sv2, FUN = g, params = c(params$g_int,
                                               params$g_slope,
                                               params$sd_g),
                 r_effect = params$g_r_5)

g_1 <- truncated_distributions(g_1, n_mesh_p = 500)
g_2 <- truncated_distributions(g_2, n_mesh_p = 500)
g_3 <- truncated_distributions(g_3, n_mesh_p = 500)
g_4 <- truncated_distributions(g_4, n_mesh_p = 500)
g_5 <- truncated_distributions(g_5, n_mesh_p = 500)

s_1 <- s(sv1, c(params$s_int, params$s_slope,
                params$f_r_int, params$f_r_slope),
         params$s_r_1)
s_2 <- s(sv1, c(params$s_int, params$s_slope,
                params$f_r_int, params$f_r_slope), params$s_r_2)
s_3 <- s(sv1, c(params$s_int, params$s_slope,
                params$f_r_int, params$f_r_slope), params$s_r_3)
s_4 <- s(sv1, c(params$s_int, params$s_slope,
                params$f_r_int, params$f_r_slope), params$s_r_4)
s_5 <- s(sv1, c(params$s_int, params$s_slope,
                params$f_r_int, params$f_r_slope), params$s_r_5)

P_1 <- t(s_1 * t(g_1))
P_2 <- t(s_2 * t(g_2))
P_3 <- t(s_3 * t(g_3))
P_4 <- t(s_4 * t(g_4))
P_5 <- t(s_5 * t(g_5))

F_1 <- h * outer(sv1, sv2, FUN = fec,
                 params = unlist(params[6:11]),
                 r_effect = params$f_s_r_1)

F_2 <- h * outer(sv1, sv2, FUN = fec,
                 params = unlist(params[6:11]),
                 r_effect = params$f_s_r_2)

F_3 <- h * outer(sv1, sv2, FUN = fec,
                 params = unlist(params[6:11]),
                 r_effect = params$f_s_r_3)

F_4 <- h * outer(sv1, sv2, FUN = fec,
                 params = unlist(params[6:11]),
                 r_effect = params$f_s_r_4)

F_5 <- h * outer(sv1, sv2, FUN = fec,
                 params = unlist(params[6:11]),
                 r_effect = params$f_s_r_5)


K_1 <- P_1 + F_1
K_2 <- P_2 + F_2
K_3 <- P_3 + F_3
K_4 <- P_4 + F_4
K_5 <- P_5 + F_5

lambdas <- c(Re(eigen(K_1)$values[1]),
             Re(eigen(K_2)$values[1]),
             Re(eigen(K_3)$values[1]),
             Re(eigen(K_4)$values[1]),
             Re(eigen(K_5)$values[1]))


## IPMR version
h_levels <- list(yr = 1:5)
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

monocarp_sys <- init_ipm('simple_di_stoch_kern') %>%
  define_kernel('P_yr',
                formula = t(s_yr * t(g_yr)) ,
                family = "CC",
                s_yr = inv_logit_r(ht_1, s_int,
                                   s_slope, s_r_yr) - inv_logit(ht_1, f_r_int,
                                                                f_r_slope),
                g_yr = dnorm(ht_2, mean = mu_g_yr, sd = sd_g) * cell_size_ht,
                mu_g_yr = g_int + g_slope * ht_1 + g_r_yr,
                data_list = params,
                state_list = list(c('ht')),
                uses_par_sets = TRUE,
                par_set_indices = h_levels,
                evict = TRUE,
                evict_fun = truncated_distributions(g_yr,
                                                    n_mesh_p = 500)) %>%
  define_kernel("F_yr",
                formula = f_r * f_s_yr * f_d,
                family = "CC",
                f_r = inv_logit(ht_1, f_r_int, f_r_slope),
                f_s_yr = pois_r(ht_1, f_s_int, f_s_slope, f_s_r_yr),
                f_d = dnorm(ht_2, mean = mu_fd, sd = sd_fd) * cell_size_ht,
                data_list = params,
                state_list = list(c('ht')),
                uses_par_sets = TRUE,
                par_set_indices = h_levels,
                evict = FALSE) %>%
  define_kernel('K_yr',
                formula = P_yr + F_yr,
                family = "IPM",
                data_list = params,
                state_list = list(c("ht")),
                uses_par_sets = TRUE,
                par_set_indices = h_levels,
                evict = FALSE) %>%
  define_domains(ht = c(0.2, 400, 500)) %>%
  make_ipm(return_all = TRUE,
           usr_funs = list(inv_logit = inv_logit,
                           inv_logit_r = inv_logit_r,
                           pois_r = pois_r))


ks <- monocarp_sys$iterators
l_bc <- vapply(ks, function(x) Re(eigen(x)$values[1]), numeric(1))

l_bc - lambdas
