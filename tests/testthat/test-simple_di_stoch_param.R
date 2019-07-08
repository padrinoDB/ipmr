# test simple_di_stoch_param
library(mvtnorm)
library(rlang)

context('Simple density independent stochastic parameter resampled models')

s <- function(s_int_yr, s_slope, sv1) {

  return(1/(1 + exp(-(s_int_yr + s_slope * sv1))))

}

g <- function(g_int_yr, g_slope, g_sd, sv1, sv2) {

  mu <- g_int_yr + g_slope * sv1

  return(dnorm(sv2, mean = mu, sd = g_sd))

}

f_r <- function(f_r_int_yr, f_r_slope, sv1) {
  return(1/(1 + exp(-(f_r_int_yr + f_r_slope * sv1))))
}

f_s <- function(f_s_int_yr, f_s_slope, sv1) {
  return(exp(f_s_int_yr + f_s_int_slope * sv1))
}

f_d <- function(f_d_mu, f_d_sd, sv2) {

  return(dnorm(sv2, mean = f_d_mu, sd = f_d_sd))

}

data_list <- list(s_slope = 0.9,
                  g_slope = 0.99,
                  f_r_slope = 0.003,
                  f_s_slope = 0.001,
                  f_d_mu = 0.9,
                  f_d_sd = 0.1)

b <- seq(0, 10, length.out = 101)
d1 <- d2 <- (b[2:101] + b[1:100]) * 0.5
h <- d1[3] - d1[2]

r_means <- c(s_int_yr = 0.8,
             g_int_yr = 0.1,
             f_r_int_yr = 1.2,
             f_s_int_yr = 1.9)

set.seed(5000)
r_sigma <- runif(16)
r_sigma <- matrix(r_sigma, nrow = 4)

r_sigma <- r_sigma %*% t(r_sigma)


t_1_ests <- rmvnorm(1, mean = r_means, sigma = r_sigma)





