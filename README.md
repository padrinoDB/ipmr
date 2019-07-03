
[![Travis build
status](https://travis-ci.org/levisc8/ipmr.svg?branch=master)](https://travis-ci.org/levisc8/ipmr)

# ipmr

Simple, density-independent, deterministic and kernel-resampled
stochastic models (`simple_di_det`, `simple_di_stoch_kern`) are now
functional. However, expect changes as more complicated methods are
implemented\! See below for an example of how to implement an IPM in
this framework.

Next on the implementation to-do list are simple, stochastic IPMs using
parameter
resampling.

``` r
# Example of the setup for a simple IPM without density dependence or environmental
# stochasticity

library(ipmr)

data_list = list(s_int = 2.2,
                 s_slope = 0.25,
                 g_int = 0.2,
                 g_slope = 1.02,
                 sd_g = 0.7,
                 f_r_int = 0.03,
                 f_r_slope = 0.015,
                 f_s_int = 1.3,
                 f_s_slope = 0.075,
                 mu_fd = 0.5,
                 sd_fd = 0.2)

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

b <- seq(0, 50, length.out = 101)
d1 <- d2 <- (b[2:101] + b[1:100]) * 0.5
h <- d1[3] - d1[2]


G <- h * outer(d1, d2, FUN = g, params = c(data_list$g_int,
                                           data_list$g_slope,
                                           data_list$sd_g))
G2 <- G/matrix(as.vector(apply(G, 2, sum)),
               nrow = length(d1),
               ncol = length(d1),
               byrow = TRUE)

S <- s(d1, c(data_list$s_int, data_list$s_slope))

P <- t(S * t(G2))

Fm <- h * outer(d1, d2, FUN = fec, params = unlist(data_list[6:11]))

K <- P + Fm

lambda_usr <- Re(eigen(K)$values[1])
w_usr <- Re(eigen(K)$vectors[ , 1])


# User specified functions can be passed to make_ipm(usr_funs = list(my_fun = my_fun)).
# inv_logit is a simple example, but more complicated ones can be specified as well. 

inv_logit <- function(int, slope, sv) {
  return(1/(1 + exp(-(int + slope * sv))))
}

state_list <- list(c('dbh'))

x <- init_ipm('simple_di_det') %>%
  define_kernel(
    # Name of the kernel
    name = "P",
    # The type of transition it describes (e.g. continuous - continuous, discrete - continuous)
    family = "CC",
    # The formula for the kernel. don't forget to transpose the growth matrix when multiplying survival!
    formula = t(s * t(g)),
    # A named set of expressions for the vital rates it includes. ipmr automatically
    # computes the bin width creates a variable called "cell_size_stateVariable"
    # for usage in expressions that require integrating a density function. NOTE - the
    # cell_size_variable will not be required much longer as ipmr should be able
    # to detect when it is needed and insert it automatically!
    
    s = inv_logit(s_int, s_slope, dbh_1), # note the use of user-specified function here
    g = cell_size_dbh * dnorm(dbh_2, mu_g, sd_g),
    mu_g = g_int + g_slope * dbh_1,
    data_list = list(s_int = 2.2,
                     s_slope = 0.25,
                     g_int = 0.2,
                     g_slope = 1.02,
                     sd_g = 0.7),
    # The implementation details will soon be split out into another function
    # that matches integration rules and domains to kernels
    state_list = state_list,
    has_hier_effs = FALSE,
    # The evict_fun argument can take any function. ipmr provides built in
    # truncated density functions, but a user specified one will work as well, 
    # provided all parameters are provided in the data_list argument.
    evict = TRUE,
    evict_fun = truncated_distributions(g,
                                        n_mesh_p = 100)) %>%
  
  # Define the fecundity kernel
  define_kernel(
    name = 'F',
    formula = f_r * f_s * f_d,
    family = 'CC',
    f_r = inv_logit(f_r_int, f_r_slope, dbh_1),
    f_s = exp(f_s_int + f_s_slope * dbh_1),
    f_d = cell_size_dbh * dnorm(dbh_2, mu_fd, sd_fd),
    data_list = list(f_r_int = 0.03,
                     f_r_slope = 0.015,
                     f_s_int = 1.3,
                     f_s_slope = 0.075,
                     mu_fd = 0.5,
                     sd_fd = 0.2),
    state_list = state_list,
    evict = FALSE
  ) %>%
  
  # K kernels will likely get their own define_k() function, so don't expect
  # this interface to remain consistent
  define_k(
    name = "K",
    formula = P + F,
    family = "IPM",
    state_list = state_list,
    evict = FALSE # this is dealt with in the sub-kernels, so no need to do so again.
  ) %>%
  define_impl(
    kernel_impl_list =  make_impl_args_list(
      kernel_names = c("K", "P", "F"),
      int_rule = rep('midpoint', 3),
      dom_start = rep('dbh', 3),
      dom_end = rep('dbh', 3)
    )
  ) %>% 
  define_domains(
    dbh = c(0, # the first entry is the lower bound of the domain.
            50, # the second entry is the upper bound of the domain.
            100 # third entry is the number of meshpoints for the domain.
    )
  ) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit))

lambda_ipmr <- Re(eigen(x$iterators$K)$values[1])
w_ipmr <- Re(eigen(x$iterators$K)$vectors[ , 1])

lambda_ipmr - lambda_usr
```

## Simple, density independent, stochastic kernel resampling models

These models are typically the result of vital rate models that are fit
in a mixed effects framework (e.g. multiple sites or multiple years of
data). They have a special syntax that mirrors the mathematical notation
of these models (and has the side effect of saving you a considerable
amount of copying/pasting/typing in general).

The syntax uses a `name_hierarchicalVariable` notation. These names are
automatically expanded to include the multiple levels of the
hierarchical variable. For example, the P kernel for a model with 5
years of data could be `name`’d `P_yr`, and the underlying vital rates
would be denoted `vitalRate_yr`. The parameter values in the `data_list`
are the exception to this - they must be labeled with the actual values
that the hierarchical variable takes. See the example
below.

``` r
# rlang is a useful shortcut for splicing named values into lists. purrr is used to manipulate said lists.
# This is intended to simulate a monocarpic perennial life history where flowering is always fatal.
# Note that this means the survival function also includes the probability of reproduction function.

library(rlang)
library(ipmr)
library(purrr)

# define functions for target ipm

# Survival - logistic regression
s <- function(sv1, params, r_effect) {
  1/(1 + exp(-(params[1] + params[2] * sv1 + r_effect))) *
    (1 - f_r(sv1, params[3:4]))
}

# Growth
g <- function(sv1, sv2, params, r_effect) {
  mu <- params[1] + params[2] * sv1 + r_effect
  dnorm(sv2, mean = mu, sd = params[3])
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
f_d <- function(sv2, params) {
  dnorm(sv2, mean = params[1], sd = params[2])
}

# constructor function for the F kernel
fec <- function(sv1, sv2, params, r_effect) {
  f_r(sv1, params[1:2]) * f_s(sv1, params[3:4], r_effect) * f_d(sv2, params[5:6])
}


set.seed(50127)


# Define some fixed parameters
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

# Now, simulate some random intercepts for growth, survival, and offspring production
g_r_int <- rnorm(5, 0, 0.3)
s_r_int <- rnorm(5, 0, 0.7)
f_s_r_int <- rnorm(5, 0, 0.2)

nms <- paste("r_", 1:5, sep = "")

names(g_r_int) <- paste('g_', nms, sep = "")
names(s_r_int) <- paste('s_', nms, sep = "")
names(f_s_r_int) <- paste('f_s_', nms, sep = "")

# The !!! operator used inside of list2 from rlang takes the named vector
# and converts it to a named list. This can be spliced into the data list 
# to rapidly make a parameter set suitable for usage in the data_list argument
# of define_kernel

g_params <- list2(!!! g_r_int)
s_params <- list2(!!! s_r_int)
f_s_params <- list2(!!! f_s_r_int)

params <- splice(data_list, g_params, s_params, f_s_params)


b <- seq(0.2, 400, length.out = 501)
sv1 <- sv2 <- (b[2:501] + b[1:500]) * 0.5
h <- sv1[2] - sv1[1]

# repetitive to demonstrate the typical kernel construction process.

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
                 r_effect = params$g_r_4)

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
                params$f_r_int, params$f_r_slope), params$s_r_1)
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

# These are not corrected for eviction, but they probably should be
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


## ipmr version

# define the levels of the hierarchical variable and save them in a named
# list that corresponds to the suffix in the kernel notation
hier_levels <- list(yr = 1:5)

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


monocarp_sys <- init_ipm('simple_di_stoch_kern') %>%
  define_kernel(
    # The yr suffix is appended to the kernel name and the parameter names
    # within each vital rate expression. ipmr substitutes in the hier_levels
    # for each suffix occurrence, thus changing P_yr in P_1, P_2, P_3, P_4, P_5,
    # s_yr in s_1, s_2, s_3, s_4, and s_5. s_r_yr is converted to s_r_1, s_r_2, 
    # etc. In the case of s_r_yr, provided that the names in the data_list match
    # the expanded names, all will go well!
    
    name = 'P_yr',
    formula = t(s_yr * t(g_yr)) ,
    family = "CC",
    s_yr = inv_logit_r(ht_1, s_int, s_slope, s_r_yr) * 
      (1 - inv_logit(ht_1, f_r_int, f_r_slope)),
    g_yr = dnorm(ht_2, mean = mu_g_yr, sd = sd_g) * cell_size_ht,
    mu_g_yr = g_int + g_slope * ht_1 + g_r_yr,
    data_list = params,
    state_list = list(c('ht')),
    has_hier_effs = TRUE,
    levels_hier_effs = hier_levels,
    evict = TRUE,
    # Note that the suffix is appended here since the growth kernel also has a random intercept.
    evict_fun = truncated_distributions(g_yr,
                                        n_mesh_p = 500)
  ) %>%
  define_kernel(
    "F_yr",
    formula = f_r * f_s_yr * f_d,
    family = "CC",
    f_r = inv_logit(ht_1, f_r_int, f_r_slope),
    f_s_yr = pois_r(ht_1, f_s_int, f_s_slope, f_s_r_yr),
    f_d = dnorm(ht_2, mean = mu_fd, sd = sd_fd) * cell_size_ht,
    data_list = params,
    state_list = list(c('ht')),
    has_hier_effs = TRUE,
    levels_hier_effs = hier_levels,
    evict = FALSE) %>%
  define_k(
    'K_yr',
    formula = P_yr + F_yr,
    family = "IPM",
    data_list = params,
    state_list = list(c("ht")),
    has_hier_effs = TRUE,
    levels_hier_effs = hier_levels,
    evict = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("K_yr", "P_yr", "F_yr"),
      int_rule = rep("midpoint", 3),
      dom_start = rep("ht", 3),
      dom_end = rep("ht", 3)
    )
  ) %>%
  define_domains(ht = c(0.2, 400, 500)) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit,
                           inv_logit_r = inv_logit_r,
                           pois_r = pois_r))


ks <- monocarp_sys$iterators
l_bc <- vapply(ks, function(x) Re(eigen(x)$values[1]), numeric(1))

l_bc - lambdas
```
