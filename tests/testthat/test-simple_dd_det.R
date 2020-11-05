
data_list = list(s_int = 2.2,
                 s_slope = 0.25,
                 dd_slope = -0.022,
                 g_int = 0.2,
                 g_slope = 1.02,
                 sd_g = 0.7,
                 f_r_int = 0.03,
                 f_r_slope = 0.015,
                 f_s_int = 1.3,
                 f_s_slope = 0.075,
                 mu_fd = 0.5,
                 sd_fd = 0.2)

state_list <- list(c("dbh"))

# set.seed(1231241)

x <- init_ipm('simple_dd_det') %>%
  define_kernel("P",
                formula = s * g, # make helper for double transposes
                family = "CC",
                s = plogis(s_int + s_slope * dbh_1 + dd_slope * sum(n_dbh_t)),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * dbh_1,
                data_list = data_list,

                # split out implementation details into a separate
                # function - named lists??
                states = state_list,
                has_hier_effs = FALSE,
                evict_cor = TRUE,
                evict_fun = truncated_distributions("norm", "g")) %>%
  define_kernel('F',
                formula = f_r * f_s * f_d,
                family = 'CC',
                f_r = plogis(f_r_int + f_r_slope * dbh_1),
                f_s = exp(f_s_int + f_s_slope * dbh_1 + dd_slope * sum(n_dbh_t)),
                f_d = dnorm(dbh_2, mu_fd, sd_fd),
                data_list = data_list,
                states = state_list,
                evict_cor = TRUE,
                evict_fun = truncated_distributions("norm", "f_d")) %>%
  define_k(name = 'K',
           family = "IPM",
           K = P + F,
           n_dbh_t_1 = K %*% n_dbh_t,
           data_list = list(),
           states = state_list,
           evict_cor = FALSE) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P","F", "K"),
      int_rule = rep("midpoint",3),
      dom_start = rep("dbh", 3),
      dom_end   = rep("dbh", 3)
    )
  ) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  define_pop_state(n_dbh = runif(100)) %>%
  make_ipm(iterate = TRUE,
           iterations = 100,
           normalize_pop_size = FALSE)


# Hand implementation

g_z1z <- function(z1, z, par_list, U, L) {

  mu <- par_list$g_int + par_list$g_slope * z

  ev <- pnorm(U, mu, par_list$sd_g) - pnorm(L, mu, par_list$sd_g)

  dnorm(z1, mu, par_list$sd_g) / ev

}

s_z <- function(z, par_list, pop_size) {
  plogis(par_list$s_int + par_list$s_slope * z + par_list$dd_slope * pop_size)
}

f_z1z <- function(z1, z, par_list, U, L, pop_size) {

  f_r <- plogis(par_list$f_r_int + par_list$f_r_slope * z)
  f_s <- exp(par_list$f_s_int + par_list$f_s_slope * z + par_list$dd_slope * pop_size)
  ev  <- pnorm(U, par_list$mu_fd, par_list$sd_fd) - pnorm(L,
                                                          par_list$mu_fd,
                                                          par_list$sd_fd)

  f_d <- dnorm(z1, par_list$mu_fd, par_list$sd_fd) / ev

  return(f_r * f_s * f_d)

}

k_dd <- function(z1, z, par_list, U, L, pop_size) {

  g <- outer(z1, z, g_z1z, par_list = par_list, U = U, L = L)
  s <- s_z(z, par_list, pop_size)
  f <- outer(z1, z, f_z1z, par_list = par_list, U = U, L = L, pop_size = pop_size)

  out <- t(s * t(g)) + f

  return(out * (z[2] - z[1]))

}

pop_holder <- matrix(NA_real_,
                     nrow = 100,
                     ncol = 101)

L <- 0
U <- 50
n <- 100

bounds  <- seq(L, U, length.out = n + 1)
z <- z1 <- (bounds[2:101] + bounds[1:100]) * 0.5

pop_holder[ , 1] <- x$pop_state$n_dbh[ , 1]

pop_size <- sum(pop_holder[ , 1])

for(i in 2:101) {

  k <- k_dd(z1, z, par_list = data_list, U, L, pop_size)

  pop_holder[ , i] <- k %*% pop_holder[ , (i - 1)]

  pop_size <- sum(pop_holder[ , i])

}

ipmr_lam <- x$pop_state$lambda
ipmr_pop_sizes <- colSums(x$pop_state$n_dbh)

hand_lam <- colSums(pop_holder[ , 2:101]) / colSums(pop_holder[ , 1:100])
hand_pop_sizes <- colSums(pop_holder)

test_that("asymptotic behavior is preserved at every time step", {

  expect_equal(as.vector(ipmr_lam), hand_lam)
  expect_equal(ipmr_pop_sizes, hand_pop_sizes)

})

test_that("sub-kernel names and values are generated correctly", {

  p_rngs <- vapply(x$sub_kernels[grepl("P", names(x$sub_kernels))],
                   range,
                   numeric(2L))

  expect_true(all(p_rngs >=0 & p_rngs <= 1))

  nms <- vapply(1:100,
                function(x) paste(c("P_it", "F_it"), x, sep = "_"),
                character(2L)) %>%
    as.vector()

  expect_equal(nms, names(x$sub_kernels))

})



data_list = list(s_int = 2.2,
                 s_slope = 0.25,
                 dd_slope = -0.022,
                 g_int = 0.2,
                 g_slope = 1.02,
                 sd_g = 0.7,
                 f_r_int = 0.03,
                 f_r_slope = 0.015,
                 f_s_int = 1.3,
                 f_s_slope = 0.075,
                 mu_fd = 0.5,
                 sd_fd = 0.2)

state_list <- list(c("dbh"))

# set.seed(1231241)

x_sub <- init_ipm('simple_dd_det') %>%
  define_kernel("P",
                formula = s * g, # make helper for double transposes
                family = "CC",
                s = plogis(s_int + s_slope * dbh_1 + dd_slope * sum(n_dbh_t[51:100]) * d_dbh),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * dbh_1,
                data_list = data_list,

                # split out implementation details into a separate
                # function - named lists??
                states = state_list,
                has_hier_effs = FALSE,
                evict_cor = TRUE,
                evict_fun = truncated_distributions("norm", "g")) %>%
  define_kernel('F',
                formula = f_r * f_s * f_d,
                family = 'CC',
                f_r = plogis(f_r_int + f_r_slope * dbh_1),
                f_s = exp(f_s_int + f_s_slope * dbh_1 + dd_slope * sum(n_dbh_t[51:100]) * d_dbh),
                f_d = dnorm(dbh_2, mu_fd, sd_fd),
                data_list = data_list,
                states = state_list,
                evict_cor = TRUE,
                evict_fun = truncated_distributions("norm", "f_d")) %>%
  define_k(name = 'K',
           family = "IPM",
           K = P + F,
           n_dbh_t_1 = K %*% n_dbh_t,
           data_list = list(),
           states = state_list,
           evict_cor = FALSE) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P","F", "K"),
      int_rule = rep("midpoint",3),
      dom_start = rep("dbh", 3),
      dom_end   = rep("dbh", 3)
    )
  ) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  define_pop_state(n_dbh = runif(100)) %>%
  make_ipm(iterate = TRUE,
           iterations = 100,
           normalize_pop_size = FALSE)


# Hand implementation

g_z1z <- function(z1, z, par_list, U, L) {

  mu <- par_list$g_int + par_list$g_slope * z

  ev <- pnorm(U, mu, par_list$sd_g) - pnorm(L, mu, par_list$sd_g)

  dnorm(z1, mu, par_list$sd_g) / ev

}

s_z <- function(z, par_list, pop_size) {
  plogis(par_list$s_int + par_list$s_slope * z + par_list$dd_slope * pop_size)
}

f_z1z <- function(z1, z, par_list, U, L, pop_size) {

  f_r <- plogis(par_list$f_r_int + par_list$f_r_slope * z)
  f_s <- exp(par_list$f_s_int + par_list$f_s_slope * z + par_list$dd_slope * pop_size)
  ev  <- pnorm(U, par_list$mu_fd, par_list$sd_fd) - pnorm(L,
                                                          par_list$mu_fd,
                                                          par_list$sd_fd)

  f_d <- dnorm(z1, par_list$mu_fd, par_list$sd_fd) / ev

  return(f_r * f_s * f_d)

}

k_dd <- function(z1, z, par_list, U, L, pop_size) {

  g <- outer(z1, z, g_z1z, par_list = par_list, U = U, L = L)
  s <- s_z(z, par_list, pop_size)
  f <- outer(z1, z, f_z1z, par_list = par_list, U = U, L = L, pop_size = pop_size)

  out <- t(s * t(g)) + f

  return(out * (z[2] - z[1]))

}

pop_holder <- matrix(NA_real_,
                     nrow = 100,
                     ncol = 101)

L <- 0
U <- 50
n <- 100

bounds  <- seq(L, U, length.out = n + 1)
z <- z1 <- (bounds[2:101] + bounds[1:100]) * 0.5
h <- z[2] - z[1]

pop_holder[ , 1] <- x_sub$pop_state$n_dbh[ , 1]

pop_size <- sum(pop_holder[51:100 , 1]) * h

for(i in 2:101) {

  k <- k_dd(z1, z, par_list = data_list, U, L, pop_size)

  pop_holder[ , i] <- k %*% pop_holder[ , (i - 1)]

  pop_size <- sum(pop_holder[51:100 , i]) * h

}

ipmr_lam <- x_sub$pop_state$lambda
ipmr_pop_sizes <- colSums(x_sub$pop_state$n_dbh) * h

hand_lam <- colSums(pop_holder[ , 2:101]) / colSums(pop_holder[ , 1:100])
hand_pop_sizes <- colSums(pop_holder) * h

test_that("asymptotic behavior is preserved at every time step (subsetted models)", {

  expect_equal(as.vector(ipmr_lam), hand_lam)
  expect_equal(ipmr_pop_sizes, hand_pop_sizes)

  expect_true(all(ipmr_lam > 0))
  expect_true(all(ipmr_pop_sizes > 0))

})



# Rees + Rose density dependent weevil infestation model. Here,
# we rebuild that to make sure we're reproducing that behavior.
# Consider adding to DD vignette

m_par <- c(
  ## herbivory intercept
  weevil_int  = -17
)

## Growth function, pdf of size z1 given size z
G_z1z <- function(z1, z, m_par) {

  mu         <- 0.83 + 0.69 * z
  sig        <- sqrt(0.19)
  p.den.grow <- dnorm(z1, mean = mu, sd = sig)

  return(p.den.grow)
}

## Survival function, logistic regression

s_z <- function(z, m_par) {
  linear.p <- -0.62 + 0.85 * z
  p        <- plogis(linear.p)
  return(p)
}

## Probability of flowering function, logistic regression
p_bz <- function(z, m_par) {

  linear.p <- -10.22 + 4.25 * z
  p        <- plogis(linear.p)

  return(p)
}

## Seed production function of a size z plant, taking into account
## size dependent mean weevil load and negative exponential distribution
## of weevils. See Rose et al. (2005) Appendix, Ecological Archives E086-022-A1

b_z <- function(z, m_par) {

  eps <- exp(m_par["weevil_int"] + 1.71 * z);
  N   <- exp(-0.55 + 2.02 * z) / ((1 + eps / 16) ^ 0.32)

  return(N)
}

## Recruit size distribution
c_0z1 <- function(z1, m_par) {

  pRecr <- dnorm(z1, mean = 0.75, sd = sqrt(0.17))

  return(pRecr)
}

## Define the survival kernel function
P_z1z <- function (z1, z, m_par) {

  return((1 - p_bz(z, m_par)) * s_z(z, m_par) * G_z1z(z1, z, m_par))

}

## Define the total recruitment function
B_z <- function (z, m_par) {

  return( p_bz(z, m_par) * b_z(z, m_par) )

}

mk_P <- function(m, m_par, L, U) {

  h       <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h

  P       <- h * (outer(meshpts, meshpts, P_z1z, m_par = m_par))

  return(list(meshpts = meshpts, P = P))
}

## Function to iterate the IPM using midpoint rule
## with mesh points meshpts
Iterate <- function(nt,meshpts,P,m_par) {

  h        <- meshpts[2] - meshpts[1]
  N        <- h*sum(nt * B_z(meshpts,m_par))
  recruits <- c_0z1(meshpts,m_par) * (N ^ 0.67)

  return(recruits + P %*% nt)

}


weevilLoad <- function(nt, meshpts, m.par) {

  nFlower <- h * sum(nt * p_bz(meshpts))
  eps     <- exp(m_par["weevil_int"] + 1.71 * meshpts)
  nWeevil <- h * sum(nt * p_bz(meshpts) * eps)

  return(nWeevil/nFlower)

}

int <- numeric(5L)
int[1] <- mean (c(-17.3,-17,-16.7))
int[2] <- mean (c(-3.74,-3.39,-2.81))
int[3] <- mean(c(-1.64,-0.96,-0.75))
int[4] <- mean(c(-1.24,-0.52))
int[5] <- 1.5; # beyond the range of the data


m_par["weevil_int"] <- int[1]

L <- -1.5
U <- 4.5
n <- 200

n_size <- matrix(1, n, 201)

consts  <- mk_P(n, m_par, L, U)
P       <- consts$P
mesh_pts <- consts$meshpts
h        <- mesh_pts[2] - mesh_pts[1]

for(i in seq_len(n)) {

  n_size[ ,(i + 1)] <- Iterate(n_size[ , i],
                              mesh_pts,
                              P,
                              m_par)

}

hand_pop_sizes <- colSums(n_size)
hand_lam       <- hand_pop_sizes[2:201] / hand_pop_sizes[1:200]

# IPMR constants

data_list <- list(
  g_int        = 0.83,
  g_slope      = 0.69,
  g_sd         = sqrt(0.19),
  s_int        = -0.62,
  s_slope      = 0.85,
  p_b_int      = -10.22,
  p_b_slope    = 4.25,
  weevil_slope = 1.71,
  f_s_int      = -0.55,
  f_s_slope    = 2.02,
  f_s_exp      = 0.32,
  recr_mu      = 0.75,
  recr_sd      = sqrt(0.17)
)

weevil_ints <- as.list(int) %>%
  setNames(
    paste("weevil_int_", 1:5, sep = "")
  )

data_list <- c(data_list, weevil_ints)

weevil_ipm <- init_ipm("simple_dd_det") %>%
  define_kernel(
    name      = "P",
    formula   = s * (1-p_b) * g,
    family    = "CC",
    s         = plogis(s_int + s_slope * z_1),
    g_mu      = g_int + g_slope * z_1,
    g         = dnorm(z_2, g_mu, g_sd),
    p_b       = plogis(p_b_int + p_b_slope * z_1),
    data_list = data_list,
    states    = list(c("z")),
    evict_cor = FALSE,
    integrate = TRUE
  ) %>%
  define_kernel(
    name       = "F",
    formula    = c_0 * N,
    family     = "CC",
    c_0        = dnorm(z_2, recr_mu, recr_sd),
    N          = (d_z * sum(n_z_t * B_z)) ^ 0.67,
    B_z        = p_b * b_z,
    p_b        = plogis(p_b_int + p_b_slope * z_1),
    b_z        = exp(f_s_int + f_s_slope * z_1) / ((1 + e_z/16) ^ (f_s_exp)),
    e_z        = exp(weevil_int + weevil_slope * z_1),
    weevil_int = weevil_int_1,
    data_list  = data_list,
    states     = list(c("z")),
    evict_cor  = FALSE,
    integrate  = FALSE
  ) %>%
  define_k(
    name      = "K",
    family    = "IPM",
    F_1       = F[ , 1],
    n_z_t_1   = P %*% n_z_t + F_1,
    data_list = data_list,
    states    = list(c("z")),
    evict_cor = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      c("P", "F", "K"),
      rep("midpoint",3),
      rep("z",3),
      rep("z", 3)
    )
  ) %>%
  define_domains(
    z = c(L, U, n)
  ) %>%
  define_pop_state(
    pop_vectors = list(
      n_z = n_size[ , 1]
    )
  ) %>%
  make_ipm(
    iterate = TRUE,
    iterations = 200,
    normalize_pop_size = FALSE,
    report_progress = TRUE
  )


ipmr_lam       <- as.vector(weevil_ipm$pop_state$lambda)
ipmr_pop_sizes <- colSums(weevil_ipm$pop_state$n_z)

test_that("recover population dynamics from weevil-thistle IPM", {

  expect_equal(ipmr_lam, hand_lam)
  expect_equal(ipmr_pop_sizes, hand_pop_sizes)

})
