
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
                uses_par_sets = FALSE,
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

