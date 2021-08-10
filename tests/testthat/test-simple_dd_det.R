
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

x <- init_ipm(sim_gen    = "simple",
              di_dd      = "dd",
              det_stoch  = "det") %>%
  define_kernel("P",
                formula = s * g,
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
                evict_fun = truncated_distributions("norm", "f_d"))  %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P","F"),
      int_rule = rep("midpoint",2),
      state_start = rep("dbh", 2),
      state_end   = rep("dbh", 2)
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

ipmr_lam <- lambda(x)
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

test_that("classes are correctly set", {

  test_ind <- vapply(x$sub_kernels,
                     function(x) inherits(x, "ipmr_matrix"),
                     logical(1L))

  expect_true(all(test_ind))
  expect_s3_class(x, "simple_dd_det_ipm")
  expect_s3_class(x, "ipmr_ipm")

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

x_sub <- init_ipm(sim_gen    = "simple",
                  di_dd      = "dd",
                  det_stoch  = "det") %>%
  define_kernel("P",
                formula = s * g,
                family = "CC",
                s = plogis(s_int + s_slope * dbh_1 + dd_slope * sum(n_dbh_t[51:100]) * d_dbh),
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
                f_s = exp(f_s_int + f_s_slope * dbh_1 + dd_slope * sum(n_dbh_t[51:100]) * d_dbh),
                f_d = dnorm(dbh_2, mu_fd, sd_fd),
                data_list = data_list,
                states = state_list,
                evict_cor = TRUE,
                evict_fun = truncated_distributions("norm", "f_d")) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P","F"),
      int_rule = rep("midpoint", 2),
      state_start = rep("dbh", 2),
      state_end   = rep("dbh", 2)
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

ipmr_lam <- lambda(x_sub)
ipmr_pop_sizes <- colSums(x_sub$pop_state$n_dbh) * h

hand_lam <- colSums(pop_holder[ , 2:101]) / colSums(pop_holder[ , 1:100])
hand_pop_sizes <- colSums(pop_holder) * h

test_that("asymptotic behavior is preserved at every time step (subsetted models)", {

  expect_equal(as.vector(ipmr_lam), hand_lam)
  expect_equal(ipmr_pop_sizes, hand_pop_sizes)

  expect_true(all(ipmr_lam > 0))
  expect_true(all(ipmr_pop_sizes > 0))

})



# par_setarachical density dependent models ----------

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

g_ints <- list(g_int_1 = -0.2,
               g_int_2 = 0.651,
               g_int_3 = -0.512)

data_list <- c(data_list, g_ints)

state_list <- list(c("dbh"))

# set.seed(1231241)

par_set_mod <- init_ipm(sim_gen    = "simple",
                     di_dd      = "dd",
                     det_stoch  = "det") %>%
  define_kernel("P_yr",
                formula = s_yr * g_yr,
                family = "CC",
                s_yr = plogis(s_int + s_slope * dbh_1 + dd_slope * sum(n_dbh_yr_t)),
                g_yr = dnorm(dbh_2, mu_g_yr, sd_g),
                mu_g_yr = g_int + g_int_yr + g_slope * dbh_1,
                data_list = data_list,
                states = state_list,
                uses_par_sets = TRUE,
                par_set_indices = list(yr = 1:3),
                evict_cor = TRUE,
                evict_fun = truncated_distributions("norm", "g_yr")) %>%
  define_kernel('F_yr',
                formula = f_r * f_s_yr * f_d,
                family = 'CC',
                f_r = plogis(f_r_int + f_r_slope * dbh_1),
                f_s_yr = exp(f_s_int + f_s_slope * dbh_1 + dd_slope * sum(n_dbh_yr_t)),
                f_d = dnorm(dbh_2, mu_fd, sd_fd),
                data_list = data_list,
                states = state_list,
                uses_par_sets = TRUE,
                par_set_indices = list(yr = 1:3),
                evict_cor = TRUE,
                evict_fun = truncated_distributions("norm", "f_d")) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P_yr","F_yr"),
      int_rule = rep("midpoint",2),
      state_start = rep("dbh", 2),
      state_end   = rep("dbh", 2)
    )
  ) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  define_pop_state(n_dbh_yr = runif(100)) %>%
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

pop_holder <- lapply(1:3,
                     function(x, ps){
                       temp <- matrix(NA_real_,
                                      nrow = 100,
                                      ncol = 101)
                       temp[, 1] <- ps[[x]][ , 1]
                       return(temp)
                     },
                     ps = par_set_mod$pop_state[1:3])

names(pop_holder) <- paste("n_dbh_", 1:3, sep = "")

L <- 0
U <- 50
n <- 100

bounds  <- seq(L, U, length.out = n + 1)
z <- z1 <- (bounds[2:101] + bounds[1:100]) * 0.5

pop_size <- lapply(pop_holder, function(x) sum(x[ , 1]))

for(yr in 1:3) {

  temp_data <- c(data_list[1:12], data_list[grepl(yr, names(data_list))])

  names(temp_data) <- gsub("_[0-9]", "", names(temp_data))

  for(i in 2:101) {

    k <- k_dd(z1, z, par_list = temp_data, U, L, pop_size[[yr]])

    pop_holder[[yr]][ , i] <- k %*% pop_holder[[yr]][ , (i - 1)]

    pop_size[[yr]] <- sum(pop_holder[[yr]][ , i])

  }
}


ipmr_pop_sizes <- lapply(pop_holder, colSums)
hand_pop_sizes <- lapply(pop_holder, colSums)



test_that("population behavior is correctly modeled", {

  expect_equal(as.vector(ipmr_pop_sizes$n_dbh_1), hand_pop_sizes$n_dbh_1)
  expect_equal(as.vector(ipmr_pop_sizes$n_dbh_2), hand_pop_sizes$n_dbh_2)
  expect_equal(as.vector(ipmr_pop_sizes$n_dbh_3), hand_pop_sizes$n_dbh_3)

})
