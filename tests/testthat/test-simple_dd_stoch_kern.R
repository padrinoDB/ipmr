
# Define some fixed parameters
data_list = list(
  s_int     = 1.03,
  s_slope   = 2.2,
  s_dd      = -0.7,
  g_int     = 8,
  g_slope   = 0.92,
  sd_g      = 0.9,
  f_r_int   = 0.09,
  f_r_slope = 0.05,
  f_s_int   = 0.1,
  f_s_slope = 0.005,
  f_s_dd    = -0.03,
  mu_fd     = 9,
  sd_fd     = 2
)

# Now, simulate some random intercepts for growth, survival, and offspring production

g_r_int   <- rnorm(5, 0, 0.3)
s_r_int   <- rnorm(5, 0, 0.7)
f_s_r_int <- rnorm(5, 0, 0.2)

nms <- paste("r_", 1:5, sep = "")

names(g_r_int) <- paste("g_", nms, sep = "")
names(s_r_int) <- paste("s_", nms, sep = "")
names(f_s_r_int) <- paste("f_s_", nms, sep = "")

params     <- c(data_list, g_r_int, s_r_int, f_s_r_int)

x <- init_ipm(sim_gen    = "simple",
              di_dd      = "dd",
              det_stoch  = "stoch",
              "kern") %>%
  define_kernel(
    name             = "P_yr",
    formula          = s_yr * g_yr,
    family           = "CC",
    s_yr             = plogis(s_int + s_r_yr + s_slope * size_1 + s_dd * sum(n_size_t)),
    g_yr             = dnorm(size_2, g_mu_yr, sd_g),
    g_mu_yr          = g_int + g_r_yr + g_slope * size_1,
    data_list        = params,
    states           = list(c("size")),
    uses_par_sets    = TRUE,
    par_set_indices = list(yr = 1:5),
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "g_yr")
  ) %>%
  define_kernel(
    name             = "F_yr",
    formula          = f_r * f_s_yr * f_d,
    family           = "CC",
    f_r              = plogis(f_r_int + f_r_slope * size_1),
    f_s_yr           = exp(f_s_int + f_s_r_yr + f_s_slope * size_1 + f_s_dd * sum(n_size_t)),
    f_d              = dnorm(size_2, mu_fd, sd_fd),
    data_list        = params,
    states           = list(c("size")),
    uses_par_sets    = TRUE,
    par_set_indices = list(yr = 1:5),
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "f_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P_yr", "F_yr"),
      int_rule     = rep("midpoint", 2),
      state_start    = rep("size", 2),
      state_end      = rep("size", 2)
    )
  ) %>%
  define_domains(
    size = c(0, 50, 200)
  ) %>%
  define_pop_state(
    n_size = runif(200)
  ) %>%
  make_ipm(
    iterate = TRUE,
    iterations = 50,
    kernel_seq = sample(1:5, 50, replace = TRUE)
  )


use_seq <- x$env_seq

g_z1z <- function(z1, z, par_list, L, U, yr) {

  g_int_r <- par_list[grepl(paste("g_r_", yr, sep = ""),
                            names(par_list))] %>%
    unlist()

  mu <- par_list$g_int + g_int_r + par_list$g_slope * z
  ev <- pnorm(U, mu, par_list$sd_g) - pnorm(L, mu, par_list$sd_g)

  out <- dnorm(z1, mu, par_list$sd_g) / ev

  return(out)

}

s_z <- function(z, par_list, yr, pop_size) {

  s_int_r <- par_list[grepl(paste("s_r_", yr, sep = ""),
                            names(par_list))] %>%
    unlist()

  out <- plogis(par_list$s_int + s_int_r + par_list$s_slope * z + par_list$s_dd * pop_size)

  return(out)

}

f_z1z <- function(z1, z, par_list, L, U, yr, pop_size) {

  f_s_r <- par_list[grepl(paste("f_s_r_", yr, sep = ""),
                          names(par_list))] %>%
    unlist()

  f_s <- exp(par_list$f_s_int + f_s_r + par_list$f_s_slope * z + par_list$f_s_dd * pop_size)

  f_r <- plogis(par_list$f_r_int + par_list$f_r_slope * z)

  ev <- pnorm(U, par_list$mu_fd, par_list$sd_fd) - pnorm(L, par_list$mu_fd, par_list$sd_fd)

  f_d <- dnorm(z1, par_list$mu_fd, par_list$sd_fd) / ev

  out <- f_r * f_s * f_d

  return(out)

}

k_dd <- function(z1, z, par_list, L, U, yr, pop_size) {

  g <- outer(z, z, FUN = g_z1z,
             par_list = par_list,
             L = L,
             U = U,
             yr = yr)
  s <- s_z(z, par_list, yr, pop_size)
  f <- outer(z, z, FUN = f_z1z,
             par_list = par_list,
             L = L,
             U = U,
             yr = yr,
             pop_size = pop_size)

  k <- t(s * t(g)) + f

  h <- z[2] - z[1]

  return(k * h)

}

pop_holder <- matrix(NA_real_,
                     nrow = 200,
                     ncol = 51)

L <- 0
U <- 50
n <- 200

bounds  <- seq(L, U, length.out = n + 1)
z <- z1 <- (bounds[2:201] + bounds[1:200]) * 0.5

pop_holder[ , 1] <- x$pop_state$n_size[ , 1]

pop_size <- sum(pop_holder[ , 1])

for(i in 2:51) {

  k <- k_dd(z1, z,
            par_list = params,
            L, U,
            yr = use_seq[(i - 1)],
            pop_size)

  pop_holder[ , i] <- k %*% pop_holder[ , (i - 1)]

  pop_size <- sum(pop_holder[ , i])

}


ipmr_lam <- lambda(x, type_lambda = "all") %>%
  as.vector()

ipmr_pop_sizes <- colSums(x$pop_state$n_size)

hand_lam <- colSums(pop_holder[ , 2:51]) / colSums(pop_holder[ , 1:50])
hand_pop_sizes <- colSums(pop_holder)

test_that("asymptotic behavior is preserved at every time step", {

  expect_equal(ipmr_lam, hand_lam, tolerance = 2e-2)
  expect_equal(ipmr_pop_sizes, hand_pop_sizes, tolerance = 1)

})


test_that("sub-kernel names and values are generated correctly", {

  p_rngs <- vapply(x$sub_kernels[grepl("P", names(x$sub_kernels))],
                   range,
                   numeric(2L))

  expect_true(all(p_rngs >=0 & p_rngs <= 1))

  nms <- vapply(1:5,
                function(x) paste(c("P", "F"), x, sep = "_"),
                character(2L)) %>%
    as.vector() %>%
    vapply(., function(x) paste(x,"it", 1:50, sep = "_"), character(50L)) %>%
    as.vector()

  expect_true(all(nms %in% names(x$sub_kernels)))

})

