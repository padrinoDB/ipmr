# General_dd_stoch_kern -------

# Levin et al 2019 control treatment ligustrum obtusifolium
data_list_control <- list(
  g_int     = 5.781,
  g_slope   = 0.988,
  g_sd      = 20.55699,
  s_int     = -0.352,
  s_slope   = 0.122,
  s_slope_2 = -0.000213,
  s_dd      = -1.5e-7,
  f_r_int   = -11.46,
  f_r_slope = 0.0835,
  f_s_int   = 2.6204,
  f_s_slope = 0.01256,
  f_s_dd    = - 0.002,
  f_d_mu    = 5.6655,
  f_d_sd    = 2.0734,
  e_p       = 0.15,
  g_i       = 0.5067
)

L <- 1.02
U <- 624
n <- 100

g_ints <- rnorm(3, sd = 4) %>%
  as.list() %>%
  setNames(paste("g_int_", LETTERS[1:3], sep = ""))


f_r_ints <- rnorm(3, sd = 2) %>%
  as.list() %>%
  setNames(paste("f_r_int_", LETTERS[1:3], sep = ""))

data_list_par_set <- c(data_list_control, g_ints, f_r_ints)


init_pop_vec <- runif(100)
init_b <- 20

gen_dd_stoch_co <- init_ipm(sim_gen    = "general",
                            di_dd      = "dd",
                            det_stoch  = "stoch",
                            "kern") %>%
  define_kernel(
    name             = "P_site",
    family           = "CC",
    formula          = s * g_site * d_ht,
    s                = plogis(s_int + s_slope * ht_1 + s_slope_2 * ht_1^2 + s_dd * sum(n_ht_t)),
    g_site           = dnorm(ht_2, g_mu_site, g_sd),
    g_mu_site        = g_int + g_int_site + g_slope * ht_1,
    states           = list(c("ht", "b")),
    data_list        = data_list_par_set,
    uses_par_sets    = TRUE,
    par_set_indices = list(site = LETTERS[1:3]),
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "g_site")
  ) %>%
  define_kernel(
    name             = "go_discrete_site",
    family           = "CD",
    formula          = f_r_site * f_s * g_i * d_ht,
    f_r_site         = plogis(f_r_int + f_r_int_site + f_r_slope * ht_1),
    f_s              = exp(f_s_int + f_s_slope * ht_1 + f_s_dd * sum(n_ht_t)),
    states           = list(c("ht", "b")),
    data_list        = data_list_par_set,
    uses_par_sets    = TRUE,
    par_set_indices = list(site = LETTERS[1:3]),
    evict_cor        = FALSE
  ) %>%
  define_kernel(
    name          = "leave_discrete",
    family        = "DC",
    formula       = e_p * f_d * d_ht,
    f_d           = dnorm(ht_2, f_d_mu, f_d_sd),
    states        = list(c("ht", "b")),
    data_list     = data_list_par_set,
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "f_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P_site", "go_discrete_site", "leave_discrete"),
      int_rule     = rep("midpoint", 3),
      state_start    = c("ht", "ht", "b"),
      state_end      = c("ht", "b", "ht")
    )
  ) %>%
  define_pop_state(
    pop_vectors = list(
      n_ht = init_pop_vec,
      n_b  = init_b
    )
  ) %>%
  define_domains(
    ht = c(L, U, n)
  ) %>%
  make_ipm(
    iterate    = TRUE,
    iterations = 50,
    kernel_seq = sample(LETTERS[1:3], 50, TRUE)
  )


# Hand version


kern_seq <- gen_dd_stoch_co$env_seq

pop_holder <- lapply(gen_dd_stoch_co$pop_state[1:2],
                     function(x) {
                       out <- matrix(0, nrow = dim(x)[1], ncol = 51)
                       out[ , 1] <- x[ , 1]

                       return(out)
                     })

iterate_model <- function(P, cd, dc, pop_holder, iteration) {

  pop_holder$n_ht[ , (iteration + 1)] <- P %*% pop_holder$n_ht[ , iteration] +
    dc %*% pop_holder$n_b[ , iteration]

  pop_holder$n_b[ , (iteration + 1)] <- cd %*% pop_holder$n_ht[ , iteration]

  return(pop_holder)

}


s_x <- function(int, slope1, slope2, s_dd, pop_size, sv1) {
  1/(1 + exp(-(int + slope1 * sv1 + slope2 * sv1^2 + s_dd * pop_size)))
}

f_r_x <- function(int, slope1,  sv1) {
  1/(1 + exp(-(int + slope1 * sv1)))
}


g_x <- function(sv2, sv1, int, slope, gsd, L, U) {
  mu <- int + slope * sv1
  ev <- pnorm(U, mu, gsd) - pnorm(L, mu, gsd)
  out <- dnorm(sv2, mean = mu, sd = gsd) / ev

  return(out)
}

b        <- seq(L, U, length.out = n + 1)
d1 <- d2 <- (b[2:(n + 1)] + b[1:n]) * 0.5
h        <- d1[3] - d1[2]

domains <- expand.grid(list(d2 = d1, d1 = d1))


for(i in seq_len(50)) {

  pop_size <- sum(pop_holder$n_ht[ , i])

  g_int_yr   <- g_ints[grepl(kern_seq[i],   names(g_ints))] %>% unlist()
  f_r_int_yr <- f_r_ints[grepl(kern_seq[i], names(f_r_ints))] %>% unlist()

  G_co <- h * g_x(domains$d1, domains$d2,
                  int   = data_list_control$g_int + g_int_yr,
                  slope = data_list_control$g_slope,
                  gsd   = data_list_control$g_sd,
                  L = L,
                  U = U)

  s_co <- s_x(data_list_control$s_int,
              data_list_control$s_slope,
              data_list_control$s_slope_2,
              data_list_control$s_dd,
              pop_size,
              d1)

  P_co <- s_co * G_co
  P_co <- matrix(P_co, nrow = 100, ncol = 100, byrow = TRUE)


  cd_co <- f_r_x(data_list_control$f_r_int + f_r_int_yr,
                 data_list_control$f_r_slope,
                 d1) *
    exp(data_list_control$f_s_int +
          data_list_control$f_s_slope * d1 +
          data_list_control$f_s_dd * pop_size) *
    data_list_control$g_i %>%
    matrix(ncol = n,
           nrow = 1) * h

  dc_co <- (dnorm(d2, mean = data_list_control$f_d_mu, sd = data_list_control$f_d_sd) /
              (pnorm(U,
                     data_list_control$f_d_mu,
                     data_list_control$f_d_sd) - pnorm(L,
                                                       data_list_control$f_d_mu,
                                                       data_list_control$f_d_sd))) *
    data_list_control$e_p *
    h %>%
    matrix(nrow = n, ncol = 1)

  pop_holder <- iterate_model(P_co, cd_co, dc_co, pop_holder, i)

}


ipmr_pop_sizes <- lapply(gen_dd_stoch_co$pop_state[1:2],
                         colSums) %>%
  do.call(what = `+`, args = .)

hand_pop_sizes <- lapply(pop_holder,
                         colSums) %>%
  do.call(what = `+`, args = .)

ipmr_lams <- lambda(gen_dd_stoch_co, type_lambda = "all") %>%
  as.vector()

hand_lams <- hand_pop_sizes[2:51] / hand_pop_sizes[1:50]

test_that("models compute population sizes correctly", {

  expect_equal(ipmr_pop_sizes, hand_pop_sizes)
  expect_equal(ipmr_lams, hand_lams)

})

test_that("classes are correctly set", {

  test_ind <- vapply(gen_dd_stoch_co$sub_kernels,
                     function(x) inherits(x, "ipmr_matrix"),
                     logical(1L))

  expect_true(all(test_ind))
  expect_s3_class(gen_dd_stoch_co, "general_dd_stoch_kern_ipm")
  expect_s3_class(gen_dd_stoch_co, "ipmr_ipm")

})
