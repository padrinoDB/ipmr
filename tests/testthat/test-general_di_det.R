

set.seed(2312)
init_pop_vec <- runif(500)
init_b <- 20

# Levin et al 2019 control treatment ligustrum obtusifolium
data_list_control <- list(
  g_int     = 5.781,
  g_slope   = 0.988,
  g_sd      = 20.55699,
  s_int     = -0.352,
  s_slope   = 0.122,
  s_slope_2 = -0.000213,
  f_r_int   = -11.46,
  f_r_slope = 0.0835,
  f_s_int   = 2.6204,
  f_s_slope = 0.01256,
  f_d_mu    = 5.6655,
  f_d_sd    = 2.0734,
  e_p       = 0.15,
  g_i       = 0.5067
)

# Levin et al 2019 cr treatment ligustrum obtusifolium
data_list_cr <- list(
  g_int     = 7.229,
  g_slope   = 0.988,
  g_sd      = 21.72262,
  s_int     = 0.0209,
  s_slope   = 0.0831,
  s_slope_2 = -0.00012999,
  f_r_int   = -11.46,
  f_r_slope = 0.0835,
  f_s_int   = 2.6204,
  f_s_slope = 0.01256,
  f_d_mu    = 5.6655,
  f_d_sd    = 2.0734,
  e_p       = 0.15,
  g_i       = 0.5067
)

s_x <- function(int, slope1, slope2, sv1) {
  1/(1 + exp(-(int + slope1 * sv1 + slope2 * sv1^2)))
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


L <- 1.02
U <- 624
n <- 500

b        <- seq(L, U, length.out = n + 1)
d1 <- d2 <- (b[2:(n + 1)] + b[1:n]) * 0.5
h        <- d1[3] - d1[2]

domains <- expand.grid(list(d2 = d1, d1 = d1))
G <- h * g_x(domains$d1, domains$d2,
             int   = data_list_control$g_int,
             slope = data_list_control$g_slope,
             gsd   = data_list_control$g_sd,
             L = L,
             U = U)

s <- s_x(data_list_control$s_int,
         data_list_control$s_slope,
         data_list_control$s_slope_2,
         d1)



P <- s * G

cd_co <- f_r_x(data_list_control$f_r_int,
               data_list_control$f_r_slope,
               d1) *
  exp(data_list_control$f_s_int + data_list_control$f_s_slope * d1) *
  data_list_control$g_i %>%
  matrix(ncol = n,
         nrow = 1)

dc_co <- (dnorm(d2, mean = data_list_control$f_d_mu, sd = data_list_control$f_d_sd) /
            (pnorm(U,
                   data_list_control$f_d_mu,
                   data_list_control$f_d_sd) - pnorm(L,
                                                     data_list_control$f_d_mu,
                                                     data_list_control$f_d_sd))) *
  data_list_control$e_p *
  h %>%
  matrix(nrow = n, ncol = 1)

cc_co <- matrix(P, nrow = n, ncol = n, byrow = TRUE)
dd_co <- 0

K_co <- rbind(
  cbind(dd_co, cd_co),
  cbind(dc_co, cc_co)
)

actual_co <- Re(eigen(K_co)$values[1])

# Repeat, but with CRs

G <- h * g_x(domains$d1, domains$d2,
             int   = data_list_cr$g_int,
             slope = data_list_cr$g_slope,
             gsd   = data_list_cr$g_sd,
             L = L, U = U)

s <- s_x(data_list_cr$s_int,
         data_list_cr$s_slope,
         data_list_cr$s_slope_2,
         d1)



P <- s * G


cd_cr <- f_r_x(data_list_cr$f_r_int,
               data_list_cr$f_r_slope,
               d1) *
  exp(data_list_cr$f_s_int + data_list_cr$f_s_slope * d1) *
  data_list_cr$g_i %>%
  matrix(ncol = n,
         nrow = 1)

dc_cr <- (dnorm(d2, mean = data_list_cr$f_d_mu, sd = data_list_cr$f_d_sd) /
            (pnorm(U,
                   data_list_cr$f_d_mu,
                   data_list_cr$f_d_sd) - pnorm(L,
                                                data_list_cr$f_d_mu,
                                                data_list_cr$f_d_sd))) *
  data_list_cr$e_p *
  h %>%
  matrix(nrow = n, ncol = 1)

cc_cr <- matrix(P, nrow = n, ncol = n, byrow = TRUE)
dd_cr <- 0

K_cr <- rbind(
  cbind(dd_cr, cd_cr),
  cbind(dc_cr, cc_cr)
)

actual_cr <- Re(eigen(K_cr)$values[1])

target_cr <- 1.26
target_co <- 1.22


# Verify test implementation actually does the published one

test_that('test implementation matches published quantities', {
  expect_equal(target_cr, actual_cr, tolerance = 1e-2)
  expect_equal(target_co, actual_co, tolerance = 1e-2)
})

# Now, set up simulation to compare w/ ipmr

n_runs <- 100
pop_state <- list(ht_cr  = matrix(NA_real_, nrow = n, ncol = n_runs + 1),
                  b_t_cr = matrix(NA_real_, nrow = 1, ncol = n_runs + 1),
                  ht_co  = matrix(NA_real_, nrow = n, ncol = n_runs + 1),
                  b_t_co = matrix(NA_real_, nrow = 1, ncol = n_runs + 1))

pop_state$ht_cr[ , 1]  <- pop_state$ht_co[ , 1]  <- init_pop_vec
pop_state$b_t_cr[ , 1] <- pop_state$b_t_co[ , 1] <- init_b
pop_state$lambda_cr <- pop_state$lambda_co <- matrix(0, nrow = 1, ncol = 100)

DD_cr <- K_cr[1 , 1] # dd
DD_co <- K_co[1 , 1] # dd
DC_cr <- K_cr[2:501 , 1] # dc
DC_co <- K_co[2:501 , 1] # dc
CD_cr <- K_cr[1 , 2:501] # cd
CD_co <- K_co[1 , 2:501] # cd
CC_cr <- K_cr[2:501 , 2:501] # cc
CC_co <- K_co[2:501 , 2:501] # cc


for(i in seq_len(n_runs)) {

  b_t_co <- matrix(pop_state$b_t_co[ ,i], nrow = 1, ncol = 1)
  b_t_cr <- matrix(pop_state$b_t_cr[ ,i], nrow = 1, ncol = 1)

  # transitions to continuous state
  n_t_1_cr <- CC_cr %*% pop_state$ht_cr[ , i] + DC_cr %*% b_t_cr
  n_t_1_co <- CC_co %*% pop_state$ht_co[ , i] + DC_co %*% b_t_co

  # Transitions to discrete state
  b_t_1_cr <- CD_cr %*% pop_state$ht_cr[ , i] + DD_cr %*% b_t_cr
  b_t_1_co <- CD_co %*% pop_state$ht_co[ , i] + DD_co %*% b_t_co

  pop_state$ht_cr[ , (i + 1)] <- n_t_1_cr
  pop_state$ht_co[ , (i + 1)] <- n_t_1_co
  pop_state$b_t_cr[ , (i + 1)] <- b_t_1_cr
  pop_state$b_t_co[ , (i + 1)] <- b_t_1_co

  pop_state$lambda_cr[ , i] <- sum(n_t_1_cr, b_t_1_cr) / sum(pop_state$ht_cr[ , i],
                                                             b_t_cr)
  pop_state$lambda_co[ , i] <- sum(n_t_1_co, b_t_1_co) / sum(pop_state$ht_co[ , i],
                                                             b_t_co)
}

pop_size_co <- list(pop_state = pop_state[3:4])
pop_size_cr <- list(pop_state = pop_state[1:2])

lam_s_cr <- pop_state$lambda_cr %>% as.vector()
lam_s_co <- pop_state$lambda_co %>% as.vector()

## ipmr version --------
states <- list(c('ht', 'b'))

inv_logit <- function(int, slope, sv) {
  1/(1 + exp(-(int + slope * sv)))
}

inv_logit_2 <- function(int, slope, slope_2, sv) {
  1/(1 + exp(-(int + slope * sv + slope_2 * sv ^ 2)))
}

ipmr_cr <- init_ipm(sim_gen    = "general",
                    di_dd      = "di",
                    det_stoch  = "det") %>%
  define_kernel(
    name          = "P",
    formula       = s * g * d_ht,
    family        = "CC",
    g             = dnorm(ht_2, g_mu, g_sd),
    g_mu          = g_int + g_slope * ht_1,
    s             = inv_logit_2(s_int, s_slope, s_slope_2, ht_1),
    data_list     = data_list_cr,
    states        = states,
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm',
                                            'g')
  ) %>%
  define_kernel(
    name          = "go_discrete",
    formula       = f_r * f_s * g_i,
    family        = 'CD',
    f_r           = inv_logit(f_r_int, f_r_slope, ht_1),
    f_s           = exp(f_s_int + f_s_slope * ht_1),
    data_list     = data_list_cr,
    states        = states,
    uses_par_sets = FALSE
  ) %>%
  define_kernel(
    name    = 'stay_discrete',
    formula = 0,
    family  = "DD",
    states  = states,
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name          = 'leave_discrete',
    formula       = e_p * f_d * d_ht,
    f_d           = dnorm(ht_2, f_d_mu, f_d_sd),
    family        = 'DC',
    data_list     = data_list_cr,
    states        = states,
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm',
                                            'f_d')
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "go_discrete", "stay_discrete", "leave_discrete"),
      int_rule     = c(rep("midpoint", 4)),
      state_start    = c('ht', "ht", "b", "b"),
      state_end      = c('ht', "b","b", 'ht')
    )
  ) %>%
  define_domains(
    ht = c(1.02, 624, 500)
  ) %>%
  define_pop_state(
    pop_vectors = list(
      n_ht = init_pop_vec,
      n_b  = init_b
    )
  ) %>%
  make_ipm(iterations = 100,
           usr_funs = list(inv_logit   = inv_logit,
                           inv_logit_2 = inv_logit_2),
           return_all_envs = TRUE,
           normalize_pop_size = FALSE)


ipmr_lam_cr <- lambda(ipmr_cr,
                      type_lambda = 'all') %>%
  as.vector()


test_that('ipmr version matches hand version', {

  expect_equal(ipmr_lam_cr, lam_s_cr, tolerance = 1e-10)

  # compare final population distributions
  pop_size_final_cr   <- pop_size_cr$pop_state$ht_cr[ , 101]
  pop_size_final_ipmr <- ipmr_cr$pop_state$n_ht[ , 101]
  bank_size_final_cr   <- pop_size_cr$pop_state$b_t_cr[ , 101]
  bank_size_final_ipmr <- ipmr_cr$pop_state$n_b[ , 101]

  expect_equal(pop_size_final_cr,
               pop_size_final_ipmr,
               tolerance = 1e-7)
  expect_equal(bank_size_final_cr,
               bank_size_final_ipmr,
               tolerance = 1e-7)

})

test_that('classes are correctly set', {

  sub_cls <- vapply(ipmr_cr$sub_kernels,
                    function(x) class(x)[1],
                    character(1L))

  expect_true(all(sub_cls == 'ipmr_matrix'))
  expect_s3_class(ipmr_cr, 'general_di_det_ipm')
  expect_s3_class(ipmr_cr, "ipmr_ipm")
})

test_that("kernel definition order doesn't matter", {

  test_order <- init_ipm(sim_gen    = "general",
                         di_dd      = "di",
                         det_stoch  = "det")  %>%
    define_kernel(
      name          = "go_discrete",
      formula       = f_r * f_s * g_i,
      family        = 'CD',
      f_r           = inv_logit(f_r_int, f_r_slope, ht_1),
      f_s           = exp(f_s_int + f_s_slope * ht_1),
      data_list     = data_list_cr,
      states        = states,
      uses_par_sets = FALSE
    ) %>%
    define_kernel(
      name    = 'stay_discrete',
      formula = 0,
      family  = "DD",
      states  = states,
      evict_cor = FALSE
    ) %>%
    define_kernel(
      name          = 'leave_discrete',
      formula       = e_p * f_d * d_ht,
      f_d           = dnorm(ht_2, f_d_mu, f_d_sd),
      family        = 'DC',
      data_list     = data_list_cr,
      states        = states,
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'f_d')
    ) %>%
    define_kernel(
      name          = "P",
      formula       = s * g * d_ht,
      family        = "CC",
      g             = dnorm(ht_2, g_mu, g_sd),
      g_mu          = g_int + g_slope * ht_1,
      s             = inv_logit_2(s_int, s_slope, s_slope_2, ht_1),
      data_list     = data_list_cr,
      states        = states,
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'g')
    )  %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "go_discrete", "stay_discrete", "leave_discrete"),
        int_rule     = c(rep("midpoint", 4)),
        state_start    = c('ht', "ht", "b", "b"),
        state_end      = c('ht', "b", "b", 'ht')
      )
    ) %>%
    define_domains(
      ht = c(1.02, 624, 500)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_ht = init_pop_vec,
        n_b  = init_b
      )
    ) %>%
    make_ipm(iterations = 100,
             usr_funs = list(inv_logit   = inv_logit,
                             inv_logit_2 = inv_logit_2),
             return_all_envs = TRUE,
             normalize_pop_size = FALSE)

  test_lam_cr <- lambda(test_order,
                        type_lambda = 'all') %>%
    as.vector()

  expect_equal(ipmr_lam_cr, test_lam_cr)

})


test_that('normalize pop vec works the right way', {


  ipmr_cr <- init_ipm(sim_gen    = "general",
                      di_dd      = "di",
                      det_stoch  = "det") %>%
    define_kernel(
      name          = "P",
      formula       = s * g * d_ht,
      family        = "CC",
      g             = dnorm(ht_2, g_mu, g_sd),
      g_mu          = g_int + g_slope * ht_1,
      s             = inv_logit_2(s_int, s_slope, s_slope_2, ht_1),
      data_list     = data_list_cr,
      states        = states,
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'g')
    ) %>%
    define_kernel(
      name          = "go_discrete",
      formula       = f_r * f_s * g_i,
      family        = 'CD',
      f_r           = inv_logit(f_r_int, f_r_slope, ht_1),
      f_s           = exp(f_s_int + f_s_slope * ht_1),
      data_list     = data_list_cr,
      states        = states,
      uses_par_sets = FALSE
    ) %>%
    define_kernel(
      name    = 'stay_discrete',
      formula = 0,
      family  = "DD",
      states  = states,
      evict_cor = FALSE
    ) %>%
    define_kernel(
      name          = 'leave_discrete',
      formula       = e_p * f_d * d_ht,
      f_d           = dnorm(ht_2, f_d_mu, f_d_sd),
      family        = 'DC',
      data_list     = data_list_cr,
      states        = states,
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'f_d')
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "go_discrete", "stay_discrete", "leave_discrete"),
        int_rule     = c(rep("midpoint", 4)),
        state_start    = c('ht', "ht", "b", "b"),
        state_end      = c('ht', "b", "b", 'ht')
      )
    ) %>%
    define_domains(
      ht = c(1.02, 624, 500)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_ht = init_pop_vec,
        n_b  = init_b
      )
    ) %>%
    make_ipm(iterations = 100,
             usr_funs = list(inv_logit   = inv_logit,
                             inv_logit_2 = inv_logit_2),
             return_all_envs = TRUE,
             normalize_pop_size = TRUE)


  norm_lam_ipmr <- lambda(ipmr_cr,
                           type_lambda = 'all') %>%
    as.vector()

  expect_equal(norm_lam_ipmr, ipmr_lam_cr)

  pop_sizes <- vapply(1:101,
                      function(it, pop_state) {

                        b <- pop_state$n_b[ , it]
                        ht <- pop_state$n_ht[ , it]

                        sum(b, ht)
                      },
                      numeric(1L),
                      pop_state = ipmr_cr$pop_state)

  expect_equal(pop_sizes, rep(1, 101), tolerance = 1e-15)

})



s_x <- function(int, int_r, slope1, slope2, sv1) {
  1/(1 + exp(-(int + int_r + slope1 * sv1 + slope2 * sv1^2)))
}

f_r_x <- function(int, slope1,  sv1) {
  1/(1 + exp(-(int + slope1 * sv1)))
}


g_x <- function(sv2, sv1, int, int_r, slope, gsd, L, U) {
  mu <- int + int_r + slope * sv1
  ev <- pnorm(U, mu, gsd) - pnorm(L, mu, gsd)
  out <- dnorm(sv2, mean = mu, sd = gsd) / ev

  return(out)
}

g_ints <- rnorm(3) %>% as.list()
s_ints <- rnorm(3) %>% as.list()

names(g_ints) <- paste("g_int_", 1:3, sep = "")
names(s_ints) <- paste("s_int_", 1:3, sep = "")

data_list_control <- c(data_list_control, g_ints, s_ints)

L <- 1.02
U <- 624
n <- 500

b        <- seq(L, U, length.out = n + 1)
d1 <- d2 <- (b[2:(n + 1)] + b[1:n]) * 0.5
h        <- d1[3] - d1[2]

domains <- expand.grid(list(d2 = d1, d1 = d1))

G_1 <- h * g_x(domains$d1, domains$d2,
             int   = data_list_control$g_int,
             slope = data_list_control$g_slope,
             gsd   = data_list_control$g_sd,
             int_r = data_list_control$g_int_1,
             L = L,
             U = U)

G_2 <- h * g_x(domains$d1, domains$d2,
             int   = data_list_control$g_int,
             slope = data_list_control$g_slope,
             gsd   = data_list_control$g_sd,
             int_r = data_list_control$g_int_2,
             L = L,
             U = U)

G_3 <- h * g_x(domains$d1, domains$d2,
               int   = data_list_control$g_int,
               slope = data_list_control$g_slope,
               gsd   = data_list_control$g_sd,
               int_r = data_list_control$g_int_3,
               L = L,
               U = U)

s_1 <- s_x(data_list_control$s_int,
           data_list_control$s_int_1,
           data_list_control$s_slope,
           data_list_control$s_slope_2,
           d1)

s_2 <- s_x(data_list_control$s_int,
           data_list_control$s_int_2,
           data_list_control$s_slope,
           data_list_control$s_slope_2,
           d1)

s_3 <- s_x(data_list_control$s_int,
           data_list_control$s_int_3,
           data_list_control$s_slope,
           data_list_control$s_slope_2,
           d1)

P_1 <- (s_1 * G_1) %>%
  matrix(nrow = n, ncol = n, byrow = TRUE)
P_2 <- (s_2 * G_2) %>%
  matrix(nrow = n, ncol = n, byrow = TRUE)
P_3 <- (s_3 * G_3) %>%
  matrix(nrow = n, ncol = n, byrow = TRUE)

cd_co <- f_r_x(data_list_control$f_r_int,
               data_list_control$f_r_slope,
               d1) *
  exp(data_list_control$f_s_int + data_list_control$f_s_slope * d1) *
  data_list_control$g_i %>%
  matrix(ncol = n,
         nrow = 1)

dc_co <- (dnorm(d2, mean = data_list_control$f_d_mu, sd = data_list_control$f_d_sd) /
            (pnorm(U,
                   data_list_control$f_d_mu,
                   data_list_control$f_d_sd) - pnorm(L,
                                                     data_list_control$f_d_mu,
                                                     data_list_control$f_d_sd))) *
  data_list_control$e_p *
  h %>%
  matrix(nrow = n, ncol = 1)

dd_co <- 0

K_co_1 <- rbind(
  cbind(dd_co, cd_co),
  cbind(dc_co, P_1)
)
K_co_2 <- rbind(
  cbind(dd_co, cd_co),
  cbind(dc_co, P_2)
)
K_co_3 <- rbind(
  cbind(dd_co, cd_co),
  cbind(dc_co, P_3)
)

eigs_1 <- eigen(K_co_1)
eigs_2 <- eigen(K_co_2)
eigs_3 <- eigen(K_co_3)

lambdas_hand <- c(Re(eigs_1$values[1]),
                  Re(eigs_2$values[1]),
                  Re(eigs_3$values[1]))

ws_hand      <- cbind(
  w_1 = Re(eigs_1$vectors[ , 1]),
  w_2 = Re(eigs_2$vectors[ , 1]),
  w_3 = Re(eigs_3$vectors[ , 1])
)

ws_hand <- apply(ws_hand, 2, function(x) x / sum(x))

# ipmr definition


ipmr_control <- init_ipm(sim_gen    = "general",
                         di_dd      = "di",
                         det_stoch  = "det") %>%
  define_kernel(
    name             = "P_site",
    formula          = s_site * g_site * d_ht,
    family           = "CC",
    g_site           = dnorm(ht_2, g_mu_site, g_sd),
    g_mu_site        = g_int + g_int_site + g_slope * ht_1,
    s_site           = inv_logit_2(s_int + s_int_site, s_slope, s_slope_2, ht_1),
    data_list        = data_list_control,
    states           = states,
    uses_par_sets    = TRUE,
    par_set_indices = list(site = 1:3),
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions('norm',
                                               'g_site')
  ) %>%
  define_kernel(
    name          = "go_discrete",
    formula       = f_r * f_s * g_i,
    family        = 'CD',
    f_r           = inv_logit(f_r_int, f_r_slope, ht_1),
    f_s           = exp(f_s_int + f_s_slope * ht_1),
    data_list     = data_list_control,
    states        = states,
    uses_par_sets = FALSE
  ) %>%
  define_kernel(
    name    = 'stay_discrete',
    formula = 0,
    family  = "DD",
    states  = states,
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name          = 'leave_discrete',
    formula       = e_p * f_d * d_ht,
    f_d           = dnorm(ht_2, f_d_mu, f_d_sd),
    family        = 'DC',
    data_list     = data_list_control,
    states        = states,
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm',
                                            'f_d')
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P_site", "go_discrete", "stay_discrete", "leave_discrete"),
      int_rule     = c(rep("midpoint", 4)),
      state_start    = c('ht', "ht", "b", "b"),
      state_end      = c('ht', "b", "b", 'ht')
    )
  ) %>%
  define_domains(
    ht = c(1.02, 624, 500)
  ) %>%
  define_pop_state(
    pop_vectors = list(
      n_ht_site = init_pop_vec,
      n_b_site  = init_b
    )
  ) %>%
  make_ipm(iterations = 200,
           usr_funs = list(inv_logit   = inv_logit,
                           inv_logit_2 = inv_logit_2),
           return_all_envs = TRUE,
           normalize_pop_size = TRUE)

# lambdas_ipmr <- vapply(ipmr_control$pop_state[grepl("lambda", names(ipmr_control$pop_state))],
#                       function(x) x[ , 200],
#                       numeric(1L))

lambdas_ipmr <- unname(lambda(ipmr_control, type_lambda = "last"))

ws <- list()

ipmr_control$pop_state <- ipmr_control$pop_state[-c(1:2)]

for(i in seq_len(3)) {


  ind <- c(i * 2 - 1, i * 2)

  pop <- do.call("rbind", ipmr_control$pop_state[rev(ind)])

  ws[[i]] <- pop[ , 200]

}

ws <- do.call("cbind", ws) %>%
  setNames(paste("w_", 1:3, sep = ""))

test_that("par_setarchical models get the same answers as hand generated models", {

  expect_equal(lambdas_ipmr, lambdas_hand, tolerance = 1e-10)
  expect_equal(ws[ , 1], ws_hand[ , 1], tolerance = 1e-10)
  expect_equal(ws[ , 2], ws_hand[ , 2], tolerance = 1e-10)
  expect_equal(ws[ , 3], ws_hand[ , 3], tolerance = 1e-10)

})

test_that("make_iter_kernel can handle arithmetic in expressions", {

  f <- matrix(runif(500 * 500), 500, 500)

  k_1 <- rbind(
    cbind(dd_co, cd_co),
    cbind(dc_co, P_1 + f)
  )

  k_2 <- rbind(
    cbind(dd_co, cd_co),
    cbind(dc_co, P_2 + f)
  )
  k_3 <- rbind(
    cbind(dd_co, cd_co),
    cbind(dc_co, P_3 + f)
  )

  ipmr_control$sub_kernels <- c(ipmr_control$sub_kernels, list(f = f))

  test_ks <- make_iter_kernel(ipmr_control,
                              mega_mat = c(stay_discrete,  go_discrete,
                                           leave_discrete, P_site + f))


  expect_equal(k_1, test_ks$mega_matrix_1, ignore_attr = c("dimnames", "class"))
  expect_equal(k_2, test_ks$mega_matrix_2, ignore_attr = c("dimnames", "class"))
  expect_equal(k_3, test_ks$mega_matrix_3, ignore_attr = c("dimnames", "class"))

  expect_s3_class(test_ks$mega_matrix_1, "ipmr_matrix")
  expect_s3_class(test_ks$mega_matrix_2, "ipmr_matrix")
  expect_s3_class(test_ks$mega_matrix_3, "ipmr_matrix")

})



test_that("DC/CD transitions without size-dependence work", {

  data_list <- list(
    g_int     = 0.6801982,
    g_slope   = 0.8186528,
    g_sd      = 0.1691976,
    s_notH    = 0.99,
    m_notH    = 0,
    d_mu    = 6.33,
    d_sd    = 2.06,
    a       = 0.92,
    b       = 0.057,
    d       = 1
  )

  general_ipm <- init_ipm(sim_gen    = "general",
                          di_dd      = "di",
                          det_stoch  = "det") %>%
    define_kernel(
      name          = "P",
      formula       = s * g * (1-m) * d_size,
      family        = "CC",
      g             = dnorm(size_2, g_mu, g_sd),
      g_mu          = g_int + g_slope * size_1,
      s             = s_notH,
      m             = m_notH,
      data_list     = data_list,
      states        = list(c('size')),
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'g')
    ) %>%

    define_kernel(
      name          = "go_sprout",
      formula       = s * m * d_size,
      family        = 'CD',
      s             = s_notH,
      m             = m_notH,
      data_list     = data_list,
      states        = list(c('size', "sprout")),
      uses_par_sets = FALSE
    ) %>%


    define_kernel(
      name    = 'stay_sprout',
      formula = a * (1-b),
      family  = "DD",
      data_list     = data_list, states  = list(c('size', "sprout")),
      evict_cor = FALSE
    ) %>%

    define_kernel(
      name          = 'leave_sprout',
      formula       = a * b * f_d * d_size,
      f_d           = dnorm(size_2, d_mu, d_sd),
      family        = 'DC',
      data_list     = data_list,
      states        = list(c('size', "sprout")),
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm',
                                              'f_d')
    )

  general_ipm <- general_ipm %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "go_sprout", "stay_sprout", "leave_sprout"),
        int_rule     = c(rep("midpoint", 4)),
        state_start    = c('size', "size", "sprout", "sprout"),
        state_end      = c('size',"sprout", "sprout", 'size')
      )
    )

  # The lower and upper bounds for the continuous state variable and the number
  # of meshpoints for the midpoint rule integration. We'll also create the initial
  # population vector from a random uniform distribution
  L <- 1.02
  U <- 4.3
  n <- 500

  init_pop_vec  <- rpois(500, 2)
  init_sprout   <- 30

  ipm <- general_ipm %>%
    define_domains(
      size = c(L, U, n)
    ) %>%
    define_pop_state(
      pop_vectors = list(
        n_size = init_pop_vec,
        n_sprout  = init_sprout
      )
    ) %>%
    make_ipm(iterations = 20,
             usr_funs = list(inv_logit   = inv_logit),
             normalize_pop_size = FALSE)


  l_pop <- lambda(ipm)

  # Lambda must always be less than 1 here

  expect_true(l_pop < 1)

  imputed_kern <- ipm$sub_kernels$go_sprout

  expect_true(all(imputed_kern == 0))

})
