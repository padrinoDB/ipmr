context('general density independent deterministic models')


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


g_x <- function(sv2, sv1, int, slope, gsd) {
  mu <- int + slope * sv1
  out <- dnorm(sv2, mean = mu, sd = gsd)

  return(out)
}


L <- 1.02
U <- 624
n <- 500

b        <- seq(L, U, length.out = n + 1)
d1 <- d2 <- (b[2:(n + 1)] + b[1:n]) * 0.5
h        <- d1[3] - d1[2]

G <- h * outer(d1, d2, FUN = g_x,
               int   = data_list_control$g_int,
               slope = data_list_control$g_slope,
               gsd   = data_list_control$g_sd)

s <- s_x(data_list_control$s_int,
         data_list_control$s_slope,
         data_list_control$s_slope_2,
         d1)

G <- truncated_distributions(G,
                             n)

P <- s_g_mult(s, G)

cd_co <- f_r_x(data_list_control$f_r_int,
            data_list_control$f_r_slope,
            d1) *
  exp(data_list_control$f_s_int + data_list_control$f_s_slope * d1) *
  data_list_control$g_i %>%
  matrix(ncol = n,
         nrow = 1)

dc_co <- dnorm(d2, mean = data_list_control$f_d_mu, sd = data_list_control$f_d_sd) *
  data_list_control$e_p *
  h %>%
  matrix(nrow = n, ncol = 1)

cc_co <- P
dd_co <- 0

K_co <- rbind(
  cbind(dd_co, cd_co),
  cbind(dc_co, cc_co)
)

actual_co <- Re(eigen(K_co)$values[1])

# Repeat, but with CRs

G <- h * outer(d1, d2, FUN = g_x,
               int   = data_list_cr$g_int,
               slope = data_list_cr$g_slope,
               gsd   = data_list_cr$g_sd)

s <- s_x(data_list_cr$s_int,
         data_list_cr$s_slope,
         data_list_cr$s_slope_2,
         d1)

G <- truncated_distributions(G,
                             n)

P <- s_g_mult(s, G)

cd_cr <- f_r_x(data_list_cr$f_r_int,
               data_list_cr$f_r_slope,
               d1) *
  exp(data_list_cr$f_s_int + data_list_cr$f_s_slope * d1) *
  data_list_cr$g_i %>%
  matrix(ncol = n,
         nrow = 1)

dc_cr <- dnorm(d2, mean = data_list_cr$f_d_mu, sd = data_list_cr$f_d_sd) *
  data_list_cr$e_p *
  h %>%
  matrix(nrow = n, ncol = 1)

cc_cr <- P
dd_cr <- 0

K_cr <- rbind(
  cbind(dd_cr, cd_cr),
  cbind(dc_cr, cc_cr)
)

actual_cr <- Re(eigen(K_cr)$values[1])

target_cr <- 1.2578
target_co <- 1.2208


# Verify test implementation actually does the published one

test_that('test implementation matches published quantities', {
  expect_equal(target_cr, actual_cr, tolerance = 1e-3)
  expect_equal(target_co, actual_co, tolerance = 1e-3)
})

# Now, set up simulation to compare w/ ipmr

n_runs <- 100
pop_state <- list(ht_cr  = matrix(NA_real_, nrow = n, ncol = n_runs + 1),
                  b_t_cr = matrix(NA_real_, nrow = 1, ncol = n_runs + 1),
                  ht_co  = matrix(NA_real_, nrow = n, ncol = n_runs + 1),
                  b_t_co = matrix(NA_real_, nrow = 1, ncol = n_runs + 1))

pop_state$ht_cr[ , 1]  <- pop_state$ht_co[ , 1]  <- init_pop_vec
pop_state$b_t_cr[ , 1] <- pop_state$b_t_co[ , 1] <- init_b

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
}

pop_size_co <- list(pop_state = pop_state[3:4])
pop_size_cr <- list(pop_state = pop_state[1:2])

lam_s_cr <- ipmr:::.stoch_lambda_pop_size(pop_size_cr)
lam_s_co <- ipmr:::.stoch_lambda_pop_size(pop_size_co)

## ipmr version --------
states <- list(c('ht', 'sb'))

inv_logit <- function(int, slope, sv) {
  1/(1 + exp(-(int + slope * sv)))
}

inv_logit_2 <- function(int, slope, slope_2, sv) {
  1/(1 + exp(-(int + slope * sv + slope_2 * sv^2)))
}

ipmr_cr <- init_ipm("general_di_det") %>%
  define_kernel(
    name          = "P",
    formula       = s_g_mult(s, g),
    family        = "CC",
    g             = dnorm(ht_2, g_mu, g_sd),
    g_mu          = g_int + g_slope * ht_1,
    s             = inv_logit_2(s_int, s_slope, s_slope_2, ht_1),
    states        = states,
    has_hier_effs = FALSE,
    evict         = TRUE,
    evict_fun     = truncated_distributions(g,
                                            500)
  ) %>%
  define_kernel(
    name          = "go_discrete",
    formula       = f_r * f_s * g_i,
    family        = 'CD',
    f_r           = inv_logit(f_r_int, f_r_slope, ht_1),
    f_s           = exp(f_s_int + f_s_slope * ht_1),
    data_list     = data_list_cr,
    states        = states,
    has_hier_effs = FALSE,
    evict         = FALSE
  ) %>%
  define_kernel(
    name    = 'stay_discrete',
    formula = 0,
    family  = "DD",
    states  = states,
    evict = FALSE
  ) %>%
  define_kernel(
    name          = 'leave_discrete',
    formula       = e_p * dnorm(ht_2, f_d_mu, f_d_sd),
    family        = 'DC',
    data_list     = data_list_cr,
    states        = states,
    has_hier_effs = FALSE,
    evict         = FALSE,
  ) %>%
  define_k(
    name     = "K",
    family   = "IPM",

    n_b_t_1  = right_mult(stay_discrete, n_b_t)  + right_mult(go_discrete,
                                                              n_ht_t),
    n_ht_t_1 = right_mult(leave_discrete, n_b_t) + right_mult(P, n_ht_t),


    data_list     = data_list_cr,
    states        = states,
    has_hier_effs = FALSE,
    evict         = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "go_discrete", "stay_discrete", "leave_discrete", "K"),
      int_rule     = c(rep("midpoint", 5)),
      dom_start    = c('ht', "ht", NA_character_, NA_character_, "ht"),
      dom_end      = c('ht', NA_character_, NA_character_, 'ht', 'ht')
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
  ) #%>%
  # make_ipm(iterations = 100,
           # usr_funs = list(inv_logit   = inv_logit,
           #                 inv_logit_2 = inv_logit_2))

