# Testing out ipmr reproducibility with Lonicera maackii from Levin, Crandall & Knight 2019

library(ipmr)

data_list <- list(
  s_int = -2.830987,
  s_slope = 0.403509,
  s_slope_2 = -0.000421,
  g_int = 16.884,
  g_slope = 0.9972,
  sd_g = 32.74774,
  f_r_int = -10.4478,
  f_r_slope = 0.0485,
  f_s_int = 3.391,
  f_s_slope = 0.0105,
  f_d_mu = 3.118,
  f_d_sd = 1.215,
  e_p = 0.003563
)

s_x <- function(int, s1, s2, sv) {
  1/(1 + exp(-(int + s1 * sv + s2 * sv^2)))
}

inv_logit <- function(int, s, sv) {
  1/(1 + exp(-(int + s * sv)))
}

l_maackii <- init_ipm('simple_di_det') %>%
  define_kernel(
    'P',
    formula = s_g_mult(s, g),
    family = 'CC',
    g = dnorm(ht_2, mean = mu_g, sd = sd_g),
    mu_g = g_int + g_slope * ht_1,
    s = s_x(s_int, s_slope, s_slope_2, ht_1),
    data_list = data_list,
    states = list(c('ht')),
    has_hier_effs = FALSE,
    evict = TRUE,
    evict_fun = truncated_distributions(g,
                                        500)
  ) %>%
  define_kernel(
    "F",
    family = "CC",
    formula = e_p * f_r * f_s * f_d,
    f_r = inv_logit(f_r_int, f_r_slope, ht_1),
    f_s = exp(f_s_int + f_s_slope * ht_1),
    f_d = dnorm(ht_2, f_d_mu, f_d_sd),
    data_list = data_list,
    states = list(c('ht')),
    has_hier_effs = FALSE,
    evict = FALSE
  ) %>%
  define_k(
    'K',
    family = "IPM",
    K = P + F,
    data_list = data_list,
    states = list(c('ht')),
    evict = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      c('K', 'P', "F"),
      rep("midpoint", 3),
      rep('ht', 3),
      rep('ht', 3)
    )
  ) %>%
  define_domains(
    ht = c(1.26, 482.9, 500)
  ) %>%
  make_ipm(usr_funs = list(s_x = s_x,
                           inv_logit = inv_logit))


cr_data_list <- list(
  s_int = -2.784776,
  s_slope = 0.272694,
  s_slope_2 = -0.000175,
  g_int = 14.6068,
  g_slope = 0.9964,
  sd_g = 28.01974,
  f_r_int = -10.4478,
  f_r_slope = 0.0485,
  f_s_int = 3.391,
  f_s_slope = 0.0105,
  f_d_mu = 3.118,
  f_d_sd = 1.215,
  e_p = 0.003563
)

l_maackii_cr <- init_ipm('simple_di_det') %>%
  define_kernel(
    'P',
    formula = s_g_mult(s, g),
    family = 'CC',
    g = dnorm(ht_2, mean = mu_g, sd = sd_g),
    mu_g = g_int + g_slope * ht_1,
    s = s_x(s_int, s_slope, s_slope_2, ht_1),
    data_list = cr_data_list,
    states = list(c('ht')),
    has_hier_effs = FALSE,
    evict = TRUE,
    evict_fun = truncated_distributions(g,
                                        500)
  ) %>%
  define_kernel(
    "F",
    family = "CC",
    formula = e_p * f_r * f_s * f_d,
    f_r = inv_logit(f_r_int, f_r_slope, ht_1),
    f_s = exp(f_s_int + f_s_slope * ht_1),
    f_d = dnorm(ht_2, f_d_mu, f_d_sd),
    data_list = cr_data_list,
    states = list(c('ht')),
    has_hier_effs = FALSE,
    evict = FALSE
  ) %>%
  define_k(
    'K',
    family = "IPM",
    K = P + F,
    data_list = cr_data_list,
    states = list(c('ht')),
    evict = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      c('K', 'P', "F"),
      rep("midpoint", 3),
      rep('ht', 3),
      rep('ht', 3)
    )
  ) %>%
  define_domains(
    ht = c(1.26, 482.9, 500)
  ) %>%
  make_ipm(usr_funs = list(s_x = s_x,
                           inv_logit = inv_logit))



# Should be ~1.2848
target <- 1.2848
ipmr_cont <- Re(eigen(l_maackii$iterators$K)$values[1])


# Should be ~1.2305
target_cr <- 1.2305
ipmr_cr <- Re(eigen(l_maackii_cr$iterators$K)$values[1])

# Control
(ipmr_cont - target)/target
#competitor removal
(ipmr_cr - target_cr)/target_cr
