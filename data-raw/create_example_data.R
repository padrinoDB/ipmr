# Create a simple deterministic IPM

devtools::load_all()

data_list = list(s_int = 2.2,
                 s_slope = 0.25,
                 g_int = 0.2,
                 g_slope = 1.02,
                 sd_g = 0.7,
                 f_r_int = 0.003,
                 f_r_slope = 0.015,
                 f_s_int = 1.3,
                 f_s_slope = 0.075,
                 mu_fd = 2,
                 sd_fd = 0.3)


inv_logit <- function(int, slope, sv) {
  return(1/(1 + exp(-(int + slope * sv))))
}

impl_args <- make_impl_args_list(c('P', 'F'),
                                 int_rule = rep('midpoint', 2),
                                 state_start = rep('dbh', 2),
                                 state_end = rep('dbh', 2))

states <- c('dbh', 'dbh')

sim_di_det_ex <- init_ipm('simple', 'di', 'det') %>%
  define_kernel("P",
                formula = s * g,
                family = "CC",
                s = inv_logit(s_int, s_slope, dbh_1),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * dbh_1,
                data_list = data_list,
                states = states,
                evict_cor = TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'g')
  ) %>%
  define_kernel('F',
                formula = f_r * f_s * f_d,
                family = 'CC',
                f_r = inv_logit(f_r_int, f_r_slope, dbh_1),
                f_s = exp(f_s_int + f_s_slope * dbh_1),
                f_d = dnorm(dbh_2, mu_fd, sd_fd),
                data_list = data_list,
                states = states,
                evict_cor = TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'f_d')
  ) %>%
  define_impl(impl_args) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  define_pop_state(n_dbh = runif(100)) %>%
  make_ipm(usr_funs = list(inv_logit = inv_logit),
           iterate = TRUE,
           iterations = 100,
           return_all_envs = FALSE,
           return_main_env = FALSE)


# Create a general deterministic ipm


init_pop_vec   <- runif(200)
init_seed_bank <- 20

data_list <- list(
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

# The lower and upper bounds for the continuous state variable and the number
# of meshpoints for the midpoint rule integration.

L <- 1.02
U <- 624
n <- 200

# Initialize the state list and add some helper functions. The survival function
# in this model is a quadratic function.

states <- list(c('ht', 'b'))

inv_logit <- function(int, slope, sv) {
  1/(1 + exp(-(int + slope * sv)))
}

inv_logit_2 <- function(int, slope, slope_2, sv) {
  1/(1 + exp(-(int + slope * sv + slope_2 * sv ^ 2)))
}

gen_di_det_ex <- init_ipm("general", "di", "det") %>%
  define_kernel(
    name          = "P",
    formula       = s * g * d_ht,
    family        = "CC",
    g             = dnorm(ht_2, g_mu, g_sd),
    g_mu          = g_int + g_slope * ht_1,
    s             = inv_logit_2(s_int, s_slope, s_slope_2, ht_1),
    data_list     = data_list,
    states        = states,
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm',
                                            'g')
  ) %>%
  define_kernel(
    name          = "go_discrete",
    formula       = f_r * f_s * g_i * d_ht,
    family        = 'CD',
    f_r           = inv_logit(f_r_int, f_r_slope, ht_1),
    f_s           = exp(f_s_int + f_s_slope * ht_1),
    data_list     = data_list,
    states        = states,
    uses_par_sets = FALSE
  ) %>%
  define_kernel(
    name      = 'stay_discrete',
    formula   = 0,
    family    = "DD",
    states    = states,
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name          = 'leave_discrete',
    formula       = e_p * f_d * d_ht,
    f_d           = dnorm(ht_2, f_d_mu, f_d_sd),
    family        = 'DC',
    data_list     = data_list,
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
    ht = c(L, U, n)
  ) %>%
  define_pop_state(
    pop_vectors = list(
      n_ht = init_pop_vec,
      n_b  = init_seed_bank
    )
  ) %>%
  make_ipm(iterations = 100,
           return_all_envs = FALSE,
           return_main_env = FALSE,
           usr_funs = list(inv_logit   = inv_logit,
                           inv_logit_2 = inv_logit_2))

data_list <- list(
  s_i  = -0.65, # Intercept of the survival model (Logistic regression)
  s_z  = 0.75,  # slope of the survival model
  f_i  = -18,   # Intercept of the probability of flowering model (Logistic regression)
  f_z  = 6.9,   # Slope of the probability of flowering model
  G_i  = 0.96,  # Intercept of the growth model (Gaussian regression)
  G_z  = 0.59,  # Slope of the growth model
  sd_G = 0.67,  # Standard deviation of residuals of growth model
  mu_r = -0.08, # mean of the recruit size distribution
  sd_r = 0.76,  # Standard deviation of the recruit size distribution
  p_i  = 1,     # Intercept of the seed production model (Poisson regression)
  p_z  = 2.2,   # Slope of the seed production model
  e_p  = 0.007  # Establishment probability
)

monocarp_proto <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det") %>%
  define_kernel(
    name      = "P",
    formula   = (1 - f_p_z) * s * G,
    f_p_z     = plogis(f_i + f_z  * z_1),
    s         = plogis(s_i + s_z * z_1),
    G         = dnorm(z_2, mu_G, sd_G),
    mu_G      = G_i + G_z * z_1,
    data_list = data_list,
    states    = list(c("z"))
  ) %>%
  define_kernel(
    name      = "F",
    formula   = f_p_z * p_n_z * r_d * e_p,
    f_p_z     = plogis(f_i + f_z  * z_1),
    p_n_z     = exp(p_i + p_z * z_1),
    r_d       = dnorm(z_2, mu_r, sd_r),
    data_list = data_list,
    states    = list(c("z"))
  ) %>%
  define_impl(
    list(
      P = list(int_rule = "midpoint", state_start = "z", state_end = "z"),
      F = list(int_rule = "midpoint", state_start = "z", state_end = "z")
    )
  ) %>%
  define_domains(
    z = c(-2.65, 4.5, 250) # Vector with c(L, U, m), where m is number of meshpoints
  ) %>%
  define_pop_state(
    n_z = rep(1/250, 250)
  )

usethis::use_data(monocarp_proto, overwrite = TRUE)
usethis::use_data(sim_di_det_ex, overwrite = TRUE)
usethis::use_data(gen_di_det_ex, overwrite = TRUE)

iceplant_ex <- read.csv("data-raw/iceplant_ex.csv",
                        stringsAsFactors = FALSE)

usethis::use_data(iceplant_ex, overwrite = TRUE)

