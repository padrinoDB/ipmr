
tab_legend <- "Math Formula corresponds to the mathematical notation for the model. R Formula corresponds to how one could write the model using R's model formula syntax (think `lm`, `glm`, `lmer`). ipmr shows an equivalent way to write this model in kernel formula or vital rate expression. Model Number indicates which rows are grouped together, as some models have too many terms to cleanly fit on a line, or have multiple equivalent representations."

styles <- c("style='width:150%;'")

knitr::kable(
  data.frame(

    Math    = c("$\\mu_G^{yr} = \\alpha_G + \\alpha_G^{yr} * \\beta_G * z$",
                         "$G(z',z) = f_G(z', \\mu_g^{yr}, \\sigma_G)$",
                "$Logit(\\mu_{s}^{plot,yr}) = \\alpha_s + \\alpha_{s}^{plot,yr} + \\beta * z$",
                "",
                "",
                "$log(\\mu_{r}^{yr}) = \\alpha_{r} + \\alpha_r^{yr} + (\\beta_r + \\beta_r^{yr}) * z$"),


    R = c("`size_2 ~ size_1 + (1 | year), family = gaussian()`",
          "",
          "`surv ~ size_1 + (1 | year) + (1 | plot), family = binomial()`",
          "`surv ~ size_1 + (1 | year) + (1 | plot), family = binomial()`",
          "",
          "`fec ~ size_1 + (size_1 | year), family = poisson()`"),


    ipmr = c("mu_g_**yr** = g_int + g_int_**yr** + g_slope * z",
             "g = dnorm(z_2, mu_g_**yr**, sd_g)",
             "s_**pl**_**yr** = s_int + s_int_**pl**_**yr** + s_slope * z",
             "s_**pl**_**yr** = s_int + s_int_**pl** + s_int_**yr** + s_slope * z",
             "s = inv_logit(s_**pl**_**yr**)",
             "r = exp(r_int + r_int_**yr** + (r_slope + r_slope_**yr**) * z)"),
    Model = c(1, 1, 2, 2, 2, 3)
  ),
  format    = "html",
  escape    = FALSE,
  col.names = c("Math Formula", "R Formula", "ipmr", "Model Number"),
  caption   = tab_legend,
  table.attr = styles
)

library(ipmr)

all_fixed_params <- list(
  alpha_g   = 0.2,
  beta_g    = 1.01,
  sigma_g   = 0.2,
  alpha_s   = -0.1,
  beta_s    = 0.3,
  alpha_r_p = -1,
  beta_r_p  = 0.7,
  alpha_r_r = 0.2,
  beta_r_r  = 0.2,
  mu_r_d    = 0.4,
  sigma_r_d = 0.1
)

g_alpha_list <- as.list(rnorm(6)) %>%
  setNames(
    paste('alpha_g_', 2001:2006, sep = "")
  )
s_alpha_list <- as.list(rnorm(6)) %>%
  setNames(
    paste('alpha_s_', 2001:2006, sep = "")
  )

all_params <- c(all_fixed_params, g_alpha_list, s_alpha_list)



 ex_ipm <-init_ipm(sim_gen    = "simple",
                   di_dd      = "di",
                   det_stoch  = "stoch",
                   kern_param = "kern") %>%
   define_kernel(

     # The _yr index is appended as a suffix with an underscore. Notice how the formula
     # from Eq 2 is translated into the formula argument and the kernel name

     name             = "P_yr",
     family           = "CC",
     formula          = s_yr * g_yr,

     # We also modify each parameter name is indexed as well.
     # Here, I've split out the linear predictor and the inverse logit
     # transformation into separate steps to avoid over cluttering

     s_lin_p          = alpha_s + alpha_s_yr + beta_s * z_1,
     s_yr             = 1 / (1 + exp(- s_lin_p)),

     # We do the same with growth as we did with survival and the P_yr formula

     g_yr             = dnorm(z_2, mu_g_yr, sigma_g),
     mu_g_yr          = alpha_g + alpha_g_yr + beta_g * z_1,

     data_list        = all_params,
     states           = list(c("z")),

     # We signal that the kernel has parameter sets

     uses_par_sets    = TRUE,

     # And provide the values that each index can take as a named list.
     # The name(s) in this list MUST match the index used in the expressions.

     par_set_indices = list(yr = 2001:2006),
     evict_cor        = TRUE,
     evict_fun        = truncated_distributions("norm", "g_yr")
   ) %>%

   # This kernel doesn't get anything special because there are no
   # time-varying parameter values.

   define_kernel(
     name          = "F",
     family        = "CC",
     formula       = r_p * r_r * r_d,
     r_p_lin_p     = alpha_r_p + beta_r_p * z_1,
     r_p           = 1 / (1 + exp( -r_p_lin_p)),
     r_r           = exp(alpha_r_r + beta_r_r * z_1),
     r_d           = dnorm(z_2, mu_r_d, sigma_r_d),
     data_list     = all_params,
     states        = list(c("z")),
     uses_par_sets = FALSE,
     evict_cor     = TRUE,
     evict_fun     = truncated_distributions("norm", "r_d")
   ) %>%
   define_impl(
     make_impl_args_list(
       kernel_names = c("P_yr", "F"),
       int_rule     = rep("midpoint", 2),
       state_start    = rep("z", 2),
       state_end      = rep("z", 2)
     )
   ) %>%
   define_domains(
     z = c(0.1, 10, 50)
   ) %>%
   define_pop_state(
     n_z = runif(50)
   ) %>%
   make_ipm(
     iterate    = TRUE,
     iterations = 100,
     kernel_seq = sample(2001:2006, size = 100, replace = TRUE),
     normalize_pop_size = TRUE,
     return_sub_kernels = TRUE
   )



 library(ipmr)

 set.seed(51)

 # Set up the sites, and the parameters

 sites <- LETTERS[1:6]

 all_pars <- list(
   r_i       = 0.6,
   r_sb1     = 0.2,
   r_sb2     = 0.3,
   s_sb1     = 0.4,
   s_sb2     = 0.1,
   mu_c_d    = 4,
   sigma_c_d = 0.9,
   sigma_g   = 3,
   beta_s_z    = 2.5,
   beta_s_prec = 0.3,
   beta_s_temp = -0.5,
   beta_g_z    = 0.99,
   beta_g_prec = 0.01,
   beta_g_temp = 0.1,
   beta_c_p_z    = 0.4,
   beta_c_p_prec = 0.6,
   beta_c_p_temp = -0.3,
   beta_c_r_z    = 0.005,
   beta_c_r_prec = -0.3,
   beta_c_r_temp = 0.9
 )

 g_alphs <- rnorm(6, mean = 4.5, sd = 1) %>%
   as.list() %>%
   setNames(
     paste('alpha_g_', sites, sep = "")
   )

 s_alphs <- rnorm(6, mean = -7, sd = 1.5) %>%
   as.list() %>%
   setNames(
     paste('alpha_s_', sites, sep = "")
   )

 c_p_alphs <- rnorm(6, mean = -10, sd = 2)  %>%
   as.list() %>%
   setNames(
     paste('alpha_c_p_', sites, sep = "")
   )

 c_r_alphs <- runif(6, 0.01, 0.1) %>%
   as.list() %>%
   setNames(
     paste('alpha_c_r_', sites, sep = "")
   )


 all_pars <- c(all_pars,
               g_alphs,
               s_alphs,
               c_p_alphs,
               c_r_alphs)

 # This list contains parameters that define the distributions for environmental
 # covariates

 env_params <- list(
   temp_mu = 12,
   temp_sd = 2,
   prec_shape = 1000,
   prec_rate  = 2
 )

 # This function samples the environmental covariate distributions and returns
 # the values in a named list. We can reference the names in the list in our
 # vital rate expressions.

 sample_fun <- function(env_params) {

   temp <- rnorm(1, env_params$temp_mu, env_params$temp_sd)

   prec <- rgamma(1, env_params$prec_shape, env_params$prec_rate)

   out  <- list(temp = temp,
                prec = prec)

   return(out)

 }



 ex_ipm <- init_ipm(sim_gen    = "general",
                    di_dd      = "di",
                    det_stoch  = "stoch",
                    kern_param = "param") %>%
   define_kernel(

     name = "P_site",
     family = "CC",

     # site is appended so that we don't have to write out s_A, s_B, etc.

     formula          = s_site * g_site * d_z,

     s_site_lin       = alpha_s_site +
                        beta_s_z * z_1 +
                        beta_s_prec * prec +
                        beta_s_temp * temp,
     s_site           = 1 / (1 + exp(-s_site_lin)),
     mu_g_site        = alpha_g_site +
                        beta_g_z * z_1 +
                        beta_g_prec * prec +
                        beta_g_temp * temp,
     g_site           = dnorm(z_2, mu_g_site, sigma_g),
     data_list        = all_pars,
     states           = list(c("z")),
     uses_par_sets    = TRUE,
     par_set_indices = list(site = LETTERS[1:6]),
     evict_cor        = TRUE,
     evict_fun        = truncated_distributions("norm", "g_site")
   ) %>%
   define_kernel(
     name             = "F_site",
     family           = "CC",
     formula          = c_p_site * c_r_site * c_d * d_z,
     c_p_site_lin     = alpha_c_p_site +
                        beta_c_p_z * z_1 +
                        beta_c_p_prec * prec +
                        beta_c_p_temp * temp,
     c_p_site         = 1 / (1 + exp(-c_p_site_lin)),
     c_r_site         = exp(alpha_c_r_site +
                            beta_c_r_z * z_1 +
                            beta_c_r_prec * prec +
                            beta_c_r_temp * temp),
     c_d              = dnorm(z_2, mu_c_d, sigma_c_d),
     data_list        = all_pars,
     states           = list(c("z", "sb1", "sb2")),
     uses_par_sets    = TRUE,
     par_set_indices = list(site = LETTERS[1:6]),
     evict_cor        = TRUE,
     evict_fun        = truncated_distributions("norm", "c_d")
     ) %>%
   define_kernel(
     name             = "go_sb_1_site",
     family           = "CD",
     formula          = (1 - r_i) * c_p_site * c_r_site * d_z,
     c_p_site_lin     = alpha_c_p_site +
                        beta_c_p_z * z_1 +
                        beta_c_p_prec * prec +
                        beta_c_p_temp * temp,
     c_p_site         = 1 / (1 + exp(-c_p_site_lin)),
     c_r_site         = exp(alpha_c_r_site +
                            beta_c_r_z * z_1 +
                            beta_c_r_prec * prec +
                            beta_c_r_temp * temp),
     data_list        = all_pars,
     states           = list(c("z", "sb1")),
     uses_par_sets    = TRUE,
     par_set_indices = list(site = LETTERS[1:6]),
     evict_cor        = FALSE
   ) %>%
   define_kernel(
     name          = "go_sb_2",
     family        = "DD",
     formula       = s_sb1 * (1 - r_sb1),
     data_list     = all_pars,
     states        = list(c("sb1", "sb2")),
     uses_par_sets = FALSE,
     evict_cor     = FALSE
   ) %>%
   define_kernel(
     name          = "leave_sb_1",
     family        = "DC",
     formula       = s_sb1 * r_sb1 * c_d * d_z,
     c_d           = dnorm(z_2, mu_c_d, sigma_c_d),
     data_list     = all_pars,
     states        = list(c("z", "sb1")),
     uses_par_sets = FALSE,
     evict_cor     = TRUE,
     evict_fun     = truncated_distributions("norm", "c_d")
   ) %>%
   define_kernel(
     name          = "leave_sb_2",
     family        = "DC",
     formula       = s_sb2 * r_sb2 * c_d * d_z,
     c_d           = dnorm(z_2, mu_c_d, sigma_c_d),
     data_list     = all_pars,
     states        = list(c("z", "sb2")),
     uses_par_sets = FALSE,
     evict_cor     = TRUE,
     evict_fun     = truncated_distributions("norm", "c_d")
   ) %>%
   define_impl(
     make_impl_args_list(
       kernel_names = c("P_site",
                        "F_site",
                        "go_sb_1_site",
                        "go_sb_2",
                        "leave_sb_1",
                        "leave_sb_2"),
       int_rule     = rep("midpoint", 6),
       state_start    = c("z", "z", "z", "sb1", "sb1", "sb2"),
       state_end      = c("z", "z", "sb1", "sb2", "z", "z")
     )
   ) %>%
   define_domains(
     z = c(0.1, 100, 200)
   ) %>%
   define_pop_state(
     n_z   = runif(200),
     n_sb1 = 20,
     n_sb2 = 20
   ) %>%
   define_env_state(
     env_state = sample_fun(env_params),
     data_list = list(env_params = env_params,
                      sample_fun = sample_fun)
   ) %>%
   make_ipm(
     iterations         = 20,
     kernel_seq         = sample(LETTERS[1:6], 20, replace = TRUE),
     normalize_pop_size = TRUE,
     return_sub_kernels = TRUE
   )

 lambda(ex_ipm)
