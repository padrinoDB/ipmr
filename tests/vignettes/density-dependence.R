
 library(ipmr)

 data_list = list(
   s_int     = 1.03,
   s_slope   = 2.2,
   s_dd      = -0.7,
   g_int     = 8,
   g_slope   = 0.92,
   sd_g      = 0.9,
   r_r_int   = 0.09,
   r_r_slope = 0.05,
   r_s_int   = 0.1,
   r_s_slope = 0.005,
   r_s_dd    = -0.03,
   mu_rd     = 9,
   sd_rd     = 2
 )

 # Now, simulate some random intercepts for growth, survival, and offspring production

 g_r_int   <- rnorm(5, 0, 0.3)
 s_r_int   <- rnorm(5, 0, 0.7)
 r_s_r_int <- rnorm(5, 0, 0.2)

 nms <- paste("r_", 1:5, sep = "")

 names(g_r_int) <- paste("g_", nms, sep = "")
 names(s_r_int) <- paste("s_", nms, sep = "")
 names(r_s_r_int) <- paste("r_s_", nms, sep = "")

 params     <- c(data_list, g_r_int, s_r_int, r_s_r_int)



 dd_ipm <- init_ipm(sim_gen = "simple",
                    di_dd = "dd",
                    det_stoch = "stoch",
                    kern_param = "kern")

 dd_ipm <- define_kernel(
   proto_ipm        = dd_ipm,
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
 )



 dd_ipm <- define_kernel(
   proto_ipm        = dd_ipm,
   name             = "F_yr",
   formula          = r_r * r_s_yr * r_d,
   family           = "CC",
   r_r              = plogis(r_r_int + r_r_slope * size_1),
   r_s_yr           = exp(r_s_int + r_s_r_yr + r_s_slope * size_1 + r_s_dd * sum(n_size_t)),
   r_d              = dnorm(size_2, mu_rd, sd_rd),
   data_list        = params,
   states           = list(c("size")),
   uses_par_sets    = TRUE,
   par_set_indices = list(yr = 1:5),
   evict_cor        = TRUE,
   evict_fun        = truncated_distributions("norm", "r_d")
   )


  dd_ipm <-  dd_ipm %>%
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



 time_step_lams <- lambda(dd_ipm, type_lambda = "all")
 stoch_lam      <- lambda(dd_ipm, type_lambda = "stochastic", burn_in = 0.15)

 pop_sizes <- colSums(dd_ipm$pop_state$n_size)

 plot(pop_sizes, type = "l")

