
 library(ipmr)

 # Set up the initial population conditions and parameters

 data_list <- list(
   g_int     = 5.781,
   g_slope   = 0.988,
   g_sd      = 20.55699,
   s_int     = -0.352,
   s_slope   = 0.122,
   s_slope_2 = -0.000213,
   r_r_int   = -11.46,
   r_r_slope = 0.0835,
   r_s_int   = 2.6204,
   r_s_slope = 0.01256,
   r_d_mu    = 5.6655,
   r_d_sd    = 2.0734,
   e_p       = 0.15,
   g_i       = 0.5067
 )




 # We'll set up some helper functions. The survival function
 # in this model is a quadratic function, so we use an additional inverse logit function
 # that can handle the quadratic term.

 inv_logit <- function(int, slope, sv) {
   1/(1 + exp(-(int + slope * sv)))
 }

 inv_logit_2 <- function(int, slope, slope_2, sv) {
   1/(1 + exp(-(int + slope * sv + slope_2 * sv ^ 2)))
 }


 general_ipm <- init_ipm(sim_gen = "general", di_dd = "di", det_stoch = "det") %>%
   define_kernel(
     name          = "P",

     # We add d_ht to formula to make sure integration is handled correctly.
     # This variable is generated internally by make_ipm(), so we don't need
     # to do anything else.

     formula       = s * g * d_ht,

     # The family argument tells ipmr what kind of transition this kernel describes.
     # it can be "CC" for continuous -> continuous, "DC" for discrete -> continuous
     # "CD" for continuous -> discrete, or "DD" for discrete -> discrete.

     family        = "CC",

     # The rest of the arguments are exactly the same as in the simple models

     g             = dnorm(ht_2, g_mu, g_sd),
     g_mu          = g_int + g_slope * ht_1,
     s             = inv_logit_2(s_int, s_slope, s_slope_2, ht_1),
     data_list     = data_list,
     states        = list(c('ht')),
     uses_par_sets = FALSE,
     evict_cor     = TRUE,
     evict_fun     = truncated_distributions('norm',
                                             'g')
   ) %>%
   define_kernel(
     name          = "go_discrete",
     formula       = r_r * r_s * g_i * d_ht,

     # Note that now, family = "CD" because it denotes a continuous -> discrete transition

     family        = 'CD',
     r_r           = inv_logit(r_r_int, r_r_slope, ht_1),
     r_s           = exp(r_s_int + r_s_slope * ht_1),
     data_list     = data_list,

     # Note that here, we add "b" to our list in states, because this kernel
     # "creates" seeds entering the seedbank

     states        = list(c('ht', "b")),
     uses_par_sets = FALSE
   ) %>%
   define_kernel(
     name    = 'stay_discrete',

     # In this case, seeds in the seed bank either germinate or die, but they
     # do not remain for multiple time steps. This can be adjusted as needed.

     formula = 0,

     # Note that now, family = "DD" becuase it denotes a discrete -> discrete transition

     family  = "DD",

     # The only state variable this operates on is "b", so we can leave "ht"
     # out of the states list

     states  = list(c('b')),
     evict_cor = FALSE
   ) %>%
   define_kernel(

     # Here, the family changes to "DC" because it is the discrete -> continuous
     # transition

     name          = 'leave_discrete',

     # We append d_ht here as well, because we need to integrate over the
     # the recruit size distribution.

     formula       = e_p * r_d * d_ht,
     r_d           = dnorm(ht_2, r_d_mu, r_d_sd),
     family        = 'DC',
     data_list     = data_list,

     # Again, we need to add "b" to the states list

     states        = list(c('ht', "b")),
     uses_par_sets = FALSE,
     evict_cor     = TRUE,
     evict_fun     = truncated_distributions('norm',
                                             'r_d')
   )



 general_ipm <- general_ipm %>%
   define_impl(
     list(
       P              = list(int_rule    = "midpoint",
                             state_start = "ht",
                             state_end   = "ht"),
       go_discrete    = list(int_rule    = "midpoint",
                             state_start = "ht",
                             state_end   = "b"),
       stay_discrete  = list(int_rule    = "midpoint",
                             state_start = "b",
                             state_end   = "b"),
       leave_discrete = list(int_rule    = "midpoint",
                             state_start = "b",
                             state_end   = "ht")
     )
   )




 general_ipm <- general_ipm %>%
   define_impl(
     make_impl_args_list(
       kernel_names = c("P", "go_discrete", "stay_discrete", "leave_discrete"),
       int_rule     = c(rep("midpoint", 4)),
       state_start    = c('ht', "ht", "b", "b"),
       state_end      = c('ht', "b", "b", 'ht')
     )
   )



 # The lower and upper bounds for the continuous state variable and the number
 # of meshpoints for the midpoint rule integration. We'll also create the initial
 # population vector from a random uniform distribution

 L <- 1.02
 U <- 624
 n <- 500

 set.seed(2312)

 init_pop_vec   <- runif(500)
 init_seed_bank <- 20


 general_ipm <- general_ipm %>%
   define_domains(

     # We can pass the variables we created above into define_domains

     ht = c(L, U, n)

   ) %>%
   define_pop_state(

     # We can also pass them into define_pop_state

     pop_vectors = list(
       n_ht = init_pop_vec,
       n_b  = init_seed_bank
     )
   ) %>%
   make_ipm(iterations = 100,
            usr_funs = list(inv_logit   = inv_logit,
                            inv_logit_2 = inv_logit_2))


 # lambda is a generic function to compute per-capita growth rates. It has a number
 # of different options depending on the type of model

 lambda(general_ipm)

 # If we are worried about whether or not the model converged to stable
 # dynamics, we can use the exported utility is_conv_to_asymptotic. The default
 # tolerance for convergence is 1e-10, but can be changed with the 'tol' argument.

 is_conv_to_asymptotic(general_ipm, tol = 1e-10)

 w <- right_ev(general_ipm)
 v <- left_ev(general_ipm)





 make_N <- function(ipm) {

   P     <- ipm$sub_kernels$P

   I     <- diag(nrow(P))
   N     <- solve(I - P)

   return(N)
 }

 eta_bar <- function(ipm) {

   N     <- make_N(ipm)
   out   <- colSums(N)

   return(as.vector(out))

 }

 sigma_eta <- function(ipm) {

   N     <- make_N(ipm)

   out <- colSums(2 * (N %^% 2L) - N) - colSums(N) ^ 2

   return(as.vector(out))
 }

 mean_l <- eta_bar(general_ipm)
 var_l  <- sigma_eta(general_ipm)

 mesh_ps <- int_mesh(general_ipm)$ht_1 %>%
   unique()

 par(mfrow = c(1,2))

 plot(mesh_ps, mean_l, type = "l", xlab = expression( "Initial size z"[0]))
 plot(mesh_ps, var_l, type = "l", xlab = expression( "Initial size z"[0]))



 library(ipmr)

 # Set up the initial population conditions and parameters
 # Here, we are simulating random intercepts for growth
 # and seed production, converting them to a list,
 # and adding them into the list of constants. Equivalent code
 # to produce the output for the output from lmer/glmer
 # is in the comments next to each line

 all_g_int   <- as.list(rnorm(5, mean = 5.781, sd = 0.9)) # as.list(unlist(ranef(my_growth_model)))
 all_r_s_int <- as.list(rnorm(5, mean = 2.6204, sd = 0.3)) # as.list(unlist(ranef(my_seed_model)))


 names(all_g_int)   <- paste("g_int_", 1:5, sep = "")
 names(all_r_s_int) <- paste("r_s_int_", 1:5, sep = "")

 constant_list <- list(
   g_slope   = 0.988,
   g_sd      = 20.55699,
   s_int     = -0.352,
   s_slope   = 0.122,
   s_slope_2 = -0.000213,
   r_r_int   = -11.46,
   r_r_slope = 0.0835,
   r_s_slope = 0.01256,
   r_d_mu    = 5.6655,
   r_d_sd    = 2.0734,
   e_p       = 0.15,
   g_i       = 0.5067
 )

 all_params <- c(constant_list, all_g_int, all_r_s_int)

 # The lower and upper bounds for the continuous state variable and the number
 # of meshpoints for the midpoint rule integration.

 L <- 1.02
 U <- 624
 n <- 500

 set.seed(2312)

 init_pop_vec   <- runif(500)
 init_seed_bank <- 20

 # add some helper functions. The survival function
 # in this model is a quadratic function, so we use an additional inverse logit function
 # that can handle the quadratic term.

 inv_logit <- function(int, slope, sv) {
   1/(1 + exp(-(int + slope * sv)))
 }

 inv_logit_2 <- function(int, slope, slope_2, sv) {
   1/(1 + exp(-(int + slope * sv + slope_2 * sv ^ 2)))
 }




 general_stoch_kern_ipm <- init_ipm(sim_gen    = "general",
                                    di_dd      = "di",
                                    det_stoch  = "stoch",
                                    kern_param = "kern") %>%
   define_kernel(

     # The kernel name gets indexed by _year to denote that there
     # are multiple possible kernels we can build with our parameter set.
     # The _year gets substituted by the values in "par_set_indices" in the
     # output, so in this example we will have P_1, P_2, P_3, P_4, and P_5

     name          = "P_year",

     # We also add _year to "g" to signify that it is going to vary across kernels.

     formula       = s * g_year * d_ht,
     family        = "CC",

     # Here, we add the suffixes again, ensuring they are expanded and replaced
     # during model building by the parameter names

     g_year           = dnorm(ht_2, g_mu_year, g_sd),
     g_mu_year        = g_int_year + g_slope * ht_1,
     s                = inv_logit_2(s_int, s_slope, s_slope_2, ht_1),
     data_list        = all_params,
     states           = list(c('ht')),

     # We set uses_par_sets to TRUE, signalling that we want to expand these
     # expressions across all levels of par_set_indices.

     uses_par_sets    = TRUE,
     par_set_indices = list(year = 1:5),
     evict_cor        = TRUE,

     # we also add the suffix to `target` here, because the value modified by
     # truncated_distributions is time-varying.

     evict_fun        = truncated_distributions('norm',
                                                target = 'g_year')
   ) %>%
   define_kernel(

     # again, we append the index to the kernel name, vital rate expressions,
     # and in the model formula.

     name          = "go_discrete_year",
     formula       = r_r * r_s_year * g_i * d_ht,
     family        = 'CD',
     r_r           = inv_logit(r_r_int, r_r_slope, ht_1),

     # Again, we modify the left and right hand side of this expression to
     # show that there is a time-varying component

     r_s_year      = exp(r_s_int_year + r_s_slope * ht_1),
     data_list     = all_params,
     states        = list(c('ht', "b")),
     uses_par_sets    = TRUE,
     par_set_indices = list(year = 1:5)
   ) %>%
   define_kernel(

     # This kernel has no time-varying parameters, and so is not indexed.

     name    = 'stay_discrete',

     # In this case, seeds in the seed bank either germinate or die, but they
     # do not remain for multipe time steps. This can be adjusted as needed.

     formula = 0,

     # Note that now, family = "DD" becuase it denotes a discrete -> discrete transition

     family  = "DD",
     states  = list(c('b')),

     # This kernel has no time-varying parameters, so we don't need to designate
     # it as such.

     uses_par_sets = FALSE,
     evict_cor = FALSE
   ) %>%
   define_kernel(

     # This kernel also doesn't get a index, because there are no varying parameters.

     name          = 'leave_discrete',
     formula       = e_p * r_d * d_ht,
     r_d           = dnorm(ht_2, r_d_mu, r_d_sd),
     family        = 'DC',
     data_list     = all_params,
     states        = list(c('ht', "b")),
     uses_par_sets = FALSE,
     evict_cor     = TRUE,
     evict_fun     = truncated_distributions('norm',
                                             'r_d')
   )  %>%
   define_impl(
     make_impl_args_list(

       # We add suffixes to the kernel names here to make sure they match the names
       # we specified above.

       kernel_names = c("P_year", "go_discrete_year", "stay_discrete", "leave_discrete"),
       int_rule     = c(rep("midpoint", 4)),
       state_start    = c('ht', "ht", "b", "b"),
       state_end      = c('ht', "b", "b", 'ht')
     )
   ) %>%
   define_domains(
     ht = c(L, U, n)
   ) %>%
   define_pop_state(
       n_ht = init_pop_vec,
       n_b  = init_seed_bank
   ) %>%
   make_ipm(
     iterations = 100,

     # We can specify a sequence of kernels to select for the simulation.
     # This helps others to reproduce what we did,
     # and lets us keep track of the consequences of different selection
     # sequences for population dynamics.

     kernel_seq = sample(1:5, size = 100, replace = TRUE),
     usr_funs = list(inv_logit   = inv_logit,
                     inv_logit_2 = inv_logit_2)
   )



 mean_kernels <- mean_kernel(general_stoch_kern_ipm)
 lam_s        <- lambda(general_stoch_kern_ipm, burn_in = 0.15) # Remove first 15% of iterations



 library(ipmr)

 # Define the fixed parameters in a list

 constant_params <- list(
   s_int     = -5,
   s_slope   = 2.2,
   s_precip  = 0.0002,
   s_temp    = -0.003,
   g_int     = 0.2,
   g_slope   = 1.01,
   g_sd      = 1.2,
   g_temp    = -0.002,
   g_precip  = 0.004,
   c_r_int   = 0.3,
   c_r_slope = 0.03,
   c_s_int   = 0.4,
   c_s_slope = 0.01,
   c_d_mu    = 1.1,
   c_d_sd    = 0.1,
   r_e       = 0.3,
   r_d       = 0.3,
   r_r       = 0.2,
   r_s       = 0.2
 )




 # Now, we create a set of environmental covariates. In this example, we use
 # a normal distribution for temperature and a Gamma for precipitation.

 env_params <- list(
   temp_mu = 8.9,
   temp_sd = 1.2,
   precip_shape = 1000,
   precip_rate  = 2
 )

 # We define a wrapper function that samples from these distributions

 sample_env <- function(env_params) {

   # We generate one value for each covariate per iteration, and return it
   # as a named list.

   temp_now <- rnorm(1,
                     env_params$temp_mu,
                     env_params$temp_sd)

   precip_now <- rgamma(1,
                        shape = env_params$precip_shape,
                        rate = env_params$precip_rate)

   # The vital rate expressions can now use the names "temp" and "precip"
   # as if they were in the data_list.

   out        <- list(temp = temp_now, precip = precip_now)

   return(out)

 }

 # Again, we can define our own functions and pass them into calls to make_ipm. This
 # isn't strictly necessary, but can make the model code more readable/less error prone.

 inv_logit <- function(lin_term) {
   1/(1 + exp(-lin_term))
 }




 general_stoch_param_model <- init_ipm(sim_gen    = "general",
                                       di_dd      = "di",
                                       det_stoch  = "stoch",
                                       kern_param = "param") %>%
   define_kernel(
     name       = "P_stoch",
     family     = "CC",

     # As in the examples above, we have to add the d_surf_area
     # to ensure the integration of the functions is done.

     formula    = s * g * d_surf_area,

     # We can reference continuously varying parameters by name
     # in the vital rate expressions just as before, even though
     # they are passed in define_env_state() as opposed to the kernel's
     # data_list

     g_mu    = g_int + g_slope * surf_area_1 + g_temp * temp + g_precip * precip,
     s_lin_p = s_int + s_slope * surf_area_1 + s_temp * temp + s_precip * precip,
     s       = inv_logit(s_lin_p),
     g       = dnorm(surf_area_2, g_mu, g_sd),

     data_list     = constant_params,
     states        = list(c("surf_area")),
     uses_par_sets = FALSE,
     evict_cor     = TRUE,
     evict_fun     = truncated_distributions("norm", "g")
   ) %>%
   define_kernel(
     name          = "F",
     family        = "CC",
     formula       = c_r * c_s * c_d * (1 - r_e) * d_surf_area,

     c_r_lin_p     = c_r_int + c_r_slope * surf_area_1,
     c_r           = inv_logit(c_r_lin_p),
     c_s           = exp(c_s_int + c_s_slope * surf_area_1),
     c_d           = dnorm(surf_area_2, c_d_mu, c_d_sd),
     data_list     = constant_params,
     states        = list(c("surf_area")),
     uses_par_sets = FALSE,
     evict_cor     = TRUE,
     evict_fun     = truncated_distributions("norm", "c_d")
   ) %>%
   define_kernel(

     # Name can be anything, but it helps to make sure they're descriptive

     name = "go_discrete",

     # Family is now "CD" because it is a continuous -> discrete transition

     family        = "CD",
     formula       = r_e * r_s * c_r * c_s * d_surf_area,
     c_r_lin_p     = c_r_int + c_r_slope * surf_area_1,
     c_r           = inv_logit(c_r_lin_p),
     c_s           = exp(c_s_int + c_s_slope * surf_area_1),
     data_list     = constant_params,
     states        = list(c("surf_area", "sb")),
     uses_par_sets = FALSE,

     # There is not eviction to correct here, so we can set this to false

     evict_cor     = FALSE
   ) %>%
   define_kernel(
     name          = "stay_discrete",
     family        = "DD",
     formula       = r_s * r_r,
     data_list     = constant_params,
     states        = list("sb"),
     uses_par_sets = FALSE,
     evict_cor     = FALSE
   ) %>%
   define_kernel(
     name          = "leave_discrete",
     family        = "DC",
     formula       = r_d * r_s * c_d * d_surf_area,
     c_d           = dnorm(surf_area_2, c_d_mu, c_d_sd),
     data_list     = constant_params,
     states        = list(c("surf_area", "sb")),
     uses_par_sets = FALSE,
     evict_cor     = TRUE,
     evict_fun     = truncated_distributions("norm", "c_d")
   ) %>%
   define_impl(
     make_impl_args_list(
       kernel_names = c("P_stoch",
                        "F",
                        "go_discrete",
                        "stay_discrete",
                        "leave_discrete"),
       int_rule     = rep("midpoint", 5),
       state_start  = c("surf_area",
                        "surf_area",
                        "surf_area",
                        "sb",
                        "sb"),
       state_end    = c("surf_area",
                        "surf_area",
                        "sb",
                        "sb",
                        "surf_area")
     )
   ) %>%
   define_domains(
     surf_area = c(0, 10, 100)
   )



 # In the first version, sample_env is provided in the data_list of
 # define_env_state.

 general_stoch_param_ipm <-  define_env_state(
   proto_ipm  = general_stoch_param_model,
   env_covs   = sample_env(env_params),
   data_list  = list(env_params = env_params,
                     sample_env = sample_env)
 ) %>%
   define_pop_state(
     n_surf_area = runif(100),
     n_sb        = rpois(1, 20)
   ) %>%
   make_ipm(usr_funs = list(inv_logit  = inv_logit),
            iterate = TRUE,
            iterations = 100)


 # in the second version, sample_env is provided in the usr_funs list of
 # make_ipm(). These two versions are equivalent.

 general_stoch_param_ipm <-  define_env_state(
   proto_ipm  = general_stoch_param_model,
   env_covs   = sample_env(env_params),
   data_list  = list(env_params = env_params)
 ) %>%
   define_pop_state(
     n_surf_area = runif(100),
     n_sb        = rpois(1, 20)
   ) %>%
   make_ipm(usr_funs = list(inv_logit  = inv_logit,
                            sample_env = sample_env),
            iterate = TRUE,
            iterations = 100)




 env_draws  <- general_stoch_param_ipm$env_seq

 mean_kerns <- mean_kernel(general_stoch_param_ipm)

 lam_s      <- lambda(general_stoch_param_ipm)

 all.equal(mean_kerns$mean_F,
           general_stoch_param_ipm$sub_kernels$F_it_1,
           check.attributes = FALSE)



 data_list <- list(
   g_int     = 5.781,
   g_slope   = 0.988,
   g_sd      = 20.55699,
   s_int     = -0.352,
   s_slope   = 0.122,
   s_slope_2 = -0.000213,
   r_r_int   = -11.46,
   r_r_slope = 0.0835,
   r_s_int   = 2.6204,
   r_s_slope = 0.01256,
   r_d_mu    = 5.6655,
   r_d_sd    = 2.0734,
   e_p       = 0.15,
   g_i       = 0.5067
 )


 L <- 1.02
 U <- 624
 n <- 500

 set.seed(2312)

 init_pop_vec   <- runif(500)
 init_seed_bank <- 20

 # Initialize the state list and add some helper functions. The survival function
 # in this model is a quadratic function.

 inv_logit <- function(int, slope, sv) {
   1/(1 + exp(-(int + slope * sv)))
 }

 inv_logit_2 <- function(int, slope, slope_2, sv) {
   1/(1 + exp(-(int + slope * sv + slope_2 * sv ^ 2)))
 }

 general_ipm <- init_ipm(sim_gen    = "general",
                         di_dd      = "di",
                         det_stoch  = "det") %>%
   define_kernel(
     name          = "P",
     formula       = s * g * d_ht,
     family        = "CC",
     g             = dnorm(ht_2, g_mu, g_sd),
     g_mu          = g_int + g_slope * ht_1,
     s             = inv_logit_2(s_int, s_slope, s_slope_2, ht_1),
     data_list     = data_list,
     states        = list(c('ht')),
     uses_par_sets = FALSE,
     evict_cor     = TRUE,
     evict_fun     = truncated_distributions('norm',
                                             'g')
   ) %>%
   define_kernel(
     name          = "go_discrete",
     formula       = r_r * r_s * g_i,
     family        = 'CD',
     r_r           = inv_logit(r_r_int, r_r_slope, ht_1),
     r_s           = exp(r_s_int + r_s_slope * ht_1),
     data_list     = data_list,
     states        = list(c('ht', "b")),
     uses_par_sets = FALSE
   ) %>%
   define_kernel(
     name      = 'stay_discrete',
     formula   = 0,
     family    = "DD",
     states    = list(c('ht', "b")),
     evict_cor = FALSE
   ) %>%
   define_kernel(
     name          = 'leave_discrete',
     formula       = e_p * r_d * d_ht,
     r_d           = dnorm(ht_2, r_d_mu, r_d_sd),
     family        = 'DC',
     data_list     = data_list,
     states        = list(c('ht', "b")),
     uses_par_sets = FALSE,
     evict_cor     = TRUE,
     evict_fun     = truncated_distributions('norm',
                                             'r_d')
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
            usr_funs = list(inv_logit   = inv_logit,
                            inv_logit_2 = inv_logit_2))


 mega_mat <- make_iter_kernel(ipm = general_ipm,
                              mega_mat = c(
                                stay_discrete, go_discrete,
                                leave_discrete, P
                              ))

 # These values should be almost identical, so this should ~0

 Re(eigen(mega_mat[[1]])$values[1]) - lambda(general_ipm)




 # Get the names of each sub_kernel
 sub_k_nms     <- names(general_ipm$sub_kernels)

 mega_mat_text <- c(sub_k_nms[3], sub_k_nms[2], sub_k_nms[4], sub_k_nms[1])

 mega_mat_2 <- make_iter_kernel(general_ipm,
                                mega_mat = mega_mat_text)

 # Should be TRUE
 identical(mega_mat, mega_mat_2)



 mega_mat <- make_iter_kernel(general_ipm,
                                mega_mat = c(P, 0,
                                             I, P))






 all_g_int   <- as.list(rnorm(5, mean = 5.781, sd = 0.9))
 all_f_s_int <- as.list(rnorm(5, mean = 2.6204, sd = 0.3))

 names(all_g_int)   <- paste("g_int_", 1:5, sep = "")
 names(all_f_s_int) <- paste("f_s_int_", 1:5, sep = "")

 constant_list <- list(
   g_slope   = 0.988,
   g_sd      = 20.55699,
   s_int     = -0.352,
   s_slope   = 0.122,
   s_slope_2 = -0.000213,
   f_r_int   = -11.46,
   f_r_slope = 0.0835,
   f_s_slope = 0.01256,
   f_d_mu    = 5.6655,
   f_d_sd    = 2.0734,
   e_p       = 0.15,
   g_i       = 0.5067
 )

 all_params <- c(constant_list, all_g_int, all_f_s_int)

 L <- 1.02
 U <- 624
 n <- 500

 set.seed(2312)

 init_pop_vec   <- runif(500)
 init_seed_bank <- 20

 inv_logit <- function(int, slope, sv) {
   1/(1 + exp(-(int + slope * sv)))
 }

 inv_logit_2 <- function(int, slope, slope_2, sv) {
   1/(1 + exp(-(int + slope * sv + slope_2 * sv ^ 2)))
 }

 general_stoch_kern_ipm <- init_ipm(sim_gen    = "general",
                                    di_dd      = "di",
                                    det_stoch  = "stoch",
                                    kern_param = "kern") %>%
   define_kernel(
     name          = "P_year",
     formula       = s * g_year * d_ht,
     family        = "CC",
     g_year           = dnorm(ht_2, g_mu_year, g_sd),
     g_mu_year        = g_int_year + g_slope * ht_1,
     s                = inv_logit_2(s_int, s_slope, s_slope_2, ht_1),
     data_list        = all_params,
     states           = list(c('ht')),
     uses_par_sets    = TRUE,
     par_set_indices = list(year = 1:5),
     evict_cor        = TRUE,
     evict_fun        = truncated_distributions('norm',
                                                'g_year')
   ) %>%
   define_kernel(
     name          = "go_discrete_year",
     formula       = f_r * f_s_year * g_i * d_ht,
     family        = 'CD',
     f_r           = inv_logit(f_r_int, f_r_slope, ht_1),
     f_s_year      = exp(f_s_int_year + f_s_slope * ht_1),
     data_list     = all_params,
     states        = list(c('ht', "b")),
     uses_par_sets    = TRUE,
     par_set_indices = list(year = 1:5)
   ) %>%
   define_kernel(
     name    = 'stay_discrete',
     formula = 0,
     family  = "DD",
     states  = list(c('b')),
     uses_par_sets = FALSE,
     evict_cor = FALSE
   ) %>%
   define_kernel(
     name          = 'leave_discrete',
     formula       = e_p * f_d * d_ht,
     f_d           = dnorm(ht_2, f_d_mu, f_d_sd),
     family        = 'DC',
     data_list     = all_params,
     states        = list(c('ht', "b")),
     uses_par_sets = FALSE,
     evict_cor     = TRUE,
     evict_fun     = truncated_distributions('norm',
                                             'f_d')
   )  %>%
   define_impl(
     make_impl_args_list(
       kernel_names = c("P_year", "go_discrete_year", "stay_discrete", "leave_discrete"),
       int_rule     = c(rep("midpoint", 4)),
       state_start    = c('ht', "ht", "b", "b"),
       state_end      = c('ht', "b", "b", 'ht')
     )
   ) %>%
   define_domains(
     ht = c(L, U, n)
   ) %>%
   define_pop_state(
       n_ht = init_pop_vec,
       n_b  = init_seed_bank
   ) %>%
   make_ipm(
     iterations = 10,
     kernel_seq = sample(1:5, size = 10, replace = TRUE),
     usr_funs = list(inv_logit   = inv_logit,
                     inv_logit_2 = inv_logit_2)
   )



 block_list <- make_iter_kernel(general_stoch_kern_ipm,
                                  mega_mat = c(stay_discrete, go_discrete_year,
                                               leave_discrete, P_year))


