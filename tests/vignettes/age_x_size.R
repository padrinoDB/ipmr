
 library(ipmr)

 param_list <- list(
   surv_int = -17,
   surv_z = 6.68,
   surv_a = -0.334,
   grow_int = 1.27,
   grow_z = 0.612,
   grow_a = -0.00724,
   grow_sd = 0.0787,
   repr_int = -7.88,
   repr_z = 3.11,
   repr_a = -0.078,
   recr_int = 1.11,
   recr_a = 0.184,
   rcsz_int = 0.362,
   rcsz_z = 0.709,
   rcsz_sd = 0.159
 )


 inv_logit <- function(x) {

   return( 1 / (1 + exp(-x)) )
 }

 f_fun <- function(age, s_age, pb_age, pr_age, recr) {

   if(age == 0) return(0)

   s_age * pb_age * pr_age * recr * 0.5

 }




 age_size_ipm <- init_ipm(sim_gen    = "general",
                          di_dd      = "di",
                          det_stoch  = "det",
                          uses_age    = TRUE) %>%
   define_kernel(
     name          = "P_age",
     family        = "CC",
     formula       = s_age * g_age * d_wt,
     s_age         = inv_logit(surv_int + surv_z * wt_1 + surv_a * age),
     g_age         = dnorm(wt_2, mu_g_age, grow_sd),
     mu_g_age      = grow_int + grow_z * wt_1 + grow_a * age,
     data_list     = param_list,
     states        = list(c("wt")),
     uses_par_sets = FALSE,
     age_indices   = list(age = c(0:20), max_age = 21),
     evict_cor     = FALSE
   )



 age_size_ipm <-
   define_kernel(
     proto_ipm     = age_size_ipm,
     name          = "F_age",
     family        = "CC",
     formula       = f_fun(age, s_age, pb_age, pr_age, rcsz) * d_wt,
     s_age         = inv_logit(surv_int + surv_z * wt_1 + surv_a * age),
     pb_age        = inv_logit(repr_int + repr_z * wt_1 + repr_a * age),
     pr_age        = inv_logit(recr_int + recr_a * age),
     rcsz          = dnorm(wt_2, rcsz_mu, rcsz_sd),
     rcsz_mu       = rcsz_int + rcsz_z * wt_1,
     data_list     = param_list,
     states        = list(c("wt")),
     uses_par_sets = FALSE,
     age_indices   = list(age = c(0:20), max_age = 21),
     evict_cor     = FALSE
   )



 age_size_ipm <- age_size_ipm %>%
   define_impl(
     make_impl_args_list(
       kernel_names = c("P_age", "F_age"),
       int_rule     = rep("midpoint", 2),
       state_start    = c("wt_age", "wt_age"),
       state_end      = c("wt_age", "wt_0")
     )
   )



 age_size_ipm <- age_size_ipm %>%
   define_domains(
     wt = c(1.6, 3.7, 100)
   ) %>%
   define_pop_state(
     n_wt_age = runif(100)
   ) %>%
   make_ipm(
     usr_funs = list(inv_logit = inv_logit,
                     f_fun     = f_fun),
     iterate  = TRUE,
     iterations = 100
   )

 lam <- lambda(age_size_ipm)
 lam



 v_a <- left_ev(age_size_ipm, iterations = 100)
 w_a <- right_ev(age_size_ipm, iterations = 100)

 par(mfrow = c(1, 2))

 plot(1:100, seq(0, max(unlist(w_a)), length.out = 100), type = "n",
      ylab = expression(paste("w"[a],"(z)")),
      xlab = "Size bin")

 for(i in seq_along(w_a)) {

   lines(w_a[[i]])

 }

 plot(1:100,
      seq(0, max(unlist(v_a)), length.out = 100),
      type = "n",
      ylab = expression(paste("v"[a], "(z)")),
      xlab = "Size bin")

 for(i in seq_along(v_a)) {

  lines(v_a[[i]])

 }

