
library(ipmr)

# Set parameter values and names

param_list <- list(
  ## Survival
  surv_int  = -1.70e+1,
  surv_z    =  6.68e+0,
  surv_a    = -3.34e-1,
  ## growth
  grow_int  =  1.27e+0,
  grow_z    =  6.12e-1,
  grow_a    = -7.24e-3,
  grow_sd   =  7.87e-2,
  ## reproduce or not
  repr_int  = -7.88e+0,
  repr_z    =  3.11e+0,
  repr_a    = -7.80e-2,
  ## recruit or not
  recr_int  =  1.11e+0,
  recr_a    =  1.84e-1,
  ## recruit size
  rcsz_int  =  3.62e-1,
  rcsz_z    =  7.09e-1,
  rcsz_sd   =  1.59e-1
)

# define a custom function to handle the F kernels. We could write a rather
# verbose if(age == 0) {0} else {other_math} in the define_kernel(), but that
# might look ugly. Note that we CANNOT use ifelse(), as its output is the same
# same length as its input (in this case, it would return 1 number, not 10000
# numbers).

r_fun <- function(age, s_age, pb_age, pr_age, recr) {

  if(age == 0) return(0)

  s_age * pb_age * pr_age * recr * 0.5

}



age_size_ipm <- init_ipm(sim_gen = "general",
                         di_dd = "di",
                         det_stoch = "det",
                         uses_age = TRUE) %>%
  define_kernel(
    name          = "P_age",
    family        = "CC",
    formula       = s_age * g_age * d_z,
    s_age         = plogis(surv_int + surv_z * z_1 + surv_a * age),
    g_age         = dnorm(z_2, mu_g_age, grow_sd),
    mu_g_age      = grow_int + grow_z * z_1 + grow_a * age,
    data_list     = param_list,
    states        = list(c("z")),
    uses_par_sets = FALSE,
    age_indices   = list(age = c(0:20), max_age = 21),
    evict_cor     = FALSE
  )


age_size_ipm <-   define_kernel(
  proto_ipm     = age_size_ipm,
  name          = "F_age",
  family        = "CC",
  formula       = r_fun(age, s_age, pb_age, pr_age, recr) * d_z,
  s_age         = plogis(surv_int + surv_z * z_1 + surv_a * age),
  pb_age        = plogis(repr_int + repr_z * z_1 + repr_a * age),
  pr_age        = plogis(recr_int + recr_a * age),
  recr          = dnorm(z_2, rcsz_mu, rcsz_sd),
  rcsz_mu       = rcsz_int + rcsz_z * z_1,
  data_list     = param_list,
  states        = list(c("z")),
  uses_par_sets = FALSE,
  age_indices   = list(age = c(0:20), max_age = 21),
  evict_cor     = FALSE
)




age_size_ipm <-  define_impl(
  proto_ipm = age_size_ipm,
  make_impl_args_list(
    kernel_names = c("P_age", "F_age"),
    int_rule     = rep("midpoint", 2),
    state_start    = c("z_age", "z_age"),
    state_end      = c("z_age", "z_0")
  )
)




age_size_ipm <- age_size_ipm %>%
  define_domains(
    z = c(1.6, 3.7, 100)
  )



age_size_ipm <-  define_pop_state(
  proto_ipm = age_size_ipm,
  n_z_age = rep(1/100, 100)
  ) %>%
  make_ipm(
    usr_funs = list(r_fun = r_fun),
    iterate  = TRUE,
    iterations = 100,
    return_all_envs = TRUE,
    return_sub_kernels = TRUE
  )




lamb <- lambda(age_size_ipm)
lamb

is_conv_to_asymptotic(age_size_ipm)



vr_funs <- vital_rate_funs(age_size_ipm)

# Age 0 survival and growth vital rate functions

vr_funs$P_0

# Age 12 fecundity functions

vr_funs$F_12





new_proto <- age_size_ipm$proto_ipm

vital_rate_exprs(new_proto,
                 kernel = "F_age",
                 vital_rate = "pr_age") <-
  new_fun_form(plogis(recr_int + recr_z * z_1 + recr_a * age))

parameters(new_proto) <- list(recr_z = 0.05)

new_ipm <- make_ipm(new_proto,
                    return_all_envs = TRUE,
                    return_sub_kernels = TRUE)

lambda(new_ipm)



options(width = 95)



d_z <- int_mesh(age_size_ipm)$d_z

stable_dists <- right_ev(age_size_ipm)

w_plot <- lapply(stable_dists, function(x, d_z) x / d_z,
         d_z = d_z)

repro_values <- left_ev(age_size_ipm)

v_plot <- lapply(repro_values, function(x, d_z) x / d_z,
                 d_z = d_z)


par(mfrow = c(1, 2))
plot(w_plot[[1]], type = 'l',
     ylab = expression(paste("w"[a],"(z)")),
     xlab = "Size bin")

x <- lapply(w_plot[2:22], function(x) lines(x))

plot(v_plot[[1]], type = 'l',
     ylab = expression(paste("v"[a],"(z)")),
     xlab = "Size bin",
     ylim = c(0, 0.2))

x <- lapply(v_plot[2:22], function(x) lines(x))



# Initialize a cohort and some vectors to hold our quantities of interest.

init_pop <- stable_dists[[1]] / sum(stable_dists[[1]])
n_yrs    <- 100L

l_a <- f_a <- numeric(n_yrs)

P_kerns <- age_size_ipm$sub_kernels[grepl("P", names(age_size_ipm$sub_kernels))]
F_kerns <- age_size_ipm$sub_kernels[grepl("F", names(age_size_ipm$sub_kernels))]

# We have to bump max_age because R indexes from 1, and our minimum age is 0.

max_age <- 22

P_a <- diag(length(init_pop))

for(yr in seq_len(n_yrs)) {

  # When we start, we want to use age-specific kernels until we reach max_age.
  # after that, all survivors have entered the "greybeard" class.

  if(yr < max_age) {

    P_now <- P_kerns[[yr]]
    F_now <- F_kerns[[yr]]

  } else {

    P_now <- P_kerns[[max_age]]
    F_now <- F_kerns[[max_age]]

  }

  l_a[yr] <- sum(colSums(P_a) * init_pop)
  f_a[yr] <- sum(colSums(F_now %*% P_a) * init_pop)

  P_a <- P_now %*% P_a
}

f_a <- f_a / l_a

# Looks like most are dead at after 25 years, so we'll restrict our
# plot range to that time span

par(mfrow = c(1, 2))

plot(l_a[1:25], type = 'l',
     ylab = expression(paste("l"[a])),
     xlab = "Age")
plot(f_a[1:25], type = 'l',
     ylab = expression(paste("f"[a])),
     xlab = "Age")




p_a <- l_a[2:26] / l_a[1:25]

plot(p_a, type = 'l',
     ylab = expression(paste("p"[a])),
     xlab = "Age")

