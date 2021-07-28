
library(ipmr)

my_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det")




my_ipm <- define_kernel(
  proto_ipm = my_ipm,
  name      = "P",

  # The formula describes how the vital rates generate the sub-kernel. This is
  # equivalent to equation 3 above.

  formula   = s * g,

  # Next, we define the vital rate expressions for the P kernel (Eqs 4-6).
  # Here, we create an expression for the inverse of the link function from
  # our GLM for survival. In this case, it's the inverse logit (Eq 4).

  s         = 1/(1 + exp(-(s_int + s_slope * dbh_1))),

  # Growth is defined by Eqs 5-6 above. There is no inverse link function required
  # because we use a Gaussian GLM with an identity link function (i.e. no
  # transformation).

  g         = dnorm(dbh_2, g_mu, g_sd), # Eq 5
  g_mu      = g_int + g_slope * dbh_1,  # Eq 6

  # The family describes the type of transition that kernel produces. See below.

  family    = "CC",

  data_list = list(
    s_int     = 0.2,   # coef(my_surv_mod)[1]
    s_slope   = 0.5,   # coef(my_surv_mod)[2]
    g_int     = 0.1,   # coef(my_grow_mod)[1]
    g_slope   = 1.033, # coef(my_grow_mod)[2]
    g_sd      = 2.2    # sd(resid(my_grow_mod))
  ),

  # states should be a list of the state variables that the kernel operates on

  states        = list(c('dbh')),
  uses_par_sets = FALSE,
  evict_cor     = TRUE,
  evict_fun     = truncated_distributions("norm", "g")
)



# the %>% takes the result of the first operation and passes it as the first
# argument to the second function.

my_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det") %>%
  define_kernel(
    name      = "P",
    formula   =  s * g,
    family    = "CC",
    s         = 1/(1 + exp(-(s_int + s_slope * dbh_1))),
    g         = dnorm(dbh_2, g_mu, g_sd),
    g_mu      = g_int + g_slope * dbh_1,

    data_list = list(
      s_int     = -3,
      s_slope   = 0.5,
      g_int     = 0.1,
      g_slope   = 1.033,
      g_sd      = 2.2
    ),
    states        = list(c('dbh')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun   = "norm",
                                            target =  "g")
  )



my_ipm <- my_ipm %>%
  define_kernel(
    name      = "F",
    formula   = r_r * r_s * r_d, # Eq 7
    family    = "CC",

    # Eq 8. Again, we use the inverse logit transformation to compute pr(flowering)
    r_r       = 1/(1 + exp(-(r_r_int + r_r_slope * dbh_1))),

    # Eq 9. We exponentiate this because of the log link in our seed production model
    r_s       = exp(r_s_int + r_s_slope * dbh_1),

    # Eq 10. In this case, both the mean and standard deviation are constants

    r_d       = dnorm(dbh_2, r_d_mu, r_d_sd),
    data_list = list(
      r_r_int   = -5,   # coef(my_flower_mod)[1]
      r_r_slope = 0.1,   # coef(my_flower_mod)[2]
      r_s_int   = -3,   # coef(my_seed_mod)[1]
      r_s_slope = 0.03,  # coef(my_seed_mod)[2]
      r_d_mu    = 1.2,   # mean(my_recr_data$size_2, na.rm = TRUE)
      r_d_sd    = 0.7    # sd(my_recr_data$size_2, na.rm = TRUE)
    ),
    states        = list(c('dbh')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "r_d")
  )



my_ipm <- my_ipm %>%

  define_impl(

    # Create one list of lists with all the info

    list(

      # The components of the big list should be the sub-kernels. Each sub-kernel
      # requires 3 pieces of information: the numerical integration rule, the
      # state variable it modifies (state_start), and the state variable it
      # generates (state_end).

      P = list(int_rule  = "midpoint",
               state_start = "dbh",
               state_end   = "dbh"),
      F = list(int_rule  = "midpoint",
               state_start = "dbh",
               state_end   = "dbh")
    )
  )



# Alternative 1 - call make_impl_args_list() before beginning the IPM creation pipe

impl_args <- make_impl_args_list(
  kernel_names = c("P", "F"),
  int_rule     = rep("midpoint", 2),
  state_start    = rep("dbh", 2),
  state_end      = rep("dbh", 2)
)

my_ipm <- my_ipm %>%
  define_impl(impl_args)

my_ipm <- my_ipm %>%

  # Alternative 2, put the call to make_impl_args_list() inside of define_impl().

  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F"),
      int_rule     = rep("midpoint", 2),
      state_start    = rep("dbh", 2),
      state_end      = rep("dbh", 2)
    )
  )




my_ipm <- my_ipm %>%
  define_domains(
    dbh = c(        # the name of the state variable.
      1,            # the lower bound for the domain. L from Eq 1
      30,           # the upper bound for the domain. U from Eq 1
      200           # the number of mesh points to use for integration
    )
  )



my_ipm <- my_ipm %>%
  define_pop_state(n_dbh = rep(1/200, 200))



my_ipm <- make_ipm(my_ipm,
                   iterations = 100)

lambda_ipmr <- lambda(my_ipm)
repro_value <- left_ev(my_ipm, iterations = 200)
stable_dist <- right_ev(my_ipm, iterations = 100)


data(iceplant_ex)

# growth model

grow_mod <- lm(log_size_next ~ log_size, data = iceplant_ex)
grow_sd  <- sd(resid(grow_mod))


# survival model

surv_mod <- glm(survival ~ log_size, data = iceplant_ex, family = binomial())

# Pr(flowering) model

repr_mod <- glm(repro ~ log_size, data = iceplant_ex, family = binomial())

# Number of flowers per plant model

seed_mod <- glm(flower_n ~ log_size, data = iceplant_ex, family = poisson())

# New recruits have no size(t), but do have size(t + 1)

recr_data <- subset(iceplant_ex, is.na(log_size))

recr_mu  <- mean(recr_data$log_size_next)
recr_sd  <- sd(recr_data$log_size_next)

# This data set doesn't include information on germination and establishment.
# Thus, we'll compute the realized recruitment parameter as the number
# of observed recruits divided by the number of flowers produced in the prior
# year.

recr_n   <- length(recr_data$log_size_next)

flow_n   <- sum(iceplant_ex$flower_n, na.rm = TRUE)

# Now, we bundle everything into a list as we did before. We can
# pass the model objects directly into the list, and do not need to extract
# coefficients.

params <- list(recr_mu = recr_mu,
               recr_sd = recr_sd,
               grow_sd = grow_sd,
               surv_mod = surv_mod,
               grow_mod = grow_mod,
               repr_mod = repr_mod,
               seed_mod = seed_mod,
               recr_n   = recr_n,
               flow_n   = flow_n)

# The lower and upper bounds for integration. Adding 20% on either end to minimize
# eviction. The minimum value is negative, so we multiply by 1.2, rather than 0.8.
# If it were positive, we would need to change that to 0.8.

L <- min(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2
U <- max(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2



pred_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det") %>%
  define_kernel(
    name          = "P",
    family        = "CC",
    formula       = s * g,

    # Instead of the inverse logit transformation, we use predict() here.
    # We have to be sure that the "newdata" argument of predict is correctly specified.
    # This means matching the names used in the model itself (log_size) to the names
    # we give the domains (sa_1). In this case, we use "sa" (short for surface area) as
    # the new data to generate predictions for.

    s             = predict(surv_mod,
                            newdata = data.frame(log_size = sa_1),
                            type    = 'response'),
    g_mu          = predict(grow_mod,
                            newdata = data.frame(log_size = sa_1),
                            type    = 'response'),

    # We specify the rest of the kernel the same way.

    g             = dnorm(sa_2, g_mu, grow_sd),
    states        = list(c("sa")),
    data_list     = params,
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun   = "norm",
                                            target =  "g")
  ) %>%
  define_kernel(
    name          = "F",
    family        = "CC",
    formula       = r_p * r_s * r_d * r_r,

    # As above, we use predict(model_object). We make sure the names of the "newdata"
    # match the names in the vital rate model formulas, and the values match the
    # names of domains they use.

    r_p           = predict(repr_mod,
                            newdata = data.frame(log_size = sa_1),
                            type    = 'response'),
    r_s           = predict(seed_mod,
                            newdata = data.frame(log_size = sa_1),
                            type    = 'response'),
    r_r           = recr_n / flow_n,
    r_d           = dnorm(sa_2, recr_mu, recr_sd),
    states        = list(c("sa")),
    data_list     = params,
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun   = "norm",
                                            target =  "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F"),
      int_rule     = rep('midpoint', 2),
      state_start    = rep("sa", 2),
      state_end      = rep("sa", 2)
    )
  ) %>%
  define_domains(
    sa = c(L,
           U,
           100)
  ) %>%
  define_pop_state(n_sa = runif(100)) %>%
  make_ipm(iterate = TRUE,
           iterations = 200)



params <- list(recr_mu = recr_mu,
               recr_sd = recr_sd,
               grow_sd = grow_sd,
               surv_mod = use_vr_model(surv_mod),  # wrap the model
               grow_mod = use_vr_model(grow_mod),  # wrap the model
               repr_mod = use_vr_model(repr_mod),  # wrap the model
               seed_mod = use_vr_model(seed_mod),  # wrap the model
               recr_n   = recr_n,
               flow_n   = flow_n)


pred_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det") %>%
  define_kernel(
    name          = "P",
    family        = "CC",
    formula       = s * g,
    s             = predict(surv_mod,
                            newdata = data.frame(log_size = sa_1),
                            type = 'response'),
    g_mu          = predict(grow_mod,
                            newdata = data.frame(log_size = sa_1),
                            type = 'response'),
    g             = dnorm(sa_2, g_mu, grow_sd),
    states        = list(c("sa")),
    data_list     = params,
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "g")
  ) %>%
  define_kernel(
    name          = "F",
    family        = "CC",
    formula       = r_p * r_s * r_d * r_r,
    r_p           = predict(repr_mod,
                            newdata = data.frame(log_size = sa_1),
                            type = 'response'),
    r_s           = predict(seed_mod,
                            newdata = data.frame(log_size = sa_1),
                            type = 'response'),
    r_r           = recr_n / flow_n,
    r_d           = dnorm(sa_2, recr_mu, recr_sd),
    states        = list(c("sa")),
    data_list     = params,
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "r_d")
  )  %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F"),
      int_rule     = rep('midpoint', 2),
      state_start    = rep("sa", 2),
      state_end      = rep("sa", 2)
    )
  ) %>%
  define_domains(
    sa = c(L,
           U,
           100)
  ) %>%
  define_pop_state(n_sa = runif(100)) %>%
  make_ipm(iterate = TRUE,
           iterations = 200)



R_0 <- function(ipm) {

  Fm <- ipm$sub_kernels$F
  Pm <- ipm$sub_kernels$P

  I  <- diag(nrow(Pm))

  N  <- solve(I - Pm)

  R  <- Fm %*% N

  out <- Re(eigen(R)$values[1])

  return(out)

}

gen_T <- function(ipm) {

  R0  <- R_0(ipm)
  lam <- lambda(ipm) %>% as.vector()

  out <- log(R0) / log(lam)

  return(out)

}

R_0(pred_ipm)
gen_T(pred_ipm)



library(ipmr)

# Define some fixed parameters

fixed_list <- list(
  s_int     = 1.03,   # fixef(my_surv_mod)[1] - uses fixef because we now have a model with random effects
  s_slope   = 2.2,    # fixef(my_surv_mod)[2]
  g_int     = 3.7,    # fixef(my_grow_mod)[1]
  g_slope   = 0.92,   # fixef(my_grow_mod)[2]
  sd_g      = 0.9,    # sd(resid(my_grow_mod))
  r_r_int   = 0.09,   # coef(my_repro_mod)[1] - uses coef because there are no random effects in this model
  r_r_slope = 0.05,   # coef(my_repro_mod)[2]
  r_s_int   = 0.1,    # fixef(my_flower_mod)[1]
  r_s_slope = 0.005,  # fixef(my_flower_mod)[2]
  mu_rd     = 9,      # mean(my_recr_data$size_2, na.rm = TRUE)
  sd_rd     = 2       # sd(my_recr_data$size_2, na.rm = TRUE)
)



# Now, simulate some random intercepts for growth (g_), survival (s_),
# and offspring production (r_s_). This part is for the purpose of the example.

# First, we create vector of values that each random component can take.

g_r_int   <- rnorm(5, 0, 0.3) # unlist(ranef(my_grow_mod)) for an lme4 output
s_r_int   <- rnorm(5, 0, 0.7) # unlist(ranef(my_surv_mod)) for an lme4 output
r_s_r_int <- rnorm(5, 0, 0.2) # unlist(ranef(my_flower_mod)) for an lme4 output

# We'll call our hypothetical sites 1, 2, 3, 4, and 5. The "r" prefix is to
# remind us that these are random quantities.

nms <- paste("r_", 1:5, sep = "")

names(g_r_int)   <- paste('g_', nms, sep = "")
names(s_r_int)   <- paste('s_', nms, sep = "")
names(r_s_r_int) <- paste('r_s_', nms, sep = "")

# Each set of parameters is converted to a named list. The names should match
# the variables referenced in each define_kernel() call.

g_params   <- as.list(g_r_int)
s_params   <- as.list(s_r_int)
r_s_params <- as.list(r_s_r_int)

# add them all together using c()

all_params_list <- c(fixed_list, g_params, s_params, r_s_params)



my_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det") %>%
  define_kernel(

    # Our P kernels will vary from site to site, so we index it with "_site"

    name             = 'P_site',

    # Similarly, our survival and growth functions will vary from site to site
    # so these are also indexed

    formula          = s_site * g_site,
    family           = "CC",

    # The linear predictor for the survival function can be split out
    # into its own expression as well. This might help keep track of things.
    # Survival is indexed by site as well.

    s_lin_site       = s_int + s_r_site + s_slope * ht_1,
    s_site           = 1 / (1 + exp(-s_lin_site)),

    # Again, we modify the vital rate expression to include "_site".

    g_site           = dnorm(ht_2, mean = mu_g_site, sd = sd_g),
    mu_g_site        = g_int + g_slope * ht_1 + g_r_site,

    data_list        = all_params_list,
    states           = list(c('ht')),

    # Here, we tell ipmr that the model has some parameter sets, and
    # provide a list describing the values the index can take. The values in
    # par_set_indices are substituted for "site" everywhere in the model, except
    # for the data list. This is why we had to make sure that the names there
    # matched the levels we supply here.

    uses_par_sets    = TRUE,
    par_set_indices  = list(site = 1:5),

    # We must also index the variables in the eviction function

    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "g_site")

  ) %>%
  define_kernel(

    # The F kernel also varies from site to site

    name             = "F_site",
    formula          = r_r * r_s_site * r_d,
    family           = "CC",

    # In this example, we didn't include a site level effect for probability
    # of flowering, only seed production. Thus, this expression is NOT indexed.

    r_r_lin          = r_r_int + r_r_slope * ht_1,
    r_r              = 1 / (1 + exp(- r_r_lin)),

    # We index the seed production expression with the site effect

    r_s_site         = exp(r_s_int + r_s_r_site + r_s_slope * ht_1),
    r_d              = dnorm(ht_2, mean = mu_rd, sd = sd_rd),
    data_list        = all_params_list,
    states           = list(c('ht')),

    # As in the P kernel, we specify the values the index can have.

    uses_par_sets    = TRUE,
    par_set_indices  = list(site = 1:5),
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(

      # The impl_args are also modified with the index

      kernel_names = c("P_site", "F_site"),
      int_rule     = rep("midpoint", 2),
      state_start    = rep("ht", 2),
      state_end      = rep("ht", 2)
    )
  ) %>%
  define_domains(ht = c(0.2, 40, 100)) %>%

  # We also append the suffix in define_pop_state(). THis will create a deterministic
  # simulation for every "site"

  define_pop_state(n_ht_site = runif(100)) %>%
  make_ipm(iterate  = TRUE,
           iterations = 100)

lambda(my_ipm)




library(ipmr)

# Define some fixed parameters

fixed_list <- list(
  s_int     = 1.03,   # fixef(my_surv_mod)[1] - uses fixef because we now have a model with random effects
  s_slope   = 2.2,    # fixef(my_surv_mod)[2]
  g_int     = 3.7,    # fixef(my_grow_mod)[1]
  g_slope   = 0.92,   # fixef(my_grow_mod)[2]
  sd_g      = 0.9,    # sd(resid(my_grow_mod))
  r_r_int   = 0.09,   # coef(my_repro_mod)[1] - uses coef because there are no random effects in this model
  r_r_slope = 0.05,   # coef(my_repro_mod)[2]
  r_s_int   = 0.1,    # fixef(my_flower_mod)[1]
  r_s_slope = 0.005,  # fixef(my_flower_mod)[2]
  mu_fd     = 9,      # mean(my_recr_data$size_2, na.rm = TRUE)
  sd_fd     = 2       # sd(my_recr_data$size_2, na.rm = TRUE)
)


# Now, simulate some random intercepts for growth (g_), survival (s_),
# and offspring production (r_s_). This part is for the purpose of the example.

# First, we create vector of values corresponding to

g_r_int   <- rnorm(5, 0, 0.3) # unlist(ranef(my_grow_mod)) for an lme4 output
s_r_int   <- rnorm(5, 0, 0.7) # unlist(ranef(my_surv_mod)) for an lme4 output
r_s_r_int <- rnorm(5, 0, 0.2) # unlist(ranef(my_flower_mod)) for an lme4 output

nms <- paste("r_", 1:5, sep = "")

names(g_r_int)   <- paste('g_', nms, sep = "")
names(s_r_int)   <- paste('s_', nms, sep = "")
names(r_s_r_int) <- paste('r_s_', nms, sep = "")

# Each set of parameters is converted to a named list. The names should match
# the variables referenced in each define_kernel() call.

g_params   <- as.list(g_r_int)
s_params   <- as.list(s_r_int)
r_s_params <- as.list(r_s_r_int)

# add them all together using c()

all_params_list <- c(fixed_list, g_params, s_params, r_s_params)



inv_logit <- function(sv, int, slope) {
  return(
    1/(1 + exp(-(int + slope * sv)))
  )
}

# same as above, but handles the extra term from the random effect we simulated.

inv_logit_r <- function(sv, int, slope, r_eff) {
  return(
    1/(1 + exp(-(int + slope * sv + r_eff)))
  )
}

pois_r <- function(sv, int, slope, r_eff) {
  return(
    exp(
      int + slope * sv + r_eff
    )
  )
}

my_funs <- list(inv_logit   = inv_logit,
                inv_logit_r = inv_logit_r,
                pois_r      = pois_r)


my_ipm <- init_ipm(sim_gen    = "simple",
                   di_dd      = "di",
                   det_stoch  = "stoch",
                   kern_param = "kern") %>%
  define_kernel(
    name             = 'P_yr',         # P becomes P_yr
    formula          = s_yr * g_yr,    # g and s become g_yr and s_yr, respectively
    family           = "CC",

    # Note the usage of the inv_logit_r, which we defined in the block above.
    # it is passed to make_ipm()

    s_yr             = inv_logit_r(ht_1, s_int, s_slope, s_r_yr),
    g_yr             = dnorm(ht_2, mean = mu_g_yr, sd = sd_g),
    mu_g_yr          = g_int + g_slope * ht_1 + g_r_yr,

    # all_params_list contains the named parameters g_r_1, g_r_2, s_r_1, s_r_2, etc.
    # This is the only level where the user is required to fully expand the name
    # X par_set_indices combinations.

    data_list        = all_params_list,
    states           = list(c('ht')),
    uses_par_sets    = TRUE,
    par_set_indices  = list(yr = 1:5),
    evict_cor        = TRUE,

    # reference to g_yr in evict_fun is also updated

    evict_fun        = truncated_distributions("norm", "g_yr")

  ) %>%
  define_kernel(
    name             = "F_yr",             # Update the names as we did for the P kernel
    formula          = r_r * r_s_yr * r_d,
    family           = "CC",
    r_r              = inv_logit(ht_1, r_r_int, r_r_slope),
    r_s_yr           = pois_r(ht_1, r_s_int, r_s_slope, r_s_r_yr),
    r_d              = dnorm(ht_2, mean = mu_fd, sd = sd_fd),
    data_list        = all_params_list,
    states           = list(c('ht')),
    uses_par_sets    = TRUE,
    par_set_indices = list(yr = 1:5),
    evict_cor        = TRUE,
    evict_fun        = truncated_distributions("norm", "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P_yr", "F_yr"),
      int_rule     = rep("midpoint", 2),
      state_start    = rep("ht", 2),
      state_end      = rep("ht", 2)
    )
  ) %>%
  define_domains(ht = c(0.2, 40, 100)) %>%
  define_pop_state(n_ht = runif(100)) %>%
  make_ipm(usr_funs   = my_funs,
           kernel_seq = sample(1:5, 100, replace = TRUE),
           iterate    = TRUE,
           iterations = 100)


# The default for stochastic lambda is to compute mean(log(ipm$pop_state$lambda)).
# It also removes the first 10% of iterations to handle "burn in". The amount
# removed can be adjusted with the burn_in parameter.

stoch_lambda <- lambda(my_ipm)

# lambda(comp_method = 'pop_size', type = 'all') will compute the population
# growth rate for every time step as the sum(n_ht_t_1) / sum(n_ht_t).

iteration_lambdas <- lambda(my_ipm, type_lambda = 'all')




library(ipmr)

# Define the fixed parameters in a list. The temperature and precipitation
# coefficients are defined as s_temp, s_precip, g_temp, and g_precip.

constant_params <- list(
  s_int     = -10,
  s_slope   = 1.5,
  s_precip  = 0.00001,
  s_temp    = -0.003,
  g_int     = 0.2,
  g_slope   = 1.01,
  g_sd      = 1.2,
  g_temp    = -0.002,
  g_precip  = 0.004,
  r_r_int   = -3.2,
  r_r_slope = 0.55,
  r_s_int   = -0.4,
  r_s_slope = 0.5,
  r_d_mu    = 1.1,
  r_d_sd    = 0.1
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
  # as a named list. We can reference the names in this list in vital rate
  # expressions.

  temp_now   <- rnorm(1,
                      env_params$temp_mu,
                      env_params$temp_sd)

  precip_now <- rgamma(1,
                       shape = env_params$precip_shape,
                       rate  = env_params$precip_rate)

  out        <- list(temp = temp_now, precip = precip_now)

  return(out)

}


# Again, we can define our own functions and pass them into calls to make_ipm. This
# isn't strictly necessary, but can make the model code more readable/less error prone.

inv_logit <- function(lin_term) {
  1/(1 + exp(-lin_term))
}




init_pop_vec <- runif(100)

param_resamp_model <- init_ipm(sim_gen    = "simple",
                               di_dd      = "di",
                               det_stoch  = "stoch",
                               kern_param = "param") %>%
  define_kernel(
    name    = 'P',
    formula = s * g,
    family  = 'CC',

    # Parameters created by define_env_state() can be referenced by name just like
    # any other parameter in the model.

    g_mu    = g_int + g_slope * surf_area_1 + g_temp * temp + g_precip * precip,
    s_lin_p = s_int + s_slope * surf_area_1 + s_temp * temp + s_precip * precip,
    s       = inv_logit(s_lin_p),
    g       = dnorm(surf_area_2, g_mu, g_sd),


    data_list = constant_params,
    states    = list(c('surf_area')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "g")
  ) %>%
  define_kernel(
    name          = 'F',
    formula       = r_r * r_s * r_d,
    family        = 'CC',
    r_r_lin_p     = r_r_int + r_r_slope * surf_area_1,
    r_r           = inv_logit(r_r_lin_p),
    r_s           = exp(r_s_int + r_s_slope * surf_area_1),
    r_d           = dnorm(surf_area_2, r_d_mu, r_d_sd),
    data_list     = constant_params,
    states        = list(c('surf_area')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F"),
      int_rule     = rep('midpoint', 2),
      state_start    = rep('surf_area', 2),
      state_end      = rep('surf_area', 2)
    )
  ) %>%
  define_domains(surf_area = c(0, 10, 100))



param_resamp_model <- param_resamp_model %>%

  define_env_state(
    env_covs   = sample_env(env_params),
    data_list  = list(env_params = env_params,
                      sample_env = sample_env)
  ) %>%
  define_pop_state(
    pop_vectors = list(
      n_surf_area = init_pop_vec
    )
  ) %>%
  make_ipm(usr_funs   = list(inv_logit  = inv_logit),
           iterate    = TRUE,
           iterations = 10)

# By default, lambda computes stochastic lambda with stochastic models

lambda(param_resamp_model)

# We can get lambdas for each time step by adding type_lambda = "all"

lambda(param_resamp_model, type_lambda = 'all')

# If we want to see the actual draws that were used at each step of the
# model iteration, we can access these using the output's $env_seq slot.

param_resamp_model$env_seq



param_resamp_model <- init_ipm(sim_gen    = "simple",
                               di_dd      = "di",
                               det_stoch  = "stoch",
                               kern_param = "param") %>%
  define_kernel(
    name    = 'P',
    formula = s * g,
    family  = 'CC',

    # Parameters created by define_env_state() can be referenced by name just like
    # any other parameter in the model.

    g_mu    = g_int + g_slope * surf_area_1 + g_temp * temp + g_precip * precip,
    s_lin_p = s_int + s_slope * surf_area_1 + s_temp * temp + s_precip * precip,
    s       = inv_logit(s_lin_p),
    g       = dnorm(surf_area_2, g_mu, g_sd),


    data_list = constant_params,
    states    = list(c('surf_area')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "g")
  ) %>%
  define_kernel(
    name          = 'F',
    formula       = r_r * r_s * r_d,
    family        = 'CC',
    r_r_lin_p     = r_r_int + r_r_slope * surf_area_1,
    r_r           = inv_logit(r_r_lin_p),
    r_s           = exp(r_s_int + r_s_slope * surf_area_1),
    r_d           = dnorm(surf_area_2, r_d_mu, r_d_sd),
    data_list     = constant_params,
    states        = list(c('surf_area')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F"),
      int_rule     = rep('midpoint', 2),
      state_start    = rep('surf_area', 2),
      state_end      = rep('surf_area', 2)
    )
  ) %>%
  define_domains(surf_area = c(0, 10, 100))





env_states <- data.frame(precip = rgamma(10, shape = 1000, rate = 2),
                         temp   = rnorm(10, mean = 8.9, sd = 1.2))



sample_env <- function(env_states, iteration) {

  out <- as.list(env_states[iteration, ])
  names(out) <- names(env_states)

  return(out)

}



param_resamp_model <- param_resamp_model %>%
  define_env_state(env_params = sample_env(env_states,
                                           iteration = t), # "t" indexes the current model iteration
                   data_list = list(
                     env_states = env_states,
                     sample_env = sample_env
                   )) %>%
  define_pop_state(
    pop_vectors = list(
      n_surf_area = init_pop_vec
    )
  ) %>%
  make_ipm(usr_funs   = list(inv_logit  = inv_logit),
           iterate    = TRUE,
           iterations = 10)



param_resamp_model <- init_ipm(sim_gen    = "simple",
                               di_dd      = "di",
                               det_stoch  = "stoch",
                               kern_param = "param") %>%
  define_kernel(
    name    = 'P',
    formula = s * g,
    family  = 'CC',

    # Parameters created by define_env_state() can be referenced by name just like
    # any other parameter in the model.

    g_mu    = g_int + g_slope * surf_area_1 + g_temp * temp + g_precip * precip,
    s_lin_p = s_int + s_slope * surf_area_1 + s_temp * temp + s_precip * precip,
    s       = inv_logit(s_lin_p),
    g       = dnorm(surf_area_2, g_mu, g_sd),


    data_list = constant_params,
    states    = list(c('surf_area')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "g")
  ) %>%
  define_kernel(
    name          = 'F',
    formula       = r_r * r_s * r_d,
    family        = 'CC',
    r_r_lin_p     = r_r_int + r_r_slope * surf_area_1,
    r_r           = inv_logit(r_r_lin_p),
    r_s           = exp(r_s_int + r_s_slope * surf_area_1),
    r_d           = dnorm(surf_area_2, r_d_mu, r_d_sd),
    data_list     = constant_params,
    states        = list(c('surf_area')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F"),
      int_rule     = rep('midpoint', 2),
      state_start    = rep('surf_area', 2),
      state_end      = rep('surf_area', 2)
    )
  ) %>%
  define_domains(surf_area = c(0, 10, 100))





env_sampler <- function(time, init_temp = 10, init_precip = 1500) {

  t_est <- init_temp + 0.2 * time + rnorm(1, mean = 0, sd = 0.2)

  p_est <- init_precip + -3.3 * time + rnorm(1, mean = 0, sd = 100)

  if(p_est <= 0) p_est <- 1

  out <- list(temp   = t_est,
              precip = p_est)

  return(out)

}




param_resamp_ipm <- param_resamp_model %>%
  define_env_state(
    env_params = env_sampler(time = t, init_temp = 10, init_precip = 1500),
    data_list  = list()
  ) %>%
  define_pop_state(
    n_surf_area = init_pop_vec
  ) %>%
  make_ipm(
    iterations = 100,
    usr_funs = list(env_sampler = env_sampler,
                    inv_logit = inv_logit)
  )

lambda(param_resamp_ipm)



library(ipmr)

data(iceplant_ex)

grow_mod <- lm(log_size_next ~ log_size, data = iceplant_ex)
grow_sd  <- sd(resid(grow_mod))

surv_mod <- glm(survival ~ log_size, data = iceplant_ex, family = binomial())
repr_mod <- glm(repro ~ log_size, data = iceplant_ex, family = binomial())
seed_mod <- glm(flower_n ~ log_size, data = iceplant_ex, family = poisson())
recr_data <- subset(iceplant_ex, is.na(log_size))

recr_mu  <- mean(recr_data$log_size_next)
recr_sd  <- sd(recr_data$log_size_next)
recr_n   <- length(recr_data$log_size_next)
flow_n   <- sum(iceplant_ex$flower_n, na.rm = TRUE)

params <- list(recr_mu = recr_mu,
               recr_sd = recr_sd,
               grow_sd = grow_sd,
               surv_mod = surv_mod,
               grow_mod = grow_mod,
               repr_mod = repr_mod,
               seed_mod = seed_mod,
               recr_n   = recr_n,
               flow_n   = flow_n)

L <- min(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2
U <- max(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2

obs_ipm <- init_ipm(sim_gen    = "simple",
                    di_dd      = "di",
                    det_stoch  = "det") %>%
  define_kernel(
    name          = "P",
    family        = "CC",
    formula       = s * g,

    # Instead of the inverse logit transformation, we use predict() here.
    # We have to be sure that the "newdata" argument of predict is correctly specified.
    # This means matching the names used in the model itself (log_size) to the names
    # we give the domains. In this case, we use "sa" (short for surface area) as
    # the new data to generate predictions for.

    s             = predict(surv_mod,
                            newdata = data.frame(log_size = sa_1),
                            type    = 'response'),
    g_mu          = predict(grow_mod,
                            newdata = data.frame(log_size = sa_1),
                            type    = 'response'),

    # We specify the rest of the kernel the same way.

    g             = dnorm(sa_2, g_mu, grow_sd),
    states        = list(c("sa")),
    data_list     = params,
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun   = "norm",
                                            target =  "g")
  ) %>%
  define_kernel(
    name          = "F",
    family        = "CC",
    formula       = r_p * r_s * r_d * r_r,

    # As above, we use predict(model_object). We make sure the names of the "newdata"
    # match the names in the vital rate model formulas, and the values match the
    # names of domains they use.

    r_p           = predict(repr_mod,
                            newdata = data.frame(log_size = sa_1),
                            type    = 'response'),
    r_s           = predict(seed_mod,
                            newdata = data.frame(log_size = sa_1),
                            type    = 'response'),
    r_r           = recr_n / flow_n,
    r_d           = dnorm(sa_2, recr_mu, recr_sd),
    states        = list(c("sa")),
    data_list     = params,
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun   = "norm",
                                            target =  "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F"),
      int_rule     = rep('midpoint', 2),
      state_start    = rep("sa", 2),
      state_end      = rep("sa", 2)
    )
  ) %>%
  define_domains(
    sa = c(L,
           U,
           100)
  ) %>%
  define_pop_state(n_sa = runif(100)) %>%
  make_ipm(iterate = TRUE,
           iterations = 100)

lambda_obs <- lambda(obs_ipm)



all_lambdas <- numeric(50L)

recr_data <- subset(iceplant_ex, is.na(log_size))

adults <- subset(iceplant_ex, !is.na(log_size))

use_proto <- obs_ipm$proto_ipm




for(i in 1:50) {

  sample_ind <- seq(1, nrow(adults), by = 1)

  boot_ind   <- sample(sample_ind, size = nrow(adults), replace = TRUE)

  boot_data  <- rbind(adults[boot_ind, ],
                      recr_data)

  grow_mod <- lm(log_size_next ~ log_size, data = boot_data)
  grow_sd  <- sd(resid(grow_mod))

  surv_mod <- glm(survival ~ log_size, data = boot_data, family = binomial())
  repr_mod <- glm(repro ~ log_size, data = boot_data, family = binomial())
  seed_mod <- glm(flower_n ~ log_size, data = boot_data, family = poisson())

  recr_mu  <- mean(recr_data$log_size_next)
  recr_sd  <- sd(recr_data$log_size_next)
  recr_n   <- length(recr_data$log_size_next)
  flow_n   <- sum(boot_data$flower_n, na.rm = TRUE)

  params <- list(recr_mu = recr_mu,
                 recr_sd = recr_sd,
                 grow_sd = grow_sd,
                 surv_mod = surv_mod,
                 grow_mod = grow_mod,
                 repr_mod = repr_mod,
                 seed_mod = seed_mod,
                 recr_n   = recr_n,
                 flow_n   = flow_n)

  # Insert the new vital rate models into the proto_ipm, then rebuild the IPM.

  parameters(use_proto) <- params

  boot_ipm <- use_proto %>%
    make_ipm(iterate = TRUE,
             iterations = 100)

  all_lambdas[i] <- lambda(boot_ipm)

}

# Plot the results

hist(all_lambdas)
abline(v = lambda_obs, col = 'red', lwd = 2, lty = 2)


