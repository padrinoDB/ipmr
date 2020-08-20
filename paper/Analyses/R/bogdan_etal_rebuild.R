# Building a simple IPM using ipmr's internal data set. These data are
# from a Caprobrotus spp. This genus of spreading succulents is native to
# South Africa, but has invaded Mediterranean-like regions all over the world.
# Data were collected in April 2018 and 2019 using drones to take aerial images
# of plants at low altitude. Plants were marked with polygons and flowers were
# counted at each sampling.

library(ipmr)

data(iceplant_ex)

# growth model

grow_mod <- lm(log_size_next ~ log_size, data = iceplant_ex)
grow_sd  <- sd(resid(grow_mod))

# survival model

surv_mod <- glm(survival ~ log_size, data = iceplant_ex, family = binomial())

# Pr(flowering) model

repr_mod <- glm(repro ~ log_size, data = iceplant_ex, family = binomial())

# Number of flowers per plant model

flow_mod <- glm(flower_n ~ log_size, data = iceplant_ex, family = poisson())

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

# Now, we put all parameters into a list. This case study shows how to use
# the mathematical notation, as well as how to use predict() methods

params <- list(
  surv_int = coef(surv_mod)[1],
  surv_slo = coef(surv_mod)[2],
  grow_int = coef(grow_mod)[1],
  grow_slo = coef(grow_mod)[2],
  grow_sdv = sd(resid(grow_mod)),
  repr_int = coef(repr_mod)[1],
  repr_slo = coef(repr_mod)[2],
  flow_int = coef(flow_mod)[1],
  flow_slo = coef(flow_mod)[2],
  recr_n   = recr_n,
  flow_n   = flow_n,
  recr_mu  = recr_mu,
  recr_sd  = recr_sd
)

# The lower and upper bounds for integration. Adding 20% on either end to minimize
# eviction

L <- min(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 0.8
U <- max(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2


math_ipm <- init_ipm("simple_di_det") %>%
  define_kernel(
    name          = "P",
    family        = "CC",
    formula       = s * g,

    # plogis is the inverse logit transformation. We need this to convert the
    # linear predictor back into a survival probability.

    s             = plogis(surv_int + surv_slo * sa_1),

    # The growth model is a simple linear model, which predicts the mean.
    # We use the residual variance to compute the actual transition probabilities
    # by integrating over the probability density function defined by the predicted
    # mean and the variance of the linear model.

    g_mu          = grow_int + grow_slo * sa_1,

    g             = dnorm(sa_2, g_mu, grow_sdv),
    states        = list(c("sa")),
    data_list     = params,
    has_hier_effs = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun   = "norm",
                                            param = "g")
  ) %>%
  define_kernel(
    name          = "F",
    family        = "CC",
    formula       = f_p * f_n * f_d * p_r,

    f_p           = plogis(repr_int + repr_slo * sa_1),
    f_n           = exp(flow_int + flow_slo * sa_1),
    p_r           = recr_n / flow_n,
    f_d           = dnorm(sa_2, recr_mu, recr_sd),
    states        = list(c("sa")),
    data_list     = params,
    has_hier_effs = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun   = "norm",
                                            param = "f_d")
  ) %>%
  define_k(
    name          = "K",
    family        = "IPM",
    K             = P + F,
    states        = list(c("sa")),
    data_list     = list(),
    has_hier_effs = FALSE,
    evict_cor     = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F", "K"),
      int_rule     = rep('midpoint', 3),
      dom_start    = rep("sa", 3),
      dom_end      = rep("sa", 3)
    )
  ) %>%
  define_domains(
    sa = c(L,
           U,
           100)
  ) %>%
  make_ipm(iterate = FALSE)

lambda_math <- lambda(math_ipm, comp_method = 'eigen')


# Now, an alternative implementation that uses predict() methods for each
# vital rate model. These models run a bit slower because of the internal machinery
# in each predict method. However, they can save a quite a bit of time when
# vital rate models are more complicated, and you don't feel like working out
# the exact functional form.

models <- list(recr_mu = recr_mu,
               recr_sd = recr_sd,
               grow_sd = grow_sd,
               surv_mod = surv_mod,
               grow_mod = grow_mod,
               repr_mod = repr_mod,
               flow_mod = flow_mod,
               recr_n   = recr_n,
               flow_n   = flow_n)

pred_ipm <- init_ipm("simple_di_det") %>%
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
    data_list     = models,
    has_hier_effs = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun   = "norm",
                                            param = "g")
  ) %>%
  define_kernel(
    name          = "F",
    family        = "CC",
    formula       = f_p * f_n * f_d * p_r,

    # As above, we use predict(model_object). We make sure the names of the "newdata"
    # match the names in the vital rate model formulas, and the values match the
    # names of domains they use.

    f_p           = predict(repr_mod,
                            newdata = data.frame(log_size = sa_1),
                            type    = 'response'),
    f_n           = predict(flow_mod,
                            newdata = data.frame(log_size = sa_1),
                            type    = 'response'),
    p_r           = recr_n / flow_n,
    f_d           = dnorm(sa_2, recr_mu, recr_sd),
    states        = list(c("sa")),
    data_list     = models,
    has_hier_effs = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun   = "norm",
                                            param = "f_d")
  ) %>%
  define_k(
    name          = "K",
    family        = "IPM",
    K             = P + F,
    states        = list(c("sa")),
    data_list     = list(),
    has_hier_effs = FALSE,
    evict_cor     = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F", "K"),
      int_rule     = rep('midpoint', 3),
      dom_start    = rep("sa", 3),
      dom_end      = rep("sa", 3)
    )
  ) %>%
  define_domains(
    sa = c(L,
           U,
           100)
  ) %>%
  make_ipm(iterate = FALSE)

lambda_pred <- lambda(pred_ipm, comp_method = 'eigen')

# Brief sanity check - should return nothing

stopifnot(all.equal(lambda_math, lambda_pred, tolerance = 1e-13))
