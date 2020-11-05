
data(iceplant_ex)

grow_mod <- lm(log_size_next ~ log_size, data = iceplant_ex)
surv_mod <- glm(survival ~ log_size, data = iceplant_ex, family = binomial())
repr_mod <- glm(repro ~ log_size, data = iceplant_ex, family = binomial())
seed_mod <- glm(flower_n ~ log_size, data = iceplant_ex, family = poisson())

recr_mu  <- mean(iceplant_ex$log_size_next[is.na(iceplant_ex$log_size)])
recr_sd  <- sd(iceplant_ex$log_size_next[is.na(iceplant_ex$log_size)])
recr_n   <- length(iceplant_ex$log_size_next[is.na(iceplant_ex$log_size)])
flow_n   <- sum(iceplant_ex$flower_n, na.rm = TRUE)

grow_sd <- sd(resid(grow_mod))

params <- list(recr_mu = recr_mu,
               recr_sd = recr_sd,
               grow_sd = grow_sd,
               surv_mod = surv_mod,
               grow_mod = grow_mod,
               repr_mod = repr_mod,
               seed_mod = seed_mod,
               recr_n   = recr_n,
               flow_n   = flow_n)


g_z1z <- function(grow_mod, d1, d2, L, U) {

  mu <- predict(grow_mod, data.frame(log_size = d1), type = 'response')
  g_sd <- sd(resid(grow_mod))

  ev <- pnorm(U, mu, g_sd) - pnorm(L, mu, g_sd)

  out <- dnorm(d2, mu, g_sd) / ev

  return(out)

}

s_z <- function(surv_mod, d1) {

  predict(surv_mod, data.frame(log_size = d1), type = 'response')

}

f_z1z <- function(f_mod, s_mod, mu, r_sd, r_n, f_n, d1, d2, L, U) {

  ev <- pnorm(U, mu, r_sd) - pnorm(L, mu, r_sd)

  s_tot <- predict(s_mod, data.frame(log_size = d1), type = 'response')
  f_r   <- r_n / f_n

  p_r   <- predict(f_mod, data.frame(log_size = d1), type = 'response')


  f_d <- dnorm(d2, mu, r_sd) / ev

  out <- s_tot * f_r * p_r * f_d

  return(out)


}

n <- 100
L <- min(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2
U <- max(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2

d <- seq(L, U, length.out = (n + 1))

sv <- (d[2:101] + d[1:100]) * 0.5

doms <- expand.grid(d1 = sv, d2 = sv)

h <- sv[2] - sv[1]


g_kern <- g_z1z(grow_mod, doms$d1, doms$d2, L, U)
s_kern <- s_z(surv_mod, doms$d1)
f_kern <- f_z1z(repr_mod, seed_mod, recr_mu, recr_sd, recr_n, flow_n,
                doms$d1, doms$d2, L, U) * h

p_kern <- g_kern * s_kern * h

p_kern <- matrix(p_kern, nrow = 100, ncol = 100, byrow = TRUE)
f_kern <- matrix(f_kern, nrow = 100, ncol = 100, byrow = TRUE)

k_kern <- p_kern + f_kern

lambda_hand <- Re(eigen(k_kern)$values[1])


pred_ipm <- init_ipm("simple_di_det") %>%
  define_kernel(
    name          = "P",
    family        = "CC",
    formula       = s * g,
    s             = predict(surv_mod, data.frame(log_size = sa_1), type = 'response'),
    g_mu          = predict(grow_mod, data.frame(log_size = sa_1), type = 'response'),
    g             = dnorm(sa_2, g_mu, grow_sd),
    states        = list(c("sa")),
    data_list     = params,
    has_hier_effs = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "g")
  ) %>%
  define_kernel(
    name          = "F",
    family        = "CC",
    formula       = f_p * f_s * f_d * f_r,
    f_p           = predict(repr_mod, data.frame(log_size = sa_1), type = 'response'),
    f_s           = predict(seed_mod, data.frame(log_size = sa_1), type = 'response'),
    f_r           = recr_n / flow_n,
    f_d           = dnorm(sa_2, recr_mu, recr_sd),
    states        = list(c("sa")),
    data_list     = params,
    has_hier_effs = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "f_d")
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
    sa = c(min(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2,
           max(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2,
           100)
  ) %>%
  make_ipm(iterate = FALSE)

lambda_ipmr <- lambda(pred_ipm, comp_method = 'eigen') %>%
  as.vector()

test_that("Kernels w/ predict() are the same as hand implemented ones", {

  expect_equal(lambda_hand, lambda_ipmr, tolerance = 1e-10)

})

protected_params <- list(recr_mu = recr_mu,
                         recr_sd = recr_sd,
                         grow_sd = grow_sd,
                         surv_mod = use_vr_model(surv_mod),
                         grow_mod = use_vr_model(grow_mod),
                         repr_mod = use_vr_model(repr_mod),
                         seed_mod = use_vr_model(seed_mod),
                         recr_n   = recr_n,
                         flow_n   = flow_n)



prot_ipm <- init_ipm("simple_di_det") %>%
  define_kernel(
    name          = "P",
    family        = "CC",
    formula       = s * g,
    s             = predict(surv_mod, data.frame(log_size = sa_1), type = 'response'),
    g_mu          = predict(grow_mod, data.frame(log_size = sa_1), type = 'response'),
    g             = dnorm(sa_2, g_mu, grow_sd),
    states        = list(c("sa")),
    data_list     = params,
    has_hier_effs = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "g")
  ) %>%
  define_kernel(
    name          = "F",
    family        = "CC",
    formula       = f_p * f_s * f_d * f_r,
    f_p           = predict(repr_mod, data.frame(log_size = sa_1), type = 'response'),
    f_s           = predict(seed_mod, data.frame(log_size = sa_1), type = 'response'),
    f_r           = recr_n / flow_n,
    f_d           = dnorm(sa_2, recr_mu, recr_sd),
    states        = list(c("sa")),
    data_list     = params,
    has_hier_effs = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "f_d")
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
    sa = c(min(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2,
           max(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2,
           100)
  ) %>%
  make_ipm(iterate = FALSE)

lambda_protected <- lambda(prot_ipm, comp_method = 'eigen') %>%
  as.vector()

test_that("use_vr_model works as expected", {

  expect_true(attr(protected_params$grow_mod, "flat_protect"))
  expect_true(attr(protected_params$surv_mod, "flat_protect"))
  expect_null(attr(protected_params$grow_sd, "flat_protect"))

  expect_equal(lambda_protected, lambda_ipmr, tolerance = 1e-10)
})





f_z1z <- function(f_mod, s_mod, mu, r_sd, r_n, d1, d2, L, U) {

  ev <- pnorm(U, mu, r_sd) - pnorm(L, mu, r_sd)

  s_tot <- predict(s_mod, data.frame(log_size = d1), type = 'response')
  f_n   <- sum(predict(s_mod, data.frame(log_size = unique(d1)), type = 'response'))
  f_r   <- r_n / f_n

  p_r   <- predict(f_mod, data.frame(log_size = d1), type = 'response')


  f_d <- dnorm(d2, mu, r_sd) / ev

  out <- s_tot * f_r * p_r * f_d

  return(out)


}

g_kern <- g_z1z(grow_mod, doms$d1, doms$d2, L, U)
s_kern <- s_z(surv_mod, doms$d1)
f_kern <- f_z1z(repr_mod, seed_mod, recr_mu, recr_sd, recr_n,
                doms$d1, doms$d2, L, U) * h

p_kern <- g_kern * s_kern * h

p_kern <- matrix(p_kern, nrow = 100, ncol = 100, byrow = TRUE)
f_kern <- matrix(f_kern, nrow = 100, ncol = 100, byrow = TRUE)

k_kern <- p_kern + f_kern

lambda_hand <- Re(eigen(k_kern)$values[1])

# params$flow_n <- sum(predict(seed_mod, data.frame(log_size = unique(doms$d1)), type = 'response'))
params$flow_n <- NULL

sum_ipm <- init_ipm("simple_di_det") %>%
  define_kernel(
    name          = "P",
    family        = "CC",
    formula       = s * g,
    s             = predict(surv_mod, data.frame(log_size = sa_1), type = 'response'),
    g_mu          = predict(grow_mod, data.frame(log_size = sa_1), type = 'response'),
    g             = dnorm(sa_2, g_mu, grow_sd),
    states        = list(c("sa")),
    data_list     = params,
    has_hier_effs = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "g")
  ) %>%
  define_kernel(
    name          = "F",
    family        = "CC",
    formula       = f_p * f_s * f_d * f_r,
    f_p           = predict(repr_mod, data.frame(log_size = sa_1), type = 'response'),
    f_s           = predict(seed_mod, data.frame(log_size = sa_1), type = 'response'),
    f_r           = recr_n / flow_n,
    flow_n        = sum(f_s),
    f_d           = dnorm(sa_2, recr_mu, recr_sd),
    states        = list(c("sa")),
    data_list     = params,
    has_hier_effs = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "f_d")
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
    sa = c(min(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2,
           max(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2,
           100)
  ) %>%
  make_ipm(iterate = FALSE)

lambda_ipmr <- lambda(sum_ipm, comp_method = 'eigen') %>%
  as.vector()

test_that("Kernels w/ predict() are the same as hand implemented ones", {

  expect_equal(lambda_hand, lambda_ipmr, tolerance = 1e-10)

  p_mat <- sum_ipm$sub_kernels$P
  class(p_mat) <- NULL

  f_mat <- sum_ipm$sub_kernels$F
  class(f_mat) <- NULL

  expect_equal(p_mat, p_kern)
  expect_equal(f_mat, f_kern)

  k_mat <- sum_ipm$iterators$K
  class(k_mat) <- NULL

  expect_equal(k_mat, k_kern)

})


