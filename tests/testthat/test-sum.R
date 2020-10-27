context("test the behavior of sum() *_det vr_exprs")

build_eval_exprs <- function(proto, kernel, dd = FALSE) {

  ind <- which(proto$kernel_id == kernel)

  proto <- .initialize_kernels.default(proto, TRUE)
  proto <- proto$others[ind, ]

  if(dd) proto <- .prep_dd_vr_exprs(proto)

  main_env <- .make_main_env(proto$domain, usr_funs = list())
  kern_env <- rlang::child_env(.parent = main_env)
  rlang::env_bind(kern_env,
                  !!! proto$params[[1]]$params)

  kern_quos <- .parse_vr_formulae(proto$params[[1]]$vr_text,
                                  kern_env,
                                  proto,
                                  main_env)

  return(kern_quos)
}

# IPMR constants

data_list <- list(
  g_int        = 0.83,
  g_slope      = 0.69,
  g_sd         = sqrt(0.19),
  s_int        = -0.62,
  s_slope      = 0.85,
  p_b_int      = -10.22,
  p_b_slope    = 4.25,
  weevil_slope = 1.71,
  f_s_int      = -0.55,
  f_s_slope    = 2.02,
  f_s_exp      = 0.32,
  recr_mu      = 0.75,
  recr_sd      = sqrt(0.17)
)

int <- numeric(5L)
int[1] <- mean (c(-17.3,-17,-16.7))
int[2] <- mean (c(-3.74,-3.39,-2.81))
int[3] <- mean(c(-1.64,-0.96,-0.75))
int[4] <- mean(c(-1.24,-0.52))
int[5] <- 1.5

weevil_ints <- as.list(int) %>%
  setNames(
    paste("weevil_int_", 1:5, sep = "")
  )


data_list <- c(data_list, weevil_ints)

weevil_ipm <- init_ipm("simple_dd_det") %>%
  define_kernel(
    name      = "P",
    formula   = s * (1-p_b) * g,
    family    = "CC",
    s         = plogis(s_int + s_slope * z_1),
    g_mu      = g_int + g_slope * z_1,
    g         = dnorm(z_2, g_mu, g_sd),
    p_b       = plogis(p_b_int + p_b_slope * z_1),
    data_list = data_list,
    states    = list(c("z")),
    evict_cor = FALSE,
    integrate = TRUE
  ) %>%
  define_kernel(
    name       = "F",
    formula    = c_0 * N,
    family     = "CC",
    c_0        = dnorm(z_2, recr_mu, recr_sd),
    N          = (d_z * sum(n_z_t * B_z)) ^ 0.67,
    B_z        = p_b * b_z,
    p_b        = plogis(p_b_int + p_b_slope * z_1),
    b_z        = exp(f_s_int + f_s_slope * z_1) / ((1 + e_z/16) ^ (f_s_exp)),
    b_d        = f_s_int + f_s_slope * z_1,
    b_2        = sum(n_z_t) * b_d,
    e_z        = exp(weevil_int + weevil_slope * z_1),
    weevil_int = weevil_int_1,
    data_list  = data_list,
    states     = list(c("z")),
    evict_cor  = FALSE,
    integrate  = FALSE
  ) %>%
  define_k(
    name      = "K",
    family    = "IPM",
    F_1       = F[ , 1],
    n_z_t_1   = P %*% n_z_t + F_1,
    data_list = data_list,
    states    = list(c("z")),
    evict_cor = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      c("P", "F", "K"),
      rep("midpoint",3),
      rep("z",3),
      rep("z", 3)
    )
  ) %>%
  define_domains(
    z = c(0, 10, 100)
  ) %>%
  define_pop_state(
    pop_vectors = list(
      n_z = runif(100)
    )
  )

kern_quos <- build_eval_exprs(weevil_ipm, "F", TRUE)

N_ipmr_text <- quo_text(kern_quos$N)
N_target_text <- "(d_z * sum(as.vector(pop_state_z_t) * B_z/n_z_1))^0.67"

b_2_ipmr_text <- rlang::quo_text(kern_quos$b_2)
b_2_target_text <- "sum(as.vector(pop_state_z_t)) * b_d"

test_that("sum is translated correctly", {

  expect_equal(N_ipmr_text, N_target_text)
  expect_equal(b_2_ipmr_text, b_2_target_text)

})

# Density-independent expressions

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


n <- 100
L <- min(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2
U <- max(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2


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
    f_r           = recr_n / sum(f_s),
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
  define_pop_state(n_sa = runif(100))

f_r_target <- "recr_n/sum(f_s/n_sa_1)"
kern_quos <- build_eval_exprs(pred_ipm, "F")
f_r_ipmr  <- quo_text(kern_quos$f_r)

test_that("sum works in di IPMs", {

  expect_equal(f_r_target, f_r_ipmr)

})


# context("dd_stoch_kern sums")



# Make up some really gnarly vr_exprs to test this

# context("highly nested sum expression")
