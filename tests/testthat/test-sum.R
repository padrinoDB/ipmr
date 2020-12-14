
build_eval_exprs <- function(proto, kernel, dd = FALSE) {

  ind <- which(proto$kernel_id == kernel)

  proto <- .initialize_kernels.default(proto, TRUE, "right")
  proto <- proto$others[ind, ]

  if(dd) proto <- .prep_dd_vr_exprs(proto)

  main_env <- .make_main_env(proto$domain, usr_funs = list(), FALSE)
  kern_env <- rlang::child_env(.parent = main_env)
  rlang::env_bind(kern_env,
                  !!! proto$params[[1]]$params)

  kern_quos <- .parse_vr_formulae(proto$params[[1]]$vr_text,
                                  kern_env,
                                  proto,
                                  main_env)

  return(kern_quos)
}



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
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F"),
      int_rule     = rep('midpoint', 2),
      state_start    = rep("sa", 2),
      state_end      = rep("sa", 2)
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

