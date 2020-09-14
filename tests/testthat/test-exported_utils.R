context("Exported accessor functions")

test_that("exported utils return expected values w general IPMs", {

  data(gen_di_det_ex)

  prot <- gen_di_det_ex$proto_ipm

  dom <- domains(prot)

  vr_exprs <- vital_rates(prot)

  forms <- kernel_formulae(prot)

  vr_types <- vapply(vr_exprs,
                     rlang::is_call,
                     logical(1L))

  kern_types <- vapply(forms,
                       function(x)
                         rlang::is_call(x) || rlang::is_bare_atomic(x),
                       logical(1L))

  expect_true(all(vr_types))
  expect_true(all(kern_types))

  dom_types <- vapply(dom, is.numeric, logical(1L))

  expect_true(all(dom_types))


})

test_that("exported utils return expected values w simple IPMs", {

  data(sim_di_det_ex)

  prot <- sim_di_det_ex$proto_ipm

  dom <- domains(prot)

  vr_exprs <- vital_rates(prot)

  forms <- kernel_formulae(prot)

  vr_types <- vapply(vr_exprs,
                     rlang::is_call,
                     logical(1L))

  kern_types <- vapply(forms,
                       rlang::is_call,
                       logical(1L))

  expect_true(all(vr_types))
  expect_true(all(kern_types))

  dom_types <- vapply(dom, is.numeric, logical(1L))

  expect_true(all(dom_types))


})

context("int_mesh and other *_ipm helper funs")

# simple_di_det methods ---------
data_list = list(s_int = 2.2,
                 s_slope = 0.25,
                 g_int = 0.2,
                 g_slope = 1.02,
                 sd_g = 0.7,
                 f_r_int = 0.003,
                 f_r_slope = 0.015,
                 f_s_int = 1.3,
                 f_s_slope = 0.075,
                 mu_fd = 2,
                 sd_fd = 0.3)

impl_args <- make_impl_args_list(c('P', 'F', 'K'),
                                 int_rule = rep('midpoint', 3),
                                 dom_start = rep('dbh', 3),
                                 dom_end = rep('dbh', 3))

states <- list(c("dbh"))

sim_di_det_2 <- init_ipm('simple_di_det') %>%
  define_kernel("P",
                formula = s * g,
                family = "CC",
                s = plogis(s_int, s_slope, dbh_1),
                g = dnorm(dbh_2, mu_g, sd_g),
                mu_g = g_int + g_slope * dbh_1,
                data_list = data_list,
                states = states,
                evict_cor = TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'g')
  ) %>%
  define_kernel('F',
                formula = f_r * f_s * f_d,
                family = 'CC',
                f_r = plogis(f_r_int, f_r_slope, dbh_1),
                f_s = exp(f_s_int + f_s_slope * dbh_1),
                f_d = dnorm(dbh_2, mu_fd, sd_fd),
                data_list = data_list,
                states = states,
                evict_cor = TRUE,
                evict_fun = truncated_distributions('norm',
                                                    'f_d')
  ) %>%
  define_k('K',
           K         = P + F,
           n_dbh_t_1 = K %*% n_dbh_t,
           family    = 'IPM',
           data_list = list(),
           states    = states,
           evict_cor = FALSE) %>%
  define_impl(impl_args) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  define_pop_state(n_dbh_t = runif(100)) %>%
  make_ipm(iterate = TRUE,
           iterations = 100,
           normalize_pop_size = FALSE,
           return_all_envs = TRUE,
           return_main_env = TRUE)

test_that("int_mesh works as expected", {

  mesh_ps <- int_mesh(sim_di_det_2)

  expect_equal(names(mesh_ps), c("d_dbh", "dbh_1", "dbh_2"))

  expect_equal(length(mesh_ps$dbh_1), 10000L)
  expect_equal(length(mesh_ps$d_dbh), 1L)

})
