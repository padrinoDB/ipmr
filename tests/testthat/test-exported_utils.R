
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
  data(gen_di_det_ex)

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

  vr_ipm <- vital_rates(sim_di_det_ex)

  expect_identical(vr_exprs, vr_ipm)

  kern_ipm <- kernel_formulae(sim_di_det_ex)

  expect_identical(kern_ipm, forms)

  dom_proto <- domains(prot)
  dom_ipm   <- domains(sim_di_det_ex)

  expect_identical(dom_proto, dom_ipm)


  # pop_state shouldn't actually return the same thing here,
  # so these expectations need to be adjusted
  prot    <- gen_di_det_ex$proto_ipm

  ps_prot <- pop_state(prot)

  prot_tst <- vapply(ps_prot,
                     function(x) x == "Pre-defined population state.",
                     logical(1L))

  expect_true(all(prot_tst))
  expect_s3_class(ps_prot, "ipmr_pop_state")

  ps_ipm  <- pop_state(gen_di_det_ex)

  expect_true(all(is.array(ps_ipm$n_ht), is.array(ps_ipm$n_b)))
  expect_type(ps_ipm, "list")

})


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

impl_args <- make_impl_args_list(c('P', 'F'),
                                 int_rule = rep('midpoint', 2),
                                 state_start = rep('dbh', 2),
                                 state_end = rep('dbh', 2))

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
  )  %>%
  define_impl(impl_args) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  define_pop_state(n_dbh = runif(100)) %>%
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

test_that("parameters gets and sets correctly", {

  pars <- parameters(sim_di_det_2$proto_ipm)
  expect_s3_class(pars, "ipmr_parameters")

  pars <- unlist(pars)

  expect_equal(unlist(data_list), pars)

  new_pars <- lapply(data_list,
                     function(x) x + rnorm(1, 0, 0.2))

  new_proto <- sim_di_det_2$proto_ipm
  parameters(new_proto) <- new_pars

  expect_equal(unlist(parameters(new_proto)), unlist(new_pars))

  # Now, test subsetted assignment

  newer_pars <- new_pars[1:5]
  newer_pars <- lapply(newer_pars, function(x) x  + rnorm(1, 0, 0.1))

  parameters(new_proto) <- newer_pars

  test_pars <- parameters(new_proto)
  test_pars <- test_pars[names(new_pars)]

  expect_equal(unlist(test_pars[names(newer_pars)]), unlist(newer_pars))
  expect_equal(unlist(test_pars[6:11]), unlist(new_pars)[6:11])


  data(gen_di_det_ex)

  ipm_params <- parameters(gen_di_det_ex)

  proto_params <- parameters(gen_di_det_ex$proto_ipm)

  expect_s3_class(ipm_params, "ipmr_parameters")
  expect_s3_class(proto_params, "ipmr_parameters")

  expect_identical(ipm_params, proto_params)


})

test_that("collapse_pop_state works", {

  data("gen_di_det_ex")

  gen_di_det_ex <- gen_di_det_ex$proto_ipm %>%
    make_ipm(iterate = TRUE,
             iterations = 100,
             return_main_env = TRUE)

  temp <- collapse_pop_state(gen_di_det_ex,
                             time_step = 100,
                             seedlings = ht <= 10,
                             NRA = ht > 10 & ht <= 200,
                             RA = ht > 200) %>%
    unlist() %>%
    round(digits = 6)

  target <- c(seedlings = 0.104561,
              NRA       = 0.148198,
              RA        = 0.001998)
  expect_equal(temp, target)

})
