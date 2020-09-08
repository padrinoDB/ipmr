context("Exported accessor functions")

test_that("exported utils return expected values w general IPMs", {

  data(gen_di_det_ex)

  prot <- gen_di_det_ex$proto_ipm

  dom <- domains(prot)

  vr_exprs <- vital_rate_functions(prot)

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

  vr_exprs <- vital_rate_functions(prot)

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
