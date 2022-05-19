test_that("make_report_fp produces expected paths", {

  skip_on_cran()

  dir <- tempdir()

  dt <- gsub("-", "", Sys.Date())

  fp <- paste0(dir, "/ipmr_report_", dt, ".Rmd")

  test_fp <- .make_ipm_report_fp(dir, TRUE)
  expect_equal(fp, test_fp)

  # Hard to test the tempfile() output path - return to this

  expect_message(test_fp <-  .make_ipm_report_fp("", TRUE),
                 regex = "'keep_rmd = TRUE' and 'rmd_dest' is not specified!")

  test_fp <- .make_ipm_report_fp(paste0(dir, "/ipm_report.rmd"), TRUE)
  fp      <- paste0(dir, "/ipm_report_", dt, ".Rmd")
  expect_equal(fp, test_fp)

  unlink(dir, TRUE, TRUE)

})


data("sim_di_det_ex")
data("gen_di_det_ex")

test_that("make_rmd_header produces correct headers", {

  # Test defaults for html
  hdr <- .make_ipm_report_header("html", "test_header", FALSE)

  regs <- c("title: 'test_header'",
            "output:\\n  html_document:",
            "toc: true\\n    toc_depth: 3")

  tst <- vapply(regs, function(x, hdr) grepl(x, hdr), logical(1L), hdr = hdr)

  expect_true(all(tst))

  expect_message(.make_ipm_report_header("word", "test_header", TRUE),
                 regex = "Block equation numbering")

  hdr <- .make_ipm_report_header("pdf", "test_header", TRUE)

  regs <- c(regs[c(1, 3)],
            "output:\\n  pdf_document:",
            "header_includes:",
            "usepackage", "amsmath")

  tst <- vapply(regs, function(x, hdr) grepl(x, hdr), logical(1L), hdr = hdr)

  expect_true(all(tst))

})


test_that("internal helpers translate correctly", {

  proto_ipm <- sim_di_det_ex$proto_ipm
  unk_env   <- .unk_env(proto_ipm)
  greek_env <- .greek_env(unk_env)
  pop_env   <- .pop_env(proto_ipm, greek_env)
  math_env  <- .math_env(pop_env)
  par_env   <- .par_env(proto_ipm, math_env)


  # Don't care about testing headers for now
  iter_exprs <- .make_ipm_report_iter_exprs(proto_ipm,
                                            par_env,
                                            block_eqs = FALSE,
                                            rmd_dest = tempdir()) %>%
    .[!grepl("#", .)]

  # Test RHS and LHS of each non-integral expressions (integral expressions)
  # are a mess.
  expect_true(grepl("s \\* g", iter_exprs[2]))
  expect_true(grepl("f\\\\_r \\* f\\\\_s \\* f\\\\_d", iter_exprs[3]))
  expect_true(grepl("P\\(dbh', dbh\\)", iter_exprs[2]))
  expect_true(grepl("F\\(dbh', dbh\\)", iter_exprs[3]))
  expect_equal(length(iter_exprs), 3L)

  vr_exprs   <- .make_ipm_report_vr_exprs(proto_ipm,
                                          par_env,
                                          block_eqs = FALSE) %>%
    .[!grepl("#", .)]

  expect_true(grepl("mathrm", vr_exprs[1]))
  expect_true(grepl("Norm\\(dbh'", vr_exprs[2]))
  expect_true(grepl("g\\\\_int \\+ g\\\\_slope \\* dbh", vr_exprs[3]))
  expect_true(grepl("mathrm", vr_exprs[4]))
  expect_true(grepl("inv\\\\_logit", vr_exprs[4]))
  expect_true(grepl("e\\^", vr_exprs[5]))
  expect_true(grepl("Norm\\(dbh'", vr_exprs[6]))
  expect_equal(length(vr_exprs), 6L)


  all_params <- .make_ipm_report_params(proto_ipm, par_env) %>%
    .[!grepl("#", .)]

  expect_equal(length(all_params), 2L)

})
