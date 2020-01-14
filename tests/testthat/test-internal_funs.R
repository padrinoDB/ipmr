context('test internal functions')


test_mat   <- matrix(c(1, 2, 3,
                       2, 4, 6,
                       3, 6, 9),
                     nrow = 3, byrow = TRUE)

test_mat_2 <- matrix(c(test_mat, test_mat, test_mat),
                   nrow = 3)

test_that('conv_to_asymptotic and is_square work/fail properly', {

  expect_true(is_square(test_mat))
  expect_true(is_conv_to_asymptotic(test_mat_2))


  expect_false(is_square(test_mat_2))


})
