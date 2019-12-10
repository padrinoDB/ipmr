# Test helper functions
library(rlang)

compare_kernels <- function(ipmr_object,
                            hand_object,
                            nms_ipmr,
                            nms_hand,
                            .env = rlang::caller_env()) {

  c1 <- 'all.equal(unclass('

  c2 <- '), '

  c3 <- ')'

  ipmr_object <- paste(ipmr_object, '$sub_kernels$', nms_ipmr, sep = '')
  hand_object <- paste(hand_object, '$', nms_hand, sep = "")

  temp <- paste(c1, ipmr_object, c2, hand_object, c3, sep = "")

  out <- lapply(temp, function(x, env) rlang::new_quosure(
    rlang::parse_expr(x), env = env),
    env = .env
  )

  lapply(out, rlang::eval_tidy) %>%
    unlist()


}
