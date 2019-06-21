# simple-ndd-det IPM test
data_list = list(s_int = 2.2,
                 s_slope = 0.25,
                 g_int = 0.2,
                 g_slope = 1.02,
                 sd_g = 0.7,
                 f_r_int = 0.03,
                 f_r_slope = 0.015,
                 f_s_int = 1.3,
                 f_s_slope = 0.075,
                 mu_fd = 0.5,
                 sd_fd = 0.2)

s <- function(sv1, params) {
  1/(1 + exp(-(params[1] + params[2] * sv1)))
}


g <- function(sv1, sv2, params) {
  mu <- params[1] + params[2] * sv1
  dnorm(sv2, mean = mu, sd = params[3])
}

f_r <- function(sv1, params) {
  1/(1 + exp(-(params[1] + params[2] * sv1)))
}

f_s <- function(sv1, params) {
  exp(params[1] + params[2] * sv1)
}

f_d <- function(sv2, params) {
  dnorm(sv2, mean = params[1], sd = params[2])
}

fec <- function(sv1, sv2, params) {
  f_r(sv1, params[1:2]) * f_s(sv1, params[3:4]) * f_d(sv2, params[5:6])
}

b <- seq(0, 50, length.out = 101)
d1 <- d2 <- (b[2:101] + b[1:100]) * 0.5
h <- d1[3] - d1[2]


G <- h * outer(d1, d2, FUN = g, params = c(data_list$g_int,
                                           data_list$g_slope,
                                           data_list$sd_g))
G2 <- G/matrix(as.vector(apply(G, 2, sum)),
               nrow = length(d1),
               ncol = length(d1),
               byrow = TRUE)

S <- s(d1, c(data_list$s_int, data_list$s_slope))

P <- t(S * t(G2))

Fm <- h * outer(d1, d2, FUN = fec, params = unlist(data_list[6:11]))

K <- P + Fm

lambda <- Re(eigen(K)$values[1])
w <- Re(eigen(K)$vectors[ , 1])


state_list <- list(c('dbh'))

x <- init_ipm('simple_di_det') %>%
  add_kernel("P",
             formula = t(s * t(g)),
             family = "CC",
             s = 1/(1 + exp(-(s_int + s_slope * dbh_1))),
             g = cell_size_dbh * dnorm(dbh_2, mu_g, sd_g),
             mu_g = g_int + g_slope * dbh_1,
             data_list = list(s_int = 2.2,
                              s_slope = 0.25,
                              g_int = 0.2,
                              g_slope = 1.02,
                              sd_g = 0.7),
             state_list = state_list,
             int_rule = 'midpoint',
             dom_start = 'dbh',
             dom_end = 'dbh',
             evict = TRUE,
             evict_fun = truncated_distributions(g,
                                                 n_mesh_p = 100)) %>%
  add_kernel('F',
             formula = f_r * f_s * f_d,
             family = 'CC',
             f_r = 1/(1 + exp(-(f_r_int + f_r_slope * dbh_1))),
             f_s = exp(f_s_int + f_s_slope * dbh_1),
             f_d = cell_size_dbh * dnorm(dbh_2, mu_fd, sd_fd),
             data_list = list(f_r_int = 0.03,
                              f_r_slope = 0.015,
                              f_s_int = 1.3,
                              f_s_slope = 0.075,
                              mu_fd = 0.5,
                              sd_fd = 0.2),
             state_list = state_list,
             int_rule = 'midpoint',
             dom_start = 'dbh',
             dom_end = 'dbh',
             evict = FALSE) %>%#,
  # evict_fun = truncated_distributions(f_d,
  #                                     n_mesh_p = 100)) %>%
  add_kernel('K',
             formula = P + F,
             family = 'IPM',
             data_list = list(),
             state_list = state_list,
             int_rule = 'midpoint',
             dom_start = 'dbh',
             dom_end = 'dbh',
             evict = FALSE) %>%
  define_domains(dbh = c(0, 50, 100)) %>%
  make_ipm()

lambda_ipmr <- Re(eigen(x$iterators$K)$values[1])

lambda - lambda_ipmr

