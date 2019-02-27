library(rlang)
m <- expr(g_int + g_slope * dbh_1)
g_int <- 0.5
g_slope <- 1.02
s <- 3.5

# Not automatically normalized
grow <- function(dbh_2, dbh_1, g_int, g_slope, s, expre) {
  mu <- eval_bare(expre)
  dnorm(dbh_2, mu, s)
}

g_mat <- outer(dbh_2, dbh_1, FUN = grow,
               g_int = g_int,
               g_slope = g_slope,
               s = s,
               expre = m)
g_mat <- g_mat * (dbh_1[2] - dbh_1[1])
colSums(g_mat)

for(i in seq_len(dim(g_mat)[1])) {
  g_mat[ , i] <- g_mat[ , i]/sum(g_mat[ , i])
}

# automatically normalized

norm_grow <- function(dbh_2, dbh_1, g_int, g_slope, s, expre, cum_prob) {
  mu <- eval_bare(expre)
  dnorm(dbh_2, mu, s) / eval_bare(cum_prob)

}

cum_prob <- expr(pnorm(max(dbh_1), mu, s) - pnorm(min(dbh_1), mu, s))

g_mat_norm <- outer(dbh_2, dbh_1, FUN = norm_grow,
                    g_int = g_int,
                    g_slope = g_slope,
                    s = s,
                    expre = m,
                    cum_prob = cum_prob)

g_mat_norm <- g_mat_norm * (dbh_1[2] - dbh_1[1])
colSums(g_mat_norm)
