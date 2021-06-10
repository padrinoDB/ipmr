
# Code straight out of IPM book, chapter 6

m.par.true <- c(## survival
  surv.int  = -1.70e+1,
  surv.z    =  6.68e+0,
  surv.a    = -3.34e-1,
  ## growth
  grow.int  =  1.27e+0,
  grow.z    =  6.12e-1,
  grow.a    = -7.24e-3,
  grow.sd   =  7.87e-2,
  ## reproduce or not
  repr.int  = -7.88e+0,
  repr.z    =  3.11e+0,
  repr.a    = -7.80e-2,
  ## recruit or not
  recr.int  =  1.11e+0,
  recr.a    =  1.84e-1,
  ## recruit size
  rcsz.int  =  3.62e-1,
  rcsz.z    =  7.09e-1,
  rcsz.sd   =  1.59e-1)

##
## Define the various demographic functions, we pass a vector of parameters "m.par" to each function
## this makes it easier to compare the results of different paramter sets, say the true values and
## estimated ones. Note the addition of the 'a' argument for the new functions
##

## Growth function, given you are size z and age a now returns the pdf of size z1 next time
g_z1z <- function(z1, z, a, m.par){
  mean <- m.par["grow.int"] + m.par["grow.z"] * z + m.par["grow.a"] * a
  sd <- m.par["grow.sd"]
  p.den.grow <- dnorm(z1, mean = mean, sd = sd)
  return(p.den.grow)
}

## Survival function, logistic regression
s_z <- function(z, a, m.par){
  linear.p <- m.par["surv.int"] + m.par["surv.z"] * z + m.par["surv.a"] * a
  p <- 1/(1+exp(-linear.p))
  return(p)
}


## Reproduction function, logistic regression
pb_z <- function(z, a, m.par){
  if (a==0) {
    p <- 0
  } else {
    linear.p <- m.par["repr.int"] + m.par["repr.z"] * z + m.par["repr.a"] * a
    p <- 1/(1+exp(-linear.p))
  }
  return(p)
}

## Recruitment function, logistic regression
pr_z <- function(a, m.par) {
  linear.p <- m.par["recr.int"] + m.par["recr.a"] * a
  p <- 1/(1+exp(-linear.p))
  return(p)
}

## Recruit size function
c_z1z <- function(z1, z, m.par){
  mean <- m.par["rcsz.int"] + m.par["rcsz.z"] * z
  sd <- m.par["rcsz.sd"]
  p.den.rcsz <- dnorm(z1, mean = mean, sd = sd)
  return(p.den.rcsz)
}


init.pop.size <- 500
n.yrs <- 400
m.par <- m.par.true

pb_z_ibm <- function(z, a, m.par){
  linear.p <- m.par["repr.int"] + m.par["repr.z"] * z + m.par["repr.a"] * a
  p <- 1/(1+exp(-linear.p))
  p <- ifelse(a==0, 0, p)
  return(p)
}

## initial size distribution (assuming everyone is a new recuit)
z <- rnorm(init.pop.size, mean = m.par["rcsz.int"] +  m.par["rcsz.z"] * 3.2, sd = m.par["rcsz.sd"])
a <- rep(0, init.pop.size)

## vectors to store pop size and mean size
pop.size.t <- mean.z.t <- mean.a.t <- mean.z.repr.t <- mean.a.repr.t <- numeric(n.yrs)

## Iterate the model using the "true" parameters and store data in a data.frame
yr <- 1
while(yr < n.yrs & length(z) < 5000) {

  ## calculate current population size
  pop.size <- length(z)

  ## generate binomial random number for survival, where survival depends on your size z,
  ## this is a vector of 0's and 1's, you get a 1 if you survive
  surv <- rbinom(n=pop.size, prob=s_z(z, a, m.par), size=1)

  ## generate the size of surviving individuals next year
  i.subset <- which(surv == 1)
  z1 <- rep(NA, pop.size)
  E.z1 <- m.par["grow.int"] + m.par["grow.z"] * z[i.subset] + m.par["grow.a"] * a[i.subset]
  z1[i.subset] <- rnorm(n = length(i.subset), mean = E.z1, sd = m.par["grow.sd"])

  ## generate a binomial random number for reproduction from surviving individuals
  repr <- rep(NA, pop.size)
  repr[i.subset] <- rbinom(n = length(i.subset), prob = pb_z_ibm(z[i.subset], a[i.subset], m.par), size=1)

  ## generate a binomial random number for offspring sex (female==1)
  i.subset <- which(surv == 1 & repr == 1)
  osex <- rep(NA, pop.size)
  osex[i.subset] <- rbinom(n = length(i.subset), prob = 1/2, size=1)

  ## generate a binomial random number for offspring recruitment from surviving / reproducing individuals
  i.subset <- which(surv == 1 & repr == 1 & osex == 1)
  recr <- rep(NA, pop.size)
  recr[i.subset] <- rbinom(n = length(i.subset), prob=pr_z(a[i.subset], m.par), size=1)

  ## generate the size of new recruits
  i.subset <- which(surv == 1 & repr == 1 & osex==1 & recr == 1)
  z1.rec <- rep(NA, pop.size)
  E.rec.z1 <- m.par["rcsz.int"] + m.par["rcsz.z"] * z[i.subset]
  z1.rec[i.subset] <- rnorm(n = length(i.subset), mean = E.rec.z1, sd = m.par["rcsz.sd"])

  ## store the simulation data for the current year, we'll use this later
  sim.data <- data.frame(z, a, surv, z1, repr, osex, recr, z1.rec)

  ## store the population size and mean body mass
  pop.size.t[yr] <- length(z)
  mean.z.t[yr] <- mean(z)
  mean.a.t[yr] <- mean(a)
  mean.z.repr.t[yr] <- mean(z[repr==1], na.rm=TRUE)
  mean.a.repr.t[yr] <- mean(a[repr==1], na.rm=TRUE)

  ## create new population body size vector
  i.recr <- which(surv == 1 & repr == 1 & osex==1 & recr == 1)
  z <- c(z1.rec[i.recr], z1[which(surv == 1)])
  a <- c(rep(0, length(i.recr)), a[which(surv == 1)]+1)

  ## iterate the year
  yr <- yr+1
}

## trim the population size vector to remove the zeros at the end
pop.size.t <- pop.size.t[pop.size.t>0]

## reassign column names to match the text of chapter 1
names(sim.data) <- c("z", "a", "Surv", "z1", "Repr", "Sex", "Recr", "Rcsz")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Estimation
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## take a random sample of 3000 observations
rows.to.keep <- sample.int(nrow(sim.data), 3000)
sim.data <- sim.data[rows.to.keep,]

## sort by size and print a sample to the screen
sim.data <- sim.data[order(sim.data$z),]
round(sim.data[sample(1:nrow(sim.data),12),][-c(3,4,6,7),],2)

##
## Fit the vital rate functions
##

mod.surv.glm <- glm(Surv ~ z + a, family = binomial, data = sim.data)
summary(mod.surv.glm)

grow.data <- subset(sim.data, !is.na(z1))
mod.grow <- lm(z1 ~ z + a, data = grow.data)
summary(mod.grow)

repr.data <- subset(sim.data, Surv==1 & a>0)
mod.repr <- glm(Repr ~ z + a, family = binomial, data = repr.data)
summary(mod.repr)

recr.data <- subset(sim.data, Surv==1 & Repr==1)
mod.recr <- glm(Recr ~ a, family=binomial, data=recr.data)
summary(mod.recr)

rcsz.data <- subset(sim.data, !is.na(Rcsz))
mod.rcsz <- lm(Rcsz ~ z, data=rcsz.data)
summary(mod.rcsz)

##
## Store the estimated parameters
##

m.par.est <- c(
  ## survival
  surv      = coef(mod.surv.glm),
  ## growth
  grow      =  coef(mod.grow),
  grow.sd   =  summary(mod.grow)$sigma,
  ## reproduce or not
  repr      =  coef(mod.repr),
  ## recruit or not
  recr      =  coef(mod.recr),
  ## recruit size
  rcsz      =  coef(mod.rcsz),
  rcsz.sd   =  summary(mod.rcsz)$sigma)

names(m.par.est) <- names(m.par.true)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 3 - Function definitions: implementing the kernel and iteration
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##
## Functions to build IPM kernels P(a), F(a),
##

## Define the survival kernel
P_z1z <- function (z1, z, a, m.par) {
  return( s_z(z, a, m.par) * g_z1z(z1, z, a, m.par) )
}


## Define the reproduction kernel
F_z1z <- function (z1, z, a, m.par) {
  return( s_z(z, a, m.par) * pb_z(z, a, m.par) * (1/2) * pr_z(a, m.par) * c_z1z(z1, z, m.par) )
}

##
## Functions to implement the IPM
##

## Calculate the mesh points, mesh width and store with upper/lower bounds and max age
mk_intpar <- function(m, L, U, M) {
  h <- (U - L) / m
  meshpts  <-  L + ((1:m) - 1/2) * h
  na <- M + 2
  return( list(meshpts = meshpts, M = M, na = na, h = h, m = m) )
}

## Build the list of age/process specific kernels + store the integration parameters in the same list
mk_age_IPM <- function(i.par, m.par) {
  within(i.par, {
    F <- P <- list()
    for (ia in seq_len(na)) {
      F[[ia]] <- outer(meshpts, meshpts, F_z1z, a = ia-1, m.par = m.par) * h
      P[[ia]] <- outer(meshpts, meshpts, P_z1z, a = ia-1, m.par = m.par) * h

    }
    rm(ia)
  })
}

##
## Functions to iterate the age-size IPM
##

## iterate 'forward' one time step to project dynamics. 'x' is the list of age-specific size dists
r_iter <- function(x, na, F, P) {
  xnew <- list(0)
  for (ia in seq_len(na)) {
    xnew[[1]] <- (xnew[[1]] + F[[ia]] %*% x[[ia]])[,,drop=TRUE]
  }
  for (ia in seq_len(na-1)) {
    xnew[[ia+1]] <- (P[[ia]] %*% x[[ia]])[,,drop=TRUE]
  }
  xnew[[na]] <- xnew[[na]] + (P[[na]] %*% x[[na]])[,,drop=TRUE]
  return(xnew)
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 4 - Compare truth and estimated model predictions of lambda, and calculate w(z) and v(z), using
## iteration for 100 steps
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# integration parameters
i.par <- mk_intpar(m = 100, M = 20, L = 1.6, U = 3.7)

# make an initial state vector
nt0 <- with(i.par, lapply(seq_len(na), function(ia) rep(0, m)))
nt0[[1]] <- with(i.par, rep(1 / m, m))

# estimate lambda with the true parameters...

IPM.sys <- mk_age_IPM(i.par, m.par.true)
IPM.sim <- with(IPM.sys, {
  x <- nt0
  for (i in seq_len(100)) {
    x1 <- r_iter(x, na, F, P)
    lam <- sum(unlist(x1))
    x <- lapply(x1, function(x) x / lam)
  }
  list(lambda = lam, x = x)
})

sim_lambda <- IPM.sim$lambda


## IPMR version

param_list <- as.list(m.par.true) %>%
  setNames(gsub(pattern = "\\.", replacement = "_", x = names(.)))

inv_logit <- function(x) {

  return( 1 / (1 + exp(-x)) )
}

f_fun <- function(age, s_age, pb_age, pr_age, recr) {

  if(age == 0) return(0)

  s_age * pb_age * pr_age * recr * 0.5

}

a_s_ipm <- init_ipm(sim_gen    = "general",
                    di_dd      = "di",
                    det_stoch  = "det",
                    uses_age = TRUE) %>%
  define_kernel(
    name          = "P_age",
    family        = "CC",
    formula       = s_age * g_age * d_wt,
    s_age         = inv_logit(surv_int + surv_z * wt_1 + surv_a * age),
    g_age         = dnorm(wt_2, mu_g_age, grow_sd),
    mu_g_age      = grow_int + grow_z * wt_1 + grow_a * age,
    data_list     = param_list,
    states        = list(c("wt")),
    uses_par_sets = FALSE,
    age_indices   = list(age = c(0:20), max_age = 21),
    evict_cor     = FALSE
  ) %>%
  define_kernel(
    name          = "F_age",
    family        = "CC",
    formula       = f_fun(age, s_age, pb_age, pr_age, recr) * d_wt,
    s_age         = inv_logit(surv_int + surv_z * wt_1 + surv_a * age),
    pb_age        = inv_logit(repr_int + repr_z * wt_1 + repr_a * age),
    pr_age        = inv_logit(recr_int + recr_a * age),
    recr          = dnorm(wt_2, rcsz_mu, rcsz_sd),
    rcsz_mu       = rcsz_int + rcsz_z * wt_1,
    data_list     = param_list,
    states        = list(c("wt")),
    uses_par_sets = FALSE,
    age_indices   = list(age = c(0:20), max_age = 21),
    evict_cor     = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P_age", "F_age"),
      int_rule     = rep("midpoint", 2),
      state_start    = c("wt_age", "wt_age"),
      state_end      = c("wt_age", "wt_0")
    )
  ) %>%
  define_domains(
    wt = c(1.6, 3.7, 100)
  ) %>%
  define_pop_state(
    pop_vectors = list(
      n_wt_age = runif(100))
  ) %>%
  make_ipm(
    usr_funs = list(inv_logit = plogis,
                    f_fun     = f_fun),
    iterate  = TRUE,
    iterations = 100,
    return_all_envs = TRUE
  )

ipmr_lambda <- unname(lambda(a_s_ipm, type_lambda = 'last'))

test_that("age_x_size lambdas are recovered properly", {

  expect_equal(ipmr_lambda, sim_lambda, tolerance = 1e-10)

  pop_state_ipmr <- a_s_ipm$pop_state %>%
    .[names(.) != 'lambda']

  final_pop_comp <- vapply(1:22,
                           function(i, ipmr, sim) {
                             isTRUE(
                               all.equal(
                                 ipmr[[i]][, 101], sim[[i]]
                               )
                             )
                           },
                           FUN.VALUE = logical(1L),
                           ipmr = pop_state_ipmr,
                           sim = IPM.sim$x)

  expect_true(all(final_pop_comp))

})

test_that("format_mega_kernel.age_x_size_ipm works", {

  subs <- a_s_ipm$sub_kernels

  target_f <- do.call("cbind",
                      subs[grepl("F", names(subs))])

  target_p <- matrix(0, nrow = 2100, ncol = 2200)
  Ps       <- subs[grepl("P", names(subs))]

  for(i in seq_len(21)) {

    row_ind <- col_ind <- seq(i * 100 - 99, i * 100, by = 1)
    target_p[row_ind, col_ind] <- Ps[[i]]

  }

  target_p[2001:2100, 2101:2200] <- Ps[[22]]

  target <- rbind(target_f, target_p)
  rm(target_f, target_p)

  ipmr_meg <- format_mega_kernel(a_s_ipm,
                                 name_ps = "P",
                                 f_forms = "F")$mega_mat

  expect_equal(ipmr_meg, target, tolerance = 1e-10)

})


# test if we can implement seperate kernels for max_age/age selectively.

P_max_z1z <- function (z1, z, a, m.par) {
  return( s_max_z(z, m.par) * g_z1z(z1, z, a, m.par) )
}

s_max_z <- function(z, m.par){
  linear.p <- m.par["surv.int"] + m.par["surv.z"] * z + m.par["surv.z.2"] * z^2
  p <- 1/(1+exp(-linear.p))
  return(p)
}


mk_age_IPM <- function(i.par, m.par) {

  M <- i.par$M
  m <- i.par$m
  meshpts <- i.par$meshpts
  h <- i.par$h
  na <- i.par$na
  F <- P <- list()

  for (ia in seq_len(na)) {
    F[[ia]] <- outer(i.par$meshpts, i.par$meshpts, F_z1z, a = ia-1, m.par = m.par) * h
    if(ia != na) {
      P[[ia]] <- outer(meshpts, meshpts, P_z1z, a = ia-1, m.par = m.par) * h
    } else {

      P[[ia]] <- outer(meshpts, meshpts, P_max_z1z, a = ia-1, m.par = m.par) * h

    }
  }

  out <- c(i.par, list(P = P, F = F))

  return(out)

}

m.par.true["surv.z.2"] <- -3

IPM.sys_max <- mk_age_IPM(i.par, m.par.true)
IPM.sim_max <- with(IPM.sys_max, {
  x <- nt0
  for (i in seq_len(100)) {
    x1 <- r_iter(x, na, F, P)
    lam <- sum(unlist(x1))
    x <- lapply(x1, function(x) x / lam)
  }
  list(lambda = lam, x = x)
})

sim_lambda_max <- IPM.sim_max$lambda

param_list <- as.list(m.par.true) %>%
  setNames(gsub(pattern = "\\.", replacement = "_", x = names(.)))

a_s_ipm <- init_ipm(sim_gen    = "general",
                    di_dd      = "di",
                    det_stoch  = "det",
                    uses_age = TRUE) %>%
  define_kernel(
    name          = "P_age",
    family        = "CC",
    formula       = s_age * g_age * d_wt,
    s_age         = inv_logit(surv_int + surv_z * wt_1 + surv_a * age),
    g_age         = dnorm(wt_2, mu_g_age, grow_sd),
    mu_g_age      = grow_int + grow_z * wt_1 + grow_a * age,
    data_list     = param_list,
    states        = list(c("wt")),
    uses_par_sets = FALSE,
    age_indices   = list(age = c(0:20), max_age = 21),
    evict_cor     = FALSE
  ) %>%
  define_kernel(
    name          = "P_max_age",
    family        = "CC",
    formula       = s_max_age * g_max_age * d_wt,
    g_max_age     = dnorm(wt_2, mu_g_max_age, grow_sd),
    mu_g_max_age  = grow_int + grow_z * wt_1 + grow_a * max_age,
    s_max_age     = inv_logit(surv_int + surv_z * wt_1 + surv_z_2 * wt_1^2),
    data_list     = param_list,
    states        = list(c("wt")),
    uses_par_sets = FALSE,
    age_indices   = list(age = c(0:20), max_age = 21),
    evict_cor     = FALSE) %>%
  define_kernel(
    name          = "F_age",
    family        = "CC",
    formula       = f_fun(age, s_age, pb_age, pr_age, recr) * d_wt,
    s_age         = inv_logit(surv_int + surv_z * wt_1 + surv_a * age),
    pb_age        = inv_logit(repr_int + repr_z * wt_1 + repr_a * age),
    pr_age        = inv_logit(recr_int + recr_a * age),
    recr          = dnorm(wt_2, rcsz_mu, rcsz_sd),
    rcsz_mu       = rcsz_int + rcsz_z * wt_1,
    data_list     = param_list,
    states        = list(c("wt")),
    uses_par_sets = FALSE,
    age_indices   = list(age = c(0:20), max_age = 21),
    evict_cor     = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P_age", "F_age", "P_max_age"),
      int_rule     = rep("midpoint", 3),
      state_start    = c("wt_age", "wt_age", "wt_max_age"),
      state_end      = c("wt_age", "wt_0", "wt_max_age")
    )
  )%>%
  define_domains(
    wt = c(1.6, 3.7, 100)
  ) %>%
  define_pop_state(
    pop_vectors = list(
      n_wt_age = runif(100))
  ) %>%
  make_ipm(
    usr_funs = list(inv_logit = plogis,
                    f_fun     = f_fun),
    iterate  = TRUE,
    iterations = 100,
    return_all_envs = TRUE
  )


ipmr_lambda_max <- unname(lambda(a_s_ipm))

test_that("ipmr can handle max-age specific expressions", {

  expect_equal(sim_lambda_max, ipmr_lambda_max)

  expect_equal(a_s_ipm$sub_kernels$P_21,
               IPM.sys_max$P[[22]],
               ignore_attr = TRUE,
               tolerance = 1e-7)

})
