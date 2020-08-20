# Age-size model rebuild

# First, reimplement their model and compute times. Next, we'll modify their
# code to produce output comparable to ipmr

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


# ipmr version and timings ----------

library(ipmr)

# Set parameter values and names

param_list <- list(
  ## Survival
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
  rcsz.sd   =  1.59e-1
)  %>%
  setNames(gsub(pattern = "\\.", replacement = "_", x = names(.)))

# define a custom function to handle the F kernels. We need this because
# F_0 == 0. It may be possible to do this with a an expression in the
# define_k(...), but I haven't thought that far ahead yet.

f_fun <- function(age, s_age, pb_age, pr_age, recr) {

  if(age == 0) return(0)

  s_age * pb_age * pr_age * recr * 0.5

}

# Implement the age x size model

age_size_ipm <- init_ipm("general_di_det", has_age = TRUE) %>%
  define_kernel(
    name          = "P_age",
    family        = "CC",
    formula       = s_age * g_age * d_z,
    s_age         = plogis(surv_int + surv_z * z_1 + surv_a * age),
    g_age         = dnorm(z_2, mu_g_age, grow_sd),
    mu_g_age      = grow_int + grow_z * z_1 + grow_a * age,
    data_list     = param_list,
    states        = list(c("z")),
    has_hier_effs = FALSE,
    levels_ages   = list(age = c(0:20), max_age = 21),
    evict_cor     = FALSE
  ) %>%
  define_kernel(
    name          = "F_age",
    family        = "CC",
    formula       = f_fun(age, s_age, pb_age, pr_age, recr) * d_z,
    s_age         = plogis(surv_int + surv_z * z_1 + surv_a * age),
    pb_age        = plogis(repr_int + repr_z * z_1 + repr_a * age),
    pr_age        = plogis(recr_int + recr_a * age),
    recr          = dnorm(z_2, rcsz_mu, rcsz_sd),
    rcsz_mu       = rcsz_int + rcsz_z * z_1,
    data_list     = param_list,
    states        = list(c("z")),
    has_hier_effs = FALSE,
    levels_ages   = list(age = c(0:20), max_age = 21),
    evict_cor     = FALSE
  ) %>%
  define_k(
    name               = "K",
    family             = "IPM",
    n_z_0_t_1         = all_ages(F_age %*% n_z_age_t, "+"),
    n_z_age_t_1       = P_age_minus_1 %*% n_z_age_minus_1_t,
    n_z_max_age_t_1   = P_max_age %*% n_z_max_age_t +
      P_max_age_minus_1 %*% n_z_max_age_minus_1_t,
    data_list          = param_list,
    states             = list (c("z")),
    has_hier_effs      = FALSE,
    levels_ages        = list(age = c(0:20), max_age = 21),
    evict_cor          = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P_age", "F_age", "K"),
      int_rule     = rep("midpoint", 3),
      dom_start    = rep("z", 3),
      dom_end      = rep("z", 3)
    )
  ) %>%
  define_domains(
    z = c(1.6, 3.7, 100)
  ) %>%
  define_pop_state(
    pop_vectors = list(
      n_z_age = runif(100))
  ) %>%
  make_ipm(
    usr_funs = list(f_fun = f_fun),
    iterate  = TRUE,
    iterations = 100
  )


lamb <- lambda(a_s_ipm)

stopifnot(
  isTRUE(
    all.equal(
      IPM.sim$lambda, lamb, tolerance = 1e-13
    )
  )
)



