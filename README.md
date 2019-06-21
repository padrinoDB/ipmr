
# ipmr

Simple, density-independent, deterministic models (`simple_di_det`) are
now functional. However, expect changes as more complicated methods are
implemented\! See below for an example of how to implement an IPM in
this framework.

Next on the implementation to-do list are simple, stochastic IPMs using
either parameter resampling or kernel
resampling.

``` r
# Example of the setupd for a simple IPM without density dependence or environmental
# stochasticity

library(ipmr)

state_list <- list(c('dbh'))

my_ipm <- init_ipm(class = 'simple_di_det') %>% # the %>% operator is re-exported, so no need to load dplyr/magrittr
  add_kernel(
    # Name of the kernel
    name = "P",
    # The type of transition it describes (e.g. continuous - continuous, discrete - continuous)
    family = "CC",
    # The formula for the kernel. don't forget to transpose the growth matrix when multiplying survival!
    formula = t(s * t(g)),
    # A named set of expressions for the vital rates it includes. ipmr automatically
    # computes the bin width creates a variable called "cell_size_stateVariable"
    # for usage in expressions that require integrating a density function
    s = 1/(1 + exp(-(bs_0 + bs_1 * dbh_1))),
    g = cell_size_dbh * dnorm(dbh_2, mean_g, sd_g),
    mean_g = bg_0 + bg_1 * dbh_1,
    # A list containing named values for each constant 
    data_list = list(bs_0 = 0.2, # could also be coef(my_surv_mod)[1]
                     bs_1 = 0.978,
                     bg_0 = 0.78,
                     bg_1 = 0.967,
                     sd_g = 3.4),
    state_list = state_list, 
    # The state variable for time T
    domain_start = "dbh",
    # The state variable for T+1. for IPM_CC's, this will almost always be the same 
    # as domain_start. Scroll down to see examples of when it isn't
    domain_end = "dbh",
    int_rule = "midpoint",
    # The type of eviction correction. Set to NA if you wish to correct on your own
    evict = TRUE,
    evict_fun = truncated_distributions(g,
                                        n_mesh_p = 100)
    ) %>%
  add_kernel(
    name = "F",
    family = "CC",
    formula = r_s * r_n * r_p,
    r_s = dnorm(dbh_2, recr_mean, recr_sd),
    r_n = exp(br_0 + br_1 * dbh_1),
    r_p = 1/(1 + exp(-(bp_0 + bp_1 * dbh_1))),
    data_list = list(recr_mean = 0.67,
                     recr_sd = 0.3,
                     br_0 = 2.9,
                     br_1 = 6.5,
                     bp_0 = 0.65,
                     bp_1 = 0.72),
    state_list = state_list,
    dom_start = "dbh",
    dom_end = "dbh",
    int_rule = "midpoint", 
    evict = TRUE,
    evict_type = truncated_distributions(r_s,
                                         n_mesh_p = 100),
  ) %>%
  add_kernel(
    name = "K",
    formula = P + F,
    family = "IPM",
    state_list = state_list,
    dom_start = 'dbh',
    dom_end = 'dbh',
    evict = FALSE # this is dealt with in the sub-kernels, so no need to do so again.
  ) %>%
  define_domains(
    dbh = c(0, # the first entry is the lower bound of the domain.
            50, # the second entry is the upper bound of the domain.
            100 # third entry is the number of meshpoints for the domain.
    )
  ) %>%
  make_ipm()
```
