# ipmr 0.0.4

Contains small tweaks, bug fixes, and new feature additions. There shouldn't be any breaking changes to the API. 

## Features

- Adds `log` argument to `lambda` so users can choose which scale to return. The default is `TRUE` for stochastic models and `FALSE` for deterministic models.

- Allows named expressions in `define_env_state()` such that the names of the expressions may be used in `define_kernel()`. Before, the name of the expression didn't matter, only the names of the outputted list. Now, this will also work as well.  

```
define_env_state(proto_ipm, 
                 temp = rnorm(1, 20, 3), 
                 precip = rgamma(1, 400, 2),
                 data_list = list())

```

- More unit tests for parameter resampled models.

- Removed innocuous warning messages. Added warnings for `NA` values in a few `define_*` functions and errors when they're produced in sub-kernels.

- Added the `ipmr_ipm` class so that most functions can tell the difference between a `PADRINO` object and an `ipmr` object.

- Changes argument `tol` to `tolerance` in `is_conv_to_asymptotic()` for consistency with other function/argument names.

## Bug fixes

- corrects a bug where functions in the `define_env_state` `data_list` argument weren't recognized.

- corrects bug in `left_ev()` and `right_ev()` where parameter set indices were ignored for deterministic models .

# ipmr 0.0.3

Contains a some tweaks and bug fixes, and a few new features:

## Features

  - Implements `right_ev()` and `left_ev()` methods for stochastic models.
  
  - Adds a new function, `conv_plot()`, to graphically check for convergence to asymptotic dynamics in deterministic models. 
  
  - Adds a new function, `discretize_pop_vector()`, to compute the empirical density function for a population trait distribution given a set of observed trait values.
  
  - Adds print methods for density dependent models.
  
  - Adds `log` argument for `lambda`.
  
## Bug fixes

  - Corrects bug where `tol` argument was ignored in `is_conv_to_asymptotic()`.
  
  - Gives output from `lambda()` names. Before, outputs from deterministic models with many parameter sets became hard to follow. 

# ipmr 0.0.2

Contains a some tweaks and bug fixes. There is one major API change that renames parameters in `define_kernel`.

  - Renames function arguments `hier_effs` -> `par_sets`, `levels_ages`/`levels_hier_effs` -> `age_indices`/`par_set_indices`. The idea was to shift from thinking about these IPMs as resulting from multilevel/hierarchical regresssion models to IPMs constructed from parameter sets (which can be derived from any number of other methods). 
  
  - Corrects some bugs that caused `vital_rate_funs()` to break for `stoch_param` and density dependent models.
  
  - Updates the age X size model interface so that `max_age` kernels can be specified separately if they have different functional forms from their non-`max_age` versions.
  
  - `make_iter_kernel` can handle computations passed into `mega_mat` (e.g. `mega_mat = c(P + F, 0, I, C))`.
  
  - Makes `plot.ipmr_matrix` more flexible, which is now the recommended default `plot` method for `ipm` objects. 
  
  - Changes to internal code that won't affect user experience. 
  
# ipmr 0.0.1

This is the first version of `ipmr`. It contains methods for constructing a variety of IPMs as well as methods for basic analysis. Complete documentation is in the vignettes and on the package website.
