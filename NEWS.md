# ipmr 0.0.5

Contains a number of bug fixes and some new functionality. The latter are mostly related to PADRINO models, and shouldn't have *too much* of an effect on existing user-specified IPMs.

## Features 

- Adds `return_sub_kernels` argument to `make_ipm()` for `*_stoch_param` and all density dependent methods. This is due to the large RAM footprint that these objects can occupy, particularly for long running models. The default is `FALSE`, which will save considerable memory space. Can be set to `TRUE` for downstream analyses that require sub-kernels/iteration kernels. 

- Prettier printing method for PADRINO objects/lists of user-specified models.

- Prettier warnings for `left/right_ev()` 

- Depending on your view, this may be a bug fix: updates the `log` argument in `lambda()` so that it only changes the scale, NOT the calculation method. The prior behavior was documented, but unlikely to be intuitive, and so caused some confusion. Thanks to @aariq for pointing this out. 

- For stochastic IPMs, `is_conv_to_asymptotic()` and `conv_plot()` now check for convergence in stochastic lambda.  That is, they use a cumulative mean `log(lambda)` over iterations after discarding a burn-in. Thanks to @aariq for implementing this in #45.

- Adds experimental function `make_ipm_report()`. This function converts `proto_ipm` objects into Rmarkdown documents that, when rendered, contain the equations and parameters used to implement the model in Latex. See the function documentation for more details, and report bugs/notation quirks/preferences in the issue tracker.

## Bug fixes 

 - Fixes bug in normalization of left/right eigenvectors in simple density independent stochastic models. (thanks to @aariq)
 
 - Fixes bug where `"drop_levels"` was not recognized in some parameter set indexed models.
 
 - Fixes bug where parameter set levels were recycled in some cases. 
 
 - Fixes bug in `lambda()` where `log = TRUE` had no effect when the IPM was stochastic and `type_lambda` was `'last'` or `'all'`.
 
 - Switched from `all.equal` to absolute tolerance. @aariq in #52.
 
 - Fixes bug where the value of `lambda()` was named for determinstic models and unnamed for stochastic models.  @aariq in #57.
 
 - Fixes bug in `make_ipm()` where a user-specified sequence wasn't getting used correctly for certain model classes. @aariq in #58. 
 
 - Warnings about `NA`s in `data_list` are no longer raised for model objects. 

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
