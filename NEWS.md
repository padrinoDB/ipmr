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
