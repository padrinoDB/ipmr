[![R build
status](https://github.com/levisc8/ipmr/workflows/R-CMD-check/badge.svg)](https://github.com/levisc8/ipmr/actions)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![CRAN
status](https://www.r-pkg.org/badges/version/ipmr)](https://cran.r-project.org/package=ipmr)
[![Codecov test
coverage](https://codecov.io/gh/levisc8/ipmr/branch/master/graph/badge.svg)](https://codecov.io/gh/levisc8/ipmr?branch=master)

ipmr
====

`ipmr` is a package for implementing Integral Projection Models (IPMs)
in *R*. It relies heavily on the mathematical syntax of the models, and
does not try to abstract over the process of fitting vital rates. It is
still very much a work in progress, and so models implemented using
current code may not run in the not-so-distant future as more
complicated methods are brought online and tweaks are made. Below is a
brief overview of how `ipmr` classifies different model types followed
by examples of how to implement those types in this framework.

Installation
------------

`ipmr` is not yet on CRAN. You can install the development version with
the snippet below:

``` r
if(!require('remotes', quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("levisc8/ipmr", build_vignettes = TRUE)
```

Package scope
-------------

Note that this package **will not** help with the process of fitting
vital rate models at all! That is a sufficiently different (and vast)
question that we decided it was not within the scope of this project.
This will only help you turn those regression models into an IPM without
shooting yourself in the foot. Thus, everything that follows assumes you
have parameterized those vital rate models and are now ready to begin
implementing your model.

Below is a brief overview of the package and some examples of how to
implement models with it. A more thorough introduction is available
[here](https://levisc8.github.io/ipmr/articles/ipmr-introduction.html).

Model classes
-------------

Once all parameters are estimated, the first step of defining a model in
`ipmr` is to initialize the model using `init_ipm()`. This function has
two arguments: `model_class`, and `has_age`. We will ignore `has_age`
for now, because age-size models are a complicated topic and have their
[own vignette](https://levisc8.github.io/ipmr/articles/age_x_size.html).

The `model_class` defines the basic infrastructure that will be
available for subsequent analyses and helps make sure the kernels are
correctly implemented from the underlying vital rates. `model_class`
should be a character string with at least 3 (but possibly 4) entries
separated by underscores (`_`). Below, the are the possible entries for
each position.

-   Position 1: `"simple"`/`"general"`

-   1.  **simple**: This describes an IPM with a single continuous state
        variable and no discrete stages.

-   1.  **general**: This describes and IPM with either more than one
        continuous state variable, one or more discrete stages, or both
        of the above. Basically, anything other than an IPM with a
        single continuous state variable.

-   Position 2: `"di"`/`"dd"`

-   A. **di**: This is used to denote a density-independent IPM.

-   B. **dd**: This is used to denote a density-dependent IPM.

-   Position 3: `"det"`/`"stoch"`

-   A. **det**: This is used to denote a deterministic IPM. If this is
    used in the third position of `model_class`, there should not be a
    fourth entry.

-   B. **stoch**: This is used to denote a stochastic IPM. If this is
    used in the third position of `model_class`, there should always be
    a fourth entry. The two possibilities for the fourth are described
    next.

-   Position 4: `"kern"`/`"param"` (Complete definitions found in
    [Metcalf et
    al. 2015](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12405))

-   A. **kern**: This describes an IPM with discretely varying
    parameters such that their values are known before the model is
    specified. This is usually the case with models that estimate fixed
    and/or random year/site effects and for which defining a
    multivariate joint distribution to sample parameters from is not
    desirable/needed. These models can be a bit more computationally
    efficient than the `param` alternative because all kernels can be
    constructed before the iteration procedure begins, as opposed to
    requiring reconstruction for every single iteration.

-   B. **param**: This describes an IPM with parameters that are
    re-sampled from some distribution at each iteration of the model.
    For example, this can be a multivariate normal defined by covarying
    slopes and intercepts, or posterior distribution(s) from a Bayesian
    model. All that is required is that the parameters/values for the
    distribution are specified and that the function that generates the
    parameters at each iteration returns named lists that correspond to
    the parameter names in the model. These are covered in greater depth
    in towards the end of [this
    article](https://levisc8.github.io/ipmr/articles/ipmr-introduction.html).

With the type of model selected, the `model_class` becomes a string and
the call to `init_ipm` is composed like so:
`init_ipm(model_class = "position1_position_2_position3_position4")`.

The following possibilities are currently or will become available in
`ipmr` (bold text denotes development progress):

-   Simple, density independent models: **Completed and ready**

1.  `"simple_di_det"`

2.  `"simple_di_stoch_kern"`

3.  `"simple_di_stoch_param"`

-   Simple, density dependent models: **In progress, not yet stable**

1.  `"simple_dd_det"`

2.  `"simple_dd_stoch_kern"`

3.  `"simple_dd_stoch_param"`

-   General, density independent models: **Completed and ready**

1.  `"general_di_det"`

2.  `"general_di_stoch_kern"`

3.  `"general_di_stoch_param"`

-   General, density dependent models: **In progress, not yet stable**

1.  `"general_dd_det"`

2.  `"general_dd_stoch_kern"`

3.  `"general_dd_stoch_param"`

Simple density-independent deterministic, simple kernel-resampled
stochastic, and simple parameter resampled stochastic models
(`simple_di_det`, `simple_di_stoch_kern`, `simple_di_stoch_param`) are
functional and detailed below as well as
[here](https://levisc8.github.io/ipmr/articles/ipmr-introduction.html).
The `general_*` versions of these are also ready, and an introduction to
them is available
[here](https://levisc8.github.io/ipmr/articles/general-ipms.html).
Density dependent versions are completed for simple models, and in
progress for general ones. These are not yet stable. Below is an example
implementing a `simple_di_det` IPM.

Next on the to-do list is to write generic functions for `lambda`,
`right_ev` (right eigenvector), `left_ev` (left eigenvector), and
`plot`/`print` methods for all exported classes (`proto_ipm`,
implemented `ipm` objects).

Examples for implemented IPM types
----------------------------------

Here is a simple model implemented with `ipmr`. It will use the
following set of vital rate models:

1.  Survival (`s`): a generalized linear model w/ a logit link.

    -   Example model formula:
        `glm(surv ~ size_1, data = my_surv_data, family = binomial())`

2.  Growth (`g`): a linear model with a Normal error distribution.

    -   Example model formula:
        `lm(size_2 ~ size_1, data = my_grow_data)`

3.  Pr(flowering) (`f_r`): a generalized linear model w/ a logit link.

    -   Example model formula:
        `glm(flower ~ size_1, data = my_repro_data, family = binomial())`

4.  Seed production (`f_s`): a generalized linear model w/ log link.

    -   Example model formula:
        `glm(seeds ~ size_1, data = my_flower_data, family = poisson())`

5.  Recruit size distribution (`f_d`): a normal distribution w
    parameters `mu_fd` (mean) and `sd_fd` (standard deviation).

The example below assumes we’ve already fit our vital rate models from
the raw data. In this example, the numbers are made up, but code that
extracts the values you need from a real regression model is provided in
comments.

``` r
# Load ipmr and get the parameter values. The data_list argument for define_kernel
# should hold every regression parameter and every constant used in the model.

library(ipmr)

data_list = list(s_int     = -2.2,  # coefficients(my_surv_mod)[1]
                 s_slope   = 0.25,  # coefficients(my_surv_mod)[2]
                 g_int     = 0.2,   # coefficients(my_grow_mod)[1]
                 g_slope   = 1.02,  # coefficients(my_grow_mod)[2]
                 sd_g      = 0.7,   # sd(resid(my_grow_mod))
                 f_r_int   = 0.003, # coefficients(my_pr_flower_mod)[1]
                 f_r_slope = 0.015, # coefficients(my_pr_flower_mod)[2]
                 f_s_int   = 1.3,   # coefficients(my_seed_mod)[1]
                 f_s_slope = 0.075, # coefficients(my_seed_mod)[2]
                 mu_fd     = 2,     # mean(recruit_data$size_next)
                 sd_fd     = 0.3)   # sd(recruit_data$size_next)

my_simple_ipm <- init_ipm('simple_di_det')

my_simple_ipm <- define_kernel(
  
  proto_ipm = my_simple_ipm,
  
  # Name of the kernel
  
  name      = "P_simple",
  
  # The type of transition it describes (e.g. continuous - continuous, discrete - continuous).
  # These must be specified for all kernels!
  
  family    = "CC",
  
  # The formula for the kernel. 
  
  formula   = s * g,
  
  # A named set of expressions for the vital rates it includes. 
  # Each state variable has a stateVariable_1 and stateVariable_2 internally
  # defined for the domain associated with it. Use these to distinguish between 
  # size/weight/etc at time t vs size/weight/etc at time t+1
  
  # Here, we create an expression for the inverse of the link function from 
  # our GLM for survival. In this case, it's the inverse logit.
  # For examples on using predict(my_surv_mod,...),
  # see the articles at https://levisc8.github.io/ipmr.
  
  s         = 1 / (1 + exp(-(s_int + s_slope * dbh_1))), 
  
  # The growth model requires a function to compute the mean as a function of dbh.
  # The SD is a constant, so we don't need to define that in ... expression, 
  # just the data_list.
  
  g         = dnorm(dbh_2, mu_g, sd_g),
  mu_g      = g_int + g_slope * dbh_1,
  
  
  # Specify the constants in the model in the data_list. 
  
  data_list = data_list,
  states    = list(c('dbh')),
  
  # If you want to correct for eviction, set evict = TRUE and specify an
  # evict_fun. ipmr provides truncated_distributions() to help.
  
  evict     = TRUE,
  evict_fun = truncated_distributions('norm',
                                      'g')
  )

# Define the F kernel now

my_simple_ipm <- define_kernel(
  proto_ipm = my_simple_ipm,
  name      = 'F_simple',
  formula   = f_r * f_s * f_d,
  family    = 'CC',
  
  # Inverse logit transformation for flowering probability
  # (because we used a logistic regression)
  
  f_r       = 1 / (1 + exp( - (f_r_int + f_r_slope * dbh_1))),
  
  # Exponential function for seed progression 
  # (because we used a Poisson)
  
  f_s       = exp(f_s_int + f_s_slope * dbh_1),
  
  # The recruit size distribution has no maternal effect for size,
  # so mu_fd and sd_fd are constants. These get passed in the 
  # data_list
  
  f_d       = dnorm(dbh_2, mu_fd, sd_fd),
  data_list = data_list,
  states    = list(c('dbh')),
  
  # Again, we'll correct for eviction in new recruits by
  # truncating the normal distribution.
  
  evict     = TRUE,
  evict_fun = truncated_distributions('norm',
                                      'f_d')
)

# Next, define the K kernel

my_simple_ipm <- define_k(
  
  # K kernels get their own special define_k function. Rather than use the formula
  # parameter, it simply takes named expressions for the format of the iteration
  # kernels. It can also take expressions showing how these kernels generate 
  # population states at t+1 as a function of population state at t - see examples
  # below for how to do that.
  
  proto_ipm = my_simple_ipm,
  
  name      = 'K',
  K         = P_simple + F_simple,
  
  # This is a new family - at the moment all K's will get the
  # 'IPM' family
  
  family    = 'IPM',
  
  # This kernel has no additional parameters, so the data_list
  # is empty
  
  data_list = list(),
  states    = list(c('dbh')),
  
  # We've already corrected eviction in the sub-kernels, so there's no
  # need to do that here
  
  evict     = FALSE
) 
  
# Next, we have to define the implementation details for the model. 
# We need to tell ipmr how each kernel is integrated, what domain
# it starts on (i.e. the size/weight/etc from above), and what domain
# it ends on. In simple_* models, dom_start and dom_end will always be the same,
# because we only have a single continuous state variable. general_*
# models will be more complicated.

my_simple_ipm <- define_impl(
  proto_ipm = my_simple_ipm,
  make_impl_args_list(
    kernel_names = c("K_simple", "P_simple", "F_simple"),
    int_rule     = rep("midpoint", 3),
    dom_start    = rep("dbh", 3),
    dom_end      = rep("dbh", 3)
  )
) 

my_simple_ipm <- define_domains(
  proto_ipm = my_simple_ipm,
  dbh = c(0, # the first entry is the lower bound of the domain.
          50, # the second entry is the upper bound of the domain.
          100 # third entry is the number of meshpoints for the domain.
  ) 
) 

my_simple_ipm <- make_ipm(proto_ipm = my_simple_ipm)


lambda_ipmr <- lambda(my_simple_ipm, comp_method = 'eigen')
w_ipmr      <- right_ev(my_simple_ipm)
```

More complicated models
-----------------------

Examples of more complicated models are included in the vignettes,
accesible using either `browseVignettes('ipmr')` or by visiting the
Articles tab on [project’s webpage](https://levisc8.github.io/ipmr/).
Please file all bug reports in the Issues tab of this repository or
contact me via [email](mailto:levisc8@gmail.com) with a reproducible
example.

Code of Conduct
---------------

We welcome contributions from other developers. Please note that the
ipmr project is released with a [Contributor Code of
Conduct](https://levisc8.github.io/ipmr/CODE_OF_CONDUCT.html). By
contributing to this project, you agree to abide by its terms.
