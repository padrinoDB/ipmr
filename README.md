[![R build
status](https://github.com/levisc8/ipmr/workflows/R-CMD-check/badge.svg)](https://github.com/levisc8/ipmr/actions)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![CRAN
status](https://www.r-pkg.org/badges/version/ipmr)](https://cran.r-project.org/package=ipmr)
[![Codecov test
coverage](https://codecov.io/gh/levisc8/ipmr/branch/master/graph/badge.svg)](https://codecov.io/gh/levisc8/ipmr?branch=master)

# ipmr

`ipmr` is a package for implementing Integral Projection Models (IPMs)
in *R*. It relies heavily on the mathematical syntax of the models, and
does not try to abstract over the process of fitting vital rates. It is
now relatively stable, though tweaks may be made. Below is a brief
overview of how `ipmr` classifies different model types followed by
examples of how to implement those types in this framework.

## Installation

`ipmr` is not yet on CRAN. You can install the development version with
the snippet below:

``` r
if(!require('remotes', quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("levisc8/ipmr", build_vignettes = TRUE)
```

## Package scope

`ipmr` is intended to assist with IPM implementation and analysis. It is
important to note that this package **will not** help with the process
of fitting regression models for vital rates at all! That is a
sufficiently different (and vast) topic that we decided it was not
within the scope of this project. This will only help you turn those
regression models into an IPM without shooting yourself in the foot.
Furthermore, most of the documentation assumes a basic knowledge of IPM
theory and construction. For those that are totally new to IPMs, it is
strongly recommended to read a theoretical overview of the models first.
Some favorites are [Easterling et
al. 2000](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/0012-9658%282000%29081%5B0694%3ASSSAAN%5D2.0.CO%3B2),
[Merow et al. 2013](https://doi.org/10.1111/2041-210X.12146), and [Rees
et al. 2014](https://doi.org/10.1111/1365-2656.12178). The [Introduction
to ipmr
vignette](https://levisc8.github.io/ipmr/articles/ipmr-introduction.html)
also contains a brief overview of IPM theory as well as a far more
detailed introduction to this package. Thus, everything that follows
assumes you have basic understanding of IPM theory, parameterized vital
rate models, and are now ready to begin implementing your IPM.

Below is a brief overview of the package and some examples of how to
implement models with it. A more thorough introduction is available
[here](https://levisc8.github.io/ipmr/articles/ipmr-introduction.html).

## Model classes

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
    (multivariate joint) distribution to sample parameters from is not
    desirable/needed. These models can be a bit more computationally
    efficient than the `param` alternative because all kernels can be
    constructed before the iteration procedure begins, as opposed to
    requiring reconstruction for every single iteration.

-   B. **param**: This describes an IPM with parameters that are
    re-sampled from some distribution at each iteration of the model.
    This could, for example, be a multivariate normal defined by
    co-varying slopes and intercepts, or posterior distribution(s) from
    a Bayesian model. All that is required is that the parameters/values
    for the distribution are specified and that the function that
    generates the parameters at each iteration returns named lists that
    correspond to the parameter names in the model. These are covered in
    greater depth in towards the end of [this
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

-   Simple, density dependent models: **Completed, likely stable**

1.  `"simple_dd_det"`

2.  `"simple_dd_stoch_kern"`

3.  `"simple_dd_stoch_param"`

-   General, density independent models: **Completed and ready**

1.  `"general_di_det"`

2.  `"general_di_stoch_kern"`

3.  `"general_di_stoch_param"`

-   General, density dependent models: **Completed, likely stable**

1.  `"general_dd_det"`

2.  `"general_dd_stoch_kern"`

3.  `"general_dd_stoch_param"`

Simple density-independent deterministic, simple kernel-resampled
stochastic, and simple parameter resampled stochastic models
(`simple_di_det`, `simple_di_stoch_kern`, `simple_di_stoch_param`) are
described in detail
[here](https://levisc8.github.io/ipmr/articles/ipmr-introduction.html).
The `general_*` versions of these are described
[here](https://levisc8.github.io/ipmr/articles/general-ipms.html).
Density dependent versions are completed for simple and general models,
and are probably stable, but have not been tested enough to be certain.
A very brief, though incomplete introduction is available
[here](https://levisc8.github.io/ipmr/articles/density-dependence.html).
Below is an example implementing a `simple_di_det` IPM.

## Quick example of a simple, deterministic IPM

Here is a simple model implemented with `ipmr`. This is a hypothetical
plant species where plants can survive and grow (*P*(*z*′, *z*)), and
reproduce sexually (*F*(*z*′, *z*)). We’ll use 4 regressions: survival
(*s*(*z*)), growth (*G*(*z*′, *z*)), probability of reproducing
(*f*<sub>*r*</sub>(*z*)), and number of seeds produced conditional on
flowering (*f*<sub>*s*</sub>(*z*)). New recruits will be generated with
a Gaussian distribution, which requires calculating the mean and
standard deviation of new recruits from the data. For simplicity, we’ll
assume there’s no maternal effect on recruit size. First, we’ll write
out the functional forms for each component of the model:

1.  *n*(*z*′, *t* + 1) = ∫<sub>*L*</sub><sup>*U*</sup>*K*(*z*′, *z*)*n*(*z*, *t*)*d**z*

2.  *K*(*z*′, *z*) = *P*(*z*′, *z*) + *F*(*z*′, *z*)

3.  *P*(*z*′, *z*) = *s*(*z*) \* *G*(*z*′, *z*)

4.  *L**o**g**i**t*(*s*(*z*)) = *α*<sub>*s*</sub> + *β*<sub>*s*</sub> \* *z*

5.  *G*(*z*′, *z*) ∼ *N**o**r**m*(*μ*<sub>*g*</sub>, *σ*<sub>*g*</sub>)

6.  *m**u*<sub>*g*</sub> = *α*<sub>*g*</sub> + *β*<sub>*g*</sub> \* *z*

7.  *F*(*z*′, *z*) = *f*<sub>*r*</sub>(*z*) \* *f*<sub>*s*</sub>(*z*) \* *f*<sub>*d*</sub>(*z*′)

8.  *L**o**g**i**t*(*f*<sub>*r*</sub>(*z*)) = *α*<sub>*f*<sub>*r*</sub></sub> + *β*<sub>*f*<sub>*r*</sub></sub> \* *z*

9.  *L**o**g*(*f*<sub>*s*</sub>(*z*)) = *α*<sub>*f*<sub>*s*</sub></sub> + *β*<sub>*f*<sub>*s*</sub></sub> \* *z*

10. *f*<sub>*d*</sub>(*z*′) ∼ *N**o**r**m*(*μ*<sub>*f*<sub>*d*</sub></sub>, *σ*<sub>*f*<sub>*d*</sub></sub>)

Equation 1 describes how all the vital rates act on the initial trait
distribution to produce a new one at *t* + 1. Equations 3-6 describe how
existing individuals can survive, and if they survive, grow or shrink.
Equations 7-10 describe how existing individuals create new individuals.
In order to implement this, we need to fit regression models to our
data. The following set of generalized linear models are equivalent to
the functional forms described above:

1.  Survival (*s*(*z*) / `s`): a generalized linear model w/ a logit
    link.

    -   Example model formula:
        `glm(surv ~ size_1, data = my_surv_data, family = binomial())`

2.  Growth (*G*(*z*′, *z*) / `g`): a linear model with a Normal error
    distribution.

    -   Example model formula:
        `lm(size_2 ~ size_1, data = my_grow_data)`

3.  Pr(flowering) (*f*<sub>*r*</sub>(*z*) / `f_r`): a generalized linear
    model w/ a logit link.

    -   Example model formula:
        `glm(flower ~ size_1, data = my_repro_data, family = binomial())`

4.  Seed production (*f*<sub>*s*</sub>(*z*) / `f_s`): a generalized
    linear model w/ log link.

    -   Example model formula:
        `glm(seeds ~ size_1, data = my_flower_data, family = poisson())`

5.  Recruit size distribution (*f*<sub>*d*</sub>(*z*′) / `f_d`): a
    normal distribution w parameters `mu_fd` (mean) and `sd_fd`
    (standard deviation).

    -   Example computations:

        -   `mu_fd = mean(seedling_data$size_2)`

        -   `sd_fd = sd(seedling_data$size_2)`

The example below assumes we’ve already fit our vital rate models from
the raw data. In this example, the numbers are made up, but code that
extracts the values you need from a real regression model is provided in
comments.

``` r
# Load ipmr and get the parameter values. The data_list argument for define_kernel
# should hold every regression parameter and every constant used in the model.

library(ipmr)

data_list = list(s_int     = 2.2,   # coefficients(my_surv_mod)[1]
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
  
  # The formula for the kernel. We dont need to tack on the "z'/z"s here.  
  
  formula   = s * g,
  
  # A named set of expressions for the vital rates it includes. 
  # note the use of user-specified functions here. Additionally, each 
  # state variable has a stateVariable_1 and stateVariable_2, corresponding to
  # z and z' in the equations above. We don't need to define these variables ourselves,
  # just reference them correctly based on the way we've set up our model on paper.
  
  # Perform the inverse logit transformation to get survival probabilities
  # from your model. plogis from the "stats" package does this for us. 

  s         = plogis(s_int + s_slope * dbh_1), 
  
  # The growth model requires a function to compute the mean as a function of dbh.
  # The SD is a constant, so we don't need to define that in ... expression, 
  # just the data_list.
  
  g         = dnorm(dbh_2, mu_g, sd_g),
  mu_g      = g_int + g_slope * dbh_1,
  
  
  # Specify the constant parameters in the model in the data_list. 
  
  data_list = data_list,
  states    = list(c('dbh')),
  
  # If you want to correct for eviction, set evict_cor = TRUE and specify an
  # evict_fun. ipmr provides truncated_distributions() to help. This function
  # takes 2 arguments - the type of distribution, and the name of the parameter/
  # vital rate that it acts on.
  
  evict_cor = TRUE,
  evict_fun = truncated_distributions('norm',
                                      'g')
  ) 

my_simple_ipm <- define_kernel(
  proto_ipm = my_simple_ipm,
  name      = 'F_simple',
  formula   = f_r * f_s * f_d,
  family    = 'CC',
  
  # Inverse logit transformation for flowering probability
  # (because we used a logistic regression)
  
  f_r       = plogis(f_r_int + f_r_slope * dbh_1),
  
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
  
  evict_cor = TRUE,
  evict_fun = truncated_distributions('norm',
                                      'f_d')
) 

# Next, we have to define the implementation details for the model. 
# We need to tell ipmr how each kernel is integrated, what state
# it starts on (i.e. z from above), and what state
# it ends on (i.e. z' above). In simple_* models, state_start and state_end will 
# always be the same, because we only have a single continuous state variable. 
# General_* models will be more complicated.

my_simple_ipm <- define_impl(
  proto_ipm = my_simple_ipm,
  make_impl_args_list(
    kernel_names = c("P_simple", "F_simple"),
    int_rule     = rep("midpoint", 2),
    state_start  = rep("dbh", 2),
    state_end    = rep("dbh", 2)
  )
) 

my_simple_ipm <- define_domains(
  proto_ipm = my_simple_ipm,
  dbh = c(0, # the first entry is the lower bound of the domain.
          50, # the second entry is the upper bound of the domain.
          100 # third entry is the number of meshpoints for the domain.
  ) 
) 

# Next, we define the initial state of the population. We must do this because
# ipmr computes everything through simulation, and simulations require a 
# population state.

my_simple_ipm <- define_pop_state(
  proto_ipm = my_simple_ipm,
  n_dbh = runif(100)
)

my_simple_ipm <- make_ipm(proto_ipm = my_simple_ipm)


lambda_ipmr <- lambda(my_simple_ipm)
w_ipmr      <- right_ev(my_simple_ipm)
v_ipmr      <- left_ev(my_simple_ipm)
```

## More complicated models

Examples of more complicated models are included in the vignettes,
accessible using either `browseVignettes('ipmr')` or by visiting the
Articles tab on [project’s webpage](https://levisc8.github.io/ipmr/).
Please file all bug reports in the Issues tab of this repository or
contact me via [email](mailto:levisc8@gmail.com) with a reproducible
example.

## Code of Conduct

We welcome contributions from other developers. Please note that the
ipmr project is released with a [Contributor Code of
Conduct](https://levisc8.github.io/ipmr/CODE_OF_CONDUCT.html). By
contributing to this project, you agree to abide by its terms.
