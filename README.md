
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `R`/`uniCATE`

> Univarite Conditional Average Treatment Effect Estimation

**Author:** [Philippe Boileau](https://pboileau.ca/)

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- badges: end -->

------------------------------------------------------------------------

`uniCATE` implements an efficient semiparametric procedure for the
estimation of individual variables’ conditional average treatment
effects’ best linear approximations. This method can identify potential
predictive biomarkers in randomized studies using assumption-lean
statistical inference techniques that rely on flexible machine learning
algorithms.

## Installation

The *development version* of the package may be installed from GitHub
using [`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
# install.packages("remotes")
remotes::install_github("insightsengineering/uniCATE")
```

## Example

We simulate a randomized control trial in which there is a heterogeneous
treatment effect for biomarkers 1 and 2. `unicate()` successfully
identifies these biomarkers as effect modifierss.

``` r
# load the required libraries
library(uniCATE)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(sl3)

# set the seed for reproducibility
set.seed(514)

# simulate some randomized control data
n <- 100
data <- tibble("treatment" = rbinom(n, 1, 0.5)) %>%
  mutate(
    bio1 = rnorm(n, mean = 2, sd = 0.2),
    bio2 = rnorm(n, mean = -2, sd = 0.2),
    bio3 = rnorm(n, mean = 0, sd = 0.1),
    bio4 = rnorm(n, mean = 0, sd = 0.1),
    covar = 0.2*rbinom(n, 1, 0.4),
    response = covar + bio1 * treatment + bio2 * treatment
  )

# define the required arguments
covariates <- c("bio1", "bio2", "bio3", "bio4", "covar")
biomarkers <- c("bio1", "bio2", "bio3", "bio4")
propensity_score_ls <- list("1" = 0.5, "0" = 0.5)

# create a simple SuperLearner using a linear model and a random forest
interactions <- lapply(biomarkers, function(b) c(b, "treatment"))
lrnr_interactions <- sl3::Lrnr_define_interactions$new(interactions)
lrnr_glm <- sl3::make_learner(
  sl3::Pipeline, lrnr_interactions, sl3::Lrnr_glm$new()
)
lrnr_sl <- Lrnr_sl$new(
  learners = make_learner(
    Stack, Lrnr_ranger$new(), lrnr_glm
  ),
  metalearner = make_learner(Lrnr_nnls)
)

# apply uniCATE to the simulated data
unicate(
  data,
  outcome = "response",
  treatment = "treatment",
  covariates = covariates,
  biomarkers = biomarkers,
  propensity_score_ls = propensity_score_ls,
  super_learner = lrnr_sl,
  v_folds = 2L
)
#> # A tibble: 4 x 7
#>   biomarker  coef     se     z  p_value p_value_bh p_value_holm
#>   <chr>     <dbl>  <dbl> <dbl>    <dbl>      <dbl>        <dbl>
#> 1 bio1      0.876 0.0988 8.87  7.52e-19   3.01e-18     3.01e-18
#> 2 bio2      0.841 0.115  7.30  2.97e-13   5.94e-13     8.91e-13
#> 3 bio3      0.117 0.255  0.457 6.48e- 1   6.64e- 1     1   e+ 0
#> 4 bio4      0.113 0.260  0.434 6.64e- 1   6.64e- 1     1   e+ 0
```

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/insightsengineering/uniCATE/issues).

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/insightsengineering/uniCATE/blob/master/.github/CONTRIBUTING.md)
prior to submitting a pull request.

## License

The contents of this repository are distributed under the MIT license.
See file
[`LICENSE.md`](https://github.com/insightsengineering/uniCATE/blob/main/LICENSE.md)
for details.

*Note: This license is a placeholder until the legal team is consulted.*
