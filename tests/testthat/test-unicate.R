# integration tests
test_that("biomarkers 1 and 2 are found by unicate() (cont. outcome)", {

  library(dplyr)
  set.seed(734251)

  # simulate some data
  data <- tibble("treatment" = rbinom(100, 1, 0.5)) %>%
    mutate(
      bio1 = rnorm(100, mean = 2, sd = 0.2),
      bio2 = rnorm(100, mean = -2, sd = 0.2),
      bio3 = rnorm(100, mean = 0, sd = 0.1),
      bio4 = rnorm(100, mean = 0, sd = 0.1),
      covar = 0.2*rbinom(100, 1, 0.4),
      response = covar + bio1 * treatment + bio2 * treatment + rnorm(1, 0, 0.1)
    )

  # define the required arguments
  covariates <- c("bio1", "bio2", "bio3", "bio4", "covar")
  biomarkers <- c("bio1", "bio2", "bio3", "bio4")
  propensity_score_ls <- list("1" = 0.5, "0" = 0.5)

  # apply the procedure
  set.seed(823145)
  results <- unicate(
    data,
    outcome = "response",
    treatment = "treatment",
    covariates = covariates,
    biomarkers = biomarkers,
    propensity_score_ls = propensity_score_ls,
    v_folds = 2L
  )

  # check significance
  expect_equal(
    as.vector(results$p_value < 0.05),
    c(TRUE, TRUE, FALSE, FALSE)
  )

  # make sure that the results are identical when run in parallel
  library(future)
  plan(sequential)

  # run unicate() in parallel
  set.seed(823145)
  results_parallel <- unicate(
    data,
    outcome = "response",
    treatment = "treatment",
    covariates = covariates,
    biomarkers = biomarkers,
    propensity_score_ls = propensity_score_ls,
    v_folds = 2L,
    parallel = TRUE
  )

  expect_equal(results, results_parallel, tolerance = 10e-4)
})
