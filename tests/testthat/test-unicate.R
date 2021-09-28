# integration tests
test_that("biomarkers 1 and 2 are found by unicate() (cont. outcome)", {

  library(dplyr)
  library(sl3)
  set.seed(6132451)

  # simulate some data
  n <- 100
  error_terms <- rnorm(n, 0, 0.1)
  data <- tibble("treatment" = rbinom(n, 1, 0.5)) %>%
    mutate(
      bio1 = rnorm(n, mean = 2, sd = 0.2),
      bio2 = rnorm(n, mean = -2, sd = 0.2),
      bio3 = rnorm(n, mean = 0, sd = 0.1),
      bio4 = rnorm(n, mean = 0, sd = 0.1),
      covar = 0.2*rbinom(n, 1, 0.4),
      response = covar + bio1 * treatment + bio2 * treatment + error_terms
    )

  # define the required arguments
  covariates <- c("bio1", "bio2", "bio3", "bio4", "covar")
  biomarkers <- c("bio1", "bio2", "bio3", "bio4")
  propensity_score_ls <- list("1" = 0.5, "0" = 0.5)

  # create the super learner
  interactions <- lapply(biomarkers, function(b) c(b, "treatment"))
  lrnr_interactions <- sl3::Lrnr_define_interactions$new(interactions)
  lrnr_glm <- sl3::make_learner(
    sl3::Pipeline, lrnr_interactions, sl3::Lrnr_glm$new()
  )
  lrnr_sl <- Lrnr_sl$new(
    learners = make_learner(
      Stack, Lrnr_mean$new(), lrnr_glm
    ),
    metalearner = make_learner(Lrnr_nnls)
  )

  # apply the procedure
  set.seed(71345)
  results <- unicate(
    data,
    outcome = "response",
    treatment = "treatment",
    covariates = covariates,
    biomarkers = biomarkers,
    propensity_score_ls = propensity_score_ls,
    super_learner = lrnr_sl,
    v_folds = 2L
  )

  # check significance
  expect_equal(
    as.vector(results$p_value < 0.05),
    c(TRUE, TRUE, FALSE, FALSE)
  )

  # make sure that the results are similar when run in parallel
  library(future)
  plan(sequential)

  # run unicate() in parallel
  set.seed(71345)
  results_parallel <- unicate(
    data,
    outcome = "response",
    treatment = "treatment",
    covariates = covariates,
    biomarkers = biomarkers,
    propensity_score_ls = propensity_score_ls,
    super_learner = lrnr_sl,
    v_folds = 2L,
    parallel = TRUE
  )

  expect_equal(
    as.vector(results_parallel$p_value < 0.05),
    c(TRUE, TRUE, FALSE, FALSE)
  )
})


test_that("biomarkers 1 and 2 are found by unicate() (binary outcome)", {

  library(dplyr)
  library(purrr)
  set.seed(6132451)

  # simulate some data
  n <- 200
  error_terms <- rnorm(n, 0, 0.1)
  data <- tibble("treatment" = rbinom(n, 1, 0.5)) %>%
    mutate(
      bio1 = rnorm(n, mean = 2, sd = 0.2),
      bio2 = rnorm(n, mean = -2, sd = 0.2),
      bio3 = rnorm(n, mean = 0, sd = 0.1),
      bio4 = rnorm(n, mean = 0, sd = 0.1),
      covar = 0.2*rbinom(n, 1, 0.4),
      response_prob = plogis(
        covar + 10 * bio1 * treatment + 10 * bio2 * treatment + error_terms
        ),
      response = map_dbl(response_prob, function(x) rbinom(1, 1, x))
    )

  # define the required arguments
  covariates <- c("bio1", "bio2", "bio3", "bio4", "covar")
  biomarkers <- c("bio1", "bio2", "bio3", "bio4")
  propensity_score_ls <- list("1" = 0.5, "0" = 0.5)

  # create the super learner
  meta_learner <- sl3::make_learner(
    sl3::Lrnr_solnp,
    loss_function = sl3::loss_loglik_binomial,
    learner_function = sl3::metalearner_logistic_binomial
  )
  lrnr_sl <- Lrnr_sl$new(
    learners = make_learner(
      Stack, Lrnr_mean$new(), Lrnr_ranger$new()
    ),
    metalearner = meta_learner
  )

  # apply the procedure
  set.seed(71345)
  results <- unicate(
    data,
    outcome = "response",
    treatment = "treatment",
    covariates = covariates,
    biomarkers = biomarkers,
    propensity_score_ls = propensity_score_ls,
    super_learner = lrnr_sl,
    v_folds = 2L
  )

  # check significance
  expect_equal(
    as.vector(results$p_value < 0.05),
    c(TRUE, TRUE, FALSE, FALSE)
  )

  # make sure that the results are similar when run in parallel
  library(future)
  plan(sequential)

  # run unicate() in parallel
  set.seed(71345)
  results_parallel <- unicate(
    data,
    outcome = "response",
    treatment = "treatment",
    covariates = covariates,
    biomarkers = biomarkers,
    propensity_score_ls = propensity_score_ls,
    super_learner = lrnr_sl,
    v_folds = 2L,
    parallel = TRUE
  )

  expect_equal(
    as.vector(results_parallel$p_value < 0.05),
    c(TRUE, TRUE, FALSE, FALSE)
  )
})
