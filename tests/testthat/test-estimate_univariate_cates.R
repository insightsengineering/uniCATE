test_that("fold function returns a vector of estimated biomarker coefficients
          and a tibble of influence curve calculations for each observation in
          the validation set (continuous outcomes)", {

  library(dplyr)
  library(sl3)
  library(origami)
  set.seed(61345)

  # prepare some mock data
  data <- mtcars %>%
    dplyr::mutate(am = factor(am)) %>%
    dplyr::select(mpg, am, disp, hp, wt)
  data <- bind_rows(data, data, data, data)
  outcome <- "mpg"
  treatment <- "am"
  biomarkers <- c("hp", "wt")
  propensity_score_ls <- list("1" = 0.4, "0" = 0.6)

  # create the super learner
  interactions <- lapply(biomarkers, function(b) c(b, treatment))
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

  # create the fold to test on
  fold <- make_folds(data, fold_fun = folds_vfold, V = 5)[[1]]

  # apply the fold function
  res_ls <- hold_out_calculation(
    fold, data, outcome, treatment, biomarkers, super_learner = lrnr_sl,
    propensity_score_ls, outcome_type = "continous"
  )

  # make sure that the coefficients vector is a numeric vector
  expect_equal(class(res_ls$beta_coefs), "numeric")
  expect_length(res_ls$beta_coefs, 2)

  # make sure that the IC table is a tibble of appropriate dimensions
  expect_equal(nrow(res_ls$ic_df), length(fold$validation_set))
  expect_equal(ncol(res_ls$ic_df), 2)

  # make sure that the column means of the IC table are approximately equal to 0
  expect_equal(as.vector(colMeans(res_ls$ic_df)), c(0, 0))
})

test_that("fold function returns a vector of estimated biomarker coefficients
          and a tibble of influence curve calculations for each observation in
          the validation set (binary outcome)", {

  library(dplyr)
  library(sl3)
  library(origami)
  set.seed(61345)

  # prepare some mock data
  data <- mtcars %>%
    dplyr::mutate(vs = factor(vs)) %>%
    dplyr::select(mpg, am, vs, hp, wt)
  data <- bind_rows(data, data, data, data)
  outcome <- "am"
  treatment <- "vs"
  biomarkers <- c("hp", "wt")
  propensity_score_ls <- list("1" = 0.4, "0" = 0.6)

  # create the super learner
  interactions <- lapply(biomarkers, function(b) c(b, treatment))
  lrnr_interactions <- sl3::Lrnr_define_interactions$new(interactions)
  lrnr_glm <- sl3::make_learner(
    sl3::Pipeline, lrnr_interactions, sl3::Lrnr_glm$new()
  )
  meta_learner <- sl3::make_learner(
    sl3::Lrnr_solnp,
    loss_function = sl3::loss_loglik_binomial,
    learner_function = sl3::metalearner_logistic_binomial
  )
  lrnr_sl <- Lrnr_sl$new(
    learners = make_learner(
      Stack, Lrnr_mean$new(), lrnr_glm
    ),
    metalearner = meta_learner
  )

  # create the fold to test on
  fold <- make_folds(data, fold_fun = folds_vfold, V = 5)[[1]]

  # apply the fold function
  res_ls <- hold_out_calculation(
    fold, data, outcome, treatment, biomarkers, super_learner = lrnr_sl,
    propensity_score_ls, outcome_type = "binomial"
  )

  # make sure that the coefficients vector is a numeric vector
  expect_equal(class(res_ls$beta_coefs), "numeric")
  expect_length(res_ls$beta_coefs, 2)

  # make sure that the IC table is a tibble of appropriate dimensions
  expect_equal(nrow(res_ls$ic_df), length(fold$validation_set))
  expect_equal(ncol(res_ls$ic_df), 2)

  # make sure that the column means of the IC table are approximately equal to 0
  expect_equal(as.vector(colMeans(res_ls$ic_df)), c(0, 0))
})


test_that("estimate_univariate_cates() returns a vector with estimate lm
          coefficients and a table of influence curves",
{
  library(dplyr)
  library(sl3)
  library(origami)
  set.seed(12312)

  # prepare some mock data
  data <- mtcars %>%
    dplyr::mutate(am = factor(am)) %>%
    dplyr::select(mpg, am, disp, hp, wt)
  data <- bind_rows(data, data, data, data)
  outcome <- "mpg"
  treatment <- "am"
  biomarkers <- c("hp", "wt")
  propensity_score_ls <- list("1" = 0.4, "0" = 0.6)

  # create the super learner
  interactions <- lapply(biomarkers, function(b) c(b, treatment))
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

  # compute the holdout potential outcome difference dataset
  res_ls <- estimate_univariate_cates(
    data, outcome, treatment, biomarkers, super_learner = lrnr_sl,
    propensity_score_ls, v_folds = 2L, parallel = FALSE
  )

  # check that the beta coefs table has the correct dimensions
  expect_equal(dim(res_ls$betas_df), c(2, 2))

  # check that the IC table has the correct dimensions
  expect_equal(dim(res_ls$ic_df), c(nrow(data), 2))

  # check that the column means of the IC table are zero
  expect_equal(as.vector(colMeans(res_ls$ic_df)), c(0, 0))

})
