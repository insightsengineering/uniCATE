test_that("fold function returns a vector of estimated biomarker coefficients
          and a tibble of efficient influence function values for each
          observation in the validation set",
{

  library(dplyr)
  library(sl3)
  library(origami)
  set.seed(967471)

  # prepare some mock data

  # define baseline characteristics
  n <- 500
  treat <- rbinom(n, 1, 0.5)
  biom1 <- runif(n, min = 2, max = 6)
  biom2 <- rnorm(n, mean = 10, sd = sqrt(10))

  # define hazard functions
  cond_surv_hazard <- function(t, treat, biom1, biom2) {
    (t < 9) / (1 + exp(-(-1 - 0.75*treat + 0.3*biom1 - biom2*treat))) +
      (t == 9)
  }
  cond_cens_hazard <- function(t, treat, biom1, biom2) 0.15

  # generate the failure events for t = 1 to 9
  failure_time <- sapply(
    seq_len(n),
    function(obs) {
      failure_time <- NA
      for (t in 1:9) {
        prob <- cond_surv_hazard(t, treat[obs], biom1[obs], biom2[obs])
        status <- rbinom(1, 1, prob)
        if (status == 1) {
          failure_time <- t
          break
        }
      }
      return(failure_time)
    })

  # generate the censoring events for t = 1 to 9
  censor_time <- sapply(
    seq_len(n),
    function(obs) {
      censor_time <- NA
      for (t in 1:9) {
        prob <- cond_cens_hazard(t, treat[obs], biom1[obs], biom2[obs])
        status <- rbinom(1, 1, prob)
        if (status == 1) {
          censor_time <- t
          break
        }
      }
      if (is.na(censor_time)) censor_time <- 10
      return(censor_time)
    })

  status_df <- tibble(
    "failure_time" = as.integer(failure_time),
    "censor_time" = as.integer(censor_time)
  ) %>%
    transmute(
      time = if_else(censor_time < failure_time, censor_time, failure_time),
      censor = if_else(censor_time < failure_time, 1, 0),
      failure = if_else(censor_time >= failure_time, 1, 0)
    )

  # assemble the data
  treat <- if_else(treat == 1, "treatment", "control")
  data <- tibble(treat, biom1, biom2) %>% bind_cols(status_df)

  # transform into long data
  long_data <- prep_long_data(
    data = data,
    failure = "failure",
    censor = "censor",
    relative_time = "time",
    treatment = "treat",
    covariates = c("biom1", "biom2"),
    biomarkers = c("biom1", "biom2"),
    time_cutoff = 6
  )

  # create the super learners
  interactions <- lapply(c("biom1", "biom2"), function(b) c(b, "treat"))
  lrnr_interactions <- sl3::Lrnr_define_interactions$new(interactions)
  lrnr_glm <- sl3::make_learner(
    sl3::Pipeline, lrnr_interactions, sl3::Lrnr_glm$new()
  )
  meta_learner <- sl3::make_learner(
    sl3::Lrnr_solnp,
    loss_function = sl3::loss_loglik_binomial,
    learner_function = sl3::metalearner_logistic_binomial
  )
  failure_lrnr_sl <- Lrnr_sl$new(
    learners = make_learner(
      Stack, Lrnr_mean$new(), lrnr_glm
    ),
    metalearner = meta_learner
  )
  censor_lrnr_sl <- Lrnr_sl$new(
    learners = make_learner(
      Stack, Lrnr_mean$new(), lrnr_glm
    ),
    metalearner = meta_learner
  )

  # create the fold to test on
  fold <- make_folds(
    long_data,
    fold_fun = folds_vfold,
    V = 2,
    cluster_ids = long_data$observation_id
  )[[1]]

  # apply the fold function
  res_ls <- hold_out_calculation_survival(
    fold,
    long_data = long_data,
    failure = "failure",
    censor = "censor",
    treatment = "treat",
    biomarkers = c("biom1", "biom2"),
    cond_surv_haz_super_learner = failure_lrnr_sl,
    cond_censor_haz_super_learner = censor_lrnr_sl,
    propensity_score_ls = list("treatment" = 0.5, "control" = 0.5)
  )

  # make sure that the coefficients vector is a numeric vector
  expect_equal(class(res_ls$beta_coefs), "numeric")
  expect_length(res_ls$beta_coefs, 2)

  # make sure that the IC table is a tibble of appropriate dimensions
  num_obs_v <- long_data[fold$validation_set, ] %>%
    pull(observation_id) %>%
    unique() %>%
    length()
  expect_equal(nrow(res_ls$ic_df), num_obs_v)
  expect_equal(ncol(res_ls$ic_df), 2)

  # make sure that the column means of the IC table are approximately equal to 0
  expect_equal(as.vector(colMeans(res_ls$ic_df)), c(0, 0))

})
