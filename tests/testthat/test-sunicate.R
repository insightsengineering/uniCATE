test_that("sunicate() identifies the effect modified in a simulated example", {

  library(dplyr)
  library(sl3)
  library(origami)
  set.seed(72345)

  # prepare some mock data

  # define baseline characteristics
  n <- 500
  treat <- rbinom(n, 1, 0.5)
  cov <- rnorm(n)
  biom1 <- rnorm(n)
  biom2 <- rnorm(n)

  # define hazard functions
  cond_surv_hazard <- function(t, treat, biom1, biom2) {
    (t < 9) / (1 + exp(-(-2 - 3*treat*biom1))) +
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
  treat <- factor(treat, levels = c("treatment", "control"))
  data <- tibble(treat, cov, biom1, biom2) %>% bind_cols(status_df)

  # run sunicate
  res_tbl <- sunicate(
    data = data,
    failure = "failure",
    censor = "censor",
    relative_time = "time",
    treatment = "treat",
    covariates = c("cov", "biom1", "biom2"),
    biomarkers = c("biom1", "biom2"),
    time_cutoff = 6,
    propensity_score_ls = list("treatment" = 0.5, "control" = 0.5),
    v_folds = 2L
  )

  # make sure that only biomarker 1 is a significant effect modifier
  biom1_pval <- res_tbl %>%
    filter(biomarker == "biom1") %>%
    pull(p_value)
  expect_true(biom1_pval < 0.05)
  biom2_pval <- res_tbl %>%
    filter(biomarker == "biom2") %>%
    pull(p_value)
  expect_true(biom2_pval > 0.05)

})
