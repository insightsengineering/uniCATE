test_that("can pass in a propensity score list to apply_aiptw_transform()", {
  library(dplyr)

  # prepare some mock data
  data <- mtcars %>%
    mutate(
      am = factor(am),
      y = 1,
      y_hat_treat = 2,
      y_hat_cont = 3
    )
  outcome <- "y"
  treatment <- "am"
  propensity_score_ls <- list("1" = 0.3, "0" = 0.7) # treatment group goes first

  # apply the AIPTW transform
  transform_data <- apply_aiptw_transform(
    data, outcome, treatment, propensity_score_ls,
    outcome_type = "continuous"
  )

  # make sure that the calculations are correct
  transform_data <- transform_data %>%
    dplyr::mutate(
      confirm_treat = dplyr::if_else(
        am == "1", 2 + -1 / 0.3, 2
      ),
      confirm_treat = (confirm_treat == y_aiptw_treat),
      confirm_cont = dplyr::if_else(
        am == "0", 3 + -2 / 0.7, 3
      ),
      confirm_cont = (confirm_cont == y_aiptw_cont)
    ) %>%
    dplyr::select(confirm_treat, confirm_cont)

  # make sure that all checks return true
  expect_true(sum(transform_data) == 2 * nrow(transform_data))
})

test_that("apply_aiptw_transform() truncates predicted outcomes between 0 and 1
          when the outcome is a binary variable", {
  library(dplyr)

  # prepare some mock data
  data <- mtcars %>%
    mutate(
      am = factor(am),
      y = 1,
      y_hat_treat = 2,
      y_hat_cont = 3
    )
  outcome <- "y"
  treatment <- "am"
  propensity_score_ls <- list("1" = 0.3, "0" = 0.7) # treatment group goes first

  # apply the AIPTW transform
  transform_data <- apply_aiptw_transform(
    data, outcome, treatment, propensity_score_ls,
    outcome_type = "binomial"
  )
  expect_equal(sum(transform_data$y_aiptw_treat < 0 |
    transform_data$y_aiptw_treat > 1), 0)
  expect_equal(sum(transform_data$y_aiptw_cont < 0 |
    transform_data$y_aiptw_cont > 1), 0)
})
