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
  n <- 50
  treat <- rbinom(n, 1, 0.5)
  biom1 <- runif(n, min = 2, max = 6)
  biom2 <- rnorm(n, mean = 10, sd = sqrt(10))

  # define hazard functions
  cond_surv_hazard <- function(t, treat, biom1, biom2) {
    (t < 9) / (1 + exp(-(-8 - 0.75*treat + 0.3*biom1^2 + 0.5*biom2*treat))) +
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
  censore_time <- sapply(
    seq_len(n),
    function(obs) {
      censore_time <- NA
      for (t in 1:9) {
        prob <- cond_cens_hazard(t, treat[obs], biom1[obs], biom2[obs])
        status <- rbinom(1, 1, prob)
        if (status == 1) {
          censore_time <- t
          break
        }
      }
      return(censore_time)
    })

  # create status and relative time vector
  status_time_df <- lapply(
    seq_len(n),
    function(obs) {

      if (is.na(censore_time[obs]) & is.na(failure_time[obs])) {
        list("status" = 0, "time" = 9)
      } else if (is.na(censore_time[obs]) & !is.na(failure_time[obs])) {
        list("status" = 1, "time" = failure_time[obs])
      } else if (!is.na(censore_time[obs]) & is.na(failure_time[obs])) {
        list("status" = 0, "time" = censore_time[obs])
      } else {
        if (censore_time[obs] <= failure_time[obs]) {
          list("status" = 0, "time" = censore_time[obs])
        } else {
          list("status" = 1, "time" = failure_time[obs])
        }
      }
    }
  ) %>% bind_rows()

  # assemble the data
  data <- tibble(treat, biom1, biom2) %>% bind_cols(status_time_df)

  # transform into long data
  long_data <- prep_long_data(
    data, "status", "time", "treat", c("biom1", "biom2"), c("biom1", "biom2"),
    6
  )
})
