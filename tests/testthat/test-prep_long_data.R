test_that("A wide data.frame of 3 rows with 3 unique relative times produces a
          tibble with 6 rows", {
  # define the dummy data
  data <- data.frame(
    sick = c(1, 0, 1),
    left = c(0, 1, 0),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3)
  )

  # generate the long data format
  long_data <- prep_long_data(
    data = data,
    failure = "sick",
    censor = "left",
    relative_time = "time",
    treatment = "treat",
    covariates = c("biom1", "biom2", "biom3"),
    biomarkers = c("biom1", "biom2", "biom3")
  )

  # check that 6 rows are generated
  expect_equal(nrow(long_data), 6)

  # check that the transformed data is a tibble object
  expect_s3_class(long_data, "tbl_df")
})

test_that("The observation with latest relative time is represented at every
          unique time", {
  # define the dummy data
  data <- tibble(
    sick = c(1, 0, 1),
    left = c(0, 1, 0),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3)
  )

  # generate the long data format
  long_data <- prep_long_data(
    data = data,
    failure = "sick",
    censor = "left",
    relative_time = "time",
    treatment = "treat",
    covariates = c("biom1", "biom2", "biom3"),
    biomarkers = c("biom1", "biom2", "biom3")
  )

  # check that the observation with the latest relative time is represented at
  # each unique relative time in the dataset
  unique_times <- data %>%
    pull(time) %>%
    unique() %>%
    sort()
  obs_3_times <- long_data %>%
    filter(observation_id == 3) %>%
    pull(time)
  expect_equal(obs_3_times, unique_times)
})

test_that("Only data.frame or tibble objects are accepted by the data argument", {
  # define the dummy data
  data <- tibble(
    sick = c(1, 0, 1),
    left = c(0, 1, 0),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3)
  )

  # accepts tibbles
  expect_silent(
    prep_long_data(
      data = data,
      failure = "sick",
      censor = "left",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    )
  )

  # accepts data.frames
  expect_silent(
    prep_long_data(
      data = data.frame(data),
      failure = "sick",
      censor = "left",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    )
  )

  # refuses matrix objects
  expect_error(
    prep_long_data(
      data = as.matrix(data),
      failure = "sick",
      censor = "left",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "data should be a data.frame or tibble object"
  )
})

test_that("Errors occure when the failure, censor, relative_time, treatment,
           covariates and biomarkers arguments are not characters", {
  # define the dummy data
  data <- tibble(
    sick = c(1, 0, 1),
    left = c(0, 1, 0),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3)
  )


  # error when failure argument is a vector
  expect_error(
    prep_long_data(
      data = data,
      failure = data$sick,
      censor = "left",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "failure argument should be a character"
  )

  # error when failure argument is a vector
  expect_error(
    prep_long_data(
      data = data,
      failure = "sick",
      censor = data$left,
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "censor argument should be a character"
  )

  # error when relative_time is a vector
  expect_error(
    prep_long_data(
      data = data,
      failure = "sick",
      censor = "left",
      relative_time = data$time,
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "relative_time argument should be a character"
  )

  # error when treatment is a vector
  expect_error(
    prep_long_data(
      data = data,
      failure = "sick",
      censor = "left",
      relative_time = "time",
      treatment = data$treat,
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "treatment argument should be a character"
  )

  # error when covariates is not a character
  expect_error(
    prep_long_data(
      data = data,
      failure = "sick",
      censor = "left",
      relative_time = "time",
      treatment = "treat",
      covariates = c(1, 2),
      biomarkers = c(1, 2)
    ),
    "covariates argument should be a character vector"
  )

  # error for biomarker variables not contained in covariates
  expect_error(
    prep_long_data(
      data = data,
      failure = "sick",
      censor = "left",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3", "biom4")
    ),
    "biomarkers vector is not a subset of the covariates vector"
  )
})

test_that("Errors occur when variable arguments are missing from the data", {
  # define the dummy data
  data <- tibble(
    sick = c(1, 0, 1),
    left = c(0, 1, 0),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3)
  )

  # error for missing failure variable
  expect_error(
    prep_long_data(
      data = data,
      failure = "failure",
      censor = "left",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "failure argument's variable is missing from the data"
  )

  # error for missing failure variable
  expect_error(
    prep_long_data(
      data = data,
      failure = "sick",
      censor = "censor",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "censor argument's variable is missing from the data"
  )

  # error for missing relative time variable
  expect_error(
    prep_long_data(
      data = data,
      failure = "sick",
      censor = "left",
      relative_time = "bad_time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "relative_time argument's variable is missing from the data"
  )

  # error for missing treatment variable
  expect_error(
    prep_long_data(
      data = data,
      failure = "sick",
      censor = "left",
      relative_time = "time",
      treatment = "bad_treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "treatment argument's variable is missing from the data"
  )

  # error for missing covariates variables
  expect_error(
    prep_long_data(
      data = data,
      failure = "sick",
      censor = "left",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3", "biom4"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "some covariates are missing from the data"
  )
})

test_that("An error is reported when the failure variable is not a
           binary numeric variables in the data argument", {
  # define the dummy data
  data <- tibble(
    sick = c("1", "0", "1"),
    left = c(0, 1, 0),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3)
  )

  # refuses non-binary, non-numeric status variable
  expect_error(
    prep_long_data(
      data = data,
      failure = "sick",
      censor = "left",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "failure argument should correspond to a numeric, binary variable in the data"
  )
})

test_that("An error is reported when the censor variable is not a
           binary numeric variables in the data argument", {
  # define the dummy data
  data <- tibble(
    sick = c(1, 0, 1),
    left = c("0", 1, "0"),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3)
  )

  # refuses non-binary, non-numeric status variable
  expect_error(
    prep_long_data(
      data = data,
      failure = "sick",
      censor = "left",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "censor argument should correspond to a numeric, binary variable in the data"
  )
})

test_that("An error is thrown for observations with a failure and a censoring
          event", {
  data <- tibble(
    sick = c(1, 0, 1),
    left = c(1, 1, 0),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3)
  )

  expect_error(
    prep_long_data(
      data = data,
      failure = "sick",
      censor = "left",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "observations may not have a both failure and a censoring event"
  )
})

test_that("An error is reported when the treatment variable is not a binary
          variable in the data argument", {
  # define the dummy data
  data <- tibble(
    sick = c(1, 0, 0),
    left = c(0, 1, 0),
    time = c(1, 2, 3),
    treat = c("t", "u", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3)
  )

  # refuses non-binary treatment variable
  expect_error(
    prep_long_data(
      data = data,
      failure = "sick",
      censor = "left",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "treatment argument should correspond to a binary variable in the data"
  )
})

test_that("Treatment variable is coerced to a factor if not already", {
  # define the dummy data
  data <- tibble(
    sick = c(1, 0, 0),
    left = c(0, 1, 0),
    time = c(1, 2, 3),
    treat = c("t", "u", "t"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3)
  )

  # produce long format data
  long_data <- prep_long_data(
    data = data,
    failure = "sick",
    censor = "left",
    relative_time = "time",
    treatment = "treat",
    covariates = c("biom1", "biom2", "biom3"),
    biomarkers = c("biom1", "biom2", "biom3")
  )

  # check the class of the treatment variable
  expect_s3_class(long_data$treat, "factor")
})

test_that("An error is reported when relative_time variable is not a positive
          numeric", {
  # define the dummy data
  data <- tibble(
    sick = c(1, 0, 0),
    left = c(0, 1, 0),
    time = c(-1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3)
  )

  # refuses relative time variable with negative value
  expect_error(
    prep_long_data(
      data = data,
      failure = "sick",
      censor = "left",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "relative_time argument's corresponding variable should be a positive numeric variable"
  )
})

test_that("The time cutoff is a numeric if non-null", {
  # define the dummy data
  data <- tibble(
    sick = c(1, 0, 0),
    left = c(0, 1, 0),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3),
    unwanted_cov = seq_len(3),
  )

  # refuses time_cutoff since it doesn't fall within the provided data range
  expect_error(
    prep_long_data(
      data = data,
      failure = "sick",
      censor = "left",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3"),
      time_cutoff = "4"
    ),
    "time_cutoff should be a single numeric value when specified"
  )
})

test_that("Warn when time cutoff not within the range of relative time", {
  # define the dummy data
  data <- tibble(
    sick = c(1, 0, 0),
    left = c(0, 1, 0),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3),
    unwanted_cov = seq_len(3),
  )

  # refuses time_cutoff since it doesn't fall within the provided data range
  expect_message(
    prep_long_data(
      data = data,
      failure = "sick",
      censor = "left",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3"),
      time_cutoff = 4
    ),
    "time_cutoff is larger than or equal to the largest value in relative_time's corresponding variable: 3"
  )
})

test_that("Only relevant variables remain after data prep", {
  # define the dummy data
  data <- tibble(
    sick = c(1, 0, 0),
    left = c(0, 1, 0),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3),
    unwanted_cov = seq_len(3),
  )

  # prepare the data
  long_data <- prep_long_data(
    data = data,
    failure = "sick",
    censor = "left",
    relative_time = "time",
    treatment = "treat",
    covariates = c("biom1", "biom2", "biom3"),
    biomarkers = c("biom1", "biom2", "biom3"),
    time_cutoff = 2
  )

  # check that unwanted_cov is not in long_data
  expect_false("unwanted_cov" %in% colnames(long_data))
})

test_that("Obs 2 doesn't have a censoring event, obs 3 doesn't have any events", {
  # define the dummy data
  data <- tibble(
    sick = c(1, 0, 0),
    left = c(0, 1, 0),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3),
    unwanted_cov = seq_len(3),
  )

  # prepare the data
  long_data <- prep_long_data(
    data = data,
    failure = "sick",
    censor = "left",
    relative_time = "time",
    treatment = "treat",
    covariates = c("biom1", "biom2", "biom3"),
    biomarkers = c("biom1", "biom2", "biom3"),
    time_cutoff = 2
  )

  # check that the max time is 2
  expect_true(max(long_data$time) == 2)

  # check tha the second observation's censor status is one at t2
  sec_censor_status <- long_data %>%
    filter(observation_id == 2, time == 2) %>%
    pull(left)
  expect_equal(sec_censor_status, 1)

  # check that the third observation's failure and censor status are zero at t2
  third_failure_status <- long_data %>%
    filter(observation_id == 3, time == 2) %>%
    pull(sick)
  expect_equal(third_failure_status, 0)
  third_censor_status <- long_data %>%
    filter(observation_id == 3, time == 2) %>%
    pull(left)
  expect_equal(third_censor_status, 0)
})

test_that("Obs 3 doesn't have a censoring or failure event", {
  # define the dummy data
  data <- tibble(
    sick = c(1, 0, 0),
    left = c(0, 1, 0),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3),
    unwanted_cov = seq_len(3),
  )

  # prepare the data
  long_data <- prep_long_data(
    data = data,
    failure = "sick",
    censor = "left",
    relative_time = "time",
    treatment = "treat",
    covariates = c("biom1", "biom2", "biom3"),
    biomarkers = c("biom1", "biom2", "biom3")
  )

  # check that the third observation's failure and censor status are zero at t3
  third_failure_status <- long_data %>%
    filter(observation_id == 3, time == 3) %>%
    pull(sick)
  expect_equal(third_failure_status, 0)
  third_censor_status <- long_data %>%
    filter(observation_id == 3, time == 3) %>%
    pull(left)
  expect_equal(third_censor_status, 0)
})
