test_that("A wide data.frame of 3 rows with 3 unique relative times produces a
          tibble with 6 rows",
{
  # define the dummy data
  data <- data.frame(
    status = c(0, 1, 1),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3)
  )

  # generate the long data format
  long_data <- prep_long_data(
    data = data,
    status = "status",
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
          unique time",
{
  # define the dummy data
  data <- tibble(
    status = c(0, 1, 1),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3)
  )

  # generate the long data format
  long_data <- prep_long_data(
    data = data,
    status = "status",
    relative_time = "time",
    treatment = "treat",
    covariates = c("biom1", "biom2", "biom3"),
    biomarkers = c("biom1", "biom2", "biom3")
  )

  # check that the observation with the latest relative time is represented at
  # each unique relative time in the dataset
  unique_times <- data %>% pull(time) %>% unique() %>% sort()
  obs_3_times <- long_data %>%
    filter(observation_id == 3) %>%
    pull(time)
  expect_equal(obs_3_times, unique_times)
})

test_that("Only data.frame or tibble objects are accepted by the data argument",
{
  # define the dummy data
  data <- tibble(
    status = c(0, 1, 1),
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
      status = "status",
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
      status = "status",
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
      status = "status",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "data should be a data.frame or tibble object"
  )

})

test_that("Errors occure when the status, relative_time, treatment, covariates
          and biomarkers arguments are not characters",
{
  # define the dummy data
  data <- tibble(
    status = c(0, 1, 1),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3)
  )

  # error when status is a vector
  expect_error(
    prep_long_data(
      data = data,
      status = data$status,
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "status argument should be a character"
  )

  # error when relative_time is a vector
  expect_error(
    prep_long_data(
      data = data,
      status = "status",
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
      status = "status",
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
      status = "status",
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
      status = "status",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3", "biom4")
    ),
    "biomarkers vector is not a subset of the covariates vector"
  )

})

test_that("Errors occur when variable arguments are missing from the data",
{
  # define the dummy data
  data <- tibble(
    status = c(0, 1, 1),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3)
  )

  # error for missing status variable
  expect_error(
    prep_long_data(
      data = data,
      status = "bad_status",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "status variable is missing from the data"
  )

  # error for missing relative time variable
  expect_error(
    prep_long_data(
      data = data,
      status = "status",
      relative_time = "bad_time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "relative_time variable is missing from the data"
  )

  # error for missing treatment variable
  expect_error(
    prep_long_data(
      data = data,
      status = "status",
      relative_time = "time",
      treatment = "bad_treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "treatment variable is missing from the data"
  )

  # error for missing covariates variables
  expect_error(
    prep_long_data(
      data = data,
      status = "status",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3", "biom4"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "some covariates are missing from the data"
  )

})

test_that("An error is reported when the status variable is not a binary numeric
          variable in the data argument",
{
  # define the dummy data
  data <- tibble(
    status = c("C", "F", "F"),
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
      status = "status",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "status argument should correspond to a numeric, binary variable in the data"
  )
})

test_that("An error is reported when the treatment variable is not a binary
          variable in the data argument",
{
  # define the dummy data
  data <- tibble(
    status = c(1, 0, 0),
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
      status = "status",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "treatment argument should correspond to a binary variable in the data"
  )
})

test_that("Treatment variable is coerced to a factor if not already",
{
  # define the dummy data
  data <- tibble(
    status = c(1, 0, 0),
    time = c(1, 2, 3),
    treat = c("t", "t", "c"),
    biom1 = seq_len(3),
    biom2 = seq_len(3),
    biom3 = seq_len(3)
  )

  # produce long format data
  long_data <- prep_long_data(
    data = data,
    status = "status",
    relative_time = "time",
    treatment = "treat",
    covariates = c("biom1", "biom2", "biom3"),
    biomarkers = c("biom1", "biom2", "biom3")
  )

  # check the class of the treatment variable
  expect_s3_class(long_data$treat, "factor")
})

test_that("An error is reported when relative_time variable is not a positive
          numeric",
{
  # define the dummy data
  data <- tibble(
    status = c(1, 0, 0),
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
      status = "status",
      relative_time = "time",
      treatment = "treat",
      covariates = c("biom1", "biom2", "biom3"),
      biomarkers = c("biom1", "biom2", "biom3")
    ),
    "relative_time argument's corresponding variable should be a positive numeric variable"
  )
})

test_that("Only relevant variables remain after data prep",
{
  # define the dummy data
  data <- tibble(
    status = c(1, 0, 0),
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
    status = "status",
    relative_time = "time",
    treatment = "treat",
    covariates = c("biom1", "biom2", "biom3"),
    biomarkers = c("biom1", "biom2", "biom3")
  )

  # check that unwanted_cov is not in long_data
  expect_false("unwanted_cov" %in% colnames(long_data))
})
