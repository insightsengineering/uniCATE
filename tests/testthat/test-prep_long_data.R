
test_that("A wide table of 3 rows with 3 unique relative times produces a table
          with 6 rows",
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

  # check that 6 rows are generated
  expect_equal(nrow(long_data), 6)
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
