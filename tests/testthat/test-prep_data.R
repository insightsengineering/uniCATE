library(tibble)

test_that("only accepts data.frames and tibbles", {
  expect_error(
    prep_data(as.matrix(mtcars), "mpg", "am", "hp", "hp"),
    "data should be a data.frame or tibble object"
  )
  expect_silent(
    prep_data(mtcars, "mpg", "am", "hp", "hp")
  )
  expect_silent(
    prep_data(
      as_tibble(mtcars, .name_repair = "minimal"), "mpg", "am", "hp", "hp"
    )
  )
})

test_that("checks that outcome is a continuous or binary numeric variable", {
  library(dplyr)
  data <- mtcars %>%
    mutate(vs = factor(vs))
  expect_error(
    prep_data(data, "vs", "am", "hp", "hp"),
    "outcome should be a continuous or binary numeric variable"
  )
  data <- mtcars
  data$mpg <- as.character(data$mpg)
  expect_error(
    prep_data(data, "mpg", "am", "hp", "hp"),
    "outcome should be a continuous or binary numeric variable"
  )
  expect_silent(
    prep_data(mtcars, "vs", "am", "hp", "hp")
  )
})

test_that("checks that treatment is a binary variable", {
  expect_error(
    prep_data(mtcars, "mpg", "cyl", "hp", "hp"),
    "treatment should be a binary variable"
  )
})

test_that("outcome and treatment designators are of length 1", {
  expect_error(
    prep_data(mtcars, c("mpg", "wt"), "cyl", "hp", "hp"),
    "outcome should be a length-one character"
  )
  expect_error(
    prep_data(mtcars, "mpg", c("am", "wt"), "hp", "hp"),
    "treatment should be a length-one character"
  )
})


test_that("specified variables are contained in the data", {
  expect_error(
    prep_data(mtcars, "mpgg", "am", "hp", "hp"),
    "outcome is not contained in data"
  )
  expect_error(
    prep_data(mtcars, "mpg", "amm", "hp", "hp"),
    "treatment is not contained in data"
  )
  expect_error(
    prep_data(mtcars, "mpg", "am", c("disp", "horsepower"), "disp"),
    "not all covariates are contained in data"
  )
  expect_error(
    prep_data(mtcars, "mpg", "am", c("disp", "hp"), c("disp", "hp", "wt")),
    "biomarkers is not a subset of covariates"
  )
})

test_that("if the treatment variable is not a factor, it is turned into one", {

  # when passing in a data.frame
  mod_mtcars <- prep_data(mtcars, "mpg", "am", "hp", "hp")
  expect_true(is.factor(mod_mtcars$am))
  expect_true(all(levels(mod_mtcars$am) %in% c("1", "0")))

  # when passing in a tibble
  mod_mtcars <- prep_data(as_tibble(mtcars), "mpg", "am", "hp", "hp")
  expect_true(is.factor(mod_mtcars$am))
  expect_true(all(levels(mod_mtcars$am) %in% c("1", "0")))

})

test_that("always returns a tibble with only the outcome, treatment, and
          covariates", {
            data <- prep_data(mtcars, "mpg", "am", "hp", "hp")
            expect_true(is_tibble(data))
            expect_true(ncol(data) == 3)
          })
