test_that("column of coefs is the mean of validation fold coefs", {

  # define a fake dataset
  betas_df <- tibble(
    "bio1" = seq_len(10),
    "bio2" = -1 * seq_len(10),
  )
  ic_df <- tibble(
    "bio1" = c(rep(5, 50), rep(-5, 50)),
    "bio2" = c(rep(6, 50), rep(-6, 50)),
  )
  cv_ls <- list("betas_df" = betas_df, "ic_df" = ic_df)

  # make sure that there's a row for each biomarker, and that their estimates
  # are correctly calculated
  expect_equal(length(compute_coefs_tbl(cv_ls)$coef), ncol(betas_df))
  expect_equal(compute_coefs_tbl(cv_ls)$coef, colMeans(betas_df))
})

test_that("column of SEs is the mean of the squared of influence curves, divided
          by the number of observations in the data", {

  # define a fake dataset
  betas_df <- tibble(
    "bio1" = seq_len(10),
    "bio2" = -1 * seq_len(10),
  )
  ic_df <- tibble(
    "bio1" = c(rep(5, 50), rep(-5, 50)),
    "bio2" = c(rep(6, 50), rep(-6, 50)),
  )
  cv_ls <- list("betas_df" = betas_df, "ic_df" = ic_df)

  # make sure that there's a row for each biomarker, and that their estimates
  # are correctly calculated
  expect_equal(length(compute_coefs_tbl(cv_ls)$se), ncol(ic_df))
  expect_equal(
    compute_coefs_tbl(cv_ls)$se, sqrt(colMeans((ic_df^2)) / nrow(ic_df))
  )
})
