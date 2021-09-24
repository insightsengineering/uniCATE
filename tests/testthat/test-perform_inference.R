test_that("assert that p-values are properly calculated", {

  # create some dummy data
  biomarkers_tbl <- tibble::tibble(
    biomarker = c("bio1", "bio2", "bio3"),
    coef = c(1, 0, 1),
    se = c(0.5, 0.2, 5)
  )

  # compute the p-values
  bio1_pval <- 2 * min(pnorm(1/0.5), pnorm(1/0.5, lower.tail = FALSE))
  bio2_pval <- 2 * min(pnorm(0/0.2), pnorm(0/0.2, lower.tail = FALSE))
  bio3_pval <- 2 * min(pnorm(1/5), pnorm(1/5, lower.tail = FALSE))

  # compute the p-values and adjusted p-values of each biomarker, and rank
  res_df <- perform_inference(biomarkers_tbl)
  expect_equal(res_df$biomarker, c("bio1", "bio3", "bio2"))
  expect_equal(res_df$p_value, c(bio1_pval, bio3_pval, bio2_pval))
  expect_equal(
    res_df$p_value_bh,
    p.adjust(c(bio1_pval, bio3_pval, bio2_pval), method = "BH")
  )
  expect_equal(
    res_df$p_value_holm,
    p.adjust(c(bio1_pval, bio3_pval, bio2_pval), method = "holm")
  )

})
