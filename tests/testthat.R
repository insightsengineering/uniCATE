Sys.setenv("R_TESTS" = "")

library(testthat)
library(uniCATE)

test_check("uniCATE")
