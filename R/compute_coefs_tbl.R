#' Compute Cross-Validated Estimate and its Standard Errors
#'
#' \code{compute_coefs_tbl()} estimates the linear model coefficients of each
#'   biomarker's variable importance parameter by aggregating the estimates of
#'   the validation folds. The cross-validated standard errors are also
#'   computed.
#'
#' @param cv_ls A \code{list} containing two \code{tibbles}. The first,
#'   \code{betas_df}, contains the estimated beta coefficients for each
#'   biomarker across all validation sets. The second, \code{ic_df}, is made up
#'   of the cross-validated emipirical efficient influence functions of every
#'   observation in \code{data}.
#'
#' @return A \code{tibble} of estimated coefficients and accompanying standard
#'   errors.
#'
#' @importFrom tibble tibble
#'
#' @keywords internal
compute_coefs_tbl <- function(cv_ls) {

  # compute the beta coefficients
  betas <- colMeans(cv_ls$betas_df)

  # compute the standard errors
  ses <- sqrt(colMeans(cv_ls$ic_df^2) / nrow(cv_ls$ic_df))

  # create the tibble
  betas_tbl <- tibble::tibble(
    "biomarker" = colnames(cv_ls$betas_df),
    "coef" = betas,
    "se" = ses
  )

  return(betas_tbl)
}
