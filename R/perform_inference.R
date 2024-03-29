#' Perform Inference on Estimated Biomarker Coefficients
#'
#' \code{perform_inference()} tests whether each biomarker's variable importance
#'   parameter is significantly different from zero. The computed test
#'   statistics are asymptotically Normal with mean zero under the null
#'   hypothesis. The resulting p-values for the two-sided tests are subsequently
#'   corrected for multiple testing using the Benjamini-Hochberg method to
#'   control the False Discovery Rate, and Holm's procedure to compute control
#'   the Family-Wise Error Rate.
#'
#' @param biomarkers_tbl A \code{tibble} of estimated coefficients and their
#'   standard errors.
#'
#' @return A \code{tibble} of estimated coefficients, standard errors, test
#'   statistics, and (adjusted) p-values for each biomarker. The biomarkers
#'   are ordered by significance.
#'
#' @importFrom purrr map_dbl
#' @importFrom dplyr mutate arrange .data
#' @importFrom stats p.adjust pnorm
#'
#' @keywords internal
perform_inference <- function(biomarkers_tbl) {

  # compute the test statistic, and then compute the p-values
  biomarkers_tbl <- biomarkers_tbl %>%
    dplyr::mutate(
      z = .data$coef / .data$se,
      p_value = purrr::map_dbl(
        .data$z, function(z_stat) {
          2 * min(
            stats::pnorm(z_stat),
            stats::pnorm(z_stat, lower.tail = FALSE)
          )
        }
      ),
      p_value_bh = stats::p.adjust(.data$p_value, method = "BH"),
      p_value_holm = stats::p.adjust(.data$p_value, method = "holm")
    ) %>%
    dplyr::arrange(.data$p_value)

  return(biomarkers_tbl)
}
