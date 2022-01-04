#' Prepare Data
#'
#' \code{prep_data()} assesses whether the data passed to \code{unicate()} is
#' correctly formatted. This function also performs minor modifications to the
#' data, (e.g. transforming the treatment variable into a factor if not
#' already).
#'
#' @param data A \code{data.frame} or \code{tibble} object containing the
#'   outcome variable, treatment indicator, and covariates. Note that the
#'   biomarkers must be a subset of the covariates.
#' @param outcome A \code{character} defining the name of the outcome variable
#'   in \code{data}.
#' @param treatment A \code{character} indicating the name of the binary
#'   treatment variable in \code{data}.
#' @param covariates A \code{character} vector listing the covariates in
#'   \code{data}.
#' @param biomarkers A \code{character} vector listing the biomarkers of
#'   interest in \code{data}. \code{biomarkers} must be a subset of
#'   \code{covariates}.
#'
#' @return A \code{tibble} containing only the outcome variable, treatment
#'   indicator, and covariates. Note that the treatment variable is transformed
#'   into a factor if not already.
#'
#' @importFrom assertthat assert_that
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select all_of
#' @importFrom tibble as_tibble is_tibble
#' @importFrom rlang !! sym :=
#'
#' @keywords internal
prep_data <- function(data, outcome, treatment, covariates, biomarkers) {

  # assert that the arguments are properly defined
  assertthat::assert_that(
    identical(class(data), "data.frame") |
      tibble::is_tibble(data),
    msg = "data should be a data.frame or tibble object"
  )

  # assert that the arguments are of the appropriate length
  assertthat::assert_that(
    identical(length(outcome), 1L) & identical(class(outcome), "character"),
    msg = "outcome should be a length-one character"
  )
  assertthat::assert_that(
    identical(length(treatment), 1L) & identical(class(treatment), "character"),
    msg = "treatment should be a length-one character"
  )

  # assert that the outcome, treatment, covariates, and biomarkers are contained
  # in the data
  assertthat::assert_that(
    outcome %in% colnames(data),
    msg = "outcome is not contained in data"
  )
  assertthat::assert_that(
    treatment %in% colnames(data),
    msg = "treatment is not contained in data"
  )
  assertthat::assert_that(
    all(covariates %in% colnames(data)),
    msg = "not all covariates are contained in data"
  )
  assertthat::assert_that(
    all(biomarkers %in% covariates),
    msg = "biomarkers is not a subset of covariates"
  )

  # assert that the outcome is a numeric variable (continuous or binary)
  assertthat::assert_that(
    (class(data[[outcome]]) == "numeric" &
       length(unique(data[[outcome]])) > 2) |
      (class(data[[outcome]]) == "numeric" &
         length(unique(data[[outcome]])) == 2 &
           sum(is.element(unique(data[[outcome]]), c(0, 1))) == 2),
    msg = "outcome should be a continuous or binary numeric variable"
  )

  # assert that the treatment is binary
  assertthat::assert_that(
    identical(nrow(unique(data[treatment])), 2L),
    msg = "treatment should be a binary variable"
  )

  # if the treatment variable isn't already a factor, then transform it
  if (!is.factor(data[treatment]))
    data <- data %>%
    dplyr::mutate(
      !!rlang::sym(treatment) := factor(!!rlang::sym(treatment))
    )

  # transform the data to a tibble if it is a data.frame
  if (identical(class(data), "data.frame"))
    data <- data %>% as_tibble(.name_repair = "minimal")

  # retain only the treatment, outcome, and covariates
  data <- data %>%
    dplyr::select(
      dplyr::all_of(outcome), dplyr::all_of(treatment),
      dplyr::all_of(covariates)
    )

  # return the vetted tibble
  return(data)

}