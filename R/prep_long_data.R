#' Prepare Longitudinal Time-to-Event Data
#'
#' \code{prep_long_data()} assesses whether the data passed to \code{sunicate()}
#' is correctly formatted. This function also performs modifications to the
#' data, (e.g. transforming the treatment variable into a factor if not
#' already, and transforming it from a wide to a long format).
#'
#' @param data A wide \code{data.frame} or \code{tibble} object containing the
#'   status (event variable), relative time of the event, treatment indicator,
#'   and covariates. Note that the biomarkers must be a subset of the
#'   covariates, and that there should only be one row per observation.
#' @param failure A \code{character} defining the name of the binary variable in
#'   the \code{data} argument that indicates a failure event. Observations
#'   can have a failure or a censoring event, but not both.
#' @param censor A \code{character} defining the name of the binary variable in
#'   the \code{data} argument that indicates a right-censoring event.
#'   Observations can have a failure or a censoring event, but not both.
#' @param relative_time A \code{character} providing the name of the time
#'   variable in \code{data}.
#' @param treatment A \code{character} indicating the name of the binary
#'   treatment variable in \code{data}.
#' @param covariates A \code{character} vector listing the covariates in
#'   \code{data}.
#' @param biomarkers A \code{character} vector listing the biomarkers of
#'   interest in \code{data}. \code{biomarkers} must be a subset of
#'   \code{covariates}.
#' @param time_cutoff A \code{numeric} representing the time at which to assess
#'   the biomarkers' importance with respect to the outcome. If not specified,
#'   this value is set to the maximum value in the \code{data} argument's
#'   \code{relative_time} variable.
#'
#' @return A longitudinal \code{tibble} containing only the event indicator
#'   variables, a relative time variable, treatment indicator, covariates, and
#'   observation identifier. The treatment variable is transformed into a factor
#'   if not already, and no rows whose relative time values are larger than
#'   \code{time_cutoff} are retained when \code{time_cutoff} is non-null.
#'
#' @importFrom assertthat assert_that validate_that
#' @importFrom tibble is_tibble
#' @importFrom dplyr setequal mutate filter pull select
#' @importFrom rlang sym
#'
#' @keywords internal
prep_long_data <- function(
  data, failure, censor, relative_time, treatment, covariates, biomarkers,
  time_cutoff = NULL
) {

  # check that data is a data.frame or tibble object
  assertthat::assert_that(
    identical(class(data), "data.frame") |
      tibble::is_tibble(data),
    msg = "data should be a data.frame or tibble object"
  )

  # check that the failure, censor, relative_time, treatment, covariates, and
  # biomarkers are characters
  assertthat::assert_that(
    identical(length(failure), 1L) & identical(class(failure), "character"),
    msg = "failure argument should be a character"
  )
  assertthat::assert_that(
    identical(length(censor), 1L) & identical(class(censor), "character"),
    msg = "censor argument should be a character"
  )
  assertthat::assert_that(
    identical(length(relative_time), 1L) &
      identical(class(relative_time), "character"),
    msg = "relative_time argument should be a character"
  )
  assertthat::assert_that(
    identical(length(treatment), 1L) & identical(class(treatment), "character"),
    msg = "treatment argument should be a character"
  )
  assertthat::assert_that(
    identical(class(covariates), "character"),
    msg = "covariates argument should be a character vector"
  )
  assertthat::assert_that(
    all(biomarkers %in% covariates),
    msg = "biomarkers vector is not a subset of the covariates vector"
  )

  # assert that the failure, censor, relative_time, treatment, covariates, and
  # time_cutoff are contained in the data
  assertthat::assert_that(
    failure %in% colnames(data),
    msg = "failure argument's variable is missing from the data"
  )
  assertthat::assert_that(
    censor %in% colnames(data),
    msg = "censor argument's variable is missing from the data"
  )
  assertthat::assert_that(
    relative_time %in% colnames(data),
    msg = "relative_time argument's variable is missing from the data"
  )
  assertthat::assert_that(
    treatment %in% colnames(data),
    msg = "treatment argument's variable is missing from the data"
  )
  assertthat::assert_that(
    all(covariates %in% colnames(data)),
    msg = "some covariates are missing from the data"
  )
  if (!is.null(time_cutoff)) {
    assertthat::assert_that(
      is.numeric(time_cutoff) & length(time_cutoff),
      msg = "time_cutoff should be a single numeric value when specified"
    )
    if (max(data[relative_time]) < time_cutoff) {
      message(paste0("time_cutoff is larger than or equal to the largest ",
                   "value in relative_time's corresponding variable: ",
                   max(data[relative_time])))
    }
  }

  # transform the data to a tibble if it is a data.frame
  if (identical(class(data), "data.frame"))
    data <- data %>% tibble::as_tibble(.name_repair = "minimal")

  # check that the failur and censor arguments' variables are binary numeric
  # variables
  assertthat::assert_that(
    dplyr::setequal(unique(dplyr::pull(data, failure)), c(0, 1)) &
      is.numeric(dplyr::pull(data, failure)),
    msg = paste0("failure argument should correspond to a numeric, binary ",
                 "variable in the data")
  )
  assertthat::assert_that(
    dplyr::setequal(unique(dplyr::pull(data, censor)), c(0, 1)) &
      is.numeric(dplyr::pull(data, censor)),
    msg = paste0("censor argument should correspond to a numeric, binary ",
                 "variable in the data")
  )

  # ensure that no observation has both a censoring and a failure reported
  assertthat::assert_that(
    !(data[[failure]] == 1 && data[[censor]] == 1),
    msg = "observations may not have a both failure and a censoring event"
  )

  # check that the treatment variable is binary
  assertthat::assert_that(
    identical(nrow(unique(data[treatment])), 2L),
    msg = paste0("treatment argument should correspond to a binary variable in",
                 " the data")
  )

  # if the treatment variable isn't already a factor, then transform it
  if (!is.factor(data[treatment])) {
    data <- data %>%
      dplyr::mutate(
        !!rlang::sym(treatment) := factor(!!rlang::sym(treatment))
      )
  }

  # check that the relative time variable is a positive, continuous variable
  assertthat::assert_that(
    all(dplyr::pull(data, relative_time) > 0) &
      is.numeric(dplyr::pull(data, relative_time)),
    msg = paste0("relative_time argument's corresponding variable should be a ",
                 "positive numeric variable")
  )

  # remove unnecessary variables from data
  data <- data %>%
    dplyr::select(
      dplyr::all_of(c(failure, censor, relative_time, treatment, covariates))
    )

  # transform the data from a wide format to a longitudinal format
  times <- data %>% dplyr::pull(relative_time)
  if (!is.null(time_cutoff)) {
    if (max(times) < time_cutoff) {
      times <- c(times, time_cutoff)
    }
  }
  times <- times %>% unique() %>% sort()
  long_data <- lapply(
    seq_len(nrow(data)),
    function(idx) {

      # lengthen data
      n_times <- which(data[[idx, relative_time]] == times)
      if (n_times > 1) {
        long_obs <- replicate(n_times - 1, data[idx, ], simplify = FALSE) %>%
          dplyr::bind_rows() %>%
          dplyr::mutate(
            !!rlang::sym(failure) := 0,
            !!rlang::sym(censor) := 0
          ) %>%
          dplyr::bind_rows(data[idx, ])
      } else {
        long_obs <- data[idx, ]
      }

      # add the relative time back, along with an observation id
      long_obs <- long_obs %>%
        dplyr::mutate(
          time = times[1:n_times],
          observation_id = idx
        )
      # if the relative time variable happens to be called "time", don't filter
      if (relative_time != "time")
        long_obs <- long_obs %>% dplyr::select(-dplyr::all_of(relative_time))

      return(long_obs)

    }
  ) %>%
    dplyr::bind_rows()

  # apply time cutoff if necessary
  if (!is.null(time_cutoff))
    long_data <- long_data %>% dplyr::filter(.data$time <= time_cutoff)

  return(long_data)

}
