prep_long_data <- function(
  data, status, relative_time, treatment, covariates, biomarkers,
  time_cutoff = NULL
) {

  # check that data is a data.frame or tibble object
  assertthat::assert_that(
    identical(class(data), "data.frame") |
      tibble::is_tibble(data),
    msg = "data should be a data.frame or tibble object"
  )

  # check that the status, relative_time, treatment, covariates, and biomarkers
  # are characters
  assertthat::assert_that(
    identical(length(status), 1L) & identical(class(status), "character"),
    msg = "status argument should be a character"
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

  # assert that the status, relative_time, treatment, covariates, and
  # time_cutoff are contained in the data
  assertthat::assert_that(
    status %in% colnames(data),
    msg = "status variable is missing from the data"
  )
  assertthat::assert_that(
    relative_time %in% colnames(data),
    msg = "relative_time variable is missing from the data"
  )
  assertthat::assert_that(
    treatment %in% colnames(data),
    msg = "treatment variable is missing from the data"
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
    assertthat::assert_that(
      max(data[relative_time]) >= time_cutoff,
      msg = paste0("time_cutoff should be smaller than or equal to the largest ",
                   "value in relative_time's corresponding variable: ",
                   max(data[relative_time]))
    )
  }

  # transform the data to a tibble if it is a data.frame
  if (identical(class(data), "data.frame"))
    data <- data %>% tibble::as_tibble(.name_repair = "minimal")

  # check that the status argument's variable is a binary numeric variable
  assertthat::assert_that(
    dplyr::setequal(unique(pull(data, status)), c(0, 1)),
    msg = paste0("status argument should correspond to a numeric, binary ",
                 "variable in the data")
  )

  # check that the treatment variable is binary
  assertthat::assert_that(
    identical(nrow(unique(data[treatment])), 2L),
    msg = "treatment argument should correspond to a binary variable in the data"
  )

  # if the treatment variable isn't already a factor, then transform it
  if (!is.factor(data[treatment]))
    data <- data %>%
    dplyr::mutate(
      !!rlang::sym(treatment) := factor(!!rlang::sym(treatment))
    )

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
      dplyr::all_of(c(status, relative_time, treatment, covariates))
    )

  # transform the data from a wide format to a longitudinal format
  times <- data %>% dplyr::pull(relative_time) %>% unique() %>% sort()
  num_times <- length(times)
  long_data <- lapply(
    seq_len(nrow(data)),
    function(idx) {

      # lengthen data
      expanded_obs <- replicate(num_times, data[idx, ], simplify = FALSE) %>%
        dplyr::bind_rows()
      expanded_obs <- expanded_obs %>%
        dplyr::mutate(
          all_times = times,
          observation_id = idx
        ) %>%
        dplyr::filter(.data$time >= .data$all_times) %>%
        dplyr::mutate(time = .data$all_times) %>%
        dplyr::select(-.data$all_times)

      # ensure that failures are properly recorded
      final_status <- data[[idx, status]]
      final_status_time <- data[[idx, relative_time]]
      expanded_obs <- expanded_obs %>%
        mutate(
          status = if_else(.data$time == final_status_time, final_status, 0)
        )

      return(expanded_obs)

    }
  ) %>%
    dplyr::bind_rows()

  if (!is.null(time_cutoff))
    long_data <- long_data %>% dplyr::filter(.data$time <= time_cutoff)

  return(long_data)

}
