prep_long_data <- function(
  data, status, relative_time, treatment, covariates, biomarkers
) {

  # transform the data from a wide format to a longitudinal format
  times <- data %>% dplyr::pull(relative_time) %>% unique() %>% sort()
  num_times <- length(times)
  long_data <- lapply(
    seq_len(nrow(data)),
    function(idx) {
      expanded_obs <- replicate(num_times, data[idx, ], simplify = FALSE) %>%
        dplyr::bind_rows()
      expanded_obs %>%
        dplyr::mutate(
          all_times = times,
          observation_id = idx
        ) %>%
        dplyr::filter(.data$time >= .data$all_times) %>%
        dplyr::mutate(time = .data$all_times) %>%
        dplyr::select(-.data$all_times)
    }
  ) %>%
    dplyr::bind_rows()

  return(long_data)

}
