#' Estimate the Univariate Survival CATEs
#'
#' \code{estimate_univariate_survival_cates()} estimates the slope of the
#'   univariate conditional average treatment effect's linear approximation for
#'   right-censored time-to-event outcomes. First, the data is pre-processed
#'   and split into \code{v_folds} folds. Then, the conditional survival and
#'   conditional censoring functions are estimated using their respective
#'   hazard functions' SuperLearners. These estimates are then used to estimate
#'   the expected difference in potential outcomes across treatment conditions
#'   for each individual in their respective validation set. The variable
#'   importance parameters are then estimated in the validation sets, along
#'   with their efficient influence functions.
#'
#' @param long_data A \code{tibble} object containing the longitudinal data
#'   output by \code{prep_long_data()}.
#' @param event A \code{character} defining the name of the binary variable in
#'   the \code{data} argument that indicates whether an event occurred.
#'   Observations can have an event or be censored, but not both.
#' @param censor A \code{character} defining the name of the binary variable in
#'   the \code{data} argument that indicates a right-censoring event.
#'   Observations can have an event or be censored, but not both.
#' @param treatment A \code{character} indicating the name of the binary
#'   treatment variable in \code{data}.
#' @param biomarkers A \code{character} vector listing the biomarkers of
#'   interest in \code{data}.
#' @param cond_surv_haz_super_learner A \code{\link[sl3:Lrnr_sl]{Lrnr_sl}}
#'   object used to estimate the conditional event hazard model. If set to
#'   \code{NULL}, the default, an elastic net regression is used instead. It is
#'   best to use this default behaviour when analyzing small datasets.
#' @param cond_censor_haz_super_learner A \code{\link[sl3:Lrnr_sl]{Lrnr_sl}}
#'   object used to estimate the conditional censoring hazard model. If set to
#'   \code{NULL}, the default, an elastic net regression is used instead. It is
#'   best to use this default behaviour when analyzing small datasets.
#' @param propensity_score_ls A named \code{numeric} \code{list} providing the
#'   propensity scores for the treatment conditions. The first element of the
#'   list should correspond to the "treatment" condition, and the second to the
#'   "control" condition, whatever their names may be.
#' @param v_folds A \code{numeric} indicating the number of folds used for
#'   V-fold cross-validation. Defaults to \code{5}.
#' @param parallel A \code{logical} determining whether to use
#'   \code{\link[origami:cross_validate]{origami}}'s built-in parallelized
#'   cross-validation routines. This parallelization framework is built upon
#'   the \href{https://cran.r-project.org/package=future}{\code{future}} suite.
#'   Defaults to \code{FALSE}.
#'
#' @return A \code{list} containing two \code{tibble} objects. The first,
#'   \code{betas_df}, contains the estimated beta coefficients for each
#'   biomarker across all validation sets. The second, \code{ic_df}, is made up
#'   of the cross-validated empirical efficient influence functions of every
#'   observation in \code{data}.
#'
#' @importFrom origami make_folds folds_vfold cross_validate
#' @importFrom assertthat assert_that
#' @importFrom dplyr bind_rows as_tibble
#'
#' @keywords internal
estimate_univariate_s_cates <- function(long_data,
                                        event,
                                        censor,
                                        treatment,
                                        biomarkers,
                                        cond_surv_haz_super_learner,
                                        cond_censor_haz_super_learner,
                                        propensity_score_ls,
                                        v_folds,
                                        parallel) {

  # split the data into folds
  folds <- origami::make_folds(
    long_data,
    fold_fun = origami::folds_vfold,
    V = v_folds,
    cluster_ids = long_data$observation_id
  )

  # assert that the propensity score list's elements are named after the
  # treatment assignment variable, and that they sum to 1
  assertthat::assert_that(
    all(names(propensity_score_ls) %in%
      levels(dplyr::pull(long_data, treatment))),
    msg = "treatment group names and propensity_score_ls names do not match"
  )
  assertthat::assert_that(
    sum(unlist(propensity_score_ls)) == 1,
    msg = "propensity_score_ls values do not sum to one"
  )

  # assert that the super learners are  NULL or SL3 SuperLearner objects
  assertthat::assert_that(
    is.null(cond_surv_haz_super_learner) |
      identical(
        class(cond_surv_haz_super_learner),
        c("Lrnr_sl", "Lrnr_base", "R6")
      ),
    msg = "cond_surv_haz_super_learner must be NULL or an sl3::Lrnr_sl object"
  )
  assertthat::assert_that(
    is.null(cond_censor_haz_super_learner) |
      identical(
        class(cond_censor_haz_super_learner),
        c("Lrnr_sl", "Lrnr_base", "R6")
      ),
    msg = "cond_censor_haz_super_learner must be NULL or an sl3::Lrnr_sl object"
  )

  # compute the holdout estimated potential outcome differences
  hold_out_calculations <- origami::cross_validate(
    cv_fun = hold_out_calculation_survival,
    folds = folds,
    long_data = long_data,
    event = event,
    censor = censor,
    treatment = treatment,
    biomarkers = biomarkers,
    cond_surv_haz_super_learner = cond_surv_haz_super_learner,
    cond_censor_haz_super_learner = cond_censor_haz_super_learner,
    propensity_score_ls = propensity_score_ls,
    use_future = parallel,
    .combine = FALSE
  )

  # aggregate the beta coefficient vectors into a tibble
  betas_df <- dplyr::bind_rows(hold_out_calculations$beta_coefs)
  colnames(betas_df) <- biomarkers
  betas_df <- betas_df %>% dplyr::as_tibble()

  # aggregate the IC tibbles into a tibble
  ic_df <- hold_out_calculations$ic_df %>% dplyr::bind_rows()

  return(
    list(
      "betas_df" = betas_df,
      "ic_df" = ic_df
    )
  )
}


################################################################################

#' Compute Validation Set Objects
#'
#' \code{hold_out_calculation_survival()} computes all the objects required to
#'   estimate the univariate CATEs for survival outcomes over the validation
#'   set. It begins by estimating the conditional survival and censoring
#'   hazard functions on the training set using the
#'   \code{cond_surv_haz_super_learner} and
#'   \code{cond_censor_haz_super_learner}, respectively. Next, the difference in
#'   survival probabilities is predicted for the observations in the validation
#'   set. Finally, these predicted outcomes are used to estimate the variable
#'   impotance parameters of each biomarker. Their respective empirical
#'   efficient influence functions are computed as well.
#'
#' @param fold A \code{fold} object (from \code{\link[origami]{make_folds}()}).
#' @param long_data A \code{tibble} object containing the longitudinal data
#'   output by \code{prep_long_data()}.
#' @param event A \code{character} defining the name of the binary variable in
#'   the \code{data} argument that indicates whether an event occurred.
#'   Observations can have an event or be censored, but not both.
#' @param censor A \code{character} defining the name of the binary variable in
#'   the \code{data} argument that indicates a right-censoring event.
#'   Observations can have an event or be censored, but not both.
#' @param treatment A \code{character} indicating the name of the binary
#'   treatment variable in \code{data}.
#' @param biomarkers A \code{character} vector listing the biomarkers of
#'   interest in \code{data}.
#' @param cond_surv_haz_super_learner A \code{\link[sl3:Lrnr_sl]{Lrnr_sl}}
#'   object used to estimate the conditional event hazard model. If set to
#'   \code{NULL}, the default, an elastic net regression is used instead. It is
#'   best to use this default behaviour when analyzing small datasets.
#' @param cond_censor_haz_super_learner A \code{\link[sl3:Lrnr_sl]{Lrnr_sl}}
#'   object used to estimate the conditional censoring hazard model. If set to
#'   \code{NULL}, the default, an elastic net regression is used instead. It is
#'   best to use this default behaviour when analyzing small datasets.
#' @param propensity_score_ls A named \code{numeric} \code{list} providing the
#'   propensity scores for the treatment conditions. The first element of the
#'   list should correspond to the "treatment" condition, and the second to the
#'   "control" condition, whatever their names may be.
#'
#' @return A \code{list} made up of two objects. The first is the \code{numeric}
#'   vector of biomarker variable importance estimates. The second is the
#'   \code{tibble} of the empirical efficient influence functions for each
#'   biomarker.
#'
#' @importFrom dplyr mutate select pull .data all_of bind_cols left_join
#' @importFrom rlang !! enquo sym
#' @importFrom magrittr %>%
#' @importFrom origami training validation
#' @importFrom purrr map
#' @importFrom stats predict
#' @importFrom glmnet cv.glmnet
#' @import sl3
#'
#' @keywords internal
hold_out_calculation_survival <- function(fold,
                                          long_data,
                                          event,
                                          censor,
                                          treatment,
                                          biomarkers,
                                          cond_surv_haz_super_learner,
                                          cond_censor_haz_super_learner,
                                          propensity_score_ls) {

  # define the training and testing set
  train_data <- origami::training(long_data)
  valid_data <- origami::validation(long_data)

  # grab the covariates' column names
  covar_names <- colnames(train_data)
  rm_noncovar_regex <- paste0("(", event, "|", censor, "|observation_id)")
  covar_names <- covar_names[which(!grepl(rm_noncovar_regex, covar_names))]

  # remove event and censor indicators, and add times up to t_0 for each
  # observation in a new version of the validation data
  times <- c(valid_data$time, train_data$time) %>%
    unique() %>%
    sort()
  n_times <- length(times)
  surv_valid_data <- valid_data %>%
    dplyr::select(-dplyr::all_of(c(event, censor, "time"))) %>%
    dplyr::distinct()
  surv_valid_data <- lapply(
    seq_len(nrow(surv_valid_data)),
    function(idx) {
      replicate(n_times, surv_valid_data[idx, ], simplify = FALSE) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(time = times)
    }
  ) %>%
    dplyr::bind_rows()

  # prepare the validation datasets for fixed treatment assignments
  valid_data_treat <- surv_valid_data
  valid_data_treat[treatment] <- names(propensity_score_ls)[1]
  valid_data_treat[treatment] <- factor(
    dplyr::pull(valid_data_treat, treatment),
    levels = names(propensity_score_ls)
  )
  valid_data_cont <- surv_valid_data
  valid_data_cont[treatment] <- names(propensity_score_ls)[2]
  valid_data_cont[treatment] <- factor(
    dplyr::pull(valid_data_cont, treatment),
    levels = names(propensity_score_ls)
  )

  # estimate conditional survival function #####################################

  # create the super learner if not already defined
  if (is.null(cond_surv_haz_super_learner)) {

    # prepare the data for glmnet
    y_train <- train_data[[event]]
    x_train <- train_data %>% dplyr::select(dplyr::all_of(covar_names))
    x_train <- matrix(as.numeric(unlist(x_train)), nrow = nrow(x_train))

    # train the elastic net
    glmnet_fit <- cv.glmnet(
      x = x_train,
      y = y_train,
      nfolds = 5,
      alpha = 0.5,
      family = "binomial"
    )

    # estimate the conditional hazard under treatment and control in
    # surv_valid_data to eventually compute the survival prob
    x_data_treat <- valid_data_treat %>%
      dplyr::select(dplyr::all_of(covar_names))
    x_data_treat <- matrix(
      as.numeric(unlist(x_data_treat)),
      nrow = nrow(x_data_treat)
    )
    surv_valid_data$cond_surv_haz_treat <- stats::predict(
      glmnet_fit,
      newx = x_data_treat,
      s = "lambda.min",
      type = "response"
    )
    x_data_cont <- valid_data_cont %>% dplyr::select(dplyr::all_of(covar_names))
    x_data_cont <- matrix(
      as.numeric(unlist(x_data_cont)),
      nrow = nrow(x_data_cont)
    )
    surv_valid_data$cond_surv_haz_control <- stats::predict(
      glmnet_fit,
      newx = x_data_cont,
      s = "lambda.min",
      type = "response"
    )

    # estimate the conditional survival hazard in the validation data
    x_valid <- valid_data %>% dplyr::select(dplyr::all_of(covar_names))
    x_valid <- matrix(as.numeric(unlist(x_valid)), nrow = nrow(x_valid))
    valid_data$cond_surv_haz <- stats::predict(
      glmnet_fit,
      newx = x_valid,
      s = "lambda.min",
      type = "response"
    )
  } else {

    # define an sl3 task
    surv_haz_train_task <- sl3::make_sl3_Task(
      data = train_data,
      covariates = covar_names,
      outcome = event,
      outcome_type = "binomial"
    )

    # estimate the conditional hazard function
    cond_surv_haz_fit <- cond_surv_haz_super_learner$train(surv_haz_train_task)

    # create the treatment and control prediction tasks
    pred_task_treat <- sl3::make_sl3_Task(
      data = valid_data_treat,
      covariates = covar_names
    )
    pred_task_cont <- sl3::make_sl3_Task(
      data = valid_data_cont,
      covariates = covar_names
    )

    # estimate the conditional hazard under treatment and control in
    # surv_valid_data to eventually compute the survival prob
    surv_valid_data$cond_surv_haz_treat <- cond_surv_haz_fit$predict(
      pred_task_treat
    )
    surv_valid_data$cond_surv_haz_control <- cond_surv_haz_fit$predict(
      pred_task_cont
    )

    # estimate the conditional hazard in the validation data
    pred_task_valid <- sl3::make_sl3_Task(
      data = valid_data,
      covariates = covar_names
    )
    valid_data$cond_surv_haz <- cond_surv_haz_fit$predict(pred_task_valid)
  }

  # estimate the survival probs under treatment and control at t0
  surv_valid_data <- surv_valid_data %>%
    dplyr::group_by(.data$observation_id) %>%
    dplyr::summarise(
      surv_t0_treat = prod(1 - .data$cond_surv_haz_treat),
      surv_t0_cont = prod(1 - .data$cond_surv_haz_control)
    ) %>%
    dplyr::ungroup()


  # estimate censoring function ################################################

  # create the super learner if not already defined
  if (is.null(cond_censor_haz_super_learner)) {

    # prepare the data for glmnet
    y_train <- train_data[[censor]]
    x_train <- train_data %>% dplyr::select(dplyr::all_of(covar_names))
    x_train <- matrix(as.numeric(unlist(x_train)), nrow = nrow(x_train))

    # train the elastic net
    glmnet_fit <- cv.glmnet(
      x = x_train,
      y = y_train,
      nfolds = 5,
      alpha = 0.5,
      family = "binomial"
    )

    # predict the conditional censoring hazard at each time for each treatment
    x_valid <- valid_data %>% dplyr::select(dplyr::all_of(covar_names))
    x_valid <- matrix(as.numeric(unlist(x_valid)), nrow = nrow(x_valid))
    valid_data$cond_cens_haz <- stats::predict(
      glmnet_fit,
      newx = x_valid,
      s = "lambda.min",
      type = "response"
    )
  } else {

    # define an sl3 task
    cens_haz_train_task <- sl3::make_sl3_Task(
      data = train_data,
      covariates = covar_names,
      outcome = censor,
      outcome_type = "binomial"
    )

    # estimate the conditional hazard function
    cond_cens_haz_fit <- cond_censor_haz_super_learner$train(
      cens_haz_train_task
    )

    # predict cond hazards at each time
    pred_task_valid <- sl3::make_sl3_Task(
      data = valid_data,
      covariates = covar_names
    )
    valid_data$cond_cens_haz <- cond_cens_haz_fit$predict(pred_task_valid)
  }

  # construct cumulative portion of empirical efficient influence curve ########

  # add the survival probs at t0
  valid_data <- valid_data %>%
    dplyr::left_join(surv_valid_data, by = "observation_id")

  # compute the ajdusted differences in survival probabilities
  surv_diff <- valid_data %>%
    dplyr::group_by(.data$observation_id) %>%
    dplyr::mutate(
      surv = cumprod(1 - .data$cond_surv_haz),
      surv_cens = cumprod(1 - .data$cond_cens_haz),
      surv_cens_prev = dplyr::lag(.data$surv_cens),
      surv_cens_prev = dplyr::if_else(is.na(.data$surv_cens_prev), 1,
        .data$surv_cens_prev
      ),
      h1 = dplyr::if_else(
        !!rlang::sym(treatment) != names(propensity_score_ls)[1], 0,
        -.data$surv_t0_treat /
          (propensity_score_ls[[1]] * .data$surv_cens_prev * .data$surv)
      ),
      h0 = dplyr::if_else(
        !!rlang::sym(treatment) != names(propensity_score_ls)[2], 0,
        -.data$surv_t0_cont /
          (propensity_score_ls[[2]] * .data$surv_cens_prev * .data$surv)
      ),
      haz_prod = ((!!rlang::sym(event) == 1) - .data$cond_surv_haz),
      inner_sum_t = (.data$h1 - .data$h0) * .data$haz_prod
    ) %>%
    dplyr::summarise(
      cumult_eif = sum(.data$inner_sum_t),
      adj_surv_diff = .data$cumult_eif + .data$surv_t0_treat -
        .data$surv_t0_cont,
      .groups = "drop"
    ) %>%
    dplyr::distinct() %>%
    dplyr::pull(.data$adj_surv_diff)

  # estimate variable important parameters and emp EIFs ########################
  valid_data <- valid_data %>%
    dplyr::select(dplyr::all_of(biomarkers)) %>%
    dplyr::distinct()
  coefs_and_ic_ls <- valid_data %>%
    purrr::map2(
      colnames(valid_data),
      function(bio, bio_name) {

        # center the biomarker measurements
        bio <- as.vector(base::scale(bio, center = TRUE, scale = FALSE))
        var_bio <- mean(bio^2)

        # if the variance of the biomarker is nonzero, then proceed
        # otherwise, fail gracefully
        if (var_bio > 0) {

          # estimate the best linear approximation using the estimating equation
          # formula
          bio_coef <- mean(surv_diff * bio) / var_bio

          # compute the empirical IC of each observation
          inf_curves <- ((surv_diff - bio_coef * bio) * bio) / var_bio

        } else {
          # set the coefficient to zero, and make the influence curve enormous
          bio_coef <- 0
          inf_curves <- rep(1000000, length(surv_diff))
          message(paste("Biomarker", bio_name, "has low variability. Remove",
                        "it and repeat the analysis."))
        }
        # return the beta coefficients and the influence curves
        return(list(
          "bio_coef" = bio_coef,
          "inf_curves" = inf_curves
        ))
      }
    )

  # extract the vector of coefficient estimates
  beta_coefs <- sapply(
    seq_len(length(coefs_and_ic_ls)),
    function(idx) coefs_and_ic_ls[[idx]]$bio_coef
  )
  names(beta_coefs) <- biomarkers

  # extract the table of un-scaled influence curves
  ic_ls <- lapply(
    seq_len(length(coefs_and_ic_ls)),
    function(idx) coefs_and_ic_ls[[idx]]$inf_curves
  )
  ic_df <- do.call(cbind, ic_ls) %>% dplyr::as_tibble(.name_repair = "minimal")
  colnames(ic_df) <- biomarkers

  # return the vector of coefficients and the table of empirical IC values
  return(
    list(
      "beta_coefs" = beta_coefs,
      "ic_df" = ic_df
    )
  )
}
