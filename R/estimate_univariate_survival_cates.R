#' Estimate the Univariate Survival CATEs
#'
#' \code{estimate_univariate_survival cates()} estimates the univariate
#' Conditional Average Treatment Effect of each biomarker on the predicted
#' potential time-to-event outcomes using cross-validation. First, the data is
#' pre-processed and split into \code{v_folds}. Then, the conditional survival
#' and censoring regressions are estimated over the training sets using
#' \code{failure_super_learner} and \code{censoring_super_learner},
#' respectively. The estimated outcome regressions are then used to estimate the
#' expected difference in potential outcomes across treatment groups for each
#' individual in their respective validation sets. Univariate linear regressions
#' are then regressed over the predicted potential outcomes in the validation
#' set, and the influence curves of all observations is computed as well.
#' \code{tibbles} of biomarker coefficients and influence curves are returned
#' for each validation set.
#'
#' @param long_data A \code{tibble} object containing the longitudinal data
#'   output by \code{prep_long_data()}.
#' @param failure A \code{character} defining the name of the binary variable in
#'   the \code{data} argument that indicates a failure event. Observations
#'   can have a failure or a censoring event, but not both.
#' @param censor A \code{character} defining the name of the binary variable in
#'   the \code{data} argument that indicates a right-censoring event.
#'   Observations can have a failure or a censoring event, but not both.
#' @param treatment A \code{character} indicating the name of the binary
#'   treatment variable in \code{data}.
#' @param biomarkers A \code{character} vector listing the biomarkers of
#'   interest in \code{data}.
#' @param failure_super_learner A \code{\link[sl3:Lrnr_sl]{SuperLearner}} object
#'   used to estimate the conditional failure model. If set to \code{NULL}, the
#'   default SuperLearner is used. The default's library consists of a linear
#'   model, penalized linear models (LASSO and elasticnet), XGBoost, a
#'   Random Forest, and the mean model.
#' @param censoring_super_learner A \code{\link[sl3:Lrnr_sl]{SuperLearner}}
#'   object used to estimate the conditional censoring model. If set to
#'   \code{NULL}, the default SuperLearner is used. The default's library
#'   consists of a linear model, penalized linear models (LASSO and elasticnet),
#'   XGBoost, a Random Forest, and the mean model.
#' @param propensity_score_ls A named \code{numeric} \code{list} providing the
#'   propensity scores for the treatment levels. The first element of the list
#'   should correspond to the "treatment" group, and the second to the "control"
#'   group, whatever their names may be.
#' @param v_folds A \code{numeric} indicating the number of folds to use for
#'   V-fold cross-validation. Defaults to \code{5}.
#' @param parallel A \code{logical} determining whether to use
#'   \code{\link[origami:cross_validate]{origami}}'s built-in parallelization
#'   capabilities. Defaults to \code{FALSE}.
#'
#' @return A \code{list} containing two \code{tibbles}. The first,
#'   \code{betas_df}, contains the estimated beta coefficients for each
#'   biomarker across all validation sets. The second, \code{ic_df}, is made up
#'   of the cross-validated influence curves of every observation in
#'   \code{data}.
#'
#' @importFrom origami make_folds folds_vfold cross_validate
#' @importFrom assertthat assert_that
#' @importFrom dplyr bind_rows as_tibble
#'
#' @keywords internal
estimate_univariate_survival_cates <- function(
  long_data,
  failure,
  censor,
  treatment,
  biomarkers,
  failure_super_learner,
  censoring_super_learner,
  propensity_score_ls,
  v_folds,
  parallel
) {

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
    is.null(failure_super_learner) |
      identical(class(failure_super_learner), c("Lrnr_sl", "Lrnr_base", "R6")),
    msg = "failure_super_learner must be NULL or an sl3::Lrnr_sl object"
  )
  assertthat::assert_that(
    is.null(censoring_super_learner) |
      identical(class(censoring_super_learner), c("Lrnr_sl", "Lrnr_base", "R6")),
    msg = "censoring_super_learner must be NULL or an sl3::Lrnr_sl object"
  )

  # compute the holdout estimated potential outcome differences
  hold_out_calculations <- origami::cross_validate(
    cv_fun = hold_out_calculation_survival,
    folds = folds,
    long_data = long_data,
    treatment = treatment,
    biomarkers = biomarkers,
    failure_super_learner = failure_super_learner,
    censoring_super_learner = censoring_super_learner,
    propensity_score_ls = propensity_score_ls,
    use_future = parallel,
    .combine = FALSE
  )

  # aggregate the beta coefficient vectors into a tibble
  betas_df <- hold_out_calculations$beta_coefs %>%
    unlist() %>%
    matrix(nrow = v_folds, byrow = TRUE)
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
#' \code{hold_out_calculation_survival} computes all the objects required to
#' estimate the univariate CATEs for survival outcomes over the validation set.
#' It begins by estimating the conditional survival and censoring survival
#' functions on the training set using the \code{failure_super_learner} and
#' \code{censoring_super_learner}, respectively. Next, the difference in
#' survival probabilities is predicted for the observations in the validation
#' set. Finally, these predicted outcomes are used to estimate the unvariate
#' regression coefficients of each biomarker. Their respective empirical
#' efficient influence curves are computed as well.
#'
#' @param fold A \code{fold} object (from \code{\link[origami]{make_folds}()}).
#' @param long_data A \code{tibble} object containing the longitudinal data
#'   output by \code{prep_long_data()}.
#' @param failure A \code{character} defining the name of the binary variable in
#'   the \code{data} argument that indicates a failure event. Observations
#'   can have a failure or a censoring event, but not both.
#' @param censor A \code{character} defining the name of the binary variable in
#'   the \code{data} argument that indicates a right-censoring event.
#'   Observations can have a failure or a censoring event, but not both.
#' @param treatment A \code{character} indicating the name of the binary
#'   treatment variable in \code{data}.
#' @param biomarkers A \code{character} vector listing the biomarkers of
#'   interest in \code{data}.
#' @param cond_surv_haz_super_learner A \code{\link[sl3:Lrnr_sl]{SuperLearner}}
#'   object used to estimate the conditional hazard model. If set to
#'   \code{NULL}, the default SuperLearner is used. The default's library
#'   consists of a linear model, penalized linear models (LASSO and elasticnet),
#'   a Random Forest, and the mean model.
#' @param cond_censor_haz_super_learner A
#'   \code{\link[sl3:Lrnr_sl]{SuperLearner}} object used to estimate the
#'   conditional hazard model. If set to \code{NULL}, the default SuperLearner
#'   is used. The default's library consists of a linear model, penalized linear
#'   models (LASSO and elasticnet), a Random Forest, and the mean model.
#' @param propensity_score_ls A named \code{numeric} \code{list} providing the
#'   propensity scores for the treatment levels. The first element of the list
#'   should correspond to the "treatment" group, and the second to the "control"
#'   group, whatever their names may be.
#'
#' @return A \code{list} made up of two objects. The first is the \code{numeric}
#'   vector of biomarker linear model coefficient estimates. The second is the
#'   \code{tibble} of the empirical influence curves for each biomarker.
#'
#' @importFrom dplyr mutate select pull .data all_of bind_cols
#' @importFrom rlang !! enquo
#' @importFrom magrittr %>%
#' @importFrom origami training validation
#' @importFrom stats var coef lm
#' @importFrom purrr map
#' @importFrom matrixStats colVars
#' @import sl3
#'
#' @keywords internal
hold_out_calculation_survival <- function(
  fold, long_data, failure, censor, treatment, biomarkers,
  cond_surv_haz_super_learner, cond_censor_haz_super_learner,
  propensity_score_ls, outcome_type
) {

  # define the training and testing set
  train_data <- origami::training(long_data)
  valid_data <- origami::validation(long_data)

  # grab the covariates' column names
  covar_names <- colnames(train_data)
  rm_noncovar_regex <- "(failure|censor|observation_id)"
  covar_names <- covar_names[which(!grepl(rm_noncovar_regex, covar_names))]


  # estimate conditional survival function #####################################

  # define an sl3 task
  surv_haz_train_task <- sl3::make_sl3_Task(
    data = train_data,
    covariates = covar_names,
    outcome = failure,
    outcome_type = "binomial"
  )

  # create the super learner if not already defined
  if (is.null(cond_surv_haz_super_learner)) {

    # define the treatment-covariate interaction for certain learners
    covars <- covar_names[which(!grepl(treatment, covar_names))]
    interactions <- lapply(covars, function(w) c(w, treatment))
    lrnr_interactions <- sl3::Lrnr_define_interactions$new(interactions)

    # define the base learners for binary outcomes
    lrnr_glm <- sl3::make_learner(
      sl3::Pipeline, lrnr_interactions, sl3::Lrnr_glm_fast$new()
    )
    lrnr_lasso <- sl3::make_learner(
      sl3::Pipeline, lrnr_interactions, sl3::Lrnr_glmnet$new()
    )
    lrnr_enet <- sl3::make_learner(
      sl3::Pipeline, lrnr_interactions, sl3::Lrnr_glmnet$new(alpha = 0.5)
    )
    lrnr_rf <- sl3::Lrnr_ranger$new()
    lrnr_mean <- sl3::Lrnr_mean$new()

    # assemble learners
    learner_library <- sl3::make_learner(
      sl3::Stack, lrnr_glm, lrnr_lasso, lrnr_enet, lrnr_xgboost, lrnr_rf,
      lrnr_mean
    )

    # define the metalearner
    meta_learner <- sl3::make_learner(
      sl3::Lrnr_solnp,
      loss_function = sl3::loss_loglik_binomial,
      learner_function = sl3::metalearner_logistic_binomial
    )

    # initialize the SuperLearner
    cond_surv_haz_super_learner <- sl3::Lrnr_sl$new(
      learners = learner_library,
      metalearner = meta_learner
    )

  }

  # estimate the conditional hazard function
  cond_surv_haz_fit <- cond_surv_haz_super_learner$train(surv_haz_train_task)



  # estimate censoring function ################################################
  # define an sl3 task
  cens_haz_train_task <- sl3::make_sl3_Task(
    data = train_data,
    covariates = covar_names,
    outcome = censor,
    outcome_type = "binomial"
  )

  # create the super learner if not already defined
  if (is.null(cond_censor_haz_super_learner)) {

    # define the treatment-covariate interaction for certain learners
    covars <- covar_names[which(!grepl(treatment, covar_names))]
    interactions <- lapply(covars, function(w) c(w, treatment))
    lrnr_interactions <- sl3::Lrnr_define_interactions$new(interactions)

    # define the base learners for binary outcomes
    lrnr_glm <- sl3::make_learner(
      sl3::Pipeline, lrnr_interactions, sl3::Lrnr_glm_fast$new()
    )
    lrnr_lasso <- sl3::make_learner(
      sl3::Pipeline, lrnr_interactions, sl3::Lrnr_glmnet$new()
    )
    lrnr_enet <- sl3::make_learner(
      sl3::Pipeline, lrnr_interactions, sl3::Lrnr_glmnet$new(alpha = 0.5)
    )
    lrnr_rf <- sl3::Lrnr_ranger$new()
    lrnr_mean <- sl3::Lrnr_mean$new()

    # assemble learners
    learner_library <- sl3::make_learner(
      sl3::Stack, lrnr_glm, lrnr_lasso, lrnr_enet, lrnr_xgboost, lrnr_rf,
      lrnr_mean
    )

    # define the metalearner
    meta_learner <- sl3::make_learner(
      sl3::Lrnr_solnp,
      loss_function = sl3::loss_loglik_binomial,
      learner_function = sl3::metalearner_logistic_binomial
    )

    # initialize the SuperLearner
    cond_censor_haz_super_learner <- sl3::Lrnr_sl$new(
      learners = learner_library,
      metalearner = meta_learner
    )

  }

  # estimate the conditional hazard function
  cond_cens_haz_fit <- cond_censor_haz_super_learner$train(cens_haz_train_task)


  # estimate variable important parameters #####################################

  # estimate the survival probs under treatment and control at t0
  times <- c(valid_data$time, train_data$time) %>% unique %>% sort()
  n_times <- length(times)
  surv_valid_data <- valid_data %>%
    dplyr::select(-dplyr::all_of(c(failure, censor, "time"))) %>%
    dplyr::distinct()
  surv_valid_data <- lapply(
    seq_len(nrow(valid_data)),
    function(idx) {
      replicate(n_times, valid_data[idx, ], simplify = FALSE) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(time = times)
    }
  ) %>% dplyr::bind_rows()

  # create the treatment dataset and prediction task
  valid_data_treat <- surv_valid_data
  valid_data_treat[treatment] <- names(propensity_score_ls)[1]
  valid_data_treat[treatment] <- factor(
    dplyr::pull(valid_data_treat, treatment),
    levels = names(propensity_score_ls)
  )
  pred_task_treat <- sl3::make_sl3_Task(
    data = valid_data_treat,
    covariates = covar_names
  )

  # create the control dataset and prediction task
  valid_data_cont <- surv_valid_data
  valid_data_cont[treatment] <- names(propensity_score_ls)[2]
  valid_data_cont[treatment] <- factor(
    dplyr::pull(valid_data_cont, treatment), levels = names(propensity_score_ls)
  )
  pred_task_cont <- sl3::make_sl3_Task(
    data = valid_data_cont,
    covariates = covar_names
  )

  # estimate conditional survival for all at max time (one est per row)
  surv_valid_data$cond_surv_haz_control <- cond_surv_haz_fit$predict(
    pred_task_cont
  )
  surv_valid_data$cond_surv_haz_treat <- cond_surv_haz_fit$predict(
    pred_task_treat
  )
  surv_valid_data <- surv_valid_data %>%
    dplyr::group_by(.data$observation_id) %>%
    dplyr::summarise(
      surv_t0_treat = prod(1 - cond_surv_haz_treat),
      surv_t0_cont = prod(1 - cond_surv_haz_control)
    ) %>%
    dplyr::ungroup()

  # construct cumulative portion of empirical efficient influence curve

  # add the survival probs at t0
  valid_data <- valid_data %>% left_join(surv_valid_data, by = "observation_id")

  # predict cond hazards at each time
  pred_task_valid <- sl3::make_sl3_Task(
    data = valid_data,
    covariates = covar_names
  )
  valid_data$cond_surv_haz <- cond_surv_haz_fit$predict(pred_task_valid)
  valid_data$cond_cens_haz <- cond_cens_haz_fit$predict(pred_task_valid)

  # compute the ajdusted differences in survival probabilities
  surv_diff_df <- valid_data %>%
    dplyr::group_by(.data$observation_id) %>%
    dplyr::mutate(
      surv = cumprod(1 - cond_surv_haz),
      surv_cens = cumprod(1 - cond_cens_haz),
      surv_cens_prev = dplyr::lag(surv_cens),
      surv_cens_prev = dplyr::if_else(is.na(surv_cens_prev), 1, surv_cens_prev),
      h1 = dplyr::if_else(
        !!rlang::sym(treatment) != names(propensity_score_ls)[1], 0,
        -surv_t0_treat/(propensity_score_ls[[1]] * surv_cens_prev * surv)
      ),
      h0 = dplyr::if_else(
        !!rlang::sym(treatment) != names(propensity_score_ls)[2], 0,
        -surv_t0_cont/(propensity_score_ls[[2]] * surv_cens_prev * surv)
      ),
      haz_prod = ((!!rlang::sym(failure) == 1) - .data$cond_surv_haz),
      inner_sum_t = (.data$h1 - .data$h0) * .data$haz_prod
    ) %>%
    dplyr::summarise(
      cumult_eif = sum(.data$inner_sum_t),
      adj_surv_diff = .data$cumult_eif + .data$surv_t0_treat -
        .data$surv_t0_cont
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(.data$adj_surv_diff) %>%
    dplyr::distinct()

  # compute variable importance parameter for all biomarkers

  # compute empirical efficient influence curves for all biomarkers


}