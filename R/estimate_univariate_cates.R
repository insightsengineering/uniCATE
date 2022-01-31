#' Estimate the parameters and empirical efficient influence functions
#'
#' \code{estimate_univariate_cates()} estimates the slope of the univariate
#'   conditional average treatment effect's linear approximation for continuous
#'   and binary outcomes. First, the data is pre-processed and split into
#'   \code{v_folds} folds. Then, the conditional outcome regression over all the
#'   training sets are estimated using the \code{super_learner}. The estimated
#'   outcome regressions are then used to estimate the expected difference in
#'   potential outcomes across treatment conditions for each individual in their
#'   respective validation set. The variable importance parameters are then
#'   estimated in the validation sets, along with their efficient influence
#'   functions.
#'
#' @param data A \code{tibble} object containing the outcome variable, treatment
#'   indicator, and covariates. The treatment variable is a factor whose levels
#'   correspond to the names of the \code{propensity_score_ls} argument.
#' @param outcome A \code{character} defining the name of the outcome variable
#'   in \code{data}. The outcome must be a continuous or a binary factor
#'   variable.
#' @param treatment A \code{character} indicating the name of the binary
#'   treatment variable in \code{data}.
#' @param biomarkers A \code{character} vector listing the pre-treatment
#'   biomarkers variables in \code{data}. \code{biomarkers} must be a subset of
#'   \code{covariates}.
#' @param super_learner A \code{\link[sl3:Lrnr_sl]{Lrnr_sl}} object. If set
#'   to \code{NULL}, a default SuperLearner is used. If the outcome variable is
#'   continuous, the default's library of base learners is made up of a
#'   linear model, penalized linear models (LASSO and elasticnet), a spline
#'   regression, XGBoost, a Random Forest, and the mean model. When the outcome
#'   variable is binary, the base learner library consists of (penalized)
#'   logistic regression models, XGBoost, a Random Forests, and the mean model.
#'   The type of outcome is automatically detected.
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
estimate_univariate_cates <- function(data,
                                      outcome,
                                      treatment,
                                      biomarkers,
                                      super_learner,
                                      propensity_score_ls,
                                      v_folds,
                                      parallel) {

  # split the data into folds
  folds <- origami::make_folds(
    data,
    fold_fun = origami::folds_vfold,
    V = v_folds
  )

  # assert that the propensity score list's elements are named after the
  # treatment assignment variable, and that they sum to 1
  assertthat::assert_that(
    all(names(propensity_score_ls) %in% levels(dplyr::pull(data, treatment))),
    msg = "treatment group names and propensity_score_ls names do not match"
  )
  assertthat::assert_that(
    sum(unlist(propensity_score_ls)) == 1,
    msg = "propensity_score_ls values do not sum to one"
  )

  # assert that the super_learner is either a NULL or an SL3 SuperLearner object
  assertthat::assert_that(
    is.null(super_learner) |
      identical(class(super_learner), c("Lrnr_sl", "Lrnr_base", "R6")),
    msg = "super_learner must be NULL or an sl3::Lrnr_sl object"
  )

  # identify the outcome type (binary or continuous)
  if (sum(is.element(unique(data[[outcome]]), c(0, 1))) == 2) {
    outcome_type <- "binomial"
  } else {
    outcome_type <- "continuous"
  }

  # compute the holdout estimated potential outcome differences
  hold_out_calculations <- origami::cross_validate(
    cv_fun = hold_out_calculation,
    folds = folds,
    data = data,
    outcome = outcome,
    treatment = treatment,
    biomarkers = biomarkers,
    super_learner = super_learner,
    propensity_score_ls = propensity_score_ls,
    outcome_type = outcome_type,
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
#' \code{hold_out_calculation()} computes all the objects required to estimate
#'   the variable importance parameter and compute their empirical efficient
#'   influence curve over the validation set, assuming continuous or binary
#'   outcomes. It begins by estimated the conditional outcome regression on the
#'   training set using the \code{super_learner}. Next, the difference in
#'   potential outcomes is predicted for the observations in the validation set.
#'   Finally, these predicted outcomes are used to estimate the variable
#'   importance parameters and efficient influence functions.
#'
#' @param fold A \code{fold} object (from \code{\link[origami]{make_folds}()}).
#' @param data A \code{tibble} object containing the outcome variable, treatment
#'   indicator, and covariates. The treatment variable is a factor whose levels
#'   correspond to the names of the \code{propensity_score_ls} argument.
#' @param outcome A \code{character} defining the name of the outcome variable
#'   in \code{data}. The outcome must be a continuous or a binary factor
#'   variable.
#' @param treatment A \code{character} indicating the name of the binary
#'   treatment variable in \code{data}.
#' @param biomarkers A \code{character} vector listing the pre-treatment
#'   biomarkers variables in \code{data}. \code{biomarkers} must be a subset of
#'   \code{covariates}.
#' @param super_learner A \code{\link[sl3:Lrnr_sl]{Lrnr_sl}} object. If set
#'   to \code{NULL}, a default SuperLearner is used. If the outcome variable is
#'   continuous, the default's library of base learners is made up of a
#'   linear model, penalized linear models (LASSO and elasticnet), a spline
#'   regression, XGBoost, a Random Forest, and the mean model. When the outcome
#'   variable is binary, the base learner library consists of (penalized)
#'   logistic regression models, XGBoost, a Random Forests, and the mean model.
#'   The type of outcome is automatically detected.
#' @param propensity_score_ls A named \code{numeric} \code{list} providing the
#'   propensity scores for the treatment conditions. The first element of the
#'   list should correspond to the "treatment" condition, and the second to the
#'   "control" condition, whatever their names may be.
#' @param outcome_type A \code{character} indicating the type of outcome.
#'   Currently supported outcomes are \code{"continuous"} and
#'   \code{"binomial"}. Here, \code{"binomial"} is used for binary outcomes.
#'
#' @return A \code{list} made up of two objects. The first is the \code{numeric}
#'   vector of biomarker variable importance estimates. The second is the
#'   \code{tibble} of the empirical efficient influence functions for each
#'   biomarker.
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
hold_out_calculation <- function(fold, data, outcome, treatment, biomarkers, super_learner,
                                 propensity_score_ls, outcome_type) {

  # define the training and testing set
  train_data <- origami::training(data)
  valid_data <- origami::validation(data)

  # construct the SuperLearner task for the outcome regression
  covar_names <- colnames(train_data)
  covar_names <- covar_names[which(!grepl(outcome, covar_names))]
  train_task <- sl3::make_sl3_Task(
    data = train_data,
    covariates = covar_names,
    outcome = outcome,
    outcome_type = outcome_type
  )

  # should the built-in SuperLearner be used to estimate the outcome regression?
  if (is.null(super_learner)) {

    # define the treatment-covariate interaction for certain learners
    covars <- covar_names[which(!grepl(treatment, covar_names))]
    interactions <- lapply(covars, function(w) c(w, treatment))
    lrnr_interactions <- sl3::Lrnr_define_interactions$new(interactions)

    if (outcome_type == "continuous") {

      # define the base learners for continuous outcomes
      lrnr_glm <- sl3::make_learner(
        sl3::Pipeline, lrnr_interactions, sl3::Lrnr_glm_fast$new()
      )
      lrnr_lasso <- sl3::make_learner(
        sl3::Pipeline, lrnr_interactions, sl3::Lrnr_glmnet$new()
      )
      lrnr_enet <- sl3::make_learner(
        sl3::Pipeline, lrnr_interactions, sl3::Lrnr_glmnet$new(alpha = 0.5)
      )
      lrnr_spline <- sl3::make_learner(
        sl3::Pipeline, lrnr_interactions, sl3::Lrnr_polspline$new()
      )
      lrnr_xgboost <- sl3::Lrnr_xgboost$new()
      lrnr_rf <- sl3::Lrnr_ranger$new()
      lrnr_mean <- sl3::Lrnr_mean$new()

      # assemble learners
      learner_library <- sl3::make_learner(
        sl3::Stack, lrnr_glm, lrnr_spline, lrnr_lasso, lrnr_enet, lrnr_xgboost,
        lrnr_rf, lrnr_mean
      )

      # define the metalearner
      meta_learner <- sl3::make_learner(
        sl3::Lrnr_solnp,
        loss_function = sl3::loss_squared_error,
        learner_function = sl3::metalearner_linear
      )
    } else {

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
      lrnr_xgboost <- sl3::Lrnr_xgboost$new()
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
    }

    # intialize the SuperLearner
    super_learner <- sl3::Lrnr_sl$new(
      learners = learner_library,
      metalearner = meta_learner
    )
  }

  # train the SuperLearner on the training data
  sl_fit <- super_learner$train(train_task)

  # estimate the outcome of each observation for each treatment level
  valid_data_treat <- valid_data
  valid_data_treat[treatment] <- names(propensity_score_ls)[1]
  valid_data_treat[treatment] <- factor(
    dplyr::pull(valid_data_treat, treatment),
    levels = names(propensity_score_ls)
  )
  pred_task_treat <- sl3::make_sl3_Task(
    data = valid_data_treat,
    covariates = covar_names,
    outcome = outcome
  )
  valid_data$Y_hat_treat <- sl_fit$predict(pred_task_treat)

  valid_data_cont <- valid_data
  valid_data_cont[treatment] <- names(propensity_score_ls)[2]
  valid_data_cont[treatment] <- factor(
    dplyr::pull(valid_data_cont, treatment),
    levels = names(propensity_score_ls)
  )
  pred_task_cont <- sl3::make_sl3_Task(
    data = valid_data_cont,
    covariates = covar_names,
    outcome = outcome
  )
  valid_data$Y_hat_cont <- sl_fit$predict(pred_task_cont)

  # apply the doubly-robust A-IPTW transform to each observation for each
  # treatment level in valid_data
  valid_data <- apply_aiptw_transform(
    valid_data, outcome, treatment, propensity_score_ls, outcome_type
  )

  # return the difference of the AIPTW transformed treatment and control
  # estimates
  valid_data <- valid_data %>%
    dplyr::mutate(Y_diff = .data$Y_aiptw_treat - .data$Y_aiptw_cont) %>%
    dplyr::select(.data$Y_diff, dplyr::all_of(biomarkers)) %>%
    tibble::as_tibble(.name_repair = "minimal") %>%
    scale(center = TRUE, scale = FALSE) %>%
    as.data.frame()

  # compute the coef estimates, and the influence functions of each observation
  # for each biomarker
  Y_diff <- as.vector(
    base::scale(valid_data$Y_diff, center = TRUE, scale = FALSE)
  )
  valid_data <- valid_data %>% dplyr::select(-Y_diff)
  coefs_and_ic_ls <- valid_data %>%
    purrr::map(
      function(bio) {

        # center the biomarker measurements
        bio <- as.vector(base::scale(bio, center = TRUE, scale = FALSE))

        # estimate the best linear approximation using the estimating equation
        # formula
        bio_coef <- sum(Y_diff * bio) / sum(bio^2)

        # compute the unscaled empirical IC of each observation
        unsc_inf_curves <- (Y_diff - bio_coef * bio) * bio

        # return the beta coefficients and the influence curves
        return(list(
          "bio_coef" = bio_coef,
          "unsc_inf_curves" = unsc_inf_curves
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
  unsc_ic_mat <- lapply(
    seq_len(
      length(coefs_and_ic_ls)
    ),
    function(idx) coefs_and_ic_ls[[idx]]$unsc_inf_curves
  )
  unsc_ic_mat <- do.call(cbind, unsc_ic_mat)

  # scale the un-scaled influence curves by the training set's biomarkers' vars
  train_var_bio_vec <- train_data %>%
    dplyr::select(dplyr::all_of(biomarkers)) %>%
    as.matrix() %>%
    matrixStats::colVars()
  ic_df <- tcrossprod(unsc_ic_mat, diag(1 / train_var_bio_vec)) %>%
    dplyr::as_tibble(.name_repair = "minimal")
  colnames(ic_df) <- biomarkers


  # return the vector of coefficients and the table of empirical IC values
  return(
    list(
      "beta_coefs" = beta_coefs,
      "ic_df" = ic_df
    )
  )
}

################################################################################

#' Apply an AIPTW Transform
#'
#' \code{apply_aiptw_transform()} applies the augmented inverse-probability
#'   treatment wright (AIPTW) transform to the predicted potential outcomes
#'   generated by the SuperLearner over the training set.
#'
#' @param data A \code{tibble} object containing the outcome variable, treatment
#'   indicator, and covariates. The treatment variable is a factor whose levels
#'   correspond to the names of the \code{propensity_score_ls} argument.
#' @param outcome A \code{character} defining the name of the outcome variable
#'   in \code{data}. The outcome must be a continuous or a binary factor
#'   variable.
#' @param treatment A \code{character} indicating the name of the binary
#'   treatment variable in \code{data}.
#' @param propensity_score_ls A named \code{numeric} \code{list} providing the
#'   propensity scores for the treatment conditions. The first element of the
#'   list should correspond to the "treatment" condition, and the second to the
#'   "control" condition, whatever their names may be.
#' @param outcome_type A \code{character} indicating the type of outcome.
#'   Currently supported outcomes are \code{"continuous"} and
#'   \code{"binomial"}. Here, \code{"binomial"} is used for binary outcomes.
#'
#' @return A \code{tibble} containing the AIPTW transformed predicted outcomes
#'   under treatment and control conditions of all observations in \code{data}.
#'
#' @importFrom dplyr mutate if_else .data
#' @importFrom rlang !! sym
#' @importFrom magrittr %>%
#'
#' @keywords internal
apply_aiptw_transform <- function(data, outcome, treatment, propensity_score_ls, outcome_type) {
  data %>%
    dplyr::mutate(
      Y_aiptw_treat = dplyr::if_else(
        !!rlang::sym(treatment) == names(propensity_score_ls)[1],
        .data$Y_hat_treat +
          (1 / propensity_score_ls[[1]]) *
            (!!rlang::sym(outcome) - .data$Y_hat_treat),
        .data$Y_hat_treat
      ),
      Y_aiptw_treat = dplyr::if_else(
        outcome_type == "binomial" & .data$Y_aiptw_treat < 0,
        0,
        .data$Y_aiptw_treat
      ),
      Y_aiptw_treat = dplyr::if_else(
        outcome_type == "binomial" & .data$Y_aiptw_treat > 1,
        1,
        .data$Y_aiptw_treat
      ),
      Y_aiptw_cont = dplyr::if_else(
        !!rlang::sym(treatment) == names(propensity_score_ls)[2],
        .data$Y_hat_cont +
          (1 / propensity_score_ls[[2]]) *
            (!!rlang::sym(outcome) - .data$Y_hat_cont),
        .data$Y_hat_cont
      ),
      Y_aiptw_cont = dplyr::if_else(
        outcome_type == "binomial" & .data$Y_aiptw_cont < 0,
        0,
        .data$Y_aiptw_cont
      ),
      Y_aiptw_cont = dplyr::if_else(
        outcome_type == "binomial" & .data$Y_aiptw_cont > 1,
        1,
        .data$Y_aiptw_cont
      )
    )
}
