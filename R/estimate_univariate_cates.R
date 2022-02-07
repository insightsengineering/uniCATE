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
#' @param super_learner A \code{\link[sl3:Lrnr_sl]{Lrnr_sl}} object used to
#'   estimate the conditional outcome regression. If set to \code{NULL}, the
#'   default, an elastic net regression is used instead. It is best to use this
#'   default behaviour when analyzing small datasets.
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
#' @param super_learner A \code{\link[sl3:Lrnr_sl]{Lrnr_sl}} object used to
#'   estimate the conditional outcome regression. If set to \code{NULL}, the
#'   default, an elastic net regression is used instead. It is best to use this
#'   default behaviour when analyzing small datasets.
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
#' @importFrom stats as.formula model.matrix predict
#' @importFrom purrr map
#' @importFrom matrixStats colVars
#' @importFrom glmnet cv.glmnet
#' @import sl3
#'
#' @keywords internal
hold_out_calculation <- function(fold, data, outcome, treatment, biomarkers,
                                 super_learner, propensity_score_ls,
                                 outcome_type) {

  # define the training and testing set
  train_data <- origami::training(data)
  valid_data <- origami::validation(data)

  # construct the modified validation data with fixed treatment assignments
  valid_data_treat <- valid_data
  valid_data_treat[treatment] <- names(propensity_score_ls)[1]
  valid_data_treat[treatment] <- factor(
    dplyr::pull(valid_data_treat, treatment),
    levels = names(propensity_score_ls)
  )
  valid_data_cont <- valid_data
  valid_data_cont[treatment] <- names(propensity_score_ls)[2]
  valid_data_cont[treatment] <- factor(
    dplyr::pull(valid_data_cont, treatment),
    levels = names(propensity_score_ls)
  )

  # should the built-in learner be used to estimate the outcome regression?
  if (is.null(super_learner)) {

    # specify the outcome
    y_train <- train_data[[outcome]]
    covar_names <- colnames(train_data)
    covar_names <- covar_names[which(!grepl(paste0(outcome, "|", treatment),
                                            covar_names))]
    formula_string <- paste(
      "~", treatment, "+", paste(covar_names, collapse = " + "), "+",
      paste(paste(covar_names, treatment, sep = ":"), collapse = " + ")
    )
    x_train <- stats::model.matrix(
      stats::as.formula(formula_string),
      data = train_data
    )

    if (outcome_type == "continuous") {

      # train the elastic net
      glmnet_fit <- cv.glmnet(
        x = x_train,
        y = y_train,
        nfolds = 5,
        alpha = 0.5,
        family = "gaussian"
      )
    } else {

      # train the elastic net
      glmnet_fit <- cv.glmnet(
        x = x_train,
        y = y_train,
        nfolds = 5,
        alpha = 0.5,
        family = "binomial"
      )
    }

    # predict the outcomes in the validation sets with fixed treatments
    x_data_treat <- stats::model.matrix(
      stats::as.formula(formula_string),
      data = valid_data_treat
    )
    valid_data$y_hat_treat <- stats::predict(
      glmnet_fit,
      newx = x_data_treat,
      s = "lambda.min",
      type = "response"
    )
    x_data_cont <- stats::model.matrix(
      stats::as.formula(formula_string),
      data = valid_data_cont
    )
    valid_data$y_hat_cont <- stats::predict(
      glmnet_fit,
      newx = x_data_cont,
      s = "lambda.min",
      type = "response"
    )
  } else {

    # construct the SuperLearner task for the outcome regression
    covar_names <- colnames(train_data)
    covar_names <- covar_names[which(!grepl(outcome, covar_names))]
    train_task <- sl3::make_sl3_Task(
      data = train_data,
      covariates = covar_names,
      outcome = outcome,
      outcome_type = outcome_type
    )

    # train the SuperLearner on the training data
    sl_fit <- super_learner$train(train_task)

    # estimate the outcome of each observation for each treatment level
    pred_task_treat <- sl3::make_sl3_Task(
      data = valid_data_treat,
      covariates = covar_names,
      outcome = outcome
    )
    valid_data$y_hat_treat <- sl_fit$predict(pred_task_treat)
    pred_task_cont <- sl3::make_sl3_Task(
      data = valid_data_cont,
      covariates = covar_names,
      outcome = outcome
    )
    valid_data$y_hat_cont <- sl_fit$predict(pred_task_cont)
  }

  # apply the doubly-robust A-IPTW transform to each observation for each
  # treatment level in valid_data
  valid_data <- apply_aiptw_transform(
    valid_data, outcome, treatment, propensity_score_ls
  )

  # return the difference of the AIPTW transformed treatment and control
  # estimates
  valid_data <- valid_data %>%
    dplyr::mutate(y_diff = .data$y_aiptw_treat - .data$y_aiptw_cont) %>%
    dplyr::select(.data$y_diff, dplyr::all_of(biomarkers)) %>%
    tibble::as_tibble(.name_repair = "minimal") %>%
    scale(center = TRUE, scale = FALSE) %>%
    as.data.frame()

  # compute the coef estimates, and the influence functions of each observation
  # for each biomarker
  y_diff <- as.vector(
    base::scale(valid_data$y_diff, center = TRUE, scale = FALSE)
  )
  valid_data <- valid_data %>% dplyr::select(-y_diff)
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
          bio_coef <- mean(y_diff * bio) / var_bio

          # compute the empirical IC of each observation
          inf_curves <- ((y_diff - bio_coef * bio) * bio) / var_bio

        } else {

          # set the coefficient to zero, and make the influence curve enormous
          bio_coef <- 0
          inf_curves <- rep(1000000, length(y_diff))
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
#'
#' @return A \code{tibble} containing the AIPTW transformed predicted outcomes
#'   under treatment and control conditions of all observations in \code{data}.
#'
#' @importFrom dplyr mutate if_else .data
#' @importFrom rlang !! sym
#' @importFrom magrittr %>%
#'
#' @keywords internal
apply_aiptw_transform <- function(data,
                                  outcome,
                                  treatment,
                                  propensity_score_ls) {
  data %>%
    dplyr::mutate(
      y_aiptw_treat = dplyr::if_else(
        !!rlang::sym(treatment) == names(propensity_score_ls)[1],
        .data$y_hat_treat +
          (1 / propensity_score_ls[[1]]) *
            (!!rlang::sym(outcome) - .data$y_hat_treat),
        .data$y_hat_treat
      ),
      y_aiptw_cont = dplyr::if_else(
        !!rlang::sym(treatment) == names(propensity_score_ls)[2],
        .data$y_hat_cont +
          (1 / propensity_score_ls[[2]]) *
            (!!rlang::sym(outcome) - .data$y_hat_cont),
        .data$y_hat_cont
      )
    )
}
