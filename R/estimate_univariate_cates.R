#' Estimate the Univariate CATEs
#'
#' \code{estimate_univariate_cates()} estimates the univariate Conditional
#' Average Treatment Effect of each biomarker on the predicted potential
#' outcomes using cross-validation. First, the data is pre-processed data is
#' split into \code{v_folds}. Then, the conditional outcome regression over all
#' the training sets are estimated using the \code{super_learner}. The estimated
#' outcome regressions are then used to estimate the expected difference in
#' potential outcomes across treatment groups for each individual in their
#' respective validation sets. Univariate linear regressions are then regressed
#' over the predicted potential outcomes in the validation set, and the
#' influence curves of all observations is computed as well. \code{tibbles} of
#' biomarker coefficients and influence curves are returned for each validation
#' set.
#'
##' @param data A \code{tibble} object containing the outcome variable, treatment
#'   indicator, and covariates. The treatment variable is a factor whose levels
#'   correspond to the names of the \code{propensity_score_ls} argument.
#' @param outcome A \code{character} defining the name of the outcome variable
#'   in \code{data}.
#' @param treatment A \code{character} indicating the name of the binary
#'   treatment variable in \code{data}.
#' @param biomarkers A \code{character} vector listing the biomarkers of
#'   interest in \code{data}.
#' @param super_learner A \code{\link[sl3:Lrnr_sl]{SuperLearner}} object. If set
#'   to \code{NULL}, the default SuperLearner is used.
#' @param propensity_score_ls A named \code{numeric} \code{list} providing the
#'   propensity scores for the treatment levels. Note that the first element
#'   of the list should correspond to the "treatment" group, and the second the
#'   "control" group, whatever their names may be.
#' @param v_folds A \code{numeric} indicating the number of folds to use for
#'   V-fold cross-validation.
#' @param parallel A \code{logical} determining whether to use
#'   \code{\link[origami:cross_validate]{origami}}'s built-in parallelization
#'   capabilities.
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
estimate_univariate_cates <- function(
  data,
  outcome,
  treatment,
  biomarkers,
  super_learner,
  propensity_score_ls,
  v_folds,
  parallel
) {

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
#' \code{hold_out_calculation} computes all the objects required to estimated
#' the univariate CATE over the validation set. It begins by estimated the
#' conditional outcome regression on the training set using the
#' \code{super_learner}. Next, the difference in potential outcomes is
#' predicted for the observations in the validation set. Finally, these
#' predicted outcomes are used to estimate the unvariate regression coefficients
#' of each biomarker, and their respective empirical influence curves are
#' computed as well.
#'
#' @param fold A \code{fold} object (from \code{\link[origami]{make_folds}()}).
#' @param data A \code{tibble} object containing the outcome variable, treatment
#'   indicator, and covariates. The treatment variable is a factor whose levels
#'   correspond to the names of the \code{propensity_score_ls} argument.
#' @param outcome A \code{character} defining the name of the outcome variable
#'   in \code{data}.
#' @param treatment A \code{character} indicating the name of the binary
#'   treatment variable in \code{data}.
#' @param biomarkers A \code{character} vector listing the biomarkers of
#'   interest in \code{data}.
#' @param super_learner A \code{\link[sl3:Lrnr_sl]{SuperLearner}} object. If set
#'   to \code{NULL}, the default SuperLearner is used.
#' @param propensity_score_ls A named \code{numeric} \code{list} providing the
#'   propensity scores for the treatment levels. Note that the first element
#'   of the list should correspond to the "treatment" group, and the second the
#'   "control" group, whatever their names may be.
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
hold_out_calculation <- function(
  fold, data, outcome, treatment, biomarkers, super_learner, propensity_score_ls
) {

  # define the training and testing set
  train_data <- origami::training(data)
  valid_data <- origami::validation(data)

  # construct the SuperLearner task for the outcome regression
  covar_names <- colnames(train_data)
  covar_names <- covar_names[which(!grepl(outcome, covar_names))]
  train_task <- sl3::make_sl3_Task(
    data = train_data,
    covariates = covar_names,
    outcome = outcome
  )

  # should the built-in SuperLearner be used to estimate the outcome regression?
  if (is.null(super_learner)) {

    # define the treatment-covariate interaction for certain learners
    covars <- covar_names[which(!grepl(treatment, covar_names))]
    interactions <- lapply(covars, function(w) c(w, treatment))
    lrnr_interactions <- sl3::Lrnr_define_interactions$new(interactions)

    # define the base learners for contiuous outcomes
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
    dplyr::pull(valid_data_cont, treatment), levels = names(propensity_score_ls)
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
    valid_data, outcome, treatment, propensity_score_ls
  )

  # return the difference of the AIPTW transformed treatment and control
  # estimates
  valid_data <- valid_data %>%
    dplyr::mutate(Y_diff = .data$Y_aiptw_treat - .data$Y_aiptw_cont) %>%
    dplyr::select(.data$Y_diff, dplyr::all_of(biomarkers)) %>%
    tibble::as_tibble(.name_repair = "minimal") %>%
    scale(center = TRUE, scale = FALSE) %>%
    as.data.frame

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
  unsc_ic_mat <- lapply(seq_len(
    length(coefs_and_ic_ls)),
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
#' \code{apply_aiptw_transform()} applies an augmented inverse-probability
#' treatment wright (AIPTW) transform to the estimated potential outcomes
#' generated by the SuperLearner over the training set.
#'
#' @param data A \code{tibble} object containing the outcome variable, treatment
#'   indicator, and covariates. The treatment variable is a factor whose levels
#'   correspond to the names of the \code{propensity_score_ls} argument.
#' @param outcome A \code{character} defining the name of the outcome variable
#'   in \code{data}.
#' @param treatment A \code{character} indicating the name of the binary
#'   treatment variable in \code{data}.
#' @param propensity_score_ls A named \code{numeric} \code{list} providing the
#'   propensity scores for the treatment levels. Note that the first element
#'   of the list should correspond to the "treatment" group, and the second the
#'   "control" group, whatever their names may be.
#'
#' @return A \code{tibble} containing the AIPTW transformed predicted outcomes
#'   under treatment and control group assignments of all observations in
#'   \code{data}.
#'
#' @importFrom dplyr mutate if_else .data
#' @importFrom rlang !! sym
#' @importFrom magrittr %>%
#'
#' @keywords internal
apply_aiptw_transform <- function(
  data, outcome, treatment, propensity_score_ls
) {
  data %>%
    dplyr::mutate(
      Y_aiptw_treat = dplyr::if_else(
        !!rlang::sym(treatment) == names(propensity_score_ls)[1],
        .data$Y_hat_treat +
          (1 / propensity_score_ls[[1]]) *
          (!!rlang::sym(outcome) - .data$Y_hat_treat),
        .data$Y_hat_treat
      ),
      Y_aiptw_cont = dplyr::if_else(
        !!rlang::sym(treatment) == names(propensity_score_ls)[2],
        .data$Y_hat_cont +
          (1 / propensity_score_ls[[2]]) *
          (!!rlang::sym(outcome) - .data$Y_hat_cont),
        .data$Y_hat_cont
      )
    )
}
