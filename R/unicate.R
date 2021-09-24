#' Univariate Conditional Average Treatment Effect Estimation
#'
#' \code{unicate()} estimates the best linear approximation of the conditional
#'   average treatment effect for each variable specified in the
#'   \code{biomarkers} argument. The estimated coefficients are then used to
#'   test whether their corresponding biomarkers are treatment-effect modifiers.
#'   Valid statistical inference is performed by employing the biomarkers'
#'   cross-validated empirical influence curves to estimate these coefficients'
#'   standard errors.
#'
#' @param data A \code{data.frame} or \code{tibble} object containing the
#'   outcome variable, treatment indicator, and covariates. Note that the
#'   biomarkers must be a subset of the covariates.
#' @param outcome A \code{character} defining the name of the outcome variable
#'   in \code{data}. The outcome must be a continuous or a binary factor
#'   variable.
#' @param treatment A \code{character} indicating the name of the binary
#'   treatment variable in \code{data}.
#' @param covariates A \code{character} vector listing the covariates in
#'   \code{data}.
#' @param biomarkers A \code{character} vector listing the biomarkers of
#'   interest in \code{data}. \code{biomarkers} must be a subset of
#'   \code{covariates}.
#' @param super_learner A \code{\link[sl3:Lrnr_sl]{SuperLearner}} object. If set
#'   to \code{NULL}, the default SuperLearner is used. The default's library
#'   consists of a linear model, penalized linear models (LASSO and elasticnet),
#'   a spline regression, a General Additive Model, XGBoost, a Random Forest,
#'   and the mean model.
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
#' @return A \code{tibble} of the biomarkers. Each row contains the estimated
#'   univariate linear model coefficient of the biomarker regressed against the
#'   predicted difference in potential outcomes, along with the standard error,
#'   test statistic, and (adjusted) p-values. The table is ordered by
#'   significance.
#'
#' @export
unicate <- function(
  data,
  outcome,
  treatment,
  covariates,
  biomarkers,
  super_learner = NULL,
  propensity_score_ls,
  v_folds = 5L,
  parallel = FALSE
) {

  # assess the data quality and formatting, and prepare it for analysis
  data <- prep_data(data, outcome, treatment, covariates, biomarkers)

  # compute CV coefficients and CV influence curves
  cv_ls <- estimate_univariate_cates(
    data, outcome, treatment, biomarkers, super_learner, propensity_score_ls,
    v_folds, parallel
  )

  # compute the table of biomarker coefficients and standard errors
  biomarkers_tbl <- compute_coefs_tbl(cv_ls)

  # perform inference using the estimated coefficients and standard errors
  biomarkers_tbl <- perform_inference(biomarkers_tbl)

  return(biomarkers_tbl)
}
