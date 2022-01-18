#' Univariate Conditional Average Treatment Effect Estimation For Right-Censored Time-To-Event Outcomes
#'
#' \code{sunicate()} estimates a linear approximation of the conditional
#'   average treatment effect of a time-to-event outcome for each variable
#'   specified in the \code{biomarkers} argument. The estimated coefficients are
#'   then used to test whether their corresponding biomarkers are
#'   treatment-effect modifiers. Valid statistical inference is performed by
#'   employing the biomarkers' cross-validated empirical influence curves to
#'   estimate these coefficients' standard errors.
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
#' @param v_folds A \code{numeric} indicating the number of folds to use for
#'   V-fold cross-validation. Defaults to \code{5}.
#' @param parallel A \code{logical} determining whether to use
#'   \code{\link[origami:cross_validate]{origami}}'s built-in parallelization
#'   capabilities. Defaults to \code{FALSE}.
#'
#' @return A \code{tibble} of the biomarkers. Each row contains the estimated
#'   univariate linear model coefficient of the biomarker regressed against the
#'   predicted difference in adjsuted survival probabilities, along with the
#'   standard error, test statistic, and (adjusted) p-values. The table is
#'   ordered by significance.
#'
#' @export
sunicate <- function(
  data,
  failure,
  censor,
  relative_time,
  treatment,
  covariates,
  biomarkers,
  time_cutoff = NULL,
  cond_surv_haz_super_learner = NULL,
  cond_censor_haz_super_learner = NULL,
  propensity_score_ls,
  v_folds = 5L,
  parallel = FALSE
) {

  # assess the data quality and formatting, and prepare it for analysis
  long_data <- prep_long_data(
    data, failure, censor, relative_time, treatment, covariates, biomarkers,
    time_cutoff
  )

  # compute CV coefficients and CV influence curves
  cv_ls <- estimate_univariate_survival_cates(
    long_data, failure, censor, treatment, biomarkers,
    cond_surv_haz_super_learner, cond_censor_haz_super_learner,
    propensity_score_ls, v_folds, parallel
  )

  # compute the table of biomarker coefficients and standard errors
  biomarkers_tbl <- compute_coefs_tbl(cv_ls)

  # perform tests using the estimated coefficients and standard errors
  biomarkers_tbl <- perform_inference(biomarkers_tbl)

  return(biomarkers_tbl)
}
