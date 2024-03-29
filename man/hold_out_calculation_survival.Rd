% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_univariate_survival_cates.R
\name{hold_out_calculation_survival}
\alias{hold_out_calculation_survival}
\title{Compute Validation Set Objects}
\usage{
hold_out_calculation_survival(
  fold,
  long_data,
  event,
  censor,
  treatment,
  biomarkers,
  cond_surv_haz_super_learner,
  cond_censor_haz_super_learner,
  propensity_score_ls
)
}
\arguments{
\item{fold}{A \code{fold} object (from \code{\link[origami]{make_folds}()}).}

\item{long_data}{A \code{tibble} object containing the longitudinal data
output by \code{prep_long_data()}.}

\item{event}{A \code{character} defining the name of the binary variable in
the \code{data} argument that indicates whether an event occurred.
Observations can have an event or be censored, but not both.}

\item{censor}{A \code{character} defining the name of the binary variable in
the \code{data} argument that indicates a right-censoring event.
Observations can have an event or be censored, but not both.}

\item{treatment}{A \code{character} indicating the name of the binary
treatment variable in \code{data}.}

\item{biomarkers}{A \code{character} vector listing the biomarkers of
interest in \code{data}.}

\item{cond_surv_haz_super_learner}{A \code{\link[sl3:Lrnr_sl]{Lrnr_sl}}
object used to estimate the conditional event hazard model. If set to
\code{NULL}, the default, an elastic net regression is used instead. It is
best to use this default behaviour when analyzing small datasets.}

\item{cond_censor_haz_super_learner}{A \code{\link[sl3:Lrnr_sl]{Lrnr_sl}}
object used to estimate the conditional censoring hazard model. If set to
\code{NULL}, the default, an elastic net regression is used instead. It is
best to use this default behaviour when analyzing small datasets.}

\item{propensity_score_ls}{A named \code{numeric} \code{list} providing the
propensity scores for the treatment conditions. The first element of the
list should correspond to the "treatment" condition, and the second to the
"control" condition, whatever their names may be.}
}
\value{
A \code{list} made up of two objects. The first is the \code{numeric}
vector of biomarker variable importance estimates. The second is the
\code{tibble} of the empirical efficient influence functions for each
biomarker.
}
\description{
\code{hold_out_calculation_survival()} computes all the objects required to
estimate the univariate CATEs for survival outcomes over the validation
set. It begins by estimating the conditional survival and censoring
hazard functions on the training set using the
\code{cond_surv_haz_super_learner} and
\code{cond_censor_haz_super_learner}, respectively. Next, the difference in
survival probabilities is predicted for the observations in the validation
set. Finally, these predicted outcomes are used to estimate the variable
impotance parameters of each biomarker. Their respective empirical
efficient influence functions are computed as well.
}
\keyword{internal}
