---
title: "uniCATE: Univariate Conditional Average Treatment Effect Estimation for Predictive Biomarker Discovery in R"
authors:
  - name: Philippe Boileau
    orcid: 0000-0002-4850-2507
    affiliation: "1, 2"
  - name: Mark J. van der Laan
    orcid: 0000-0003-1432-5511
    affiliation: "1, 2, 3"
  - name: Sandrine Dudoit
    orcid: 0000-0002-6069-8629
    affiliation: "1, 2, 3"
  - name: Ning Leng
    affiliation: "4"
affiliations:
  - name: Division of Biostatistics, School of Public Health, University of California, Berkeley
    index: 1
  - name: Center for Computational Biology, University of California, Berkeley
    index: 2
  - name: Department of Statistics, University of California, Berkeley
    index: 3
  - name: Genentech Inc.
    index: 4

date: "\\today"
bibliography: paper.bib
---

## Summary

Predictive biomarkers are biometric measurements delineating patient
populations into groups that draw differing benefits from a given treatment
[@royston2008; @kraus2018]. These biomarkers play an important role in precision
medicine: as indicators of treatment effect modification, they can inform
patient treatment decisions. Predictive biomarkers have other applications as
well. When used in drug target discovery, they may provide insight on the
functioning of novel medications.

The identification of predictive biomarkers through statistical techniques has,
to date, largely been a byproduct of treatment rule estimation. Examples include
the methods of @tian2014 and @chen2017. This tasks consists of learning which
treatment would most benefit a patient given their characteristics. These
characteristics could include, but are not limited to, potentially predictive
biomarkers. When interpretable statistical learning methods are used to estimate
these rules, biomarkers thought to most meaningfully modify the treatment effect
might be labeled predictive. In clinical trials, for example, estimation
procedures based on regularized linear regression that include
treatment-biomarker interactions in the model define a biomarker as predictive
if its treatment interaction coefficient estimate is non-zero.

While this approach to predictive biomarker discovery works well in traditional
asymptotic settings --- assuming that the relationship between the outcome, the
treatment, and the biomarkers admits a linear form ---, it is not so when the
number of biomarkers is similar to or larger than the number of patients in a
study. Popular treatment rule estimation methods do not guarantee false positive
rate control [@tian2014; @chen2017], which can have severe consequences in
practice. In particular, the inclusion of false positive predictive biomarkers
in diagnostic assays decreases their clinical utility, and already limited
resources are wasted on biological validation experiments of duds for drug
target discovery purposes. Patient outcomes are ultimately negatively affected.

We proposed in recent work a causal variable importance parameter that directly
measures individual biomarkers' roles in treatment effect modification
[@boileau2022]. The methodology, dubbed *uniCATE*, relies on semiparametric
theory to perform assumption-lean inference about this parameter in randomized
control trials. We demonstrated that our procedure controls the rate of false
discoveries in high-dimensional clinical trials, allowing for the accurate
identification of predictive biomarkers. This methodology is implemented in the
`uniCATE` package.

## Statement of Need

Modern semiparametric methods rely on complex statistical theory that is beyond
the training afforded to many of the clinical scientists most interested in
their use. These methods also rely on modern machine learning algorithms to
provide assumption-lean inference, further increasing the amount technical
know-how required for one-off implementations. Our method is no different. The
`uniCATE` package provides a rigorously-tested, open-source implementation of
this methodology for the `R` language and environment for statistical computing
[@rlang] that is ready for use by members of the clinical research community.

## Availability

The `uniCATE` package is publicly available on GitHub at
[`insightsengineering/uniCATE`](https://github.com/insightsengineering/uniCATE).
The package's `README` provides installation instructions and a minimal example,
and comprehensive documentation and a tutorial are included in its `pkgdown`
[documentation website](https://insightsengineering.github.io/uniCATE).

## Acknowledgements

PB gratefully acknowledges the support of the Fonds de recherche du Québec -
Nature et technologies and the Natural Sciences and Engineering Research Council
of Canada.

## References
