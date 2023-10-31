# Overview
The mtp package provides functions to apply
various forms of statistical control on the
occurrence of Type I errors (incorrect rejections)
when conducting multiple statistical tests.

The package offers to use different error criteria,
such as the (generalized) familywise error rate and
the false rejection exceedance, as well as different
methods, e.g. step-down or step-up procedures.

## Installation

To install **mtp** from R:

```R
# Install/load R package devtools
install.packages("devtools")
library(devtools)

# Install/load R package beam from github
install_github("gleday/mtp")
library(mtp)
```

## Introduction
When performing a test, the probability of
committing an error when rejecting the
null hypothesis, is bounded by \eqn{\alpha},
which is conventionally set at \eqn{5\%}.

When testing multiple null hypotheses,
the risk of committing errors increases with
the number of tests and is practically guaranteed
when the number of tests is large.

Statistical control is exerted by capping
the probability of producing Type I errors,
(in the sense of a criterion).

Multiple testing procedures make possible this
control by adjusting the significance level of
each test (possibly using different methods).

Overall, criteria and methods for multiple
testing allow to:
- improve the balance between the overall amount of
incorrectly and correctly rejected hypotheses (i.e.
between overall Type I errors and statistical power).
- provide statistical guarantees on the set of
rejected null hypotheses, e.g. that it is unlikely
to exceed 5 errors.
