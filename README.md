# mtp
An R package providing multiple testing procedures for familywise error and false exceedance rates control.

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
When performing the test of a null hypothesis
against an alternative hypothesis, the probability of
committing an error when rejecting the
null hypothesis, is bounded by \eqn{\alpha},
which is conventionally set at \eqn{5\%}.

When testing multiple null hypotheses,
the risk of committing errors increases with
the number of tests and is practically guaranteed
when the number of tests is large.

Statistical control is exerted by capping
the probability of producing Type I errors
(in the sense of the chosen criterion).

Multiple testing methods implemented in this package
allow to control the overall amount of incorrect rejections,
and to improve the balance between correctly
and incorrectly rejected hypotheses.
